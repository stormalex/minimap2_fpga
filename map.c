#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "kthread.h"
#include "kvec.h"
#include "kalloc.h"
#include "sdust.h"
#include "mmpriv.h"
#include "bseq.h"
#include "khash.h"

struct mm_tbuf_s {
	void *km;
};

mm_tbuf_t *mm_tbuf_init(void)
{
	mm_tbuf_t *b;
	b = (mm_tbuf_t*)calloc(1, sizeof(mm_tbuf_t));
	if (!(mm_dbg_flag & 1)) b->km = km_init();
	return b;
}

void mm_tbuf_destroy(mm_tbuf_t *b)
{
	if (b == 0) return;
	km_destroy(b->km);
	free(b);
}

static int mm_dust_minier(void *km, int n, mm128_t *a, int l_seq, const char *seq, int sdust_thres)
{
	int n_dreg, j, k, u = 0;
	const uint64_t *dreg;
	sdust_buf_t *sdb;
	if (sdust_thres <= 0) return n;
	sdb = sdust_buf_init(km);
	dreg = sdust_core((const uint8_t*)seq, l_seq, sdust_thres, 64, &n_dreg, sdb);
	for (j = k = 0; j < n; ++j) { // squeeze out minimizers that significantly overlap with LCRs
		int32_t qpos = (uint32_t)a[j].y>>1, span = a[j].x&0xff;
		int32_t s = qpos - (span - 1), e = s + span;
		while (u < n_dreg && (uint32_t)dreg[u] <= s) ++u;
		if (u < n_dreg && dreg[u]>>32 < e) {
			int v, l = 0;
			for (v = u; v < n_dreg && dreg[v]>>32 < e; ++v) { // iterate over LCRs overlapping this minimizer
				int ss = s > dreg[v]>>32? s : dreg[v]>>32;
				int ee = e < (uint32_t)dreg[v]? e : (uint32_t)dreg[v];
				l += ee - ss;
			}
			if (l <= span>>1) a[k++] = a[j]; // keep the minimizer if less than half of it falls in masked region
		} else a[k++] = a[j];
	}
	sdust_buf_destroy(sdb);
	return k; // the new size
}

static void collect_minimizers(void *km, const mm_mapopt_t *opt, const mm_idx_t *mi, int n_segs, const int *qlens, const char **seqs, mm128_v *mv)
{
	int i, j, n, sum = 0;
	mv->n = 0;
	for (i = n = 0; i < n_segs; ++i) {
		mm_sketch(km, seqs[i], qlens[i], mi->w, mi->k, i, mi->flag&MM_I_HPC, mv);
		for (j = n; j < mv->n; ++j)
			mv->a[j].y += sum << 1;
		if (opt->sdust_thres > 0) // mask low-complexity minimizers
			mv->n = n + mm_dust_minier(km, mv->n - n, mv->a + n, qlens[i], seqs[i], opt->sdust_thres);
		sum += qlens[i], n = mv->n;
	}
}

#include "ksort.h"
#define heap_lt(a, b) ((a).x > (b).x)
KSORT_INIT(heap, mm128_t, heap_lt)

typedef struct {
	uint32_t n;
	uint32_t q_pos, q_span;
	uint32_t seg_id:31, is_tandem:1;
	const uint64_t *cr;
} mm_match_t;

static mm_match_t *collect_matches(void *km, int *_n_m, int max_occ, const mm_idx_t *mi, const mm128_v *mv, int64_t *n_a, int *rep_len, int *n_mini_pos, uint64_t **mini_pos)
{
	int i, rep_st = 0, rep_en = 0, n_m;
	mm_match_t *m;
	*n_mini_pos = 0;
	*mini_pos = (uint64_t*)kmalloc(km, mv->n * sizeof(uint64_t));
	m = (mm_match_t*)kmalloc(km, mv->n * sizeof(mm_match_t));
	for (i = n_m = 0, *rep_len = 0, *n_a = 0; i < mv->n; ++i) {
		const uint64_t *cr;
		mm128_t *p = &mv->a[i];
		uint32_t q_pos = (uint32_t)p->y, q_span = p->x & 0xff;
		int t;
		cr = mm_idx_get(mi, p->x>>8, &t);
		if (t >= max_occ) {
			int en = (q_pos >> 1) + 1, st = en - q_span;
			if (st > rep_en) {
				*rep_len += rep_en - rep_st;
				rep_st = st, rep_en = en;
			} else rep_en = en;
		} else {
			mm_match_t *q = &m[n_m++];
			q->q_pos = q_pos, q->q_span = q_span, q->cr = cr, q->n = t, q->seg_id = p->y >> 32;
			q->is_tandem = 0;
			if (i > 0 && p->x>>8 == mv->a[i - 1].x>>8) q->is_tandem = 1;
			if (i < mv->n - 1 && p->x>>8 == mv->a[i + 1].x>>8) q->is_tandem = 1;
			*n_a += q->n;
			(*mini_pos)[(*n_mini_pos)++] = (uint64_t)q_span<<32 | q_pos>>1;
		}
	}
	*rep_len += rep_en - rep_st;
	*_n_m = n_m;
	return m;
}

static inline int skip_seed(int flag, uint64_t r, const mm_match_t *q, const char *qname, unsigned int bid, int qlen, const mm_idx_t *mi, int *is_self)
{
	*is_self = 0;

	if (qname && (flag & (MM_F_NO_DIAG|MM_F_NO_DUAL))) {
		//const mm_idx_seq_t *s = &mi->seq[r>>32];
        int cmp = 0;
#if 1
        //int rankID = mi->rever_rid[r>>32];
		uint32_t rankID = (uint32_t)r&0x1FFFFF;
        int flg = (bid &0x80000000)>>31;
        int val = bid  &0x7fffffff;
        if (val > rankID) cmp = 1;
        else if (val < rankID) cmp = -1;
        else if (val == rankID) {
            if (flg) cmp = 0;
            else cmp = -1;
        }
#endif
		//int cmp2 = strcmp(qname, s->name);
        //if (cmp2 > 0) cmp2 = 1;
        //else if(cmp2 < 0) cmp2 = -1;
        //if (cmp != cmp2) fprintf(stderr, "xxxxxxxxxxxxxx,cmp %d, cmp2 %d, qname %s,rname %s\n", qname, s->name);

        //   cmp = cmp2;

		if ((flag&MM_F_NO_DIAG) && cmp == 0) {
			if ((((r>>22) & 0x1fffff)) == (q->q_pos>>1)) return 1; // avoid the diagnonal anchors
			if ((r&MM_P_STRAND)>>21 == (q->q_pos&1)) *is_self = 1; // this flag is used to avoid spurious extension on self chain
		}
		if ((flag&MM_F_NO_DUAL) && cmp > 0) // all-vs-all mode: map once
			return 1;
	}
	if (flag & (MM_F_FOR_ONLY|MM_F_REV_ONLY)) {
		if ((r&MM_P_STRAND)>>21 == (q->q_pos&1)) { // forward strand
			if (flag & MM_F_REV_ONLY) return 1;
		} else {
			if (flag & MM_F_FOR_ONLY) return 1;
		}
	}
	return 0;
}

static mm128_t *collect_seed_hits_heap(void *km, const mm_mapopt_t *opt, int max_occ, const mm_idx_t *mi, const char *qname, const mm128_v *mv, int qlen, int64_t *n_a, int *rep_len,
								  int *n_mini_pos, uint64_t **mini_pos)
{
	int i, n_m, heap_size = 0;
	int64_t j, n_for = 0, n_rev = 0;
	mm_match_t *m;
	mm128_t *a, *heap;

	m = collect_matches(km, &n_m, max_occ, mi, mv, n_a, rep_len, n_mini_pos, mini_pos);

	heap = (mm128_t*)kmalloc(km, n_m * sizeof(mm128_t));
	a = (mm128_t*)kmalloc(km, *n_a * sizeof(mm128_t));

	for (i = 0, heap_size = 0; i < n_m; ++i) {
		if (m[i].n > 0) {
			heap[heap_size].x = m[i].cr[0];
			heap[heap_size].y = (uint64_t)i<<32;
			++heap_size;
		}
	}
	ks_heapmake_heap(heap_size, heap);
	while (heap_size > 0) {
		mm_match_t *q = &m[heap->y>>32];
		mm128_t *p;
		uint64_t r = heap->x;
		int32_t is_self, rpos = (uint32_t)r >> 1;
		if (skip_seed(opt->flag, r, q, qname, 0, qlen, mi, &is_self)) continue;
		if ((r&1) == (q->q_pos&1)) { // forward strand
			p = &a[n_for++];
			p->x = (r&0xffffffff00000000ULL) | rpos;
			p->y = (uint64_t)q->q_span << 32 | q->q_pos >> 1;
		} else { // reverse strand
			p = &a[(*n_a) - (++n_rev)];
			p->x = 1ULL<<63 | (r&0xffffffff00000000ULL) | rpos;
			p->y = (uint64_t)q->q_span << 32 | (qlen - ((q->q_pos>>1) + 1 - q->q_span) - 1);
		}
		p->y |= (uint64_t)q->seg_id << MM_SEED_SEG_SHIFT;
		if (q->is_tandem) p->y |= MM_SEED_TANDEM;
		if (is_self) p->y |= MM_SEED_SELF;
		// update the heap
		if ((uint32_t)heap->y < q->n - 1) {
			++heap[0].y;
			heap[0].x = m[heap[0].y>>32].cr[(uint32_t)heap[0].y];
		} else {
			heap[0] = heap[heap_size - 1];
			--heap_size;
		}
		ks_heapdown_heap(0, heap_size, heap);
	}
	kfree(km, m);
	kfree(km, heap);

	// reverse anchors on the reverse strand, as they are in the descending order
	for (j = 0; j < n_rev>>1; ++j) {
		mm128_t t = a[(*n_a) - 1 - j];
		a[(*n_a) - 1 - j] = a[(*n_a) - (n_rev - j)];
		a[(*n_a) - (n_rev - j)] = t;
	}
	if (*n_a > n_for + n_rev) {
		memmove(a + n_for, a + (*n_a) - n_rev, n_rev * sizeof(mm128_t));
		*n_a = n_for + n_rev;
	}
	return a;
}

typedef struct __attribute__((__packed__)) _FPGAHDR {
    uint32_t magic;
    uint32_t size;
    uint16_t tid;
    uint16_t num;
    uint8_t type;
    uint8_t lat;
    uint8_t padding[50];
} FPGAHDR;
typedef struct __attribute__((__packed__)) _DPHDR {
    uint32_t gap_ref;
    uint32_t gap_qry;
    uint32_t seednum;
    uint32_t qlensum;
    uint32_t ctxpos;
    uint32_t bid;
    uint16_t n_segs;
    uint8_t b;
    uint8_t padding[37];
} DPHDR;
typedef struct __attribute__((__packed__)) _ODPHDR {
    uint32_t err_flag;
    uint32_t ctxpos;
    uint32_t subsize;
    uint32_t n_a;
    uint32_t n_minipos;
    uint32_t rep_len;
    uint8_t padding[40];
} ODPHDR;

#define BUFFSIZE (16 * 1024 * 1024)

char txtbuff[BUFFSIZE];
int ntxt;
char * gentxt(uint8_t *dat, int datlen)
{
    int i, j;
    int n = datlen / 64;
    for (i = 0; i < n; i++) {
        for (j = 0; j < 64; j++) {
            sprintf(txtbuff + 129 * i + j * 2, "%02x", dat[64 * i + (63 - j)]);
        }
        txtbuff[129 * i + 128] = '\n';
    }
    ntxt = 129 * n;
    return txtbuff;
}

FILE *fcshi, *fcsho, *fdpi, *fdpo;
FILE *fbcshi, *fbcsho, *fbdpi, *fbdpo;
FILE *fidxcall;
FILE *fbreq, *fhreq, *fvreq, *fpreq;
FILE *fback, *fhack, *fvack, *fpack;
static uint8_t gbuff[BUFFSIZE];
static int gbufflen = 0;
static uint8_t gobuff[BUFFSIZE];
static int gobufflen = 0;
static int maxdp = 8;
static int maxcsh = 8;
static int gn = 0;
static uint32_t gmagic = 0;

static mm128_t *collect_seed_hits(void *km, const mm_mapopt_t *opt, int max_occ, const mm_idx_t *mi, const char *qname, const mm128_v *mv, unsigned int bid, int qlen, int64_t *n_a, int *rep_len,
								  int *n_mini_pos, uint64_t **mini_pos)
{
	int i, k, n_m;
	mm_match_t *m;
	mm128_t *a;
	//fprintf(stderr,"collect_matches begin\n");

    // sim input
    int orglen;
    if (gn == 0) {
        memset(gbuff, 0, sizeof(FPGAHDR));
        memset(gobuff, 0, sizeof(FPGAHDR));
        FPGAHDR *a = (FPGAHDR *)gbuff;
        a->magic = gmagic;
        a->type = 3;
        a->tid = 66;
        a->lat = 77;
        gbufflen = sizeof(FPGAHDR);
        FPGAHDR *o = (FPGAHDR *)gobuff;
        o->magic = gmagic;
        o->type = 3;
        o->tid = 66;
        o->lat = 77;
        gobufflen = sizeof(FPGAHDR);
    }
    DPHDR *input = (DPHDR *)(gbuff + gbufflen);
    memset(input, 0, sizeof(DPHDR));
    input->seednum = mv->n;
    input->qlensum = qlen;
    input->ctxpos = 88;
    input->bid = bid;
    input->b = mi->b;
    gbufflen += sizeof(DPHDR);
    memcpy(gbuff + gbufflen, mv->a, mv->n * 16);
    gbufflen += mv->n * 16;
    if (gbufflen % 64 != 0) {
        gbufflen = (gbufflen >> 6 << 6) + 64;
    }
    assert(gbufflen < BUFFSIZE);

	m = collect_matches(km, &n_m, max_occ, mi, mv, n_a, rep_len, n_mini_pos, mini_pos);

	a = (mm128_t*)kmalloc(km, *n_a * sizeof(mm128_t));
	for (i = 0, *n_a = 0; i < n_m; ++i) {
		mm_match_t *q = &m[i];
		const uint64_t *r = q->cr;
		for (k = 0; k < q->n; ++k) {
			int32_t is_self;
			//modified by LQX
			int32_t rpos = (r[k]>>22) & 0x1fffff;
			//fprintf(stderr,"get rpos OK\n");
			mm128_t *p;
			if (skip_seed(opt->flag, r[k], q, qname, bid, qlen, mi, &is_self)) continue;
			//fprintf(stderr,"skip_seed OK\n");
			p = &a[(*n_a)++];
			//modified by LQX
			/*
			n=1 val ----63...43|42...22| 21   |20...0
						 refid |refpos |strand|rankid
			n>1 p[k]----63...43|42...22| 21   |20...0
						 refid |refpos |strand|rankid
			*/
			//if ((r[k]&1) == (q->q_pos&1)) { // forward strand
			if((r[k]&MM_P_STRAND)>>21 == (q->q_pos&1)){// forward strand
				//p->x = (r[k]&0xffffffff00000000ULL) | rpos;
				p->x = ((r[k]&0xfffff80000000000ULL)>>11) | rpos;
				p->y = (uint64_t)q->q_span << 32 | q->q_pos >> 1;
			} else { // reverse strand
				//p->x = 1ULL<<63 | (r[k]&0xffffffff00000000ULL) | rpos;
				p->x = 1ULL<<63 | ((r[k]&0xfffff80000000000ULL)>>11)  | rpos;
				p->y = (uint64_t)q->q_span << 32 | (qlen - ((q->q_pos>>1) + 1 - q->q_span) - 1);
			}
			//fprintf(stderr,"seed assemble OK\n");
			p->y |= (uint64_t)q->seg_id << MM_SEED_SEG_SHIFT;
			if (q->is_tandem) p->y |= MM_SEED_TANDEM;
			if (is_self) p->y |= MM_SEED_SELF;
		}
	}
	kfree(km, m);
	//added by LQX
	FILE *fp = 0;
	fp = fopen("./out0.txt", "a");
	if(fp == NULL) fprintf(stderr,"open out0 file fail!\n");
	for (int k=0; k< *n_a; k++) 
	{
	  /*fwrite(&mv.a[k], 16, 1, fp);*/
	  char v[33]={0};
	  snprintf(v,sizeof(v),"%016lx%016lx",a[k].x,a[k].y);
	  fwrite(v,sizeof(v)-1,1,fp);
	  fwrite("\r\n",2,1,fp);
	  /*fprintf(stderr, "%016lx,%016lx\n", a[k].x,a[k].y);*/
	}
	fclose(fp);
	radix_sort_128x(a, a + (*n_a));
	//fprintf(stderr,"seed num %d\n",*n_a);

    // sim output
    orglen = gobufflen;
    ODPHDR *output = (ODPHDR *)(gobuff + gobufflen);
    memset(output, 0, sizeof(ODPHDR));
    output->ctxpos = 88;
    output->n_a = *n_a;
    output->n_minipos = *n_mini_pos;
    output->rep_len = *rep_len;
    gobufflen += sizeof(ODPHDR);
    memcpy(gobuff + gobufflen, a, *n_a * 16);
    gobufflen += *n_a * 16;
    if (gobufflen % 64 != 0) {
        gobufflen = (gobufflen >> 6 << 6) + 64;
    }
    memcpy(gobuff + gobufflen, *mini_pos, *n_mini_pos * 8);
    gobufflen += *n_mini_pos * 8;
    if (gobufflen % 64 != 0) {
        gobufflen = (gobufflen >> 6 << 6) + 64;
    }
    output->subsize = gobufflen - orglen;
    assert(gobufflen < BUFFSIZE);

    gn++;
    if (gn == maxcsh) {
        assert(gbufflen % 64 == 0);
        assert(gobufflen % 64 == 0);

        //set nread,size
        FPGAHDR *a = (FPGAHDR *)gbuff;
        a->size = gbufflen;
        a->num = gn;
        FPGAHDR *o = (FPGAHDR *)gobuff;
        o->size = gobufflen;
        o->num = gn;

        //write sim
        gentxt(gbuff, gbufflen);
        fwrite(txtbuff, 1, ntxt, fcshi);
        gentxt(gobuff, gobufflen);
        fwrite(txtbuff, 1, ntxt, fcsho);

        //write
        fwrite(gbuff, 1, gbufflen, fbcshi);
        fwrite(gobuff, 1, gobufflen, fbcsho);

        gn = 0;
        gmagic++;
    }

	return a;
}

static void chain_post(const mm_mapopt_t *opt, int max_chain_gap_ref, const mm_idx_t *mi, void *km, int qlen, int n_segs, const int *qlens, int *n_regs, mm_reg1_t *regs, mm128_t *a)
{
	if (!(opt->flag & MM_F_ALL_CHAINS)) { // don't choose primary mapping(s)
		mm_set_parent(km, opt->mask_level, *n_regs, regs, opt->a * 2 + opt->b);
		if (n_segs <= 1) mm_select_sub(km, opt->pri_ratio, mi->k*2, opt->best_n, n_regs, regs);
		else mm_select_sub_multi(km, opt->pri_ratio, 0.2f, 0.7f, max_chain_gap_ref, mi->k*2, opt->best_n, n_segs, qlens, n_regs, regs);
		if (!(opt->flag & (MM_F_SPLICE|MM_F_SR|MM_F_NO_LJOIN))) // long join not working well without primary chains
			mm_join_long(km, opt, qlen, n_regs, regs, a);
	}
}

static mm_reg1_t *align_regs(const mm_mapopt_t *opt, const mm_idx_t *mi, void *km, int qlen, const char *seq, int *n_regs, mm_reg1_t *regs, mm128_t *a)
{
	if (!(opt->flag & MM_F_CIGAR)) return regs;
	regs = mm_align_skeleton(km, opt, mi, qlen, seq, n_regs, regs, a); // this calls mm_filter_regs()
	if (!(opt->flag & MM_F_ALL_CHAINS)) { // don't choose primary mapping(s)
		mm_set_parent(km, opt->mask_level, *n_regs, regs, opt->a * 2 + opt->b);
		mm_select_sub(km, opt->pri_ratio, mi->k*2, opt->best_n, n_regs, regs);
		mm_set_sam_pri(*n_regs, regs);
	}
	return regs;
}

#include "minimap.h"
unsigned int dichotomy_sort(const char *qname, rname_rid_t* ref_name, int ref_name_size)
{
    unsigned int start=0, end = ref_name_size-1, mid;
    int cmp;
    while(start < end)
        {
            mid = (start + end)>>1;
            cmp = strcmp(qname, ref_name[mid].rname);
            if(cmp == 0)
                return (mid | 1UL<<31);
            else if( cmp < 0 )
                end = mid;
            else
                start = mid+1;
        }
    if (start == end) {
        cmp = strcmp(qname, ref_name[start].rname);
        if(cmp == 0)
            return (start | 1UL<<31);
        else if (cmp > 0) return start+1;
    }
    return start;
}


void mm_map_frag(const mm_idx_t *mi, int n_segs, const int *qlens, const char **seqs, int *n_regs, mm_reg1_t **regs, mm_tbuf_t *b, const mm_mapopt_t *opt, const char *qname)
{
	int i, j, rep_len, qlen_sum, n_regs0, n_mini_pos;
	int max_chain_gap_qry, max_chain_gap_ref, is_splice = !!(opt->flag & MM_F_SPLICE), is_sr = !!(opt->flag & MM_F_SR);
	uint32_t hash;
	int64_t n_a;
	uint64_t *u, *mini_pos;
	mm128_t *a;
	mm128_v mv = {0,0,0};
	mm_reg1_t *regs0;
	km_stat_t kmst;

	for (i = 0, qlen_sum = 0; i < n_segs; ++i)
		qlen_sum += qlens[i], n_regs[i] = 0, regs[i] = 0;

	if (qlen_sum == 0 || n_segs <= 0 || n_segs > MM_MAX_SEG) return;

	hash  = qname? __ac_X31_hash_string(qname) : 0;
	hash ^= __ac_Wang_hash(qlen_sum) + __ac_Wang_hash(opt->seed);
	hash  = __ac_Wang_hash(hash);

    unsigned int bid = dichotomy_sort(qname, mi->rname_rid, mi->n_seq);

	collect_minimizers(b->km, opt, mi, n_segs, qlens, seqs, &mv);
	//if (opt->flag & MM_F_HEAP_SORT) a = collect_seed_hits_heap(b->km, opt, opt->mid_occ, mi, qname, &mv, qlen_sum, &n_a, &rep_len, &n_mini_pos, &mini_pos);
	//else
		a = collect_seed_hits(b->km, opt, opt->mid_occ, mi, qname, &mv, bid, qlen_sum, &n_a, &rep_len, &n_mini_pos, &mini_pos);
	//if(n_a == 0) {fprintf(stderr,"error!\n");}
	//fprintf(stderr,"finish one read!\n");
	FILE *fp = 0;
	fp = fopen("./out1.txt", "a");
	if(fp == NULL) fprintf(stderr,"open out1 file fail!\n");
	for (int k=0; k<n_mini_pos; k++) 
	{
	  /*fwrite(&mv.a[k], 16, 1, fp);*/
	  char v[17]={0};
	  snprintf(v,sizeof(v),"%016lx",mini_pos[k]);
	  fwrite(v,sizeof(v)-1,1,fp);
	  fwrite("\r\n",2,1,fp);
	  //mini_pos++;
	  /*fprintf(stderr, "%016lx,%016lx\n", a[k].x,a[k].y);*/
	}
	fclose(fp);
	fp = fopen("./out2.txt", "a");
	{
	  /*fwrite(&mv.a[k], 16, 1, fp);*/
	  char v[25]={0};
	  snprintf(v,sizeof(v),"%08x%08x%08x",(uint32_t)n_a,n_mini_pos,rep_len);
	  fwrite(v,sizeof(v)-1,1,fp);
	  fwrite("\r\n",2,1,fp);
	  /*fprintf(stderr, "%016lx,%016lx\n", a[k].x,a[k].y);*/
	}
	fclose(fp);
	if (mm_dbg_flag & MM_DBG_PRINT_SEED) {
		fprintf(stderr, "RS\t%d\n", rep_len);
		for (i = 0; i < n_a; ++i)
			fprintf(stderr, "SD\t%s\t%d\t%c\t%d\t%d\t%d\n", mi->seq[a[i].x<<1>>33].name, (int32_t)a[i].x, "+-"[a[i].x>>63], (int32_t)a[i].y, (int32_t)(a[i].y>>32&0xff),
					i == 0? 0 : ((int32_t)a[i].y - (int32_t)a[i-1].y) - ((int32_t)a[i].x - (int32_t)a[i-1].x));
	}

	// set max chaining gap on the query and the reference sequence
	if (is_sr)
		max_chain_gap_qry = qlen_sum > opt->max_gap? qlen_sum : opt->max_gap;
	else max_chain_gap_qry = opt->max_gap;
	if (opt->max_gap_ref > 0) {
		max_chain_gap_ref = opt->max_gap_ref; // always honor mm_mapopt_t::max_gap_ref if set
	} else if (opt->max_frag_len > 0) {
		max_chain_gap_ref = opt->max_frag_len - qlen_sum;
		if (max_chain_gap_ref < opt->max_gap) max_chain_gap_ref = opt->max_gap;
	} else max_chain_gap_ref = opt->max_gap;

	if(n_a > 65535) fprintf(stderr,"seed num over 65535 read name is %s\n",qname);
	a = mm_chain_dp(max_chain_gap_ref, max_chain_gap_qry, opt->bw, opt->max_chain_skip, opt->min_cnt, opt->min_chain_score, is_splice, n_segs, n_a, a, &n_regs0, &u, b->km);
	//fprintf(stderr,"finish cdp one read!\n");
	if (opt->max_occ > opt->mid_occ && rep_len > 0) {
		int rechain = 0;
		if (n_regs0 > 0) { // test if the best chain has all the segments
			int n_chained_segs = 1, max = 0, max_i = -1, max_off = -1, off = 0;
			for (i = 0; i < n_regs0; ++i) { // find the best chain
				if (max < u[i]>>32) max = u[i]>>32, max_i = i, max_off = off;
				off += (uint32_t)u[i];
			}
			for (i = 1; i < (uint32_t)u[max_i]; ++i) // count the number of segments in the best chain
				if ((a[max_off+i].y&MM_SEED_SEG_MASK) != (a[max_off+i-1].y&MM_SEED_SEG_MASK))
					++n_chained_segs;
			if (n_chained_segs < n_segs)
				rechain = 1;
		} else rechain = 1;
		if (rechain) { // redo chaining with a higher max_occ threshold
			kfree(b->km, a);
			kfree(b->km, u);
			kfree(b->km, mini_pos);
			//if (opt->flag & MM_F_HEAP_SORT) a = collect_seed_hits_heap(b->km, opt, opt->max_occ, mi, qname, &mv, qlen_sum, &n_a, &rep_len, &n_mini_pos, &mini_pos);
			//else
		a = collect_seed_hits(b->km, opt, opt->max_occ, mi, qname, &mv, bid, qlen_sum, &n_a, &rep_len, &n_mini_pos, &mini_pos);
			a = mm_chain_dp(max_chain_gap_ref, max_chain_gap_qry, opt->bw, opt->max_chain_skip, opt->min_cnt, opt->min_chain_score, is_splice, n_segs, n_a, a, &n_regs0, &u, b->km);
		}
	}

	regs0 = mm_gen_regs(b->km, hash, qlen_sum, n_regs0, u, a);
	//fprintf(stderr,"finish mm_gen_regs one read!\n");
	if (mm_dbg_flag & MM_DBG_PRINT_SEED)
		for (j = 0; j < n_regs0; ++j)
			for (i = regs0[j].as; i < regs0[j].as + regs0[j].cnt; ++i)
				fprintf(stderr, "CN\t%d\t%s\t%d\t%c\t%d\t%d\t%d\n", j, mi->seq[a[i].x<<1>>33].name, (int32_t)a[i].x, "+-"[a[i].x>>63], (int32_t)a[i].y, (int32_t)(a[i].y>>32&0xff),
						i == regs0[j].as? 0 : ((int32_t)a[i].y - (int32_t)a[i-1].y) - ((int32_t)a[i].x - (int32_t)a[i-1].x));

	chain_post(opt, max_chain_gap_ref, mi, b->km, qlen_sum, n_segs, qlens, &n_regs0, regs0, a);
	//fprintf(stderr,"finish chain_post one read!\n");
	if (!is_sr) mm_est_err(mi, qlen_sum, n_regs0, regs0, a, n_mini_pos, mini_pos);

	if (n_segs == 1) { // uni-segment
		regs0 = align_regs(opt, mi, b->km, qlens[0], seqs[0], &n_regs0, regs0, a);
		mm_set_mapq(b->km, n_regs0, regs0, opt->min_chain_score, opt->a, rep_len, is_sr);
		n_regs[0] = n_regs0, regs[0] = regs0;
	} else { // multi-segment
		mm_seg_t *seg;
		seg = mm_seg_gen(b->km, hash, n_segs, qlens, n_regs0, regs0, n_regs, regs, a); // split fragment chain to separate segment chains
		free(regs0);
		for (i = 0; i < n_segs; ++i) {
			mm_set_parent(b->km, opt->mask_level, n_regs[i], regs[i], opt->a * 2 + opt->b); // update mm_reg1_t::parent
			regs[i] = align_regs(opt, mi, b->km, qlens[i], seqs[i], &n_regs[i], regs[i], seg[i].a);
			mm_set_mapq(b->km, n_regs[i], regs[i], opt->min_chain_score, opt->a, rep_len, is_sr);
		}
		mm_seg_free(b->km, n_segs, seg);
		if (n_segs == 2 && opt->pe_ori >= 0 && (opt->flag&MM_F_CIGAR))
			mm_pair(b->km, max_chain_gap_ref, opt->pe_bonus, opt->a * 2 + opt->b, opt->a, qlens, n_regs, regs); // pairing
	}
	//fprintf(stderr,"finish all one read!\n");
	kfree(b->km, mv.a);
	kfree(b->km, a);
	kfree(b->km, u);
	kfree(b->km, mini_pos);

	if (b->km) {
		km_stat(b->km, &kmst);
		if (mm_dbg_flag & MM_DBG_PRINT_QNAME)
			fprintf(stderr, "QM\t%s\t%d\tcap=%ld,nCore=%ld,largest=%ld\n", qname, qlen_sum, kmst.capacity, kmst.n_cores, kmst.largest);
		assert(kmst.n_blocks == kmst.n_cores); // otherwise, there is a memory leak
		if (kmst.largest > 1U<<28) {
			km_destroy(b->km);
			b->km = km_init();
		}
	}
}

mm_reg1_t *mm_map(const mm_idx_t *mi, int qlen, const char *seq, int *n_regs, mm_tbuf_t *b, const mm_mapopt_t *opt, const char *qname)
{
	mm_reg1_t *regs;
	mm_map_frag(mi, 1, &qlen, &seq, n_regs, &regs, b, opt, qname);
	return regs;
}

/**************************
 * Multi-threaded mapping *
 **************************/

typedef struct {
	int mini_batch_size, n_processed, n_threads, n_fp;
	const mm_mapopt_t *opt;
	mm_bseq_file_t **fp;
	const mm_idx_t *mi;
	kstring_t str;
} pipeline_t;

typedef struct {
	const pipeline_t *p;
    int n_seq, n_frag;
	mm_bseq1_t *seq;
	int *n_reg, *seg_off, *n_seg;
	mm_reg1_t **reg;
	mm_tbuf_t **buf;
} step_t;

static void worker_for(void *_data, long i, int tid) // kt_for() callback
{
    step_t *s = (step_t*)_data;
	int qlens[MM_MAX_SEG], j, off = s->seg_off[i], pe_ori = s->p->opt->pe_ori;
	const char *qseqs[MM_MAX_SEG];
	mm_tbuf_t *b = s->buf[tid];
	assert(s->n_seg[i] <= MM_MAX_SEG);
	if (mm_dbg_flag & MM_DBG_PRINT_QNAME)
		fprintf(stderr, "QR\t%s\t%d\t%d\n", s->seq[off].name, tid, s->seq[off].l_seq);
	for (j = 0; j < s->n_seg[i]; ++j) {
		if (s->n_seg[i] == 2 && ((j == 0 && (pe_ori>>1&1)) || (j == 1 && (pe_ori&1))))
			mm_revcomp_bseq(&s->seq[off + j]);
		qlens[j] = s->seq[off + j].l_seq;
		qseqs[j] = s->seq[off + j].seq;
	}
	if (s->p->opt->flag & MM_F_INDEPEND_SEG) {
		for (j = 0; j < s->n_seg[i]; ++j)
			mm_map_frag(s->p->mi, 1, &qlens[j], &qseqs[j], &s->n_reg[off+j], &s->reg[off+j], b, s->p->opt, s->seq[off+j].name);
	} else {
		mm_map_frag(s->p->mi, s->n_seg[i], qlens, qseqs, &s->n_reg[off], &s->reg[off], b, s->p->opt, s->seq[off].name);
	}
	for (j = 0; j < s->n_seg[i]; ++j) // flip the query strand and coordinate to the original read strand
		if (s->n_seg[i] == 2 && ((j == 0 && (pe_ori>>1&1)) || (j == 1 && (pe_ori&1)))) {
			int k, t;
			mm_revcomp_bseq(&s->seq[off + j]);
			for (k = 0; k < s->n_reg[off + j]; ++k) {
				mm_reg1_t *r = &s->reg[off + j][k];
				t = r->qs;
				r->qs = qlens[j] - r->qe;
				r->qe = qlens[j] - t;
				r->rev = !r->rev;
			}
		}
}

static void *worker_pipeline(void *shared, int step, void *in)
{
	int i, j, k;
    pipeline_t *p = (pipeline_t*)shared;
    if (step == 0) { // step 0: read sequences
		int with_qual = (!!(p->opt->flag & MM_F_OUT_SAM) && !(p->opt->flag & MM_F_NO_QUAL));
		int with_comment = !!(p->opt->flag & MM_F_COPY_COMMENT);
		int frag_mode = (p->n_fp > 1 || !!(p->opt->flag & MM_F_FRAG_MODE));
        step_t *s;
        s = (step_t*)calloc(1, sizeof(step_t));
		if (p->n_fp > 1) s->seq = mm_bseq_read_frag2(p->n_fp, p->fp, p->mini_batch_size, with_qual, with_comment, &s->n_seq);
		else s->seq = mm_bseq_read3(p->fp[0], p->mini_batch_size, with_qual, with_comment, frag_mode, &s->n_seq);
		if (s->seq) {
			s->p = p;
			for (i = 0; i < s->n_seq; ++i)
				s->seq[i].rid = p->n_processed++;
			s->buf = (mm_tbuf_t**)calloc(p->n_threads, sizeof(mm_tbuf_t*));
			for (i = 0; i < p->n_threads; ++i)
				s->buf[i] = mm_tbuf_init();
			s->n_reg = (int*)calloc(3 * s->n_seq, sizeof(int));
			s->seg_off = s->n_reg + s->n_seq; // seg_off and n_seg are allocated together with n_reg
			s->n_seg = s->seg_off + s->n_seq;
			s->reg = (mm_reg1_t**)calloc(s->n_seq, sizeof(mm_reg1_t*));
			for (i = 1, j = 0; i <= s->n_seq; ++i)
				if (i == s->n_seq || !frag_mode || !mm_qname_same(s->seq[i-1].name, s->seq[i].name)) {
					s->n_seg[s->n_frag] = i - j;
					s->seg_off[s->n_frag++] = j;
					j = i;
				}
			return s;
		} else free(s);
    } else if (step == 1) { // step 1: map
		kt_for(p->n_threads, worker_for, in, ((step_t*)in)->n_frag);
		return in;
    } else if (step == 2) { // step 2: output
		void *km = 0;
        step_t *s = (step_t*)in;
		const mm_idx_t *mi = p->mi;
		for (i = 0; i < p->n_threads; ++i) mm_tbuf_destroy(s->buf[i]);
		free(s->buf);
		if ((p->opt->flag & MM_F_OUT_CS) && !(mm_dbg_flag & MM_DBG_NO_KALLOC)) km = km_init();
		for (k = 0; k < s->n_frag; ++k) {
			int seg_st = s->seg_off[k], seg_en = s->seg_off[k] + s->n_seg[k];
			for (i = seg_st; i < seg_en; ++i) {
				mm_bseq1_t *t = &s->seq[i];
				for (j = 0; j < s->n_reg[i]; ++j) {
					mm_reg1_t *r = &s->reg[i][j];
					assert(!r->sam_pri || r->id == r->parent);
					if ((p->opt->flag & MM_F_NO_PRINT_2ND) && r->id != r->parent)
						continue;
					if (p->opt->flag & MM_F_OUT_SAM)
						mm_write_sam2(&p->str, mi, t, i - seg_st, j, s->n_seg[k], &s->n_reg[seg_st], (const mm_reg1_t*const*)&s->reg[seg_st], km, p->opt->flag);
					else
						mm_write_paf(&p->str, mi, t, r, km, p->opt->flag);
					mm_err_puts(p->str.s);
				}
				if (s->n_reg[i] == 0 && (p->opt->flag & MM_F_OUT_SAM)) { // write an unmapped record
					mm_write_sam2(&p->str, mi, t, i - seg_st, -1, s->n_seg[k], &s->n_reg[seg_st], (const mm_reg1_t*const*)&s->reg[seg_st], km, p->opt->flag);
					mm_err_puts(p->str.s);
				}
			}
			for (i = seg_st; i < seg_en; ++i) {
				for (j = 0; j < s->n_reg[i]; ++j) free(s->reg[i][j].p);
				free(s->reg[i]);
				free(s->seq[i].seq); free(s->seq[i].name);
				if (s->seq[i].qual) free(s->seq[i].qual);
			}
		}
		free(s->reg); free(s->n_reg); free(s->seq); // seg_off and n_seg were allocated with reg; no memory leak here
		km_destroy(km);
		if (mm_verbose >= 3)
			fprintf(stderr, "[M::%s::%.3f*%.2f] mapped %d sequences\n", __func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0), s->n_seq);
		free(s);
	}
    return 0;
}

int mm_map_file_frag(const mm_idx_t *idx, int n_segs, const char **fn, const mm_mapopt_t *opt, int n_threads)
{
	int i, j, pl_threads;
	pipeline_t pl;
	if (n_segs < 1) return -1;
	memset(&pl, 0, sizeof(pipeline_t));
	pl.n_fp = n_segs;
	pl.fp = (mm_bseq_file_t**)calloc(n_segs, sizeof(mm_bseq_file_t*));
	for (i = 0; i < n_segs; ++i) {
		pl.fp[i] = mm_bseq_open(fn[i]);
		if (pl.fp[i] == 0) {
			if (mm_verbose >= 1)
				fprintf(stderr, "ERROR: failed to open file '%s'\n", fn[i]);
			for (j = 0; j < i; ++j)
				mm_bseq_close(pl.fp[j]);
			free(pl.fp);
			return -1;
		}
	}
	pl.opt = opt, pl.mi = idx;
	pl.n_threads = n_threads > 1? n_threads : 1;
	pl.mini_batch_size = opt->mini_batch_size;
	pl_threads = n_threads == 1? 1 : (opt->flag&MM_F_2_IO_THREADS)? 3 : 2;

    fprintf(stderr, "------ bw: 0x%x, is_cnda: 0x%x\n", opt->bw, !!(opt->flag & MM_F_SPLICE));
    fprintf(stderr, "------ max_skip: 0x%x, min_sc: 0x%x\n", opt->max_chain_skip, opt->min_chain_score);
    fprintf(stderr, "------ max_occ: 0x%x, flag: 0x%x\n", opt->mid_occ, opt->flag);
    if ((fcshi = fopen("cshinput.txt", "w")) == NULL) {
        perror("open data file failed!");
        exit(1);
    }
    if ((fcsho = fopen("cshoutput.txt", "w")) == NULL) {
        perror("open data file failed!");
        exit(1);
    }
    if ((fdpi = fopen("dpinput.txt", "w")) == NULL) {
        perror("open data file failed!");
        exit(1);
    }
    if ((fdpo = fopen("dpoutput.txt", "w")) == NULL) {
        perror("open data file failed!");
        exit(1);
    }

    if ((fbcshi = fopen("cshinput.dat", "wb")) == NULL) {
        perror("open data file failed!");
        exit(1);
    }
    if ((fbcsho = fopen("cshoutput.dat", "wb")) == NULL) {
        perror("open data file failed!");
        exit(1);
    }
    if ((fbdpi = fopen("dpinput.dat", "wb")) == NULL) {
        perror("open data file failed!");
        exit(1);
    }
    if ((fbdpo = fopen("dpoutput.dat", "wb")) == NULL) {
        perror("open data file failed!");
        exit(1);
    }
    if ((fidxcall = fopen("idxcall.txt", "w")) == NULL) {
        perror("open data file failed!");
        exit(1);
    }
    //BReq.hex、BAck.hex、HReq.hex、HAck.hex、VReq.hex、VAck.hex、PReq.hex、PAck.hex
    if ((fbreq = fopen("breq.txt", "w")) == NULL) {
        perror("open data file failed!");
        exit(1);
    }
    if ((fhreq = fopen("hreq.txt", "w")) == NULL) {
        perror("open data file failed!");
        exit(1);
    }
    if ((fvreq = fopen("vreq.txt", "w")) == NULL) {
        perror("open data file failed!");
        exit(1);
    }
    if ((fpreq = fopen("preq.txt", "w")) == NULL) {
        perror("open data file failed!");
        exit(1);
    }
    if ((fback = fopen("back.txt", "w")) == NULL) {
        perror("open data file failed!");
        exit(1);
    }
    if ((fhack = fopen("hack.txt", "w")) == NULL) {
        perror("open data file failed!");
        exit(1);
    }
    if ((fvack = fopen("vack.txt", "w")) == NULL) {
        perror("open data file failed!");
        exit(1);
    }
    if ((fpack = fopen("pack.txt", "w")) == NULL) {
        perror("open data file failed!");
        exit(1);
    }
    system("rm DP_IN.dat");
    system("rm DP_SEED.dat");
    system("rm REQACK.dat");

	kt_pipeline(pl_threads, worker_pipeline, &pl, 3);

    fclose(fcshi);
    fclose(fcsho);
    fclose(fdpi);
    fclose(fdpo);

    fclose(fbcshi);
    fclose(fbcsho);
    fclose(fbdpi);
    fclose(fbdpo);
    system("./revert.py");

    fclose(fidxcall);
    fclose(fbreq);
    fclose(fhreq);
    fclose(fvreq);
    fclose(fpreq);
    fclose(fback);
    fclose(fhack);
    fclose(fvack);
    fclose(fpack);

	free(pl.str.s);
	for (i = 0; i < n_segs; ++i)
		mm_bseq_close(pl.fp[i]);
	free(pl.fp);
	return 0;
}

int mm_map_file(const mm_idx_t *idx, const char *fn, const mm_mapopt_t *opt, int n_threads)
{
	return mm_map_file_frag(idx, 1, &fn, opt, n_threads);
}
