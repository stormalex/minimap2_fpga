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
#include "minimap.h"
#include "ksw2.h"

#include "kfifo.h"
#include "fpga_sim.h"
#include "context_data.h"

#include "user_common.h"

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

static inline int skip_seed(int flag, uint64_t r, const mm_match_t *q, const char *qname, int qlen, const mm_idx_t *mi, int *is_self)
{
	*is_self = 0;
	if (qname && (flag & (MM_F_NO_DIAG|MM_F_NO_DUAL))) {
		const mm_idx_seq_t *s = &mi->seq[r>>32];
		int cmp;
		cmp = strcmp(qname, s->name);
		if ((flag&MM_F_NO_DIAG) && cmp == 0 && s->len == qlen) {
			if ((uint32_t)r>>1 == (q->q_pos>>1)) return 1; // avoid the diagnonal anchors
			if ((r&1) == (q->q_pos&1)) *is_self = 1; // this flag is used to avoid spurious extension on self chain
		}
		if ((flag&MM_F_NO_DUAL) && cmp > 0) // all-vs-all mode: map once
			return 1;
	}
	if (flag & (MM_F_FOR_ONLY|MM_F_REV_ONLY)) {
		if ((r&1) == (q->q_pos&1)) { // forward strand
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
		if (skip_seed(opt->flag, r, q, qname, qlen, mi, &is_self)) continue;
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

mm128_t *collect_seed_hits(void *km, const mm_mapopt_t *opt, int max_occ, const mm_idx_t *mi, const char *qname, const mm128_v *mv, int qlen, int64_t *n_a, int *rep_len,
								  int *n_mini_pos, uint64_t **mini_pos)
{
	int i, k, n_m;
	mm_match_t *m;
	mm128_t *a;
	m = collect_matches(km, &n_m, max_occ, mi, mv, n_a, rep_len, n_mini_pos, mini_pos);
	a = (mm128_t*)kmalloc(km, *n_a * sizeof(mm128_t));
	for (i = 0, *n_a = 0; i < n_m; ++i) {
		mm_match_t *q = &m[i];
		const uint64_t *r = q->cr;
		for (k = 0; k < q->n; ++k) {
			int32_t is_self, rpos = (uint32_t)r[k] >> 1;
			mm128_t *p;
			if (skip_seed(opt->flag, r[k], q, qname, qlen, mi, &is_self)) continue;
			p = &a[(*n_a)++];
			if ((r[k]&1) == (q->q_pos&1)) { // forward strand
				p->x = (r[k]&0xffffffff00000000ULL) | rpos;
				p->y = (uint64_t)q->q_span << 32 | q->q_pos >> 1;
			} else { // reverse strand
				p->x = 1ULL<<63 | (r[k]&0xffffffff00000000ULL) | rpos;
				p->y = (uint64_t)q->q_span << 32 | (qlen - ((q->q_pos>>1) + 1 - q->q_span) - 1);
			}
			p->y |= (uint64_t)q->seg_id << MM_SEED_SEG_SHIFT;
			if (q->is_tandem) p->y |= MM_SEED_TANDEM;
			if (is_self) p->y |= MM_SEED_SELF;
		}
	}
	kfree(km, m);
	radix_sort_128x(a, a + (*n_a));
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

int dp_consumer(void *topt, int tid);

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

	collect_minimizers(b->km, opt, mi, n_segs, qlens, seqs, &mv);
	/*
	if (opt->flag & MM_F_HEAP_SORT) a = collect_seed_hits_heap(b->km, opt, opt->mid_occ, mi, qname, &mv, qlen_sum, &n_a, &rep_len, &n_mini_pos, &mini_pos);
	else a = collect_seed_hits(b->km, opt, opt->mid_occ, mi, qname, &mv, qlen_sum, &n_a, &rep_len, &n_mini_pos, &mini_pos);
	*/
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

#define USE_FPGA_SIM
#ifndef USE_FPGA_SIM
#error "hhahah"
	a = mm_chain_dp(max_chain_gap_ref, max_chain_gap_qry, opt->bw, opt->max_chain_skip, opt->min_cnt, opt->min_chain_score, is_splice, n_segs, n_a, a, &n_regs0, &u, b->km);
#else
	int ctxspos = save_seed_data(max_chain_gap_ref,max_chain_gap_qry,/*n_a,*/opt->flag,opt->mid_occ,qlen_sum,n_segs,/*a,*/b->km,qname,mi,&mv);
	save_dpctx(ctxspos,mi,n_regs,regs,b,mv.a,/*mv.n,mv.m,mini_pos,*/hash,n_segs,qlen_sum,max_chain_gap_ref,max_chain_gap_qry,/*rep_len,n_mini_pos,*/qlens,seqs,qname);

#endif

	/*if (opt->max_occ > opt->mid_occ && rep_len > 0) { //这些代码应该放到dp_consumer中
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
			if (opt->flag & MM_F_HEAP_SORT) a = collect_seed_hits_heap(b->km, opt, opt->max_occ, mi, qname, &mv, qlen_sum, &n_a, &rep_len, &n_mini_pos, &mini_pos);
			else a = collect_seed_hits(b->km, opt, opt->max_occ, mi, qname, &mv, qlen_sum, &n_a, &rep_len, &n_mini_pos, &mini_pos);
			a = mm_chain_dp(max_chain_gap_ref, max_chain_gap_qry, opt->bw, opt->max_chain_skip, opt->min_cnt, opt->min_chain_score, is_splice, n_segs, n_a, a, &n_regs0, &u, b->km);
		}
	}*/
	/*
	snd_ctrl.stop = 1;
	chaindp_task_sender(&snd_ctrl);
	*/
}//end func mm_


/*
typedef struct 
{
    volatile uint8_t enable;
    volatile uint8_t readen[DP_BATCHSIZE];//0:read disable, need free chaindp mem, 1:readable, 2:free already,do nothing, 3:read disable for one time only, skip curr reg
} donext_t ;

static volatile donext_t donext_read[MAX_DP_THRD] = {
    [ 0 ... MAX_DP_THRD-1 ] = {1, {0}}
};

static inline void set_donext_enable(uint8_t tid, uint8_t eb)
{
    donext_read[tid].enable = eb;
}


static inline void set_donext(uint8_t tid, uint8_t eb, uint8_t rd, int n)
{
    donext_read[tid].readen[rd] = n;
    set_donext_enable(tid, eb);
}

static inline uint8_t get_donext_enable(uint8_t tid)
{
    return donext_read[tid].enable;
}

static inline int get_donext_readen(uint8_t tid, uint8_t rd)
{
    return donext_read[tid].readen[rd];
}
*/

void debug_reg(int n,mm_reg1_t *r)
{
    fprintf(stderr,"reg:nreg0 %d,id %d,cnt %d,rid %d,score %d,qs %d,qe %d,rs %d,re %d,par %d,subsc %d,as %d,mlen %d,blen %d,nsub %d,score0 %d,mapq %d,split %d,rev %d,inv %d,pri %d,p %p\n",
        n,r->id,r->cnt,r->rid,r->score,r->qs,r->qe,r->rs,r->re,r->parent,r->subsc,r->as,r->mlen,r->blen,r->n_sub,r->score0,r->mapq,r->split,r->rev,r->inv,r->sam_pri,r->p);
}

static volatile int dpthrd_stop = 0;//0=running, else thread exit

void set_stat_dpthrd(int stop)
{
    dpthrd_stop = stop;
}

int get_stat_dpthrd(void)
{
    return dpthrd_stop;
}

void reverse_cigar(int n_cigar, uint32_t *cigar)
{
  uint32_t tmp;
  for (int i=0; i<n_cigar>>1; i++)
    tmp = cigar[i],cigar[i]=cigar[n_cigar-1-i],cigar[n_cigar-1-i]=tmp;
}

static int sw_consumer(void *topt, int tid)
{
    const mm_idx_t *mi;
    mm_tbuf_t *b;
    void *newkm;
    mm128_t *a;
    mm_reg1_t *regs0, **regs;
    int *n_regs, *qlens;
    int n_regs0, rep_len;
    const mm_mapopt_t *opt = (const mm_mapopt_t *)topt;
    int is_sr = !!(opt->flag & MM_F_SR);;


    swtest_rcvhdr_t *ret = 0;
    get_from_ring_buf(rcv_ctrl.sw_ring_buf, (void **)&ret);
    if (!ret) {
        if (get_stat_dpthrd()/*lastblk_comeback_sw() && (add_rsp_sw(0) == add_req_sw(0))*/) {fprintf(stderr, "....tid %d..sw consumer exit....hhahahahah.....\n",tid);return -1;}
        else return 0;
    }

    void *km = ret->km;

    add_rsp_sw(1);
    if (IS_LAST == ret->lat) set_last_flag_sw(ret->lat);//avoid rewrite by NT_LAST

    fprintf(stderr, "sw consumer begn! magic %d,factnum %d,size %u,mytid %d, dptid %d,last %d,req %u,rsp %u\n",ret->magic,ret->num,ret->size,tid,ret->tid,lastblk_comeback_sw(),add_req_sw(0),add_rsp_sw(0));

    for (int k=0; k<ret->num; k++) {
        sw_readhdr_t *dr = ret->data + k;
        sw_reghdr_t  *rr = &dr->reg;
        sw_context_t *swctx  = get_swctx_by_pos(dr->ctxpos);
        //reg_context_t*regctx = swctx->regctx;

        uint8_t *qseq0[2];
        qseq0[0] = swctx->qseq0[0];
        qseq0[1] = swctx->qseq0[1];

        mi = swctx->mi, b = swctx->b, regs0 = swctx->regs0, regs = swctx->regs;
        n_regs = swctx->n_regs, qlens = swctx->qlens, n_regs0 = swctx->n_regs0;
        rep_len = swctx->rep_len, a = swctx->a, newkm = swctx->newkm;

        //if (dr->regnum != 0) 
        { 
            sw_reghdr_t *cu = rr;
            uint32_t regpos = cu->regpos;
            //per reg has only one ez
            ksw_extz_t *ez  = (ksw_extz_t *)((uint8_t *)ret + cu->offset);
            reg_context_t *ct = &swctx->regctx[regpos];
            mm_reg1_t *r  = &swctx->regs0[regpos];
            //mm_reg1_t r02;
            mm_reg1_t *r2 = &ct->r2;memset(r2, 0, sizeof(mm_reg1_t));

            int32_t rs1 = ct->rs1,  qs1 = ct->qs1,  re1 = ct->re1,qe1 = ct->qe1;
            int32_t rs  = ct->rs,   qs  = ct->qs,   re  = ct->re, qe  = ct->qe;
            int32_t llen= ct->llen,rlen = ct->rlen,qlen = ct->qlen;
            int32_t rid = ct->rid,  rev = ct->rev;
            int32_t as1 = ct->as1,  cnt1= ct->cnt1;
            int32_t rs0 = ct->rs0,  re0 = ct->re0;

            //fprintf(stderr, "sw consumer k %d,mytid %d,dptid %d, swpos %d, enable %d, regctx %p, regpos %d, nregs0 %d,ctx qe1 %d,llen %d,rlen %d, qlen %d, midnum %x\n",k,tid,ret->tid,dr->ctxpos,swctx->enable,swctx->regctx,regpos,n_regs0,qe1,llen,rlen,qlen,cu->midnum);
            //left
            if (get_left(cu->midnum)) {
                if (ez->n_cigar > 0) {
		  if (REV_CIGAR_FPGA == ez->revcigar) reverse_cigar(ez->n_cigar, ez->cigar);
                    mm_append_cigar(r, ez->n_cigar, ez->cigar);
                    r->p->dp_score += ez->max;
                }
                rs1 = rs - (ez->reach_end? ez->mqe_t + 1 : ez->max_t + 1);
                qs1 = qs - (ez->reach_end? llen/*qs - qs0*/ : ez->max_q + 1);
                //fprintf(stderr, "lft:append cigar,n_cigar %d,dpscore %d,rs1/rs %d/%d,qs1/qs %d/%d,ez->reach_end %d,ez->mqe_t %d,  ez->max_t %d,llen %d,ez->max_q %d,ez->max %d\n",ez->n_cigar,r->p->dp_score,rs1,rs,qs1,qs,ez->reach_end,ez->mqe_t,ez->max_t,llen,ez->max_q,ez->max);
                if (FREE_CIGAR_X86 == ret->freecigar) kfree(km, ez->cigar);
                ez++;
            }

            //middle
            int8_t mat[25];
            ksw_gen_simple_mat(5, mat, opt->a, opt->b);
            uint8_t *ztseq = (uint8_t*)kmalloc(b->km, re0 - rs0);
            int i, dropped = 0;
            for (i = is_sr? cnt1 - 1 : 1; i < cnt1; ++i) { // gap filling
                if ((a[as1+i].y & (MM_SEED_IGNORE|MM_SEED_TANDEM)) && i != cnt1 - 1) continue;
                if (is_sr && !(mi->flag & MM_I_HPC)) {
                    re = (int32_t)a[as1 + i].x + 1;
                    qe = (int32_t)a[as1 + i].y + 1;
                } else mm_adjust_minier(mi, qseq0, &a[as1 + i], &re, &qe);
                re1 = re, qe1 = qe;
                if (i == cnt1 - 1 || (a[as1+i].y&MM_SEED_LONG_JOIN) || (qe - qs >= opt->min_ksw_len && re - rs >= opt->min_ksw_len)) {
                    int j, zdrop_code;
                    zdrop_context_t *zctx = &ct->zdropctx[i];
                    uint8_t *zqseq = &qseq0[zctx->rev][zctx->qs];
                    mm_idx_getseq(mi, zctx->rid, zctx->rs, zctx->re, ztseq);
		    if (REV_CIGAR_FPGA == ez->revcigar) reverse_cigar(ez->n_cigar, ez->cigar);
                    ez->m_cigar = get_power_2_big(ez->n_cigar);
                    if ((zdrop_code = mm_test_zdrop(b->km, opt, zqseq, ztseq, ez->n_cigar, ez->cigar, mat)) != 0) {

                        mm_align_pair(b->km, opt, zctx->qe-zctx->qs, zqseq, zctx->re-zctx->rs, ztseq, mat, zctx->bw1, -1, zdrop_code == 2? opt->zdrop_inv : opt->zdrop, zctx->extra_flag, ez);
                        //fprintf(stderr, "mid zdrop,qlen %d,tlen %d,bw %d,endbos %d,zcode %d, zinv %d,zdrop %d,factzdrop %d,flag %x,ez max %d,zdrp %d,maxq %d,maxt %d,mqe %d,mqet %d,score %d,ncig %d,reaend %d\n",zctx->qe-zctx->qs,zctx->re-zctx->rs,zctx->bw1,-1,zdrop_code,opt->zdrop_inv,opt->zdrop,zdrop_code == 2?opt->zdrop_inv : opt->zdrop,zctx->extra_flag,ez->max,ez->zdropped,ez->max_q,ez->max_t,ez->mqe,ez->mqe_t,ez->score,ez->n_cigar,ez->reach_end);

		    if (REV_CIGAR_FPGA == ez->revcigar) reverse_cigar(ez->n_cigar, ez->cigar);
                    }
                    // update CIGAR
                    if (ez->n_cigar > 0)
                        mm_append_cigar(r, ez->n_cigar, ez->cigar);
                    if (ez->zdropped) { // truncated by Z-drop; TODO: sometimes Z-drop kicks in because the next seed placement is wrong. This can be fixed in principle.
                        for (j = i - 1; j >= 0; --j)
                            if ((int32_t)a[as1 + j].x <= rs + ez->max_t) break;
                        dropped = 1;
                        //fprintf(stderr,"have ez->zdroped, dropped=1\n");
                        if (j < 0) j = 0;
                        r->p->dp_score += ez->max;
                        re1 = rs + (ez->max_t + 1);
                        qe1 = qs + (ez->max_q + 1);
                        if (cnt1 - (j + 1) >= opt->min_cnt) {
                            r2->cnt = 0;
                            mm_split_reg(r, r2, as1 + j + 1 - r->as, qlen, a);
                            if (zdrop_code == 2) r2->split_inv = 1;
                            //fprintf(stderr, "re1 %d,qe1 %d,rs %d,qs %d,maxt %d,maxq %d,r2cnt %d,qlen %d,code %d,split %d\n",re1,qe1,rs,qs,ez->max_t,ez->max_q,r2->cnt,qlen,zdrop_code,r2->split_inv);
                        }
                        break;
                    } else r->p->dp_score += ez->score;
                    //fprintf(stderr, "mid:append cigar,n_cigar %d,dpscore %d,re1/re %d/%d,qe1/qe %d/%d,drop %d\n",ez->n_cigar,r->p->dp_score,re1,re,qe1,qe,dropped);
                    rs = re, qs = qe;
                    if (FREE_CIGAR_X86 == ret->freecigar) kfree(km, ez->cigar);
                    ez++;
                }//if
            }//for
            kfree(b->km, ztseq);
            //right
            if (!dropped && get_right(cu->midnum)) {
                if (ez->n_cigar > 0) {
		  if (REV_CIGAR_FPGA == ez->revcigar) reverse_cigar(ez->n_cigar, ez->cigar);
                    mm_append_cigar(r, ez->n_cigar, ez->cigar);
                    r->p->dp_score += ez->max;
                }
                re1 = re + (ez->reach_end? ez->mqe_t + 1 : ez->max_t + 1);
                qe1 = qe + (ez->reach_end? rlen/*qe0 - qe*/ : ez->max_q + 1);
                //fprintf(stderr, "rht:append cigar,n_cigar %d,dpscore %d,re1/re %d/%d,qe1/qe %d/%d,qlen %d\n",ez->n_cigar,r->p->dp_score,re1,re,qe1,qe,qlen);
                if (FREE_CIGAR_X86 == ret->freecigar) kfree(km, ez->cigar);
                ez++;
            }
            assert(qe1 <= qlen);
            r->rs = rs1, r->re = re1;
            if (rev) r->qs = qlen - qe1, r->qe = qlen - qs1;
            else r->qs = qs1, r->qe = qe1;

            //fprintf(stderr, "befor update, rs %d,re %d,qs %d, qe %d,rid %d,rs1 %d,re1 %d,re0 %d,rs0 %d, det1 %d,det0 %d\n",r->rs,r->re,r->qs,r->qe,rid,rs1,re1,re0,rs0,re1 - rs1,re0 - rs0);
            assert(re1 - rs1 <= re0 - rs0);
            if (r->p) {
                //int8_t mat[25];
                //ksw_gen_simple_mat(5, mat, opt->a, opt->b);
                uint8_t *tseq = (uint8_t*)kmalloc(km, re1 - rs1);
                mm_idx_getseq(mi, rid, rs1, re1, tseq);

                mm_update_extra(r, &qseq0[r->rev][qs1], tseq, mat, opt->q, opt->e);
                kfree(km, tseq);
                if (rev && r->p->trans_strand)
                r->p->trans_strand ^= 3; // flip to the read strand
            }
            int ii = regpos;
            //int oldnreg = n_regs0;
            if (r2->cnt > 0) regs0 = mm_insert_reg(r2, ii, &n_regs0, regs0);
            if (ii > 0 && regs0[ii].split_inv) {
                ksw_extz_t ez;
                memset(&ez, 0, sizeof(ksw_extz_t));
                if (mm_align1_inv(km, opt, mi, qlens[0], qseq0, &regs0[ii-1], &regs0[ii], r2, &ez)) {
                    regs0 = mm_insert_reg(r2, ii, &n_regs0, regs0);
                    ++swctx->n_regs0_swdone;
                    kfree(km, ez.cigar);
                    //fprintf(stderr, "############zdrop sw in consumer.....r2cnt %d,qlen %d, newnregs %d,oldnreg %d,readseq %d,cur reg %d,skip %d,ptr %p\n",r2->cnt,qlens[0],n_regs0,oldnreg,k,ii,ii+1,regs0);
                }
            }
            if (n_regs0 >= swctx->n_regs0_cap) {
                swctx->n_regs0_cap <<= 1;
		fprintf(stderr, "reg realloc, regctx %p, size %d\n", swctx->regctx,swctx->n_regs0_cap*sizeof(reg_context_t));
                swctx->regctx = (reg_context_t *)krealloc(b->km, swctx->regctx, swctx->n_regs0_cap*sizeof(reg_context_t));
		fprintf(stderr, "reg realloc done!\n");
            }
            //TODO: free regs0 first???
            swctx->regs0   = regs0;
            swctx->n_regs0 = n_regs0;
            ++swctx->n_regs0_swdone;
            //kfree(km, ez->cigar);//TODO, has bug

            if (swctx->n_regs0 == swctx->n_regs0_swdone) {
                //do this after read done!
                //int oldregnum = swctx->n_regs0;
                if (opt->flag & MM_F_CIGAR) {
                    //regs0 = 
                    //*n_regs_ = n_regs;
                    //kfree(b->km, seqs[0]);//TODO after sendto fpga
                    //kfree(b->km, ez.cigar);
                    mm_filter_regs(b->km, opt, qlens[0], &swctx->n_regs0, swctx->regs0);
                    //debug_reg(swctx->n_regs0, swctx->regs0);
                    mm_hit_sort_by_dp(b->km, &swctx->n_regs0, swctx->regs0);
                    if (!(opt->flag & MM_F_ALL_CHAINS)) { // don't choose primary mapping(s)
                        mm_set_parent(b->km, opt->mask_level, swctx->n_regs0, swctx->regs0, opt->a * 2 + opt->b);
                        mm_select_sub(b->km, opt->pri_ratio, mi->k*2, opt->best_n, &swctx->n_regs0, swctx->regs0);
                        mm_set_sam_pri(swctx->n_regs0, swctx->regs0);
                    }
                }//CIGAR

                mm_set_mapq(b->km, swctx->n_regs0, swctx->regs0, opt->min_chain_score, opt->a, rep_len, is_sr);
                n_regs[0] = swctx->n_regs0, regs[0] = swctx->regs0;

                //fprintf(stderr, "sw consumer freemem, magic %d, swpos %d, dppos %d,regctx %p, qseq0 %p\n",ret->magic,dr->ctxpos,swctx->dppos,swctx->regctx,qseq0[0]);
                kfree(b->km, swctx->regctx);
                kfree(b->km, qseq0[0]);
                kfree(newkm, a);
                swctx->enable = 0;
                //move to dpconsumer thread to avoid lock between dp consumer and sw consumer in ring
                //put_free_dppos(swctx->dppos);
                //put_free_swpos(dr->ctxpos);
            }//if ==swdone ok
            //else swctx->enable = 1;
        }
    }//end for read
#if 1
    for (int k=0; k<ret->num; k++)    {
        sw_readhdr_t *dr = ret->data + k;
        sw_context_t *swctx  = get_swctx_by_pos(dr->ctxpos);
        if (0 != swctx->enable) swctx->enable = 1;
    }
#endif
    fprintf(stderr, "sw consumer done! magic %d,mytid %d,ret tid %d\n",ret->magic,tid,ret->tid);

    kfree(km, ret);

    return 0;
}

int dp_consumer(void *topt, int tid)
{
    int i,j;
    const char *qname;
    const mm_idx_t *mi;
    const mm_mapopt_t *opt = (const mm_mapopt_t *)topt;
    int max_chain_gap_qry, max_chain_gap_ref;
    int is_splice = !!(opt->flag & MM_F_SPLICE), is_sr = !!(opt->flag & MM_F_SR);
    uint32_t hash;
    uint64_t *mini_pos;
    mm128_t *a;
    mm128_v mv = {0,0,0};
    mm_reg1_t *regs0, **regs;
    km_stat_t kmst;
    mm_tbuf_t *b;
    int n_segs,rep_len,n_regs0,qlen_sum;
    int *n_regs, *qlens;
    char **seqs;
    void *newkm = 0;

    n_segs = 1;
    int unavi = 0, avi = 0;
    int unavi_pos[SW_LFT_NUM] = {0}, avi_pos[DP_BATCHSIZE] = {0};
    for (;;) {
        int *pos = (int*)0xdead;
        int idx  = get_from_ring_buf_st(left_sw_pos[tid], (void**)&pos);
        if (idx >= 0) {
            //ring have elements
            sw_context_t *curr = get_swctx_by_pos((int)pos);
            //fprintf(stderr, "dp consumer, tid %d, get swpos %d, enable %d from left ring\n", tid,(int)pos,curr->enable);
            if (0 == curr->enable) { 
                put_free_dppos(curr->dppos); 
                put_free_swpos((int)pos); 
                continue;
            }//sw done
            else if (3 == curr->enable) { unavi_pos[unavi++] = (int)pos;}//skip
            else if (1 == curr->enable) { avi_pos[avi++] = (int)pos;}
            
            if (DP_BATCHSIZE == avi) break; 
        } else break;
    }

    //put unavi back to ring
    for (int k=0; k<unavi; k++) put_to_ring_buf_st(left_sw_pos[tid], (void *)unavi_pos[k]);
    
    int bakunavi = unavi;
    sw_context_t *tmp = get_swctx_by_pos(avi_pos[0]);
    if (tmp) newkm = tmp->newkm;
    
    if (avi > 0) {

        struct timespec tv;
        sw_to_ring_t *swring = sw_ring_init(newkm,opt,0,tid,NT_LAST);
        memset(swring->data, 0x0, sizeof(swring->data));
        clock_gettime(CLOCK_MONOTONIC, &tv);
        swring->magic = tv.tv_nsec;
        unavi = 0;//reuse unavi_pos array
        for (int c=0; c<avi; c++) {
            uint32_t swpos = avi_pos[c];
            sw_context_t *curr = get_swctx_by_pos(swpos);
            if (0 == curr->n_regs0) {
                if (1 == curr->enable) {
                    //fprintf(stderr, "dp consumer freemem, magic %d,swpos %d, dppos %d,old enable %d, regctx %p, qseq0 %p\n",swring->magic,swpos,curr->dppos,curr->enable,curr->regctx,curr->qseq0[0]);
                    curr->enable = 0;
                    kfree(curr->b->km, curr->regctx);
                    kfree(curr->b->km, curr->qseq0[0]);
                    kfree(curr->newkm, curr->a);
                    put_free_dppos(curr->dppos);
                    put_free_swpos(swpos);
                }
                continue;
            }

            if (n_segs == 1) { // uni-segment
                    
                if ((opt->flag & MM_F_CIGAR)) {
                    sw_readhdr_t *read = &swring->data[swring->factnum];
                    read->ctxpos = swpos;
                    //fprintf(stderr, "dp consumer prepare avi %d,regseq %d,nregs %d,magic %u,swpos %d\n",avi,curr->n_regs0_swdone,curr->n_regs0,swring->magic,swpos);
                    uint32_t size = mm_align_skeleton_pre(curr->b->km, opt, curr->mi, curr->qseq0, curr->qlens[0], curr->qseqs[0], &curr->n_regs0, curr->regs0, curr->a, curr->regctx, read, swpos, curr->n_regs0_swdone, curr->n_a);
                    if (size>0) /*read->regnum = 1,*/ swring->factnum++, swring->total_size += size;

                    curr->enable = 3;
                    unavi_pos[unavi++] = swpos;
                }

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
            }//else

        }//for c<avi

        for (int k=0; k<unavi; k++) put_to_ring_buf_st(left_sw_pos[tid], (void *)unavi_pos[k]);

        if (swring->factnum > 0) sw_ring_submit(swring, snd_ctrl.sw_ring_buf);
        else kfree(newkm, swring);
       
        if (swring->factnum > 0) fprintf(stderr, "ring submit comp! magic %d,factnum %u,avi %d,unavi %d,size %d,tid %d,last %d from left\n",swring->magic,swring->factnum,avi,unavi,swring->total_size,swring->tid,swring->lastblk);
    }//if avi>0

    //restore unavi for thread exit condition
    unavi = bakunavi;
	dptest_rcvhdr_t *hdr = 0;
	get_from_ring_buf(rcv_ctrl.dp_ring_buf, (void **)&hdr);
	if (!hdr) {
        if ((0 == unavi + avi) && lastblk_comeback() && (add_rsp(0) == add_req(0))) return -1;
        else return 0;
	}
    
	add_rsp(hdr->num);
	if (IS_LAST == hdr->lat) set_last_flag(hdr->lat);//avoid rewrite by NT_LAST
    
	//fprintf(stderr, "new chaindp rett, magic %d,tid %d,num %u, flag %u,req %u,rsp %u,getlast %u, avi %d, unavi %d\n",hdr->magic,tid,hdr->num,hdr->lat,add_req(0),add_rsp(0),lastblk_comeback(),avi,unavi);
    
    unavi = 0;
	newkm = hdr->km;

	uint32_t preoff = 0;
	for (int c=0; c<hdr->num; c++) {
        
        int n_regs0 = 0;
        uint64_t *u = 0;
        
        dptest_rcvsubhdr_t *sub = (dptest_rcvsubhdr_t *)((uint8_t *)hdr->subhdr + preoff);
        preoff += sub->subsize;
        
        int64_t  n_a = sub->n_a;
        uint32_t n_u = ((pv_t*)sub->pvf)[n_a].p;
        int n_mini_pos = sub->n_mini_pos;
        
        rep_len = sub->rep_len;
        a = (mm128_t *)(sub->pvf + sub->aoffset);
        mini_pos = (uint64_t *)(sub->pvf + sub->moffset);
        
        pv_t *pv = (pv_t*)sub->pvf;
        int32_t *f = (int32_t *)(pv + n_a + 1);
        
        dp_context_t *ctx = get_dpctx_by_pos(sub->ctxpos);
        regs = ctx->regs;
        n_regs = ctx->n_regs;
        b = ctx->b;
        mv.a = ctx->mva;
        hash = ctx->hash;
        n_segs = ctx->n_segs;
        qlen_sum = ctx->qlen_sum;
        max_chain_gap_qry = ctx->gap_qry;
        max_chain_gap_ref = ctx->gap_ref;
        qlens = ctx->qlen;
        seqs = ctx->qseqs;
        mi = ctx->mi;
        qname = ctx->qname;
        
        a = mm_chain_dp_fpos(opt->min_cnt, opt->min_chain_score, n_u, n_a, f,pv,a, &n_regs0, &u, b->km);

        if ((opt->max_occ > opt->mid_occ) && rep_len > 0) {
            int rechain = 0;
            if (n_regs0 > 0) { // test if the best chain has all the segments
                int n_chained_segs = 1, max = 0, max_i = -1, max_off = -1, off = 0;
                for (i = 0; i < n_regs0; ++i) { // find the best chain
                    if (max < u[i]>>32) max = u[i]>>32, max_i = i, max_off = off;
                    off += (uint32_t)u[i];
                }
                for (i = 1; i < (uint32_t)u[max_i]; ++i) // count the number of segments in the best chain
                    if ((a[max_off+i].y&MM_SEED_SEG_MASK) != (a[max_off+i-1].y&MM_SEED_SEG_MASK)) ++n_chained_segs;
                if (n_chained_segs < n_segs) rechain = 1;
            } else rechain = 1;
            if (rechain) { // redo chaining with a higher max_occ threshold
                if (opt->flag & MM_F_HEAP_SORT) a = collect_seed_hits_heap(b->km, opt, opt->max_occ, mi, qname, &mv, qlen_sum, &n_a, &rep_len, &n_mini_pos, &mini_pos);
                else a = collect_seed_hits(b->km, opt, opt->max_occ, mi, qname, &mv, qlen_sum, &n_a, &rep_len, &n_mini_pos, &mini_pos);
                a = mm_chain_dp(max_chain_gap_ref, max_chain_gap_qry, opt->bw, opt->max_chain_skip, opt->min_cnt, opt->min_chain_score, is_splice, n_segs, n_a, a, &n_regs0, &u, b->km);
            }
        }

        regs0 = mm_gen_regs(b->km, hash, qlen_sum, n_regs0, u, a);

        chain_post(opt, max_chain_gap_ref, mi, b->km, qlen_sum, n_segs, qlens, &n_regs0, regs0, a);
        if (!is_sr) mm_est_err(mi, qlen_sum, n_regs0, regs0, a, n_mini_pos, mini_pos);

        if (n_segs == 1) { // uni-segment
            //get ctxpos and save ctx
            int i, n_a;
            uint8_t *qseq0[2];
            uint32_t swpos = get_free_swpos();
            sw_context_t *curr = get_swctx_by_pos(swpos);
            curr->enable = 1;
            curr->regctx = 0;
            curr->n_regs0_cap = 3*n_regs0;
            if (n_regs0 > 0) {
                curr->regctx = (reg_context_t *)kmalloc(b->km, curr->n_regs0_cap*sizeof(reg_context_t));
                memset(curr->regctx, 0xdeaddead, sizeof(curr->n_regs0_cap*sizeof(reg_context_t)));
            }

            unavi_pos[unavi++] = swpos;

            qseq0[0] = (uint8_t*)kmalloc(b->km, qlens[0] * 2);
            qseq0[1] = qseq0[0] + qlens[0];
            for (i = 0; i < qlens[0]; ++i) {
                qseq0[0][i] = seq_nt4_table[(uint8_t)seqs[0][i]];
                qseq0[1][qlens[0] - 1 - i] = qseq0[0][i] < 4? 3 - qseq0[0][i] : 4;
            }

            //fprintf(stderr, "init, tid %d, readseq-c %d,swpos %d,regnum %d,regctx %p,qseq0[0] %p,qname:[%s]\n", tid, c,swpos,n_regs0,curr->regctx,qseq0[0],ctx->qname);
	    
            n_a = mm_squeeze_a(b->km, n_regs0, regs0, a);
            save_swctx(swpos, newkm, mi, n_regs, regs, regs0, a, b, qlens, qseq0, seqs, n_regs0, rep_len, n_a, sub->ctxpos);
        }
	}//end first for()
    
    for (int k=0; k<unavi; k++) {
        int idx = put_to_ring_buf_st(left_sw_pos[tid], (void *)unavi_pos[k]);
        if (idx >= 0) continue;
        for (;;) {
            int *pos = (int*)0xdead;
            int idx  = get_from_ring_buf_st(left_sw_pos[tid], (void**)&pos);
            if (idx >= 0) {
                int realpos = (int)pos;
                sw_context_t *curr = get_swctx_by_pos(realpos);
                if (0 == curr->enable) {
                    put_free_dppos(curr->dppos); 
                    put_free_swpos((int)pos); 
                    realpos = unavi_pos[k];
                }

                  put_to_ring_buf_st(left_sw_pos[tid], (void *)realpos);
                  if (0 == curr->enable) break;
              }//if idx>=0
          }//for(;;)
      }//for k<unavi
      
	  //if (swring->factnum > 0) fprintf(stderr, "ring submit comp! magic %d,km %p,ring %p,hdrnum %d,num %u,size %d,tid %d,last %d, unavi %d\n",swring->magic,swring->km,swring,hdr->num,swring->factnum,swring->total_size,swring->tid,swring->lastblk,unavi);
#ifdef DEBUG_API
	  //debug
	  sw_task_sender(&snd_ctrl);
	  sw_consumer(opt, tid);
#endif
      
      
      //fprintf(stderr, "dp consumer done! last %d,tid %d,put %d new reads to left ring\n",hdr->lat&IS_LAST,tid,unavi);
      kfree(newkm, hdr);
	
      return 0;
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
	  for (j = 0; j < s->n_seg[i]; ++j) {
	    mm_map_frag(s->p->mi, 1, &qlens[j], &qseqs[j], &s->n_reg[off+j], &s->reg[off+j], b, s->p->opt, s->seq[off+j].name);
	    if (b->km) {
	      km_stat_t kmst;
	      km_stat(b->km, &kmst);
	      //assert(kmst.n_blocks == kmst.n_cores); // otherwise, there is a memory leak
	      
	      if (kmst.largest > 1U<<28) {
		km_destroy(b->km);
		b->km = km_init();
	      }
	    }

	  }
	} else {
		mm_map_frag(s->p->mi, s->n_seg[i], qlens, qseqs, &s->n_reg[off], &s->reg[off], b, s->p->opt, s->seq[off].name);
		if (b->km) {
		  km_stat_t kmst;
		  km_stat(b->km, &kmst);
		  
		  //assert(kmst.n_blocks == kmst.n_cores); // otherwise, there is a memory leak
		  
		  if (kmst.largest > 1U<<28) {
		    km_destroy(b->km);
		    b->km = km_init();
		  }
		}

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
        set_reqrsp_zero();
		set_reqrsp_zero_sw();
		set_stat_dpthrd(0);

		pthread_t *tid0=0;
		void *tw0=0;
		void *kf0=0;
		kt_for_nowait(p->n_threads, worker_for, in, ((step_t*)in)->n_frag, chaindp_put_last, &tid0, &tw0, &kf0);

		pthread_t *tid1=0;
		void *tw1=0;
		void *kf1=0;
		kt_for_nowait_clone(p->n_threads, dp_consumer, (void *)p->opt, &tid1, &tw1, &kf1);
#ifndef DEBUG_API
		pthread_t *tid2=0;
		void *tw2=0;
		void *kf2=0;
		kt_for_nowait_clone(p->n_threads, sw_consumer, (void *)p->opt, &tid2, &tw2, &kf2);
#endif
		if (tid0 && tw0 && kf0) {
			for (int i = 0; i < p->n_threads; ++i) pthread_join(tid0[i], 0);
			free(tid0); free(tw0); free(kf0);
		}

		if (tid1 && tw1 && kf1) {
			for (int i = 0; i < p->n_threads; ++i) pthread_join(tid1[i], 0);
			free(tid1); free(tw1); free(kf1);
		}
		set_stat_dpthrd(1);
		fprintf(stderr, "dp consumer exit!\n");
#ifndef DEBUG_API
		if (tid2 && tw2 && kf2) {
			for (int i = 0; i < p->n_threads; ++i) pthread_join(tid2[i], 0);
			free(tid2); free(tw2); free(kf2);
		}
#endif
		fprintf(stderr, "sw consumer exit!all thread exit in step1\n\n\n");
		set_last_flag(0);
		set_last_flag_sw(0);

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
		if (0)//(mm_verbose >= 3)
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
	kt_pipeline(pl_threads, worker_pipeline, &pl, 3);
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
