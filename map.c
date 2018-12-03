#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <pthread.h>
#include <errno.h>
#include <unistd.h> 
#include <fcntl.h> 

#include "kthread.h"
#include "kvec.h"
#include "kalloc.h"
#include "sdust.h"
#include "mmpriv.h"
#include "bseq.h"
#include "khash.h"


static void* sw_result_thread(void* arg);

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

static mm128_t *collect_seed_hits(void *km, const mm_mapopt_t *opt, int max_occ, const mm_idx_t *mi, const char *qname, const mm128_v *mv, int qlen, int64_t *n_a, int *rep_len,
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

static void align_regs(const mm_mapopt_t *opt, const mm_idx_t *mi, void *km, int qlen, const char *seq, int *n_regs, mm_reg1_t *regs, mm128_t *a, long read_index, context_t* context, user_params_t* params)
{
	if (!(opt->flag & MM_F_CIGAR)) return;
	mm_align_skeleton(km, opt, mi, qlen, seq, n_regs, regs, a, read_index, context, params); // this calls mm_filter_regs()
	/*if (!(opt->flag & MM_F_ALL_CHAINS)) { // don't choose primary mapping(s)
		mm_set_parent(km, opt->mask_level, *n_regs, regs, opt->a * 2 + opt->b);
		mm_select_sub(km, opt->pri_ratio, mi->k*2, opt->best_n, n_regs, regs);
		mm_set_sam_pri(*n_regs, regs);
	}*/
	return;
}

static mm_reg1_t *align_regs_ori(const mm_mapopt_t *opt, const mm_idx_t *mi, void *km, int qlen, const char *seq, int *n_regs, mm_reg1_t *regs, int n_a, mm128_t *a)
{
	if (!(opt->flag & MM_F_CIGAR)) return regs;
	regs = mm_align_skeleton_ori(km, opt, mi, qlen, seq, n_regs, regs, n_a, a); // this calls mm_filter_regs()
	if (!(opt->flag & MM_F_ALL_CHAINS)) { // don't choose primary mapping(s)
		mm_set_parent(km, opt->mask_level, *n_regs, regs, opt->a * 2 + opt->b);
		mm_select_sub(km, opt->pri_ratio, mi->k*2, opt->best_n, n_regs, regs);
		mm_set_sam_pri(*n_regs, regs);
	}
	return regs;
}


void mm_map_frag(const mm_idx_t *mi, int n_segs, const int *qlens, const char **seqs, int *n_regs, mm_reg1_t **regs, mm_tbuf_t *b, const mm_mapopt_t *opt, const char *qname, long read_index, user_params_t* params)
{
	int i, j, rep_len, qlen_sum, n_regs0, n_mini_pos;
	int max_chain_gap_qry, max_chain_gap_ref, is_splice = !!(opt->flag & MM_F_SPLICE), is_sr = !!(opt->flag & MM_F_SR);
	uint32_t hash;
	int64_t n_a;
	uint64_t *u, *mini_pos;
	mm128_t *a;
	mm128_v mv = {0,0,0};
	mm_reg1_t *regs0;
	//km_stat_t kmst;

	for (i = 0, qlen_sum = 0; i < n_segs; ++i)
		qlen_sum += qlens[i], n_regs[i] = 0, regs[i] = 0;

	if (qlen_sum == 0 || n_segs <= 0 || n_segs > MM_MAX_SEG) return;

	hash  = qname? __ac_X31_hash_string(qname) : 0;
	hash ^= __ac_Wang_hash(qlen_sum) + __ac_Wang_hash(opt->seed);
	hash  = __ac_Wang_hash(hash);

	collect_minimizers(b->km, opt, mi, n_segs, qlens, seqs, &mv);
	if (opt->flag & MM_F_HEAP_SORT) a = collect_seed_hits_heap(b->km, opt, opt->mid_occ, mi, qname, &mv, qlen_sum, &n_a, &rep_len, &n_mini_pos, &mini_pos);
	else a = collect_seed_hits(b->km, opt, opt->mid_occ, mi, qname, &mv, qlen_sum, &n_a, &rep_len, &n_mini_pos, &mini_pos);

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

	a = mm_chain_dp(max_chain_gap_ref, max_chain_gap_qry, opt->bw, opt->max_chain_skip, opt->min_cnt, opt->min_chain_score, is_splice, n_segs, n_a, a, &n_regs0, &u, b->km);

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
			if (opt->flag & MM_F_HEAP_SORT) a = collect_seed_hits_heap(b->km, opt, opt->max_occ, mi, qname, &mv, qlen_sum, &n_a, &rep_len, &n_mini_pos, &mini_pos);
			else a = collect_seed_hits(b->km, opt, opt->max_occ, mi, qname, &mv, qlen_sum, &n_a, &rep_len, &n_mini_pos, &mini_pos);
			a = mm_chain_dp(max_chain_gap_ref, max_chain_gap_qry, opt->bw, opt->max_chain_skip, opt->min_cnt, opt->min_chain_score, is_splice, n_segs, n_a, a, &n_regs0, &u, b->km);
		}
	}

	regs0 = mm_gen_regs(b->km, hash, qlen_sum, n_regs0, u, a);

	if (mm_dbg_flag & MM_DBG_PRINT_SEED)
		for (j = 0; j < n_regs0; ++j)
			for (i = regs0[j].as; i < regs0[j].as + regs0[j].cnt; ++i)
				fprintf(stderr, "CN\t%d\t%s\t%d\t%c\t%d\t%d\t%d\n", j, mi->seq[a[i].x<<1>>33].name, (int32_t)a[i].x, "+-"[a[i].x>>63], (int32_t)a[i].y, (int32_t)(a[i].y>>32&0xff),
						i == regs0[j].as? 0 : ((int32_t)a[i].y - (int32_t)a[i-1].y) - ((int32_t)a[i].x - (int32_t)a[i-1].x));

	chain_post(opt, max_chain_gap_ref, mi, b->km, qlen_sum, n_segs, qlens, &n_regs0, regs0, a);
	if (!is_sr) mm_est_err(mi, qlen_sum, n_regs0, regs0, a, n_mini_pos, mini_pos);

    kfree(b->km, mv.a);
	kfree(b->km, u);
	kfree(b->km, mini_pos);

	if (n_segs == 1) { // uni-segment
        //创建read上下文并保存在全局数组中
        if(params->read_contexts[read_index] != NULL) {
            fprintf(stderr, "%ld read context is not NULL\n", read_index);
            exit(1);
        }
        context_t* context = create_context(read_index);
        params->read_contexts[read_index] = context;
        context->b = b;
        context->opt = opt;
        context->mi = mi;
        context->qlen = qlens[0];
        context->seq = (char*)malloc(qlens[0]);
        memcpy(context->seq, seqs[0], qlens[0]);
        context->read_index = read_index;
        context->rep_len = rep_len;
        context->n_regs = n_regs;
        context->regs = regs;

		align_regs(opt, mi, b->km, qlens[0], seqs[0], &n_regs0, regs0, a, read_index, context, params);
		//mm_set_mapq(b->km, n_regs0, regs0, opt->min_chain_score, opt->a, rep_len, is_sr);
		//n_regs[0] = n_regs0, regs[0] = regs0;
	} else { // multi-segment
		/*mm_seg_t *seg;
		seg = mm_seg_gen(b->km, hash, n_segs, qlens, n_regs0, regs0, n_regs, regs, a); // split fragment chain to separate segment chains
		free(regs0);
		for (i = 0; i < n_segs; ++i) {
			mm_set_parent(b->km, opt->mask_level, n_regs[i], regs[i], opt->a * 2 + opt->b); // update mm_reg1_t::parent
			regs[i] = align_regs(opt, mi, b->km, qlens[i], seqs[i], &n_regs[i], regs[i], seg[i].a, read_index);
			mm_set_mapq(b->km, n_regs[i], regs[i], opt->min_chain_score, opt->a, rep_len, is_sr);
		}
		mm_seg_free(b->km, n_segs, seg);
		if (n_segs == 2 && opt->pe_ori >= 0 && (opt->flag&MM_F_CIGAR))
			mm_pair(b->km, max_chain_gap_ref, opt->pe_bonus, opt->a * 2 + opt->b, opt->a, qlens, n_regs, regs); // pairing
        */
    }

	//kfree(b->km, mv.a);       //可以提前释放
	//kfree(b->km, a);          //添加到上下文中，处理完结果后才能释放
	//kfree(b->km, u);          //可以提前释放
	//kfree(b->km, mini_pos);   //可以提前释放

	/*if (b->km) {      //添加到上下文，处理完结果再处理
		km_stat(b->km, &kmst);
		if (mm_dbg_flag & MM_DBG_PRINT_QNAME)
			fprintf(stderr, "QM\t%s\t%d\tcap=%ld,nCore=%ld,largest=%ld\n", qname, qlen_sum, kmst.capacity, kmst.n_cores, kmst.largest);
		assert(kmst.n_blocks == kmst.n_cores); // otherwise, there is a memory leak
		if (kmst.largest > 1U<<28) {
			km_destroy(b->km);
			b->km = km_init();
		}
	}*/
}

mm_reg1_t *mm_map(const mm_idx_t *mi, int qlen, const char *seq, int *n_regs, mm_tbuf_t *b, const mm_mapopt_t *opt, const char *qname)
{
	mm_reg1_t *regs;
	mm_map_frag(mi, 1, &qlen, &seq, n_regs, &regs, b, opt, qname, -1, NULL);
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

static void worker_for(void *_data, long i, int tid, void* params) // kt_for() callback
{
    step_t *s = (step_t*)_data;
	int qlens[MM_MAX_SEG], j, off = s->seg_off[i], pe_ori = s->p->opt->pe_ori;
	const char *qseqs[MM_MAX_SEG];
	mm_tbuf_t *b = s->buf[tid];
    user_params_t* user_params = (user_params_t*)params;
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
			mm_map_frag(s->p->mi, 1, &qlens[j], &qseqs[j], &s->n_reg[off+j], &s->reg[off+j], b, s->p->opt, s->seq[off+j].name, i, user_params);
	} else {
		mm_map_frag(s->p->mi, s->n_seg[i], qlens, qseqs, &s->n_reg[off], &s->reg[off], b, s->p->opt, s->seq[off].name, i, user_params);
	}

    //这些代码目前不会进入
	/*for (j = 0; j < s->n_seg[i]; ++j) // flip the query strand and coordinate to the original read strand
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
		}*/
}

typedef struct {
    int tid;
    user_params_t* params;
}user_args_t;

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
        user_params_t params;
        memset(&params, 0, sizeof(params));

        
        params.read_num = ((step_t*)in)->n_frag;
        long tmp_read_num = params.read_num;

        params.read_flag = (char*)malloc(params.read_num);
        memset(params.read_flag, 0, params.read_num);

        fprintf(stderr, "read_num=%ld\n", params.read_num);
        params.read_contexts = (context_t**)malloc(params.read_num * sizeof(context_t*));
        memset(params.read_contexts, 0, params.read_num * sizeof(context_t*));

        params.chain_num = (int*)malloc(params.read_num * sizeof(int));
        memset(params.chain_num, 0xff, params.read_num * sizeof(int));

        params.read_is_complete = (char*)malloc(params.read_num);
        memset(params.read_is_complete, 0, params.read_num);

        params.read_results = (read_result_t*)malloc(params.read_num * sizeof(read_result_t));
        memset(params.read_results, 0, params.read_num * sizeof(read_result_t));


        int i = 0;
        
        if(result_empty() == 0) {
            fprintf(stderr, "result is empty\n");
        }
        else {
            fprintf(stderr, "result is not empty, have %d results\n", result_empty());
        }
        //TODO 这里还需要清理残留的结果及结果的上下文数据
        init_result_array();

        //TODO 处理结果的线程
        user_args_t user_args[10];
        pthread_t sw_tid[10];
        for(i = 0; i < 10; i++) {
            user_args[i].tid = i;
            user_args[i].params = &params;
            pthread_create(&sw_tid[i], NULL, sw_result_thread, &user_args[i]);
        }

		kt_for_map(p->n_threads, worker_for, in, ((step_t*)in)->n_frag, (void *)&params);

        for(i = 0; i < 10; i++)
            pthread_join(sw_tid[i], NULL);

        if(params.read_num == 0) {
            fprintf(stderr, "all read ok!\n");
        }
        else {
            fprintf(stderr, "%ld read do not process!\n", params.read_num);
        }
        
        free(params.read_results);
        free(params.read_is_complete);
        free(params.read_contexts);
        free(params.chain_num);
        free(params.read_flag);
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

void reverse_cigar(int n_cigar, uint32_t *cigar)
{
    uint32_t tmp;
    int i = 0;
    for (i=0; i<n_cigar>>1; i++)
        tmp = cigar[i],cigar[i]=cigar[n_cigar-1-i],cigar[n_cigar-1-i]=tmp;
    return;
}
void save_result(int count, long read_id, int chain_id, sw_result_t* result)
{
    int i;
    char* buf = NULL;
    sw_result_t* p_result = NULL;
    ksw_extz_t* p_ez = NULL;
    int* id = NULL;
    uint32_t* p_cigar = NULL;
    int size = 0;
    int ret = 0;

    char* file_name = (char*)malloc(100);
    memset(file_name, 0, 100);
    sprintf(file_name, "./10_th_1/result_%d_read_%ld_chain_%d.bin", count, read_id, chain_id);
    
    if((access(file_name,F_OK))!=-1)   
    {   
        fprintf(stderr, "file %s exist\n", file_name);  
        exit(0);
    }  
    
    buf = (char*)malloc(4*1024*1024);
    if(buf == NULL) {
        fprintf(stderr, "malloc failed:%s\n", strerror(errno));
        exit(0);
    }
    id = (int*)buf;
    *id = 0x1234abcd;
    p_result = (sw_result_t*)(buf + sizeof(int));
    *p_result = *result;
    p_ez = (ksw_extz_t*)((char*)p_result + sizeof(sw_result_t));
    p_cigar = (uint32_t*)((char*)p_ez + result->result_num * sizeof(ksw_extz_t));
    
    for(i = 0; i < result->result_num; i++) {
        p_ez[i] = *result->ezs[i];
        memcpy(p_cigar, result->ezs[i]->cigar, result->ezs[i]->n_cigar * sizeof(uint32_t));
        p_cigar += result->ezs[i]->n_cigar;
    }
    
    size = (char*)p_cigar - buf;
    
    FILE* fp = fopen(file_name, "a");
    if(fp == NULL) {
        fprintf(stderr, "fopen %s failed:%s\n", file_name, strerror(errno));
        exit(0);
    }
    
    ret = fwrite(buf, 1, size, fp);
    if(ret != size) {
        fprintf(stderr, "fwrite failed:%s\n", strerror(errno));
        exit(0);
    }
    
    free(buf);
    fclose(fp);
    free(file_name);
    return;
}
static int save_read_result(long read_id, read_result_t* read_result, context_t* context, int num, char* read_flag, int* chain_num)
{
    int ret = 0;
    int k = 0;
    int chain_i = 0;
    //fprintf(stderr, "read_id=%ld\n", read_id);
    for(chain_i = 0; chain_i < num; chain_i++) {
        chain_context_t* chain_context = context->chain_contexts[chain_i];
        sw_context_t** sw_contexts = chain_context->sw_contexts;
        sw_result_t* result = read_result->chain_results[chain_i];
        
        //save_result(g_count, read_id, chain_i, result);

        int result_num = result->result_num;
        
        int32_t rs1, qs1, re1, qe1;
        int32_t dropped = 0;
        mm_reg1_t r2;
        //r2.cnt = 0;
        memset(&r2, 0, sizeof(r2));

        mm_reg1_t *r = &(context->regs0[chain_context->i]);
        mm_reg1_t *regs = context->regs0;

        //取出read上下文
        const mm_idx_t *mi = context->mi;
        //void* km = context->b->km;
        void* km = NULL;
        const mm_mapopt_t *opt = context->opt;
        mm128_t *a = context->a;
        int qlen = context->qlen;
        int8_t* mat = context->mat;
        
        int32_t rs = chain_context->rs;
        int32_t qs = chain_context->qs;
        int32_t qs0 = chain_context->qs0;
        int32_t rev = chain_context->rev;

        int is_sr = !!(opt->flag & MM_F_SR);
        
        
        
        //fprintf(stderr, "n_regs=%d, i=%d, read_index=%ld\n", context->n_regs0, chain_context->i, read_index);
        if(sw_contexts[0]->pos_flag == 0) {  //left
            result_num--;
            k = 1;  //有左扩展
            ksw_extz_t* ez = result->ezs[0];
            sw_context_t* sw_context = sw_contexts[0];
            rs = sw_context->rs;
            qs = sw_context->qs;
            qs0 = sw_context->qs0;
            //fprintf(stderr, "1.qlen=%d, tlen=%d, w=%d\n", sw_context->qlen, sw_context->tlen, sw_context->w);
            
            if (ez->n_cigar > 0) {
#if 0
                if (0xA0B1C2D3 == ez->revcigar) {
                    reverse_cigar(ez->n_cigar, ez->cigar);
                }
#endif
                mm_append_cigar(r, ez->n_cigar, ez->cigar);
                r->p->dp_score += ez->max;
            }
            rs1 = rs - (ez->reach_end? ez->mqe_t + 1 : ez->max_t + 1);
            qs1 = qs - (ez->reach_end? qs - qs0 : ez->max_q + 1);
        }
        else {
            k = 0;  //没有左扩展
            rs1 = rs, qs1 = qs;
        }
        
        re1 = rs, qe1 = qs;
        assert(qs1 >= 0 && rs1 >= 0);

        for(; k < result->result_num; k++) {    //gap
            if(sw_contexts[k]->pos_flag == 1) {
                result_num--;
                int zdrop_code;
                ksw_extz_t tmp_ez;
                memset(&tmp_ez, 0, sizeof(tmp_ez));
                ksw_extz_t* ez = result->ezs[k];
                sw_context_t* sw_context = sw_contexts[k];
                uint32_t i = sw_context->i;
                int j = 0;
                uint8_t* qseq = sw_context->query;
                uint8_t* tseq = sw_context->target;
                int32_t cnt1 = sw_context->cnt1;
                int32_t re = sw_context->re;
                int32_t qe = sw_context->qe;
                int32_t as1 = sw_context->as1;

                rs = sw_context->rs;
                qs = sw_context->qs;
                qs0 = sw_context->qs0;
                //fprintf(stderr, "2.qlen=%d, tlen=%d, w=%d\n", sw_context->qlen, sw_context->tlen, sw_context->w);

                if (ez->n_cigar > 0) {
#if 0
                    if (0xA0B1C2D3 == ez->revcigar) {
                        reverse_cigar(ez->n_cigar, ez->cigar);
                    }
#endif
                }

                if ((zdrop_code = mm_test_zdrop(km, opt, qseq, tseq, ez->n_cigar, ez->cigar, mat)) != 0) {
                    mm_align_pair(km, opt, sw_context->qlen, qseq, sw_context->tlen, tseq, mat, sw_context->w, -1, zdrop_code == 2? opt->zdrop_inv : opt->zdrop, sw_context->zdrop_flag, &tmp_ez); // second pass: lift approximate
                    ez = &tmp_ez;
                }
                if (ez->n_cigar > 0)
                    mm_append_cigar(r, ez->n_cigar, ez->cigar);
                if (ez->zdropped) { // truncated by Z-drop; TODO: sometimes Z-drop kicks in because the next seed placement is wrong. This can be fixed in principle.
                    for (j = i - 1; j >= 0; --j) {
                        //fprintf(stderr, "a=%p as1=%d j=%d a[as1 + j].x=%d rs=%d rs+ez->max_t=%d\n", a, as1, j, (int32_t)a[as1 + j].x, rs, rs + ez->max_t);
                        if ((int32_t)a[as1 + j].x <= rs + ez->max_t){
                            if(zdrop_code != 0) {
                                if(tmp_ez.cigar != NULL) {
                                    kfree(km, tmp_ez.cigar);
                                    tmp_ez.cigar = NULL;
                                }
                            }
                            break;
                        }
                    }
                    dropped = 1;

                    if (j < 0) j = 0;
                    r->p->dp_score += ez->max;
                    re1 = rs + (ez->max_t + 1);
                    qe1 = qs + (ez->max_q + 1);
                    if (cnt1 - (j + 1) >= opt->min_cnt) {
                        mm_split_reg(r, &r2, as1 + j + 1 - r->as, qlen, a);
                        if (zdrop_code == 2) r2.split_inv = 1;
                    }
                    if(zdrop_code != 0) {
                        if(tmp_ez.cigar != NULL) {
                            kfree(km, tmp_ez.cigar);
                            tmp_ez.cigar = NULL;
                        }
                    }
                    break;
                } else r->p->dp_score += ez->score;
                rs = re, qs = qe;

                if(zdrop_code != 0) {
                    if(tmp_ez.cigar != NULL) {
                        kfree(km, tmp_ez.cigar);
                        tmp_ez.cigar = NULL;
                    }
                }
            }
            else if(sw_contexts[k]->pos_flag == 2) {     //right
                break;
            }
            else {
                fprintf(stderr, "error position flag(%d), read_index:%ld\n", sw_contexts[k]->pos_flag, context->read_index);
            }
        }
        
        if(sw_contexts[k]->pos_flag == 2) {      //right
            if(!dropped) {
                result_num--;
                ksw_extz_t* ez = result->ezs[k];
                sw_context_t* sw_context = sw_contexts[k];
                int32_t re = sw_context->re;
                int32_t qe = sw_context->qe;
                int32_t qe0 = sw_context->qe0;

                rs = sw_context->rs;
                qs = sw_context->qs;
                qs0 = sw_context->qs0;
                //fprintf(stderr, "3.qlen=%d, tlen=%d, w=%d\n", sw_context->qlen, sw_context->tlen, sw_context->w);
                
                if (ez->n_cigar > 0) {
#if 0
                    if (0xA0B1C2D3 == ez->revcigar) {
                        reverse_cigar(ez->n_cigar, ez->cigar);
                    }
#endif
                    mm_append_cigar(r, ez->n_cigar, ez->cigar);
                    r->p->dp_score += ez->max;
                }
                re1 = re + (ez->reach_end? ez->mqe_t + 1 : ez->max_t + 1);
                qe1 = qe + (ez->reach_end? qe0 - qe : ez->max_q + 1);
            }
        }

        assert(qe1 <= qlen);
        r->rs = rs1, r->re = re1;
        if (rev) r->qs = qlen - qe1, r->qe = qlen - qs1;
        else r->qs = qs1, r->qe = qe1;
        
        
        assert(re1 - rs1 <= chain_context->re0 - chain_context->rs0);
        if (r->p) {
            mm_idx_getseq(mi, chain_context->rid, rs1, re1, chain_context->tseq);
            mm_update_extra(r, &chain_context->qseq0[r->rev][qs1], chain_context->tseq, mat, opt->q, opt->e);
            if (rev && r->p->trans_strand)
                r->p->trans_strand ^= 3; // flip to the read strand
        }

        if(r2.cnt > 0) {   //此时这条read需要重新做sw
            //fprintf(stderr, "1.insert chain\n");
            *read_flag = 1;
            int n_regs0 = context->n_regs0;
            context->regs0_ori = align_regs_ori(opt, mi, km, qlen, context->seq, &n_regs0, context->regs0_ori, context->n_a, context->a_ori);
            mm_set_mapq(km, n_regs0, context->regs0_ori, opt->min_chain_score, opt->a, context->rep_len, is_sr);
            *(context->n_regs) = n_regs0;           //保存结果
            context->regs[0] = context->regs0_ori;  //保存结果
            //释放该read所有上下文内容
            return 1;   //返回1表示该条read成功处理
        }
        if (chain_context->i > 0 && regs[chain_context->i].split_inv) {
            ksw_extz_t tmp_ez;
            memset(&tmp_ez, 0, sizeof(tmp_ez));
            if (mm_align1_inv(km, opt, mi, qlen, chain_context->qseq0, &regs[chain_context->i-1], &regs[chain_context->i], &r2, &tmp_ez)) {
                //此时这条read需要重新做sw
                //fprintf(stderr, "2.insert chain\n");
                *read_flag = 1;
                int n_regs0 = context->n_regs0;
                context->regs0_ori = align_regs_ori(opt, mi, km, qlen, context->seq, &n_regs0, context->regs0_ori, context->n_a, context->a_ori);
                mm_set_mapq(km, n_regs0, context->regs0_ori, opt->min_chain_score, opt->a, context->rep_len, is_sr);
                *(context->n_regs) = n_regs0;           //保存结果
                context->regs[0] = context->regs0_ori;  //保存结果
                //释放该read所有上下文内容

                kfree(km, tmp_ez.cigar);
                return 1;   //返回1表示该条read成功处理
            }
            kfree(km, tmp_ez.cigar);
        }
        else {  //判断是否整条read做完的
            *chain_num = *chain_num - 1;
            //fprintf(stderr, "1.%ld read done, shengyude read %ld\n", read_index, params->read_num);
            if(*chain_num == 0 && *read_flag == 0) {
                mm_filter_regs(km, opt, qlen, &context->n_regs0, regs);
                mm_hit_sort_by_dp(km, &context->n_regs0, regs);
                if (!(opt->flag & MM_F_ALL_CHAINS)) { // don't choose primary mapping(s)
                    mm_set_parent(km, opt->mask_level, context->n_regs0, regs, opt->a * 2 + opt->b);
                    mm_select_sub(km, opt->pri_ratio, mi->k*2, opt->best_n, &context->n_regs0, regs);
                    mm_set_sam_pri(context->n_regs0, regs);
                }
                mm_set_mapq(km, context->n_regs0, regs, opt->min_chain_score, opt->a, context->rep_len, is_sr);
                *(context->n_regs) = context->n_regs0;           //保存结果
                context->regs[0] = regs;  //保存结果
                ret = 1;
            }
        }

        destroy_results(result);
        //销毁结果数组和所有上下文
        destroy_chain_context(chain_context);
        context->chain_contexts[chain_i] = NULL;
    }
    return ret;
}
/*
void save_chain_num(long read_id, int chain_num)
{
    int ret;
    char* file_name = (char*)malloc(100);
    memset(file_name, 0, 100);
    sprintf(file_name, "./1th_chain_num/%d_read_%ld.bin", g_count, read_id);
    char buf[100];
    memset(buf, 0, 100);
    
    if((access(file_name,F_OK))!=-1)
    {   
        fprintf(stderr, "file %s exist\n", file_name);  
        exit(0);
    }
    FILE* fp = fopen(file_name, "a");
    if(fp == NULL) {
        fprintf(stderr, "fopen %s failed\n", file_name);
        exit(0);
    }
    snprintf(buf, 100, "%d", chain_num);
    int len = strlen(buf);
    ret = fwrite(&buf, 1, len, fp);
    if(ret != len) {
        fprintf(stderr, "fwrite %s failed\n", file_name);
        exit(0);
    }
    fclose(fp);
    free(file_name);
    return;
}
*/
static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
long recv_task = 0;
static void* sw_result_thread(void* arg)
{
    user_args_t* args = (user_args_t*)arg;
    user_params_t *params = args->params;
    long read_num = params->read_num;
    while(1) {
        sw_result_t* result = get_result();
        if(result == NULL) {
            if(params->exit == 1) {
                fprintf(stderr, "1.thread exit\n");
                return NULL;
            }
            continue;
        }

        pthread_mutex_lock(&mutex);
        recv_task++;
        pthread_mutex_unlock(&mutex);
        
        long read_index = result->read_id;
        int chain_id = result->chain_id;

        //将结果放入对应的位置

        if(params->read_results[read_index].chain_results[chain_id] != NULL) {
            fprintf(stderr, "[%s][%d]read_index=%ld, chain_id=%d\n", __FILE__, __LINE__, read_index, chain_id);
            exit(0);
        }
        params->read_results[read_index].chain_results[chain_id] = result;
        
        int chain_result_num = __sync_add_and_fetch(&params->read_results[read_index].chain_result_num, 1);    //加1然后返回加1后的结果
        
        //判断整条read的chain的结果都返回
        if(chain_result_num == params->chain_num[read_index]) {
            //save_chain_num(read_index, chain_result_num);
            //处理一条read的结果
            context_t* context = params->read_contexts[read_index];
            int ret = save_read_result(read_index, &params->read_results[read_index], context, params->chain_num[read_index],
                                       &params->read_flag[read_index], &params->chain_num[read_index]);
            if(ret == 1) {
                free(context);
                params->read_contexts[read_index] = NULL;
                if(params->read_is_complete[read_index] == 1) {
                    fprintf(stderr, "read(%ld) already complete\n", read_index);
                    exit(0);
                }
                params->read_is_complete[read_index] = 1;
                long num = __sync_sub_and_fetch(&params->read_num, 1);

                if(num == 0) {  //判断所有read是否做完
                    long i = 0;
                    int complete = 0;
                    for(i = 0; i < read_num; i++) {     //再检查所有read是否做完
                        if(params->read_is_complete[i] == 0) {
                            complete = 1;
                        }
                    }
                    if(complete == 1) {
                        fprintf(stderr, "some read have not complete\n");
                    }
                    else {
                        sleep(3);
                        
                        fprintf(stderr, "2.thread exit\n");
                        params->exit = 1;
                        fprintf(stderr, "sw result thread exit, zero_seed_num=%ld read num=%ld\n", params->zero_seed_num, params->read_num);
                        return NULL;
                    }
                }
            }
            else {
                fprintf(stderr, "save_read_result ret is not 1\n");
                exit(0);
            }
        }
        else if(chain_result_num > params->chain_num[read_index]) {
            fprintf(stderr, "chain_result_num=%d, params->chain_num[%ld]=%d\n", chain_result_num, read_index, params->chain_num[read_index]);
            exit(0);
        }
        else {  //继续接收结果
            continue;
        }
    }
    return NULL;
}
