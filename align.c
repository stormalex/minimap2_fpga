#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <errno.h>
#include <unistd.h>

#include "minimap.h"
#include "mmpriv.h"
#include "ksw2.h"

#include "soft_sw.h"

#include <sys/time.h>
#include<time.h>
extern double process_task_time[100];
extern double sw_task_time[100];
extern double sw_soft_time[100];
static double realtime_msec(void)
{
	/*struct timeval tp;
	struct timezone tzp;
	gettimeofday(&tp, &tzp);
	return tp.tv_sec*1000 + tp.tv_usec * 1e-3;*/

    struct timespec tp;
    clock_gettime(CLOCK_MONOTONIC, &tp);
    return tp.tv_sec*1000 + tp.tv_nsec*1e-6;
}

static void ksw_gen_simple_mat(int m, int8_t *mat, int8_t a, int8_t b)
{
	int i, j;
	a = a < 0? -a : a;
	b = b > 0? -b : b;
	for (i = 0; i < m - 1; ++i) {
		for (j = 0; j < m - 1; ++j)
			mat[i * m + j] = i == j? a : b;
		mat[i * m + m - 1] = 0;
	}
	for (j = 0; j < m; ++j)
		mat[(m - 1) * m + j] = 0;
}

static inline void mm_seq_rev(uint32_t len, uint8_t *seq)
{
	uint32_t i;
	uint8_t t;
	for (i = 0; i < len>>1; ++i)
		t = seq[i], seq[i] = seq[len - 1 - i], seq[len - 1 - i] = t;
}

static inline void update_max_zdrop(int32_t score, int i, int j, int32_t *max, int *max_i, int *max_j, int e, int *max_zdrop, int pos[2][2])
{
	if (score < *max) {
		int li = i - *max_i;
		int lj = j - *max_j;
		int diff = li > lj? li - lj : lj - li;
		int z = *max - score - diff * e;
		if (z > *max_zdrop) {
			*max_zdrop = z;
			pos[0][0] = *max_i, pos[0][1] = i + 1;
			pos[1][0] = *max_j, pos[1][1] = j + 1;
		}
	} else *max = score, *max_i = i, *max_j = j;
}

int mm_test_zdrop(void *km, const mm_mapopt_t *opt, const uint8_t *qseq, const uint8_t *tseq, uint32_t n_cigar, uint32_t *cigar, const int8_t *mat)
{
	uint32_t k;
	int32_t score = 0, max = INT32_MIN, max_i = -1, max_j = -1, i = 0, j = 0, max_zdrop = 0;
	int pos[2][2] = {{-1, -1}, {-1, -1}}, q_len, t_len;

	// find the score and the region where score drops most along diagonal
	for (k = 0, score = 0; k < n_cigar; ++k) {
		uint32_t l, op = cigar[k]&0xf, len = cigar[k]>>4;
		if (op == 0) {
			for (l = 0; l < len; ++l) {
				score += mat[tseq[i + l] * 5 + qseq[j + l]];
				update_max_zdrop(score, i+l, j+l, &max, &max_i, &max_j, opt->e, &max_zdrop, pos);
			}
			i += len, j += len;
		} else if (op == 1 || op == 2 || op == 3) {
			score -= opt->q + opt->e * len;
			if (op == 1) j += len; // insertion
			else i += len;         // deletion
			update_max_zdrop(score, i, j, &max, &max_i, &max_j, opt->e, &max_zdrop, pos);
		}
	}

	// test if there is an inversion in the most dropped region
	q_len = pos[1][1] - pos[1][0], t_len = pos[0][1] - pos[0][0];
	if (!(opt->flag&(MM_F_SPLICE|MM_F_SR|MM_F_FOR_ONLY|MM_F_REV_ONLY)) && max_zdrop > opt->zdrop_inv && q_len < opt->max_gap && t_len < opt->max_gap) {
		uint8_t *qseq2;
		void *qp;
		int q_off, t_off;
		qseq2 = (uint8_t*)kmalloc(km, q_len);
		for (i = 0; i < q_len; ++i) {
			int c = qseq[pos[1][1] - i - 1];
			qseq2[i] = c >= 4? 4 : 3 - c;
		}
		qp = ksw_ll_qinit(km, 2, q_len, qseq2, 5, mat);
		score = ksw_ll_i16(qp, t_len, tseq + pos[0][0], opt->q, opt->e, &q_off, &t_off);
		kfree(km, qseq2);
		kfree(km, qp);
		if (score >= opt->min_chain_score * opt->a && score >= opt->min_dp_max)
			return 2; // there is a potential inversion
	}
	return max_zdrop > opt->zdrop? 1 : 0;
}

static void mm_fix_cigar(mm_reg1_t *r, const uint8_t *qseq, const uint8_t *tseq, int *qshift, int *tshift)
{
	mm_extra_t *p = r->p;
	int32_t k, toff = 0, qoff = 0, to_shrink = 0;
	*qshift = *tshift = 0;
	if (p->n_cigar <= 1) return;
	for (k = 0; k < p->n_cigar; ++k) { // indel left alignment
		uint32_t op = p->cigar[k]&0xf, len = p->cigar[k]>>4;
		if (len == 0) to_shrink = 1;
		if (op == 0) {
			toff += len, qoff += len;
		} else if (op == 1 || op == 2) { // insertion or deletion
			if (k > 0 && k < p->n_cigar - 1 && (p->cigar[k-1]&0xf) == 0 && (p->cigar[k+1]&0xf) == 0) {
				int l, prev_len = p->cigar[k-1] >> 4;
				if (op == 1) {
					for (l = 0; l < prev_len; ++l)
						if (qseq[qoff - 1 - l] != qseq[qoff + len - 1 - l])
							break;
				} else {
					for (l = 0; l < prev_len; ++l)
						if (tseq[toff - 1 - l] != tseq[toff + len - 1 - l])
							break;
				}
				if (l > 0)
					p->cigar[k-1] -= l<<4, p->cigar[k+1] += l<<4, qoff -= l, toff -= l;
				if (l == prev_len) to_shrink = 1;
			}
			if (op == 1) qoff += len;
			else toff += len;
		} else if (op == 3) {
			toff += len;
		}
	}
	assert(qoff == r->qe - r->qs && toff == r->re - r->rs);
	if (to_shrink) { // squeeze out zero-length operations
		int32_t l = 0;
		for (k = 0; k < p->n_cigar; ++k) // squeeze out zero-length operations
			if (p->cigar[k]>>4 != 0)
				p->cigar[l++] = p->cigar[k];
		p->n_cigar = l;
		for (k = l = 0; k < p->n_cigar; ++k) // merge two adjacent operations if they are the same
			if (k == p->n_cigar - 1 || (p->cigar[k]&0xf) != (p->cigar[k+1]&0xf))
				p->cigar[l++] = p->cigar[k];
			else p->cigar[k+1] += p->cigar[k]>>4<<4; // add length to the next CIGAR operator
		p->n_cigar = l;
	}
	if ((p->cigar[0]&0xf) == 1 || (p->cigar[0]&0xf) == 2) { // get rid of leading I or D
		int32_t l = p->cigar[0] >> 4;
		if ((p->cigar[0]&0xf) == 1) {
			if (r->rev) r->qe -= l;
			else r->qs += l;
			*qshift = l;
		} else r->rs += l, *tshift = l;
		--p->n_cigar;
		memmove(p->cigar, p->cigar + 1, p->n_cigar * 4);
	}
}

void mm_update_extra(mm_reg1_t *r, const uint8_t *qseq, const uint8_t *tseq, const int8_t *mat, int8_t q, int8_t e)
{
	uint32_t k, l, toff = 0, qoff = 0;
	int32_t s = 0, max = 0, qshift, tshift;
	mm_extra_t *p = r->p;
	if (p == 0) return;
	mm_fix_cigar(r, qseq, tseq, &qshift, &tshift);
	qseq += qshift, tseq += tshift; // qseq and tseq may be shifted due to the removal of leading I/D
	r->blen = r->mlen = 0;
	for (k = 0; k < p->n_cigar; ++k) {
		uint32_t op = p->cigar[k]&0xf, len = p->cigar[k]>>4;
		if (op == 0) { // match/mismatch
			int n_ambi = 0, n_diff = 0;
			for (l = 0; l < len; ++l) {
				int cq = qseq[qoff + l], ct = tseq[toff + l];
				if (ct > 3 || cq > 3) ++n_ambi;
				else if (ct != cq) ++n_diff;
				s += mat[ct * 5 + cq];
				if (s < 0) s = 0;
				else max = max > s? max : s;
			}
			r->blen += len - n_ambi, r->mlen += len - (n_ambi + n_diff), p->n_ambi += n_ambi;
			toff += len, qoff += len;
		} else if (op == 1) { // insertion
			int n_ambi = 0;
			for (l = 0; l < len; ++l)
				if (qseq[qoff + l] > 3) ++n_ambi;
			r->blen += len - n_ambi, p->n_ambi += n_ambi;
			s -= q + e * len;
			if (s < 0) s = 0;
			qoff += len;
		} else if (op == 2) { // deletion
			int n_ambi = 0;
			for (l = 0; l < len; ++l)
				if (tseq[toff + l] > 3) ++n_ambi;
			r->blen += len - n_ambi, p->n_ambi += n_ambi;
			s -= q + e * len;
			if (s < 0) s = 0;
			toff += len;
		} else if (op == 3) { // intron
			toff += len;
		}
	}
	p->dp_max = max;
	assert(qoff == r->qe - r->qs && toff == r->re - r->rs);
}

void mm_append_cigar(mm_reg1_t *r, uint32_t n_cigar, uint32_t *cigar) // TODO: this calls the libc realloc()
{
	mm_extra_t *p;
	if (n_cigar == 0) return;
	if (r->p == 0) {
		uint32_t capacity = n_cigar + sizeof(mm_extra_t);
		kroundup32(capacity);
		r->p = (mm_extra_t*)calloc(capacity, 4);
		r->p->capacity = capacity;
	} else if (r->p->n_cigar + n_cigar + sizeof(mm_extra_t) > r->p->capacity) {
		r->p->capacity = r->p->n_cigar + n_cigar + sizeof(mm_extra_t);
		kroundup32(r->p->capacity);
		r->p = (mm_extra_t*)realloc(r->p, r->p->capacity * 4);
	}
	p = r->p;
	if (p->n_cigar > 0 && (p->cigar[p->n_cigar-1]&0xf) == (cigar[0]&0xf)) { // same CIGAR op at the boundary
		p->cigar[p->n_cigar-1] += cigar[0]>>4<<4;
		if (n_cigar > 1) memcpy(p->cigar + p->n_cigar, cigar + 1, (n_cigar - 1) * 4);
		p->n_cigar += n_cigar - 1;
	} else {
		memcpy(p->cigar + p->n_cigar, cigar, n_cigar * 4);
		p->n_cigar += n_cigar;
	}
}

void mm_align_pair(void *km, const mm_mapopt_t *opt, int qlen, const uint8_t *qseq, int tlen, const uint8_t *tseq, const int8_t *mat, int w, int end_bonus, int zdrop, int flag, ksw_extz_t *ez)
{
	if (mm_dbg_flag & MM_DBG_PRINT_ALN_SEQ) {
		int i;
		fprintf(stderr, "===> q=(%d,%d), e=(%d,%d), bw=%d, flag=%d, zdrop=%d <===\n", opt->q, opt->q2, opt->e, opt->e2, w, flag, opt->zdrop);
		for (i = 0; i < tlen; ++i) fputc("ACGTN"[tseq[i]], stderr);
		fputc('\n', stderr);
		for (i = 0; i < qlen; ++i) fputc("ACGTN"[qseq[i]], stderr);
		fputc('\n', stderr);
	}
	if (opt->flag & MM_F_SPLICE)
		ksw_exts2_sse(km, qlen, qseq, tlen, tseq, 5, mat, opt->q, opt->e, opt->q2, opt->noncan, zdrop, flag, ez);
	else if (opt->q == opt->q2 && opt->e == opt->e2)
		ksw_extz2_sse(km, qlen, qseq, tlen, tseq, 5, mat, opt->q, opt->e, w, zdrop, end_bonus, flag, ez);
	else
		ksw_extd2_sse(km, qlen, qseq, tlen, tseq, 5, mat, opt->q, opt->e, opt->q2, opt->e2, w, zdrop, end_bonus, flag, ez);
	if (mm_dbg_flag & MM_DBG_PRINT_ALN_SEQ) {
		int i;
		fprintf(stderr, "score=%d, cigar=", ez->score);
		for (i = 0; i < ez->n_cigar; ++i)
			fprintf(stderr, "%d%c", ez->cigar[i]>>4, "MIDN"[ez->cigar[i]&0xf]);
		fprintf(stderr, "\n");
	}
}

static inline int mm_get_hplen_back(const mm_idx_t *mi, uint32_t rid, uint32_t x)
{
	int64_t i, off0 = mi->seq[rid].offset, off = off0 + x;
	int c = mm_seq4_get(mi->S, off);
	for (i = off - 1; i >= off0; --i)
		if (mm_seq4_get(mi->S, i) != c) break;
	return (int)(off - i);
}

static inline void mm_adjust_minier(const mm_idx_t *mi, uint8_t *const qseq0[2], mm128_t *a, int32_t *r, int32_t *q)
{
	if (mi->flag & MM_I_HPC) {
		const uint8_t *qseq = qseq0[a->x>>63];
		int i, c;
		*q = (int32_t)a->y;
		for (i = *q - 1, c = qseq[*q]; i > 0; --i)
			if (qseq[i] != c) break;
		*q = i + 1;
		c = mm_get_hplen_back(mi, a->x<<1>>33, (int32_t)a->x);
		*r = (int32_t)a->x + 1 - c;
	} else {
		*r = (int32_t)a->x - (mi->k>>1);
		*q = (int32_t)a->y - (mi->k>>1);
	}
}

static void mm_filter_bad_seeds(void *km, int as1, int cnt1, mm128_t *a, int min_gap, int diff_thres, int max_ext_len, int max_ext_cnt)
{
	int max_st, max_en, n, i, k, max, *K;
	for (i = 1, n = 0; i < cnt1; ++i) { // count the number of gaps longer than min_gap
		int gap = ((int32_t)a[as1 + i].y - a[as1 + i - 1].y) - ((int32_t)a[as1 + i].x - a[as1 + i - 1].x);
		if (gap < -min_gap || gap > min_gap) ++n;
	}
	if (n <= 1) return;
	K = (int*)kmalloc(km, n * sizeof(int));
	for (i = 1, n = 0; i < cnt1; ++i) { // store the positions of long gaps
		int gap = ((int32_t)a[as1 + i].y - a[as1 + i - 1].y) - ((int32_t)a[as1 + i].x - a[as1 + i - 1].x);
		if (gap < -min_gap || gap > min_gap)
			K[n++] = i;
	}
	max = 0, max_st = max_en = -1;
	for (k = 0;; ++k) { // traverse long gaps
		int gap, l, n_ins = 0, n_del = 0, qs, rs, max_diff = 0, max_diff_l = -1;
		if (k == n || k >= max_en) {
			if (max_en > 0)
				for (i = K[max_st]; i < K[max_en]; ++i)
					a[as1 + i].y |= MM_SEED_IGNORE;
			max = 0, max_st = max_en = -1;
			if (k == n) break;
		}
		i = K[k];
		gap = ((int32_t)a[as1 + i].y - a[as1 + i - 1].y) - ((int32_t)a[as1 + i].x - a[as1 + i - 1].x);
		if (gap > 0) n_ins += gap;
		else n_del += -gap;
		qs = (int32_t)a[as1 + i - 1].y;
		rs = (int32_t)a[as1 + i - 1].x;
		for (l = k + 1; l < n && l <= k + max_ext_cnt; ++l) {
			int j = K[l], diff;
			if ((int32_t)a[as1 + j].y - qs > max_ext_len || (int32_t)a[as1 + j].x - rs > max_ext_len) break;
			gap = ((int32_t)a[as1 + j].y - (int32_t)a[as1 + j - 1].y) - (a[as1 + j].x - a[as1 + j - 1].x);
			if (gap > 0) n_ins += gap;
			else n_del += -gap;
			diff = n_ins + n_del - abs(n_ins - n_del);
			if (max_diff < diff)
				max_diff = diff, max_diff_l = l;
		}
		if (max_diff > diff_thres && max_diff > max)
			max = max_diff, max_st = k, max_en = max_diff_l;
	}
	kfree(km, K);
}

static void mm_fix_bad_ends(const mm_reg1_t *r, const mm128_t *a, int bw, int min_match, int32_t *as, int32_t *cnt)
{
	int32_t i, l, m;
	*as = r->as, *cnt = r->cnt;
	if (r->cnt < 3) return;
	m = l = a[r->as].y >> 32 & 0xff;
	for (i = r->as + 1; i < r->as + r->cnt - 1; ++i) {
		int32_t lq, lr, min, max;
		int32_t q_span = a[i].y >> 32 & 0xff;
		if (a[i].y & MM_SEED_LONG_JOIN) break;
		lr = (int32_t)a[i].x - (int32_t)a[i-1].x;
		lq = (int32_t)a[i].y - (int32_t)a[i-1].y;
		min = lr < lq? lr : lq;
		max = lr > lq? lr : lq;
		if (max - min > l >> 1) *as = i;
		l += min;
		m += min < q_span? min : q_span;
		if (l >= bw << 1 || (m >= min_match && m >= bw) || m >= r->mlen >> 1) break;
	}
	*cnt = r->as + r->cnt - *as;
	m = l = a[r->as + r->cnt - 1].y >> 32 & 0xff;
	for (i = r->as + r->cnt - 2; i > *as; --i) {
		int32_t lq, lr, min, max;
		int32_t q_span = a[i+1].y >> 32 & 0xff;
		if (a[i+1].y & MM_SEED_LONG_JOIN) break;
		lr = (int32_t)a[i+1].x - (int32_t)a[i].x;
		lq = (int32_t)a[i+1].y - (int32_t)a[i].y;
		min = lr < lq? lr : lq;
		max = lr > lq? lr : lq;
		if (max - min > l >> 1) *cnt = i + 1 - *as;
		l += min;
		m += min < q_span? min : q_span;
		if (l >= bw << 1 || (m >= min_match && m >= bw) || m >= r->mlen >> 1) break;
	}
}

static void mm_max_stretch(const mm_mapopt_t *opt, const mm_reg1_t *r, const mm128_t *a, int32_t *as, int32_t *cnt)
{
	int32_t i, score, max_score, len, max_i, max_len;

	*as = r->as, *cnt = r->cnt;
	if (r->cnt < 2) return;

	max_score = -1, max_i = -1, max_len = 0;
	score = a[r->as].y >> 32 & 0xff, len = 1;
	for (i = r->as + 1; i < r->as + r->cnt; ++i) {
		int32_t lq, lr, q_span;
		q_span = a[i].y >> 32 & 0xff;
		lr = (int32_t)a[i].x - (int32_t)a[i-1].x;
		lq = (int32_t)a[i].y - (int32_t)a[i-1].y;
		if (lq == lr) {
			score += lq < q_span? lq : q_span;
			++len;
		} else {
			if (score > max_score)
				max_score = score, max_len = len, max_i = i - len;
			score = q_span, len = 1;
		}
	}
	if (score > max_score)
		max_score = score, max_len = len, max_i = i - len;
	*as = max_i, *cnt = max_len;
}

static int mm_seed_ext_score(void *km, const mm_mapopt_t *opt, const mm_idx_t *mi, const int8_t mat[25], int qlen, uint8_t *qseq0[2], const mm128_t *a)
{
	uint8_t *qseq, *tseq;
	int q_span = a->y>>32&0xff, qs, qe, rs, re, rid, score, q_off, t_off, ext_len = opt->anchor_ext_len;
	void *qp;
	rid = a->x<<1>>33;
	re = (uint32_t)a->x + 1, rs = re - q_span;
	qe = (uint32_t)a->y + 1, qs = qe - q_span;
	rs = rs - ext_len > 0? rs - ext_len : 0;
	qs = qs - ext_len > 0? qs - ext_len : 0;
	re = re + ext_len < mi->seq[rid].len? re + ext_len : mi->seq[rid].len;
	qe = qe + ext_len < qlen? qe + ext_len : qlen;
	tseq = (uint8_t*)kmalloc(km, re - rs);
	mm_idx_getseq(mi, rid, rs, re, tseq);
	qseq = qseq0[a->x>>63] + qs;
	qp = ksw_ll_qinit(km, 2, qe - qs, qseq, 5, mat);
	score = ksw_ll_i16(qp, re - rs, tseq, opt->q, opt->e, &q_off, &t_off);
	kfree(km, tseq);
	kfree(km, qp);
	return score;
}

static void mm_fix_bad_ends_splice(void *km, const mm_mapopt_t *opt, const mm_idx_t *mi, const mm_reg1_t *r, const int8_t mat[25], int qlen, uint8_t *qseq0[2], const mm128_t *a, int *as1, int *cnt1)
{ // this assumes a very crude k-mer based mode; it is not necessary to use a good model just for filtering bounary exons
	int score;
	double log_gap;
	*as1 = r->as, *cnt1 = r->cnt;
	if (r->cnt < 3) return;
	log_gap = log((int32_t)a[r->as + 1].x - (int32_t)a[r->as].x);
	if ((a[r->as].y>>32&0xff) < log_gap + opt->anchor_ext_shift) {
		score = mm_seed_ext_score(km, opt, mi, mat, qlen, qseq0, &a[r->as]);
		if ((double)score / mat[0] < log_gap + opt->anchor_ext_shift) // a more exact format is "score < log_4(gap) + shift"
			++(*as1), --(*cnt1);
	}
	log_gap = log((int32_t)a[r->as + r->cnt - 1].x - (int32_t)a[r->as + r->cnt - 2].x);
	if ((a[r->as + r->cnt - 1].y>>32&0xff) < log_gap + opt->anchor_ext_shift) {
		score = mm_seed_ext_score(km, opt, mi, mat, qlen, qseq0, &a[r->as + r->cnt - 1]);
		if ((double)score / mat[0] < log_gap + opt->anchor_ext_shift)
			--(*cnt1);
	}
}
#include "fpga_sw.h"
static unsigned long long g_tag = 0;

extern double get_mem_time[100];
extern double package_time[100];
extern double send_time[100];

void send_to_fpga(chain_sw_task_t* tasks[], int chain_num, int data_size)
{
    char* fpga_buf = NULL;
    char* p = NULL;
    fpga_task_t* head = NULL;
    fpga_task_id_t* chain_head = NULL;
    fpga_sw_task* sw_task = NULL;
    task_t fpga_task;

/*retry:
    fpga_buf = (char*)fpga_get_writebuf(data_size + 4096, BUF_TYPE_SW);
    if(fpga_buf == NULL) {
        goto retry;
        fprintf(stderr, "fpga_get_writebuf sw error:%s\n", strerror(errno));
        exit(0);
    }*/
    fpga_buf = (char*)malloc(data_size + 4096);

    head = (fpga_task_t*)fpga_buf;        //指向buff首地址
    head->check_id = CHECK_ID;
    head->sw_num = 0;
    chain_head = (fpga_task_id_t*)(fpga_buf + sizeof(fpga_task_t));
    sw_task = (fpga_sw_task*)(fpga_buf + 4096);
    //fprintf(stderr, "fpga_buf=%p\n", fpga_buf);
    //设置头信息
    head->tag = g_tag++;
    head->chain_task_num = chain_num;

    int chain_index = 0;
    int sw_index = 0;
    for(chain_index = 0; chain_index < chain_num; chain_index++) {
        assert((char*)chain_head < (char*)(fpga_buf + 4096));

        //设置chain的信息
        //fprintf(stderr, "chain_head=%p\n", chain_head);
        chain_head->read_id = tasks[chain_index]->read_id;
        chain_head->chain_id = tasks[chain_index]->chain_id;
        chain_head->sw_num = tasks[chain_index]->sw_num;

        //设置sw任务的数据
        head->sw_num += chain_head->sw_num;     //统计一次调用sw任务数
        for(sw_index = 0; sw_index < chain_head->sw_num; sw_index++) {
            sw_task->qlen = tasks[chain_index]->sw_tasks[sw_index]->qlen;
            sw_task->tlen = tasks[chain_index]->sw_tasks[sw_index]->tlen;
            sw_task->flag = tasks[chain_index]->sw_tasks[sw_index]->flag;
            sw_task->zdrop = tasks[chain_index]->sw_tasks[sw_index]->zdrop;
            sw_task->bw = tasks[chain_index]->sw_tasks[sw_index]->w;
            sw_task->end_bonus = tasks[chain_index]->sw_tasks[sw_index]->end_bonus;
            
            uint8_t *qseq = (uint8_t*)sw_task + sizeof(fpga_sw_task);    //设置qseq的地址
            memcpy(qseq, tasks[chain_index]->sw_tasks[sw_index]->query, sw_task->qlen);
            uint8_t *tseq = qseq + ADDR_ALIGN(sw_task->qlen, 16);
            memcpy(tseq, tasks[chain_index]->sw_tasks[sw_index]->target, sw_task->tlen);
            
            p = (char*)sw_task;
            sw_task = (fpga_sw_task*)(p + sizeof(fpga_sw_task) + ADDR_ALIGN(sw_task->qlen, 16) + ADDR_ALIGN(sw_task->tlen, 16));    //指向下一个sw任务的地址
        }
        chain_head += 1;   //指向下一个chain的头
        destroy_chain_sw_task(tasks[chain_index]);     //销毁chain任务的数据
    }

    assert((data_size + 4096) == ((char*)sw_task - fpga_buf));

    
    //fpga_writebuf_submit(fpga_buf, data_size + 4096, TYPE_SW);
    fpga_task.data = (void*)fpga_buf;
    fpga_task.size = data_size + 4096;
    while(send_fpga_task(fpga_task)) {
        usleep(5000);
        ;
    }
}

void last_send(void *data, int tid)
{
    user_params_t* params = (user_params_t*)data;
    char* fpga_buf = NULL;
    char* p = NULL;
    fpga_task_t* head = NULL;
    fpga_task_id_t* chain_head = NULL;
    fpga_sw_task* sw_task = NULL;
    chain_sw_task_t** tasks = params->send_task[tid].chain_tasks;
    int chain_num = params->send_task[tid].num;
    int data_size = params->send_task[tid].data_size;
    task_t fpga_task;

    //fprintf(stderr, "last send, num=%d, tid=%d\n", chain_num, tid);
    
    if(chain_num == 0) {
        return;
    }
    
/*retry_last:
    fpga_buf = (char*)fpga_get_writebuf(data_size + 4096, BUF_TYPE_SW);
    if(fpga_buf == NULL) {
        goto retry_last;
        fprintf(stderr, "fpga_get_writebuf sw error:%s\n", strerror(errno));
        exit(0);
    }*/
    fpga_buf = (char*)malloc(data_size + 4096);

    head = (fpga_task_t*)fpga_buf;        //指向buff首地址
    head->check_id = CHECK_ID;
    head->sw_num = 0;
    chain_head = (fpga_task_id_t*)(fpga_buf + sizeof(fpga_task_t));
    sw_task = (fpga_sw_task*)(fpga_buf + 4096);
    //fprintf(stderr, "fpga_buf=%p\n", fpga_buf);
    //设置头信息
    head->tag = g_tag++;
    head->chain_task_num = chain_num;

    int chain_index = 0;
    int sw_index = 0;
    for(chain_index = 0; chain_index < chain_num; chain_index++) {
        assert((char*)chain_head < (char*)(fpga_buf + 4096));

        //设置chain的信息
        //fprintf(stderr, "chain_head=%p\n", chain_head);
        chain_head->read_id = tasks[chain_index]->read_id;
        chain_head->chain_id = tasks[chain_index]->chain_id;
        chain_head->sw_num = tasks[chain_index]->sw_num;

        //设置sw任务的数据
        head->sw_num += chain_head->sw_num;     //统计一次调用sw任务数
        for(sw_index = 0; sw_index < chain_head->sw_num; sw_index++) {
            sw_task->qlen = tasks[chain_index]->sw_tasks[sw_index]->qlen;
            sw_task->tlen = tasks[chain_index]->sw_tasks[sw_index]->tlen;
            sw_task->flag = tasks[chain_index]->sw_tasks[sw_index]->flag;
            sw_task->zdrop = tasks[chain_index]->sw_tasks[sw_index]->zdrop;
            sw_task->bw = tasks[chain_index]->sw_tasks[sw_index]->w;
            sw_task->end_bonus = tasks[chain_index]->sw_tasks[sw_index]->end_bonus;
            
            uint8_t *qseq = (uint8_t*)sw_task + sizeof(fpga_sw_task);    //设置qseq的地址
            memcpy(qseq, tasks[chain_index]->sw_tasks[sw_index]->query, sw_task->qlen);
            uint8_t *tseq = qseq + ADDR_ALIGN(sw_task->qlen, 16);
            memcpy(tseq, tasks[chain_index]->sw_tasks[sw_index]->target, sw_task->tlen);
            
            p = (char*)sw_task;
            sw_task = (fpga_sw_task*)(p + sizeof(fpga_sw_task) + ADDR_ALIGN(sw_task->qlen, 16) + ADDR_ALIGN(sw_task->tlen, 16));    //指向下一个sw任务的地址
        }
        chain_head += 1;   //指向下一个chain的头
        destroy_chain_sw_task(tasks[chain_index]);     //销毁chain任务的数据
    }

    assert((data_size + 4096) == ((char*)sw_task - fpga_buf));

    //fpga_writebuf_submit(fpga_buf, data_size + 4096, TYPE_SW);
    fpga_task.data = (void*)fpga_buf;
    fpga_task.size = data_size + 4096;
    while(send_fpga_task(fpga_task)) {
        usleep(5000);
        ;
    }
    
    params->send_task[tid].num = 0;
}

#define VSCMIN(x,y) ({ \
  typeof(x) _x = (x);     \
  typeof(y) _y = (y);     \
  (void) (&_x == &_y);            \
  _x < _y ? _x : _y; })

static pthread_mutex_t long_chain_mutex = PTHREAD_MUTEX_INITIALIZER;
long long_chain_counter = 0;
extern volatile int g_total_task_num;
extern volatile int g_process_task_num;

void process_task(send_task_t* send_task, chain_sw_task_t* task, int tid, long read_index, const long read_num)
{
    /*if(send_task->num == SEND_ARRAY_MAX) {
        send_to_fpga(send_task->chain_tasks, send_task->num, send_task->data_size, tid);
        __sync_add_and_fetch(&g_process_task_num, 1);
        send_task->num = 0;
        send_task->data_size = 0;
    }*/
    
    send_task->chain_tasks[send_task->num++] = task;
    send_task->data_size += task->data_size;
    
    if(send_task->num > 64 && send_task->data_size < 6144000) {
        send_to_fpga(send_task->chain_tasks, send_task->num, send_task->data_size);
        __sync_add_and_fetch(&g_process_task_num, 1);
        send_task->num = 0;
        send_task->data_size = 0;
    }
    else if(send_task->num > 200 && send_task->data_size >= 6144000) {
        fprintf(stderr, "too much task, num:%d, size:%d\n", send_task->num, send_task->data_size);
    }

    /*if(g_process_task_num < (g_total_task_num/2)) {     //如果任务是总任务的一半以下，则立即发送任务
        if(send_task->num > 32) {
            send_to_fpga(send_task->chain_tasks, send_task->num, send_task->data_size);
            __sync_add_and_fetch(&g_process_task_num, 1);
            send_task->num = 0;
            send_task->data_size = 0;
        }
    }
    else if(g_process_task_num >= (g_total_task_num/2)) {
        if(send_task->num < 32) {
            send_to_fpga(send_task->chain_tasks, send_task->num, send_task->data_size);
            __sync_add_and_fetch(&g_process_task_num, 1);
            send_task->num = 0;
            send_task->data_size = 0;
        }
    }*/
    return;
}

static void mm_align1(void *km, const mm_mapopt_t *opt, const mm_idx_t *mi, int qlen, uint8_t *qseq0[2], mm_reg1_t *r, mm_reg1_t *r2, int n_a, mm128_t *a, ksw_extz_t *ez, int splice_flag, long read_index, int chain_index, context_t* context, chain_context_t* chain_context, send_task_t* send_task, int tid, const long read_num)
{
    double start, end;
	int is_sr = !!(opt->flag & MM_F_SR), is_splice = !!(opt->flag & MM_F_SPLICE);
	int32_t rid = a[r->as].x<<1>>33, rev = a[r->as].x>>63, as1, cnt1;
	uint8_t *tseq, *qseq;
	int32_t i, l, bw, dropped = 0, extra_flag = 0, rs0, re0, qs0, qe0;
	int32_t rs, re, qs, qe;
	int32_t rs1, qs1, re1, qe1;
	int8_t mat[25];
    int sw_index = 0;

    start = realtime_msec();
    
	if (is_sr) assert(!(mi->flag & MM_I_HPC)); // HPC won't work with SR because with HPC we can't easily tell if there is a gap

	r2->cnt = 0;
	if (r->cnt == 0) return;
	ksw_gen_simple_mat(5, mat, opt->a, opt->b);
	bw = (int)(opt->bw * 1.5 + 1.);

	if (is_sr && !(mi->flag & MM_I_HPC)) {
		mm_max_stretch(opt, r, a, &as1, &cnt1);
		rs = (int32_t)a[as1].x + 1 - (int32_t)(a[as1].y>>32&0xff);
		qs = (int32_t)a[as1].y + 1 - (int32_t)(a[as1].y>>32&0xff);
		re = (int32_t)a[as1+cnt1-1].x + 1;
		qe = (int32_t)a[as1+cnt1-1].y + 1;
	} else {
		if (is_splice) {
			mm_fix_bad_ends_splice(km, opt, mi, r, mat, qlen, qseq0, a, &as1, &cnt1);
		} else {
			mm_fix_bad_ends(r, a, opt->bw, opt->min_chain_score * 2, &as1, &cnt1);
		}
		mm_filter_bad_seeds(km, as1, cnt1, a, 10, 40, opt->max_gap>>1, 10);
		mm_adjust_minier(mi, qseq0, &a[as1], &rs, &qs);
		mm_adjust_minier(mi, qseq0, &a[as1 + cnt1 - 1], &re, &qe);
	}
	assert(cnt1 > 0);

	if (is_splice) {
		if (splice_flag & MM_F_SPLICE_FOR) extra_flag |= rev? KSW_EZ_SPLICE_REV : KSW_EZ_SPLICE_FOR;
		if (splice_flag & MM_F_SPLICE_REV) extra_flag |= rev? KSW_EZ_SPLICE_FOR : KSW_EZ_SPLICE_REV;
		if (opt->flag & MM_F_SPLICE_FLANK) extra_flag |= KSW_EZ_SPLICE_FLANK;
	}

	/* Look for the start and end of regions to perform DP. This sounds easy
	 * but is in fact tricky. Excessively small regions lead to unnecessary
	 * clippings and lose alignable sequences. Excessively large regions
	 * occasionally lead to large overlaps between two chains and may cause
	 * loss of alignments in corner cases. */
	if (is_sr) {
		qs0 = 0, qe0 = qlen;
		l = qs;
		l += l * opt->a + opt->end_bonus > opt->q? (l * opt->a + opt->end_bonus - opt->q) / opt->e : 0;
		rs0 = rs - l > 0? rs - l : 0;
		l = qlen - qe;
		l += l * opt->a + opt->end_bonus > opt->q? (l * opt->a + opt->end_bonus - opt->q) / opt->e : 0;
		re0 = re + l < mi->seq[rid].len? re + l : mi->seq[rid].len;
	} else {
		// compute rs0 and qs0
		rs0 = (int32_t)a[r->as].x + 1 - (int32_t)(a[r->as].y>>32&0xff);
		qs0 = (int32_t)a[r->as].y + 1 - (int32_t)(a[r->as].y>>32&0xff);
		if (rs0 < 0) rs0 = 0; // this may happen when HPC is in use
		assert(qs0 >= 0); // this should never happen, or it is logic error
		rs1 = qs1 = 0;
		for (i = r->as - 1, l = 0; i >= 0 && a[i].x>>32 == a[r->as].x>>32; --i) { // inspect nearby seeds
			int32_t x = (int32_t)a[i].x + 1 - (int32_t)(a[i].y>>32&0xff);
			int32_t y = (int32_t)a[i].y + 1 - (int32_t)(a[i].y>>32&0xff);
			if (x < rs0 && y < qs0) {
				if (++l > opt->min_cnt) {
					l = rs0 - x > qs0 - y? rs0 - x : qs0 - y;
					rs1 = rs0 - l, qs1 = qs0 - l;
					break;
				}
			}
		}
		if (qs > 0 && rs > 0) {
			l = qs < opt->max_gap? qs : opt->max_gap;
			qs1 = qs1 > qs - l? qs1 : qs - l;
			qs0 = qs0 < qs1? qs0 : qs1; // at least include qs0
			l += l * opt->a > opt->q? (l * opt->a - opt->q) / opt->e : 0;
			l = l < opt->max_gap? l : opt->max_gap;
			l = l < rs? l : rs;
			rs1 = rs1 > rs - l? rs1 : rs - l;
			rs0 = rs0 < rs1? rs0 : rs1;
		} else rs0 = rs, qs0 = qs;
		// compute re0 and qe0
		re0 = (int32_t)a[r->as + r->cnt - 1].x + 1;
		qe0 = (int32_t)a[r->as + r->cnt - 1].y + 1;
		re1 = mi->seq[rid].len, qe1 = qlen;
		for (i = r->as + r->cnt, l = 0; i < n_a && a[i].x>>32 == a[r->as].x>>32; ++i) { // inspect nearby seeds
			int32_t x = (int32_t)a[i].x + 1;
			int32_t y = (int32_t)a[i].y + 1;
			if (x > re0 && y > qe0) {
				if (++l > opt->min_cnt) {
					l = x - re0 > y - qe0? x - re0 : y - qe0;
					re1 = re0 + l, qe1 = qe0 + l;
					break;
				}
			}
		}
		if (qe < qlen && re < mi->seq[rid].len) {
			l = qlen - qe < opt->max_gap? qlen - qe : opt->max_gap;
			qe1 = qe1 < qe + l? qe1 : qe + l;
			qe0 = qe0 > qe1? qe0 : qe1; // at least include qe0
			l += l * opt->a > opt->q? (l * opt->a - opt->q) / opt->e : 0;
			l = l < opt->max_gap? l : opt->max_gap;
			l = l < mi->seq[rid].len - re? l : mi->seq[rid].len - re;
			re1 = re1 < re + l? re1 : re + l;
			re0 = re0 > re1? re0 : re1;
		} else re0 = re, qe0 = qe;
	}
	if (a[r->as].y & MM_SEED_SELF) {
		int max_ext = r->qs > r->rs? r->qs - r->rs : r->rs - r->qs;
		if (r->rs - rs0 > max_ext) rs0 = r->rs - max_ext;
		if (r->qs - qs0 > max_ext) qs0 = r->qs - max_ext;
		max_ext = r->qe > r->re? r->qe - r->re : r->re - r->qe;
		if (re0 - r->re > max_ext) re0 = r->re + max_ext;
		if (qe0 - r->qe > max_ext) qe0 = r->qe + max_ext;
	}

	assert(re0 > rs0);
    int tseq_len = re0 - rs0;
	tseq = (uint8_t*)kmalloc(km, re0 - rs0);

    //创建一个chain的sw任务
    int small;
    chain_sw_task_t* chain_task = create_chain_sw_task(read_index, chain_index);
    chain_sw_task_t* soft_sw_task = create_chain_sw_task(read_index, chain_index);
    //context->tseq = tseq;

    memcpy(context->mat, mat, sizeof(context->mat));

    chain_context->rs = rs;
    chain_context->qs = qs;
    chain_context->qs0 = qs0;
    chain_context->rid = rid;
    chain_context->rev = rev;

    sw_task_t* sw_task = NULL;
    sw_context_t* sw_context = NULL;
	if (qs > 0 && rs > 0) { // left extension
		qseq = &qseq0[rev][qs0];
		mm_idx_getseq(mi, rid, rs0, rs, tseq);
		mm_seq_rev(qs - qs0, qseq);
		mm_seq_rev(rs - rs0, tseq);

        sw_context = (sw_context_t*)malloc(sizeof(sw_context_t));
        add_sw_context(chain_context, sw_context);
        //TODO 保存上下文
        sw_context->rs = rs;
        sw_context->qs = qs;
        sw_context->qs0 = qs0;
        sw_context->pos_flag = 0;
        
        sw_context->query = (uint8_t*)malloc((qs - qs0) * sizeof(uint8_t));
        sw_context->target = (uint8_t*)malloc((rs - rs0) * sizeof(uint8_t));
        memcpy(sw_context->query, qseq, (qs - qs0) * sizeof(uint8_t));
        memcpy(sw_context->target, tseq, (rs - rs0) * sizeof(uint8_t));
        sw_context->qlen = qs - qs0;
        sw_context->tlen = rs - rs0;
        sw_context->w = bw;
        sw_context->end_bonus = opt->end_bonus;
        sw_context->zdrop = r->split_inv? opt->zdrop_inv : opt->zdrop;
        sw_context->flag = extra_flag|KSW_EZ_EXTZ_ONLY|KSW_EZ_RIGHT|KSW_EZ_REV_CIGAR;
        
        //fprintf(stderr, "1.qlen=%d, tlen=%d, w=%d\n", qs - qs0, rs - rs0, bw);
        sw_task = create_sw_task(qs - qs0, sw_context->query, rs - rs0, sw_context->target,
                        mat, opt->q, opt->e, opt->q2, opt->e2, bw, 
                        r->split_inv? opt->zdrop_inv : opt->zdrop,
                        opt->end_bonus, extra_flag|KSW_EZ_EXTZ_ONLY|KSW_EZ_RIGHT|KSW_EZ_REV_CIGAR, 
                        0);

        small = VSCMIN(sw_task->qlen, sw_task->tlen);
        small = VSCMIN(small, sw_task->w);
        if(sw_task->qlen >= 16300 || sw_task->tlen >= 16300 || small>=1024) {
            chain_context->soft_sw_num++;
            
            //造一个小的激励让fpga做，保持chain task数据和结果的完整性
            sw_task_t* tmp_sw_task = create_sw_task(100, sw_context->query, 100, sw_context->target,
                                                    mat, opt->q, opt->e, opt->q2, opt->e2, 75, 
                                                    r->split_inv? opt->zdrop_inv : opt->zdrop,
                                                    opt->end_bonus, extra_flag|KSW_EZ_EXTZ_ONLY|KSW_EZ_RIGHT|KSW_EZ_REV_CIGAR, 
                                                    0);
            add_sw_task(chain_task, tmp_sw_task);
            
            //设置该sw的index值，用来进行结果替换
            sw_task->sw_index = sw_index;
            add_sw_task(soft_sw_task, sw_task);
        }
        else {
            add_sw_task(chain_task, sw_task);
        }
        sw_index++;
		/*mm_align_pair(km, opt, qs - qs0, qseq, rs - rs0, tseq, mat, bw, opt->end_bonus, r->split_inv? opt->zdrop_inv : opt->zdrop, extra_flag|KSW_EZ_EXTZ_ONLY|KSW_EZ_RIGHT|KSW_EZ_REV_CIGAR, ez);
        if (ez->n_cigar > 0) {
			mm_append_cigar(r, ez->n_cigar, ez->cigar);
			r->p->dp_score += ez->max;
		}
		rs1 = rs - (ez->reach_end? ez->mqe_t + 1 : ez->max_t + 1);
		qs1 = qs - (ez->reach_end? qs - qs0 : ez->max_q + 1);*/
		mm_seq_rev(qs - qs0, qseq);
	}
    else {
        rs1 = rs;
        qs1 = qs;
    }
	re1 = rs, qe1 = qs;
	//assert(qs1 >= 0 && rs1 >= 0);

	for (i = is_sr? cnt1 - 1 : 1; i < cnt1; ++i) { // gap filling
		if ((a[as1+i].y & (MM_SEED_IGNORE|MM_SEED_TANDEM)) && i != cnt1 - 1) continue;
		if (is_sr && !(mi->flag & MM_I_HPC)) {
			re = (int32_t)a[as1 + i].x + 1;
			qe = (int32_t)a[as1 + i].y + 1;
		} else mm_adjust_minier(mi, qseq0, &a[as1 + i], &re, &qe);
		re1 = re, qe1 = qe;
		if (i == cnt1 - 1 || (a[as1+i].y&MM_SEED_LONG_JOIN) || (qe - qs >= opt->min_ksw_len && re - rs >= opt->min_ksw_len)) {
			//int j, bw1 = bw, zdrop_code;
            int j, bw1 = bw;
			if (a[as1+i].y & MM_SEED_LONG_JOIN)
				bw1 = qe - qs > re - rs? qe - qs : re - rs;
			// perform alignment
			qseq = &qseq0[rev][qs];
			mm_idx_getseq(mi, rid, rs, re, tseq);
			if (is_sr) { // perform ungapped alignment
				assert(qe - qs == re - rs);
				ksw_reset_extz(ez);
				for (j = 0, ez->score = 0; j < qe - qs; ++j) {
					if (qseq[j] >= 4 || tseq[j] >= 4) ez->score += opt->e2;
					else ez->score += qseq[j] == tseq[j]? opt->a : -opt->b;
				}
				ez->cigar = ksw_push_cigar(km, &ez->n_cigar, &ez->m_cigar, ez->cigar, 0, qe - qs);
			} else { // perform normal gapped alignment
                //fprintf(stderr, "2.qlen=%d, tlen=%d, w=%d\n", qe - qs, re - rs, bw1);
				
                sw_context = (sw_context_t*)malloc(sizeof(sw_context_t));
                add_sw_context(chain_context, sw_context);
                //TODO 保存上下文
                sw_context->i = i;
                sw_context->pos_flag = 1;
                sw_context->cnt1 = cnt1;
                sw_context->re = re;
                sw_context->qe = qe;

                sw_context->rs = rs;
                sw_context->qs = qs;
                sw_context->qs0 = qs0;
                
                sw_context->query = (uint8_t*)malloc((qe - qs) * sizeof(uint8_t));
                sw_context->target = (uint8_t*)malloc((re - rs) * sizeof(uint8_t));
                memcpy(sw_context->query, qseq, (qe - qs) * sizeof(uint8_t));
                memcpy(sw_context->target, tseq, (re - rs) * sizeof(uint8_t));
                sw_context->qlen = qe - qs;
                sw_context->tlen = re - rs;
                sw_context->w = bw1;
                sw_context->end_bonus = -1;
                sw_context->zdrop = opt->zdrop;
                sw_context->flag = extra_flag|KSW_EZ_APPROX_MAX;
                sw_context->zdrop_flag = extra_flag;
                sw_context->as1 = as1;

                sw_task = create_sw_task(qe - qs, sw_context->query, re - rs, sw_context->target,
                        mat, opt->q, opt->e, opt->q2, opt->e2, bw1, 
                        opt->zdrop,
                        -1, extra_flag|KSW_EZ_APPROX_MAX, 
                        1);
                small = VSCMIN(sw_task->qlen, sw_task->tlen);
                small = VSCMIN(small, sw_task->w);
                if(sw_task->qlen >= 16300 || sw_task->tlen >= 16300 || small>=1024) {
                    chain_context->soft_sw_num++;
                    //造一个小的激励让fpga做，保持chain task数据和结果的完整性
                    sw_task_t* tmp_sw_task = create_sw_task(100, sw_context->query, 100, sw_context->target,
                                                            mat, opt->q, opt->e, opt->q2, opt->e2, 75, 
                                                            r->split_inv? opt->zdrop_inv : opt->zdrop,
                                                            opt->end_bonus, extra_flag|KSW_EZ_EXTZ_ONLY|KSW_EZ_RIGHT|KSW_EZ_REV_CIGAR, 
                                                            0);
                    add_sw_task(chain_task, tmp_sw_task);
                    //设置该sw的index值，用来进行结果替换
                    sw_task->sw_index = sw_index;
                    add_sw_task(soft_sw_task, sw_task);
                }
                else {
                    add_sw_task(chain_task, sw_task);
                }
                sw_index++;
                //mm_align_pair(km, opt, qe - qs, qseq, re - rs, tseq, mat, bw1, -1, opt->zdrop, extra_flag|KSW_EZ_APPROX_MAX, ez); // first pass: with approximate Z-drop
            }
			// test Z-drop and inversion Z-drop
			//if ((zdrop_code = mm_test_zdrop(km, opt, qseq, tseq, ez->n_cigar, ez->cigar, mat)) != 0)
			//	mm_align_pair(km, opt, qe - qs, qseq, re - rs, tseq, mat, bw1, -1, zdrop_code == 2? opt->zdrop_inv : opt->zdrop, extra_flag, ez); // second pass: lift approximate
			// update CIGAR
			/*if (ez->n_cigar > 0)
				mm_append_cigar(r, ez->n_cigar, ez->cigar);
			if (ez->zdropped) { // truncated by Z-drop; TODO: sometimes Z-drop kicks in because the next seed placement is wrong. This can be fixed in principle.
				for (j = i - 1; j >= 0; --j)
					if ((int32_t)a[as1 + j].x <= rs + ez->max_t){
						fprintf(stderr, "1.break\n");
                        break;
                    }
				dropped = 1;
				if (j < 0) j = 0;
				r->p->dp_score += ez->max;
				re1 = rs + (ez->max_t + 1);
				qe1 = qs + (ez->max_q + 1);
				if (cnt1 - (j + 1) >= opt->min_cnt) {
					mm_split_reg(r, r2, as1 + j + 1 - r->as, qlen, a);
					if (zdrop_code == 2) r2->split_inv = 1;
				}
                fprintf(stderr, "2.break\n");
				break;
			} else r->p->dp_score += ez->score;*/
			rs = re, qs = qe;
		}
	}

	if (!dropped && qe < qe0 && re < re0) { // right extension
		qseq = &qseq0[rev][qe];
		mm_idx_getseq(mi, rid, re, re0, tseq);
        //fprintf(stderr, "3.qlen=%d, tlen=%d, w=%d\n", qe0 - qe, re0 - re, bw);
        
        sw_context = (sw_context_t*)malloc(sizeof(sw_context_t));
        add_sw_context(chain_context, sw_context);
        //TODO 保存上下文
        sw_context->w = bw;
        sw_context->pos_flag = 2;
        sw_context->re = re;
        sw_context->qe = qe;
        sw_context->qe0 = qe0;

        sw_context->rs = rs;
        sw_context->qs = qs;
        sw_context->qs0 = qs0;
        
        sw_context->query = (uint8_t*)malloc((qe0 - qe) * sizeof(uint8_t));
        sw_context->target = (uint8_t*)malloc((re0 - re) * sizeof(uint8_t));
        memcpy(sw_context->query, qseq, (qe0 - qe) * sizeof(uint8_t));
        memcpy(sw_context->target, tseq, (re0 - re) * sizeof(uint8_t));
        sw_context->qlen = qe0 - qe;
        sw_context->tlen = re0 - re;
        sw_context->w = bw;
        sw_context->end_bonus = opt->end_bonus;
        sw_context->zdrop = opt->zdrop;
        sw_context->flag = extra_flag|KSW_EZ_EXTZ_ONLY;

        sw_task = create_sw_task(qe0 - qe, sw_context->query, re0 - re, sw_context->target,
                        mat, opt->q, opt->e, opt->q2, opt->e2, bw, 
                        opt->zdrop,
                        opt->end_bonus, extra_flag|KSW_EZ_EXTZ_ONLY, 
                        2);

        small = VSCMIN(sw_task->qlen, sw_task->tlen);
        small = VSCMIN(small, sw_task->w);
        if(sw_task->qlen >= 16300 || sw_task->tlen >= 16300 || small>=1024) {
            chain_context->soft_sw_num++;
            //造一个小的激励让fpga做，保持chain task数据和结果的完整性
            sw_task_t* tmp_sw_task = create_sw_task(64, sw_context->query, 64, sw_context->target,
                                                    mat, opt->q, opt->e, opt->q2, opt->e2, 75, 
                                                    r->split_inv? opt->zdrop_inv : opt->zdrop,
                                                    opt->end_bonus, extra_flag|KSW_EZ_EXTZ_ONLY|KSW_EZ_RIGHT|KSW_EZ_REV_CIGAR, 
                                                    0);
            add_sw_task(chain_task, tmp_sw_task);
            //设置该sw的index值，用来进行结果替换
            sw_task->sw_index = sw_index;
            add_sw_task(soft_sw_task, sw_task);
        }
        else {
            add_sw_task(chain_task, sw_task);
        }
        sw_index++;
		/*mm_align_pair(km, opt, qe0 - qe, qseq, re0 - re, tseq, mat, bw, opt->end_bonus, opt->zdrop, extra_flag|KSW_EZ_EXTZ_ONLY, ez);
		if (ez->n_cigar > 0) {
			mm_append_cigar(r, ez->n_cigar, ez->cigar);
			r->p->dp_score += ez->max;
		}
		re1 = re + (ez->reach_end? ez->mqe_t + 1 : ez->max_t + 1);
		qe1 = qe + (ez->reach_end? qe0 - qe : ez->max_q + 1);*/
	}
	//assert(qe1 <= qlen);

    chain_context->rid = rid;
    chain_context->re0 = re0;
    chain_context->rs0 = rs0;
    chain_context->tseq = (uint8_t*)malloc(tseq_len);
    memcpy(chain_context->tseq, tseq, tseq_len);
    chain_context->qseq0[0] = (uint8_t*)malloc(qlen * 2);
    chain_context->qseq0[1] = chain_context->qseq0[0] + qlen;
    memcpy(chain_context->qseq0[0], qseq0[0], qlen * 2);

    end = realtime_msec();
    sw_task_time[tid] += (end - start);
    
	/*r->rs = rs1, r->re = re1;
	if (rev) r->qs = qlen - qe1, r->qe = qlen - qs1;
	else r->qs = qs1, r->qe = qe1;

	assert(re1 - rs1 <= re0 - rs0);
	if (r->p) {
		mm_idx_getseq(mi, rid, rs1, re1, tseq);
		mm_update_extra(r, &qseq0[r->rev][qs1], tseq, mat, opt->q, opt->e);
		if (rev && r->p->trans_strand)
			r->p->trans_strand ^= 3; // flip to the read strand
	}*/
    
    //先做软件的sw任务，保证硬件处理完成后软件已经有了结果
    if(soft_sw_task->sw_num > 0) {
        int index = 0;
        start = realtime_msec();
        sw_result_t* result = create_result();
        for(index = 0; index < soft_sw_task->sw_num; index++) {
            ksw_extz_t* ez = (ksw_extz_t*)malloc(sizeof(ksw_extz_t));
            memset(ez, 0, sizeof(ksw_extz_t));
            sw_task_t* task = soft_sw_task->sw_tasks[index];
            ksw_extd2_sse(NULL, task->qlen, task->query, task->tlen, task->target, 5, task->mat, task->q, task->e, task->q2, task->e2, task->w, task->zdrop, task->end_bonus, task->flag, ez);

            add_result(result, ez, soft_sw_task->sw_tasks[index]->sw_index);
        }
        result->read_id = soft_sw_task->read_id;
        result->chain_id = soft_sw_task->chain_id;
        
        destroy_chain_sw_task(soft_sw_task);
        pthread_mutex_lock(&long_chain_mutex);
        long_chain_counter++;
        pthread_mutex_unlock(&long_chain_mutex);
        chain_context->soft_sw_result = result;     //将软件处理的sw结果挂在chain context上
        end = realtime_msec();
        sw_soft_time[tid] += (end - start);
    }
    if(chain_task->sw_num > 0) {
        double start, end;
        start = realtime_msec();
        process_task(send_task, chain_task, tid, read_index, read_num);  //将任务放到待发送队列
        end = realtime_msec();
        process_task_time[tid] += (end - start);
    }
    
	kfree(km, tseq);
}

int mm_align1_inv(void *km, const mm_mapopt_t *opt, const mm_idx_t *mi, int qlen, uint8_t *qseq0[2], const mm_reg1_t *r1, const mm_reg1_t *r2, mm_reg1_t *r_inv, ksw_extz_t *ez)
{
	int tl, ql, score, ret = 0, q_off, t_off;
	uint8_t *tseq, *qseq;
	int8_t mat[25];
	void *qp;

	memset(r_inv, 0, sizeof(mm_reg1_t));
	if (!(r1->split&1) || !(r2->split&2)) return 0;
	if (r1->id != r1->parent && r1->parent != MM_PARENT_TMP_PRI) return 0;
	if (r2->id != r2->parent && r2->parent != MM_PARENT_TMP_PRI) return 0;
	if (r1->rid != r2->rid || r1->rev != r2->rev) return 0;
	ql = r1->rev? r1->qs - r2->qe : r2->qs - r1->qe;
	tl = r2->rs - r1->re;
	if (ql < opt->min_chain_score || ql > opt->max_gap) return 0;
	if (tl < opt->min_chain_score || tl > opt->max_gap) return 0;

	ksw_gen_simple_mat(5, mat, opt->a, opt->b);
	tseq = (uint8_t*)kmalloc(km, tl);
	mm_idx_getseq(mi, r1->rid, r1->re, r2->rs, tseq);
	qseq = r1->rev? &qseq0[0][r2->qe] : &qseq0[1][qlen - r2->qs];

	mm_seq_rev(ql, qseq);
	mm_seq_rev(tl, tseq);
	qp = ksw_ll_qinit(km, 2, ql, qseq, 5, mat);
	score = ksw_ll_i16(qp, tl, tseq, opt->q, opt->e, &q_off, &t_off);
	kfree(km, qp);
	mm_seq_rev(ql, qseq);
	mm_seq_rev(tl, tseq);
	if (score < opt->min_dp_max) goto end_align1_inv;
	q_off = ql - (q_off + 1), t_off = tl - (t_off + 1);
	mm_align_pair(km, opt, ql - q_off, qseq + q_off, tl - t_off, tseq + t_off, mat, (int)(opt->bw * 1.5), -1, opt->zdrop, KSW_EZ_EXTZ_ONLY, ez);
	if (ez->n_cigar == 0) goto end_align1_inv; // should never be here
	mm_append_cigar(r_inv, ez->n_cigar, ez->cigar);
	r_inv->p->dp_score = ez->max;
	r_inv->id = -1;
	r_inv->parent = MM_PARENT_UNSET;
	r_inv->inv = 1;
	r_inv->rev = !r1->rev;
	r_inv->rid = r1->rid;
	r_inv->div = -1.0f;
	if (r_inv->rev == 0) {
		r_inv->qs = r2->qe + q_off;
		r_inv->qe = r_inv->qs + ez->max_q + 1;
	} else {
		r_inv->qe = r2->qs - q_off;
		r_inv->qs = r_inv->qe - (ez->max_q + 1);
	}
	r_inv->rs = r1->re + t_off;
	r_inv->re = r_inv->rs + ez->max_t + 1;
	mm_update_extra(r_inv, &qseq[q_off], &tseq[t_off], mat, opt->q, opt->e);
	ret = 1;
end_align1_inv:
	kfree(km, tseq);
	return ret;
}

mm_reg1_t *mm_insert_reg(const mm_reg1_t *r, int i, int *n_regs, mm_reg1_t *regs)
{
	regs = (mm_reg1_t*)realloc(regs, (*n_regs + 1) * sizeof(mm_reg1_t));
	if (i + 1 != *n_regs)
		memmove(&regs[i + 2], &regs[i + 1], sizeof(mm_reg1_t) * (*n_regs - i - 1));
	regs[i + 1] = *r;
	++*n_regs;
	return regs;
}

void mm_align_skeleton(void *km, const mm_mapopt_t *opt, const mm_idx_t *mi, int qlen, const char *qstr, int *n_regs_, mm_reg1_t *regs, mm128_t *a, long read_index, context_t* context, user_params_t* params, int tid)
{
	extern unsigned char seq_nt4_table[256];
	int32_t i, n_regs = *n_regs_, n_a;
	uint8_t *qseq0[2];
	//ksw_extz_t ez;

	// encode the query sequence
	qseq0[0] = (uint8_t*)kmalloc(km, qlen * 2);
	qseq0[1] = qseq0[0] + qlen;
	for (i = 0; i < qlen; ++i) {
		qseq0[0][i] = seq_nt4_table[(uint8_t)qstr[i]];
		qseq0[1][qlen - 1 - i] = qseq0[0][i] < 4? 3 - qseq0[0][i] : 4;
	}

    // align through seed hits
	n_a = mm_squeeze_a(km, n_regs, regs, a);

    //保存regs和a的原始数据
    context->regs0_ori = (mm_reg1_t*)malloc(n_regs * sizeof(mm_reg1_t));
    memcpy(context->regs0_ori, regs, n_regs * sizeof(mm_reg1_t));
    context->regs0 = regs;
    context->a_ori = (mm128_t*)malloc(n_a * sizeof(mm128_t));
    memcpy(context->a_ori, a, n_a * sizeof(mm128_t));
    context->qseq0[0] = (uint8_t*)malloc(qlen * 2);
    context->qseq0[1] = context->qseq0[0] + qlen;
    memcpy(context->qseq0[0], qseq0[0], qlen * 2);
    context->qlen = qlen;
    context->a = a;
    context->n_regs0 = n_regs;
    context->n_a = n_a;
    context->read_index = read_index;
    context->next_chain = 0;        //表示在进入save_read_result函数后从第0条chain开始处理

	//memset(&ez, 0, sizeof(ksw_extz_t));
    assert(params->chain_num[read_index] == 0xffffffff);
    params->chain_num[read_index] = n_regs;
    params->read_results[read_index].chain_result_num = 0;
    params->read_results[read_index].chain_results = (sw_result_t**)malloc(n_regs * sizeof(sw_result_t*));
    memset(params->read_results[read_index].chain_results, 0, n_regs * sizeof(sw_result_t*));
    if(n_regs == 0) {
        params->read_is_complete[read_index] = 1;
        __sync_add_and_fetch(&params->zero_seed_num, 1);
        __sync_sub_and_fetch(&params->read_num, 1);
    }

    assert(tid >= 0);
    
	for (i = 0; i < n_regs; ++i) {
		mm_reg1_t r2;
        chain_context_t* chain_context = create_chain_context();
        add_chain_context(context, chain_context);
        chain_context->i = i;

        //fprintf(stderr, "n_regs=%d, i=%d, read_index=%ld\n", n_regs, i, read_index);
		if ((opt->flag&MM_F_SPLICE) && (opt->flag&MM_F_SPLICE_FOR) && (opt->flag&MM_F_SPLICE_REV)) { // then do two rounds of alignments for both strands
			/*mm_reg1_t s[2], s2[2];
			int which, trans_strand;
			s[0] = s[1] = regs[i];
			mm_align1(km, opt, mi, qlen, qseq0, &s[0], &s2[0], n_a, a, NULL, MM_F_SPLICE_FOR, read_index);
			mm_align1(km, opt, mi, qlen, qseq0, &s[1], &s2[1], n_a, a, NULL, MM_F_SPLICE_REV, read_index);
			if (s[0].p->dp_score > s[1].p->dp_score) which = 0, trans_strand = 1;
			else if (s[0].p->dp_score < s[1].p->dp_score) which = 1, trans_strand = 2;
			else trans_strand = 3, which = (qlen + s[0].p->dp_score) & 1; // randomly choose a strand, effectively
			if (which == 0) {
				regs[i] = s[0], r2 = s2[0];
				free(s[1].p);
			} else {
				regs[i] = s[1], r2 = s2[1];
				free(s[0].p);
			}
			regs[i].p->trans_strand = trans_strand;*/
		} else { // one round of alignment
			mm_align1(km, opt, mi, qlen, qseq0, &regs[i], &r2, n_a, a, NULL, opt->flag, read_index, i, context, chain_context, &params->send_task[tid], tid, params->num);
			if (opt->flag&MM_F_SPLICE)
				regs[i].p->trans_strand = opt->flag&MM_F_SPLICE_FOR? 1 : 2;
		}
		/*if (r2.cnt > 0) regs = mm_insert_reg(&r2, i, &n_regs, regs);
		if (i > 0 && regs[i].split_inv) {
			if (mm_align1_inv(km, opt, mi, qlen, qseq0, &regs[i-1], &regs[i], &r2, &ez)) {
				regs = mm_insert_reg(&r2, i, &n_regs, regs);
				++i; // skip the inserted INV alignment
			}
		}*/
	}
	//*n_regs_ = n_regs;
	//kfree(km, qseq0[0]);
	//kfree(km, ez.cigar);
	//mm_filter_regs(km, opt, qlen, n_regs_, regs);
	//mm_hit_sort_by_dp(km, n_regs_, regs);
	return;
}


static void mm_align1_ori(void *km, const mm_mapopt_t *opt, const mm_idx_t *mi, int qlen, uint8_t *qseq0[2], mm_reg1_t *r, mm_reg1_t *r2, int n_a, mm128_t *a, ksw_extz_t *ez, int splice_flag)
{
	int is_sr = !!(opt->flag & MM_F_SR), is_splice = !!(opt->flag & MM_F_SPLICE);
	int32_t rid = a[r->as].x<<1>>33, rev = a[r->as].x>>63, as1, cnt1;
	uint8_t *tseq, *qseq;
	int32_t i, l, bw, dropped = 0, extra_flag = 0, rs0, re0, qs0, qe0;
	int32_t rs, re, qs, qe;
	int32_t rs1, qs1, re1, qe1;
	int8_t mat[25];

	if (is_sr) assert(!(mi->flag & MM_I_HPC)); // HPC won't work with SR because with HPC we can't easily tell if there is a gap

	r2->cnt = 0;
	if (r->cnt == 0) return;
	ksw_gen_simple_mat(5, mat, opt->a, opt->b);
	bw = (int)(opt->bw * 1.5 + 1.);

	if (is_sr && !(mi->flag & MM_I_HPC)) {
		mm_max_stretch(opt, r, a, &as1, &cnt1);
		rs = (int32_t)a[as1].x + 1 - (int32_t)(a[as1].y>>32&0xff);
		qs = (int32_t)a[as1].y + 1 - (int32_t)(a[as1].y>>32&0xff);
		re = (int32_t)a[as1+cnt1-1].x + 1;
		qe = (int32_t)a[as1+cnt1-1].y + 1;
	} else {
		if (is_splice) {
			mm_fix_bad_ends_splice(km, opt, mi, r, mat, qlen, qseq0, a, &as1, &cnt1);
		} else {
			mm_fix_bad_ends(r, a, opt->bw, opt->min_chain_score * 2, &as1, &cnt1);
		}
		mm_filter_bad_seeds(km, as1, cnt1, a, 10, 40, opt->max_gap>>1, 10);
		mm_adjust_minier(mi, qseq0, &a[as1], &rs, &qs);
		mm_adjust_minier(mi, qseq0, &a[as1 + cnt1 - 1], &re, &qe);
	}
	assert(cnt1 > 0);

	if (is_splice) {
		if (splice_flag & MM_F_SPLICE_FOR) extra_flag |= rev? KSW_EZ_SPLICE_REV : KSW_EZ_SPLICE_FOR;
		if (splice_flag & MM_F_SPLICE_REV) extra_flag |= rev? KSW_EZ_SPLICE_FOR : KSW_EZ_SPLICE_REV;
		if (opt->flag & MM_F_SPLICE_FLANK) extra_flag |= KSW_EZ_SPLICE_FLANK;
	}

	/* Look for the start and end of regions to perform DP. This sounds easy
	 * but is in fact tricky. Excessively small regions lead to unnecessary
	 * clippings and lose alignable sequences. Excessively large regions
	 * occasionally lead to large overlaps between two chains and may cause
	 * loss of alignments in corner cases. */
	if (is_sr) {
		qs0 = 0, qe0 = qlen;
		l = qs;
		l += l * opt->a + opt->end_bonus > opt->q? (l * opt->a + opt->end_bonus - opt->q) / opt->e : 0;
		rs0 = rs - l > 0? rs - l : 0;
		l = qlen - qe;
		l += l * opt->a + opt->end_bonus > opt->q? (l * opt->a + opt->end_bonus - opt->q) / opt->e : 0;
		re0 = re + l < mi->seq[rid].len? re + l : mi->seq[rid].len;
	} else {
		// compute rs0 and qs0
		rs0 = (int32_t)a[r->as].x + 1 - (int32_t)(a[r->as].y>>32&0xff);
		qs0 = (int32_t)a[r->as].y + 1 - (int32_t)(a[r->as].y>>32&0xff);
		if (rs0 < 0) rs0 = 0; // this may happen when HPC is in use
		assert(qs0 >= 0); // this should never happen, or it is logic error
		rs1 = qs1 = 0;
		for (i = r->as - 1, l = 0; i >= 0 && a[i].x>>32 == a[r->as].x>>32; --i) { // inspect nearby seeds
			int32_t x = (int32_t)a[i].x + 1 - (int32_t)(a[i].y>>32&0xff);
			int32_t y = (int32_t)a[i].y + 1 - (int32_t)(a[i].y>>32&0xff);
			if (x < rs0 && y < qs0) {
				if (++l > opt->min_cnt) {
					l = rs0 - x > qs0 - y? rs0 - x : qs0 - y;
					rs1 = rs0 - l, qs1 = qs0 - l;
					break;
				}
			}
		}
		if (qs > 0 && rs > 0) {
			l = qs < opt->max_gap? qs : opt->max_gap;
			qs1 = qs1 > qs - l? qs1 : qs - l;
			qs0 = qs0 < qs1? qs0 : qs1; // at least include qs0
			l += l * opt->a > opt->q? (l * opt->a - opt->q) / opt->e : 0;
			l = l < opt->max_gap? l : opt->max_gap;
			l = l < rs? l : rs;
			rs1 = rs1 > rs - l? rs1 : rs - l;
			rs0 = rs0 < rs1? rs0 : rs1;
		} else rs0 = rs, qs0 = qs;
		// compute re0 and qe0
		re0 = (int32_t)a[r->as + r->cnt - 1].x + 1;
		qe0 = (int32_t)a[r->as + r->cnt - 1].y + 1;
		re1 = mi->seq[rid].len, qe1 = qlen;
		for (i = r->as + r->cnt, l = 0; i < n_a && a[i].x>>32 == a[r->as].x>>32; ++i) { // inspect nearby seeds
			int32_t x = (int32_t)a[i].x + 1;
			int32_t y = (int32_t)a[i].y + 1;
			if (x > re0 && y > qe0) {
				if (++l > opt->min_cnt) {
					l = x - re0 > y - qe0? x - re0 : y - qe0;
					re1 = re0 + l, qe1 = qe0 + l;
					break;
				}
			}
		}
		if (qe < qlen && re < mi->seq[rid].len) {
			l = qlen - qe < opt->max_gap? qlen - qe : opt->max_gap;
			qe1 = qe1 < qe + l? qe1 : qe + l;
			qe0 = qe0 > qe1? qe0 : qe1; // at least include qe0
			l += l * opt->a > opt->q? (l * opt->a - opt->q) / opt->e : 0;
			l = l < opt->max_gap? l : opt->max_gap;
			l = l < mi->seq[rid].len - re? l : mi->seq[rid].len - re;
			re1 = re1 < re + l? re1 : re + l;
			re0 = re0 > re1? re0 : re1;
		} else re0 = re, qe0 = qe;
	}
	if (a[r->as].y & MM_SEED_SELF) {
		int max_ext = r->qs > r->rs? r->qs - r->rs : r->rs - r->qs;
		if (r->rs - rs0 > max_ext) rs0 = r->rs - max_ext;
		if (r->qs - qs0 > max_ext) qs0 = r->qs - max_ext;
		max_ext = r->qe > r->re? r->qe - r->re : r->re - r->qe;
		if (re0 - r->re > max_ext) re0 = r->re + max_ext;
		if (qe0 - r->qe > max_ext) qe0 = r->qe + max_ext;
	}

	assert(re0 > rs0);
	tseq = (uint8_t*)kmalloc(km, re0 - rs0);

	if (qs > 0 && rs > 0) { // left extension 左边一条做sw，lvjingbang
		qseq = &qseq0[rev][qs0];
		mm_idx_getseq(mi, rid, rs0, rs, tseq);
		mm_seq_rev(qs - qs0, qseq);
		mm_seq_rev(rs - rs0, tseq);
		mm_align_pair(km, opt, qs - qs0, qseq, rs - rs0, tseq, mat, bw, opt->end_bonus, r->split_inv? opt->zdrop_inv : opt->zdrop, extra_flag|KSW_EZ_EXTZ_ONLY|KSW_EZ_RIGHT|KSW_EZ_REV_CIGAR, ez);
		if (ez->n_cigar > 0) {
			mm_append_cigar(r, ez->n_cigar, ez->cigar);
			r->p->dp_score += ez->max;
		}
		rs1 = rs - (ez->reach_end? ez->mqe_t + 1 : ez->max_t + 1);
		qs1 = qs - (ez->reach_end? qs - qs0 : ez->max_q + 1);
		mm_seq_rev(qs - qs0, qseq);
	} else rs1 = rs, qs1 = qs;
	re1 = rs, qe1 = qs;
	assert(qs1 >= 0 && rs1 >= 0);

	for (i = is_sr? cnt1 - 1 : 1; i < cnt1; ++i) { // gap filling，中间n条做sw，lvjingbang
		if ((a[as1+i].y & (MM_SEED_IGNORE|MM_SEED_TANDEM)) && i != cnt1 - 1) continue;
		if (is_sr && !(mi->flag & MM_I_HPC)) {
			re = (int32_t)a[as1 + i].x + 1;
			qe = (int32_t)a[as1 + i].y + 1;
		} else mm_adjust_minier(mi, qseq0, &a[as1 + i], &re, &qe);
		re1 = re, qe1 = qe;
		if (i == cnt1 - 1 || (a[as1+i].y&MM_SEED_LONG_JOIN) || (qe - qs >= opt->min_ksw_len && re - rs >= opt->min_ksw_len)) {
			int j, bw1 = bw, zdrop_code;
			if (a[as1+i].y & MM_SEED_LONG_JOIN)
				bw1 = qe - qs > re - rs? qe - qs : re - rs;		//需要判断这个bw1是否大于bw，如果大于了，要在软件做sw，lvjingbang
			// perform alignment
			qseq = &qseq0[rev][qs];
			mm_idx_getseq(mi, rid, rs, re, tseq);
			if (is_sr) { // perform ungapped alignment
				assert(qe - qs == re - rs);
				ksw_reset_extz(ez);
				for (j = 0, ez->score = 0; j < qe - qs; ++j) {
					if (qseq[j] >= 4 || tseq[j] >= 4) ez->score += opt->e2;
					else ez->score += qseq[j] == tseq[j]? opt->a : -opt->b;
				}
				ez->cigar = ksw_push_cigar(km, &ez->n_cigar, &ez->m_cigar, ez->cigar, 0, qe - qs);
			} else { // perform normal gapped alignment
				mm_align_pair(km, opt, qe - qs, qseq, re - rs, tseq, mat, bw1, -1, opt->zdrop, extra_flag|KSW_EZ_APPROX_MAX, ez); // first pass: with approximate Z-drop
			}
			// test Z-drop and inversion Z-drop
			if ((zdrop_code = mm_test_zdrop(km, opt, qseq, tseq, ez->n_cigar, ez->cigar, mat)) != 0)
				mm_align_pair(km, opt, qe - qs, qseq, re - rs, tseq, mat, bw1, -1, zdrop_code == 2? opt->zdrop_inv : opt->zdrop, extra_flag, ez); // second pass: lift approximate
			// update CIGAR
			if (ez->n_cigar > 0)
				mm_append_cigar(r, ez->n_cigar, ez->cigar);
			if (ez->zdropped) { // truncated by Z-drop; TODO: sometimes Z-drop kicks in because the next seed placement is wrong. This can be fixed in principle.
				for (j = i - 1; j >= 0; --j)
					if ((int32_t)a[as1 + j].x <= rs + ez->max_t)
						break;
				dropped = 1;
				if (j < 0) j = 0;
				r->p->dp_score += ez->max;
				re1 = rs + (ez->max_t + 1);
				qe1 = qs + (ez->max_q + 1);
				if (cnt1 - (j + 1) >= opt->min_cnt) {
					mm_split_reg(r, r2, as1 + j + 1 - r->as, qlen, a);
					if (zdrop_code == 2) r2->split_inv = 1;
				}
				break;
			} else r->p->dp_score += ez->score;
			rs = re, qs = qe;
		}
	}

	if (!dropped && qe < qe0 && re < re0) { // right extension  右边一条做sw，lvjingbang
		qseq = &qseq0[rev][qe];
		mm_idx_getseq(mi, rid, re, re0, tseq);
		mm_align_pair(km, opt, qe0 - qe, qseq, re0 - re, tseq, mat, bw, opt->end_bonus, opt->zdrop, extra_flag|KSW_EZ_EXTZ_ONLY, ez);
		if (ez->n_cigar > 0) {
			mm_append_cigar(r, ez->n_cigar, ez->cigar);
			r->p->dp_score += ez->max;
		}
		re1 = re + (ez->reach_end? ez->mqe_t + 1 : ez->max_t + 1);
		qe1 = qe + (ez->reach_end? qe0 - qe : ez->max_q + 1);
	}
	assert(qe1 <= qlen);

	r->rs = rs1, r->re = re1;
	if (rev) r->qs = qlen - qe1, r->qe = qlen - qs1;
	else r->qs = qs1, r->qe = qe1;

	assert(re1 - rs1 <= re0 - rs0);
	if (r->p) {
		mm_idx_getseq(mi, rid, rs1, re1, tseq);
		mm_update_extra(r, &qseq0[r->rev][qs1], tseq, mat, opt->q, opt->e);
		if (rev && r->p->trans_strand)
			r->p->trans_strand ^= 3; // flip to the read strand
	}

	kfree(km, tseq);
}

mm_reg1_t * mm_align_skeleton_ori(void *km, const mm_mapopt_t *opt, const mm_idx_t *mi, int qlen, const char *qstr, int *n_regs_, mm_reg1_t *regs, int n_a, mm128_t *a)
{
	extern unsigned char seq_nt4_table[256];
	int32_t i, n_regs = *n_regs_;
	uint8_t *qseq0[2];
	ksw_extz_t ez;

	// encode the query sequence
	qseq0[0] = (uint8_t*)kmalloc(km, qlen * 2);
	qseq0[1] = qseq0[0] + qlen;
	for (i = 0; i < qlen; ++i) {
		qseq0[0][i] = seq_nt4_table[(uint8_t)qstr[i]];
		qseq0[1][qlen - 1 - i] = qseq0[0][i] < 4? 3 - qseq0[0][i] : 4;
	}

	memset(&ez, 0, sizeof(ksw_extz_t));
	for (i = 0; i < n_regs; ++i) {
		mm_reg1_t r2;

        //fprintf(stderr, "n_regs=%d, i=%d, read_index=%ld\n", n_regs, i, read_index);
		if ((opt->flag&MM_F_SPLICE) && (opt->flag&MM_F_SPLICE_FOR) && (opt->flag&MM_F_SPLICE_REV)) { // then do two rounds of alignments for both strands
			mm_reg1_t s[2], s2[2];
			int which, trans_strand;
			s[0] = s[1] = regs[i];
			mm_align1_ori(km, opt, mi, qlen, qseq0, &s[0], &s2[0], n_a, a, &ez, MM_F_SPLICE_FOR);
			mm_align1_ori(km, opt, mi, qlen, qseq0, &s[1], &s2[1], n_a, a, &ez, MM_F_SPLICE_REV);
			if (s[0].p->dp_score > s[1].p->dp_score) which = 0, trans_strand = 1;
			else if (s[0].p->dp_score < s[1].p->dp_score) which = 1, trans_strand = 2;
			else trans_strand = 3, which = (qlen + s[0].p->dp_score) & 1; // randomly choose a strand, effectively
			if (which == 0) {
				regs[i] = s[0], r2 = s2[0];
				free(s[1].p);
			} else {
				regs[i] = s[1], r2 = s2[1];
				free(s[0].p);
			}
			regs[i].p->trans_strand = trans_strand;
		} else { // one round of alignment
			mm_align1_ori(km, opt, mi, qlen, qseq0, &regs[i], &r2, n_a, a, &ez, opt->flag);
			if (opt->flag&MM_F_SPLICE)
				regs[i].p->trans_strand = opt->flag&MM_F_SPLICE_FOR? 1 : 2;
		}
		if (r2.cnt > 0) regs = mm_insert_reg(&r2, i, &n_regs, regs);
		if (i > 0 && regs[i].split_inv) {
			if (mm_align1_inv(km, opt, mi, qlen, qseq0, &regs[i-1], &regs[i], &r2, &ez)) {
				regs = mm_insert_reg(&r2, i, &n_regs, regs);
				++i; // skip the inserted INV alignment
			}
		}
	}
    *n_regs_ = n_regs;
    kfree(km, qseq0[0]);
    kfree(km, ez.cigar);
    mm_filter_regs(km, opt, qlen, n_regs_, regs);
    mm_hit_sort_by_dp(km, n_regs_, regs);
    return regs;
}



void mm_align2(void *km, const mm_mapopt_t *opt, const mm_idx_t *mi, int qlen, uint8_t *qseq0[2], mm_reg1_t *r, mm_reg1_t *r2, int n_a, mm128_t *a, ksw_extz_t *ez, int splice_flag, long read_index, int chain_index, chain_context_t* chain_context)
{
	int is_sr = !!(opt->flag & MM_F_SR), is_splice = !!(opt->flag & MM_F_SPLICE);
	int32_t rid = a[r->as].x<<1>>33, rev = a[r->as].x>>63, as1, cnt1;
	uint8_t *tseq, *qseq;
	int32_t i, l, bw, dropped = 0, extra_flag = 0, rs0, re0, qs0, qe0;
	int32_t rs, re, qs, qe;
	int32_t rs1, qs1, re1, qe1;
	int8_t mat[25];
    int sw_index = 0;

	if (is_sr) assert(!(mi->flag & MM_I_HPC)); // HPC won't work with SR because with HPC we can't easily tell if there is a gap

	r2->cnt = 0;
	if (r->cnt == 0) return;
	ksw_gen_simple_mat(5, mat, opt->a, opt->b);
	bw = (int)(opt->bw * 1.5 + 1.);

	if (is_sr && !(mi->flag & MM_I_HPC)) {
		mm_max_stretch(opt, r, a, &as1, &cnt1);
		rs = (int32_t)a[as1].x + 1 - (int32_t)(a[as1].y>>32&0xff);
		qs = (int32_t)a[as1].y + 1 - (int32_t)(a[as1].y>>32&0xff);
		re = (int32_t)a[as1+cnt1-1].x + 1;
		qe = (int32_t)a[as1+cnt1-1].y + 1;
	} else {
		if (is_splice) {
			mm_fix_bad_ends_splice(km, opt, mi, r, mat, qlen, qseq0, a, &as1, &cnt1);
		} else {
			mm_fix_bad_ends(r, a, opt->bw, opt->min_chain_score * 2, &as1, &cnt1);
		}
		mm_filter_bad_seeds(km, as1, cnt1, a, 10, 40, opt->max_gap>>1, 10);
		mm_adjust_minier(mi, qseq0, &a[as1], &rs, &qs);
		mm_adjust_minier(mi, qseq0, &a[as1 + cnt1 - 1], &re, &qe);
	}
	assert(cnt1 > 0);

	if (is_splice) {
		if (splice_flag & MM_F_SPLICE_FOR) extra_flag |= rev? KSW_EZ_SPLICE_REV : KSW_EZ_SPLICE_FOR;
		if (splice_flag & MM_F_SPLICE_REV) extra_flag |= rev? KSW_EZ_SPLICE_FOR : KSW_EZ_SPLICE_REV;
		if (opt->flag & MM_F_SPLICE_FLANK) extra_flag |= KSW_EZ_SPLICE_FLANK;
	}

	/* Look for the start and end of regions to perform DP. This sounds easy
	 * but is in fact tricky. Excessively small regions lead to unnecessary
	 * clippings and lose alignable sequences. Excessively large regions
	 * occasionally lead to large overlaps between two chains and may cause
	 * loss of alignments in corner cases. */
	if (is_sr) {
		qs0 = 0, qe0 = qlen;
		l = qs;
		l += l * opt->a + opt->end_bonus > opt->q? (l * opt->a + opt->end_bonus - opt->q) / opt->e : 0;
		rs0 = rs - l > 0? rs - l : 0;
		l = qlen - qe;
		l += l * opt->a + opt->end_bonus > opt->q? (l * opt->a + opt->end_bonus - opt->q) / opt->e : 0;
		re0 = re + l < mi->seq[rid].len? re + l : mi->seq[rid].len;
	} else {
		// compute rs0 and qs0
		rs0 = (int32_t)a[r->as].x + 1 - (int32_t)(a[r->as].y>>32&0xff);
		qs0 = (int32_t)a[r->as].y + 1 - (int32_t)(a[r->as].y>>32&0xff);
		if (rs0 < 0) rs0 = 0; // this may happen when HPC is in use
		assert(qs0 >= 0); // this should never happen, or it is logic error
		rs1 = qs1 = 0;
		for (i = r->as - 1, l = 0; i >= 0 && a[i].x>>32 == a[r->as].x>>32; --i) { // inspect nearby seeds
			int32_t x = (int32_t)a[i].x + 1 - (int32_t)(a[i].y>>32&0xff);
			int32_t y = (int32_t)a[i].y + 1 - (int32_t)(a[i].y>>32&0xff);
			if (x < rs0 && y < qs0) {
				if (++l > opt->min_cnt) {
					l = rs0 - x > qs0 - y? rs0 - x : qs0 - y;
					rs1 = rs0 - l, qs1 = qs0 - l;
					break;
				}
			}
		}
		if (qs > 0 && rs > 0) {
			l = qs < opt->max_gap? qs : opt->max_gap;
			qs1 = qs1 > qs - l? qs1 : qs - l;
			qs0 = qs0 < qs1? qs0 : qs1; // at least include qs0
			l += l * opt->a > opt->q? (l * opt->a - opt->q) / opt->e : 0;
			l = l < opt->max_gap? l : opt->max_gap;
			l = l < rs? l : rs;
			rs1 = rs1 > rs - l? rs1 : rs - l;
			rs0 = rs0 < rs1? rs0 : rs1;
		} else rs0 = rs, qs0 = qs;
		// compute re0 and qe0
		re0 = (int32_t)a[r->as + r->cnt - 1].x + 1;
		qe0 = (int32_t)a[r->as + r->cnt - 1].y + 1;
		re1 = mi->seq[rid].len, qe1 = qlen;
		for (i = r->as + r->cnt, l = 0; i < n_a && a[i].x>>32 == a[r->as].x>>32; ++i) { // inspect nearby seeds
			int32_t x = (int32_t)a[i].x + 1;
			int32_t y = (int32_t)a[i].y + 1;
			if (x > re0 && y > qe0) {
				if (++l > opt->min_cnt) {
					l = x - re0 > y - qe0? x - re0 : y - qe0;
					re1 = re0 + l, qe1 = qe0 + l;
					break;
				}
			}
		}
		if (qe < qlen && re < mi->seq[rid].len) {
			l = qlen - qe < opt->max_gap? qlen - qe : opt->max_gap;
			qe1 = qe1 < qe + l? qe1 : qe + l;
			qe0 = qe0 > qe1? qe0 : qe1; // at least include qe0
			l += l * opt->a > opt->q? (l * opt->a - opt->q) / opt->e : 0;
			l = l < opt->max_gap? l : opt->max_gap;
			l = l < mi->seq[rid].len - re? l : mi->seq[rid].len - re;
			re1 = re1 < re + l? re1 : re + l;
			re0 = re0 > re1? re0 : re1;
		} else re0 = re, qe0 = qe;
	}
	if (a[r->as].y & MM_SEED_SELF) {
		int max_ext = r->qs > r->rs? r->qs - r->rs : r->rs - r->qs;
		if (r->rs - rs0 > max_ext) rs0 = r->rs - max_ext;
		if (r->qs - qs0 > max_ext) qs0 = r->qs - max_ext;
		max_ext = r->qe > r->re? r->qe - r->re : r->re - r->qe;
		if (re0 - r->re > max_ext) re0 = r->re + max_ext;
		if (qe0 - r->qe > max_ext) qe0 = r->qe + max_ext;
	}

	assert(re0 > rs0);
    int tseq_len = re0 - rs0;
	tseq = (uint8_t*)kmalloc(km, re0 - rs0);

    //创建一个chain的sw任务
    int small;
    chain_sw_task_t* chain_task = create_chain_sw_task(read_index, chain_index);
    chain_sw_task_t* soft_sw_task = create_chain_sw_task(read_index, chain_index);
    
    chain_context->rs = rs;
    chain_context->qs = qs;
    chain_context->qs0 = qs0;
    chain_context->rid = rid;
    chain_context->rev = rev;

    sw_task_t* sw_task = NULL;
    sw_context_t* sw_context = NULL;
	if (qs > 0 && rs > 0) { // left extension
		qseq = &qseq0[rev][qs0];
		mm_idx_getseq(mi, rid, rs0, rs, tseq);
		mm_seq_rev(qs - qs0, qseq);
		mm_seq_rev(rs - rs0, tseq);

        sw_context = (sw_context_t*)malloc(sizeof(sw_context_t));
        add_sw_context(chain_context, sw_context);
        //TODO 保存上下文
        sw_context->rs = rs;
        sw_context->qs = qs;
        sw_context->qs0 = qs0;
        sw_context->pos_flag = 0;
        
        sw_context->query = (uint8_t*)malloc((qs - qs0) * sizeof(uint8_t));
        sw_context->target = (uint8_t*)malloc((rs - rs0) * sizeof(uint8_t));
        memcpy(sw_context->query, qseq, (qs - qs0) * sizeof(uint8_t));
        memcpy(sw_context->target, tseq, (rs - rs0) * sizeof(uint8_t));
        sw_context->qlen = qs - qs0;
        sw_context->tlen = rs - rs0;
        sw_context->w = bw;
        sw_context->end_bonus = opt->end_bonus;
        sw_context->zdrop = r->split_inv? opt->zdrop_inv : opt->zdrop;
        sw_context->flag = extra_flag|KSW_EZ_EXTZ_ONLY|KSW_EZ_RIGHT|KSW_EZ_REV_CIGAR;
        
        //fprintf(stderr, "1.qlen=%d, tlen=%d, w=%d\n", qs - qs0, rs - rs0, bw);
        sw_task = create_sw_task(qs - qs0, sw_context->query, rs - rs0, sw_context->target,
                        mat, opt->q, opt->e, opt->q2, opt->e2, bw, 
                        r->split_inv? opt->zdrop_inv : opt->zdrop,
                        opt->end_bonus, extra_flag|KSW_EZ_EXTZ_ONLY|KSW_EZ_RIGHT|KSW_EZ_REV_CIGAR, 
                        0);

        small = VSCMIN(sw_task->qlen, sw_task->tlen);
        small = VSCMIN(small, sw_task->w);
        if(sw_task->qlen >= 16300 || sw_task->tlen >= 16300 || small>=1024) {
            chain_context->soft_sw_num++;
            //造一个小的激励让fpga做，保持chain task数据和结果的完整性
            sw_task_t* tmp_sw_task = create_sw_task(100, sw_context->query, 100, sw_context->target,
                                                    mat, opt->q, opt->e, opt->q2, opt->e2, 75, 
                                                    r->split_inv? opt->zdrop_inv : opt->zdrop,
                                                    opt->end_bonus, extra_flag|KSW_EZ_EXTZ_ONLY|KSW_EZ_RIGHT|KSW_EZ_REV_CIGAR, 
                                                    0);
            add_sw_task(chain_task, tmp_sw_task);
            //设置该sw的index值，用来进行结果替换
            sw_task->sw_index = sw_index;
            add_sw_task(soft_sw_task, sw_task);
        }
        else {
            add_sw_task(chain_task, sw_task);
        }
        sw_index++;
		/*mm_align_pair(km, opt, qs - qs0, qseq, rs - rs0, tseq, mat, bw, opt->end_bonus, r->split_inv? opt->zdrop_inv : opt->zdrop, extra_flag|KSW_EZ_EXTZ_ONLY|KSW_EZ_RIGHT|KSW_EZ_REV_CIGAR, ez);
        if (ez->n_cigar > 0) {
			mm_append_cigar(r, ez->n_cigar, ez->cigar);
			r->p->dp_score += ez->max;
		}
		rs1 = rs - (ez->reach_end? ez->mqe_t + 1 : ez->max_t + 1);
		qs1 = qs - (ez->reach_end? qs - qs0 : ez->max_q + 1);*/
		mm_seq_rev(qs - qs0, qseq);
	}
    else {
        rs1 = rs;
        qs1 = qs;
    }
	re1 = rs, qe1 = qs;
	//assert(qs1 >= 0 && rs1 >= 0);

	for (i = is_sr? cnt1 - 1 : 1; i < cnt1; ++i) { // gap filling
		if ((a[as1+i].y & (MM_SEED_IGNORE|MM_SEED_TANDEM)) && i != cnt1 - 1) continue;
		if (is_sr && !(mi->flag & MM_I_HPC)) {
			re = (int32_t)a[as1 + i].x + 1;
			qe = (int32_t)a[as1 + i].y + 1;
		} else mm_adjust_minier(mi, qseq0, &a[as1 + i], &re, &qe);
		re1 = re, qe1 = qe;
		if (i == cnt1 - 1 || (a[as1+i].y&MM_SEED_LONG_JOIN) || (qe - qs >= opt->min_ksw_len && re - rs >= opt->min_ksw_len)) {
			//int j, bw1 = bw, zdrop_code;
            int j, bw1 = bw;
			if (a[as1+i].y & MM_SEED_LONG_JOIN)
				bw1 = qe - qs > re - rs? qe - qs : re - rs;
			// perform alignment
			qseq = &qseq0[rev][qs];
			mm_idx_getseq(mi, rid, rs, re, tseq);
			if (is_sr) { // perform ungapped alignment
				assert(qe - qs == re - rs);
				ksw_reset_extz(ez);
				for (j = 0, ez->score = 0; j < qe - qs; ++j) {
					if (qseq[j] >= 4 || tseq[j] >= 4) ez->score += opt->e2;
					else ez->score += qseq[j] == tseq[j]? opt->a : -opt->b;
				}
				ez->cigar = ksw_push_cigar(km, &ez->n_cigar, &ez->m_cigar, ez->cigar, 0, qe - qs);
			} else { // perform normal gapped alignment
                //fprintf(stderr, "2.qlen=%d, tlen=%d, w=%d\n", qe - qs, re - rs, bw1);
				
                sw_context = (sw_context_t*)malloc(sizeof(sw_context_t));
                add_sw_context(chain_context, sw_context);
                //TODO 保存上下文
                sw_context->i = i;
                sw_context->pos_flag = 1;
                sw_context->cnt1 = cnt1;
                sw_context->re = re;
                sw_context->qe = qe;

                sw_context->rs = rs;
                sw_context->qs = qs;
                sw_context->qs0 = qs0;
                
                sw_context->query = (uint8_t*)malloc((qe - qs) * sizeof(uint8_t));
                sw_context->target = (uint8_t*)malloc((re - rs) * sizeof(uint8_t));
                memcpy(sw_context->query, qseq, (qe - qs) * sizeof(uint8_t));
                memcpy(sw_context->target, tseq, (re - rs) * sizeof(uint8_t));
                sw_context->qlen = qe - qs;
                sw_context->tlen = re - rs;
                sw_context->w = bw1;
                sw_context->end_bonus = -1;
                sw_context->zdrop = opt->zdrop;
                sw_context->flag = extra_flag|KSW_EZ_APPROX_MAX;
                sw_context->zdrop_flag = extra_flag;
                sw_context->as1 = as1;

                sw_task = create_sw_task(qe - qs, sw_context->query, re - rs, sw_context->target,
                        mat, opt->q, opt->e, opt->q2, opt->e2, bw1, 
                        opt->zdrop,
                        -1, extra_flag|KSW_EZ_APPROX_MAX, 
                        1);
                small = VSCMIN(sw_task->qlen, sw_task->tlen);
                small = VSCMIN(small, sw_task->w);
                if(sw_task->qlen >= 16300 || sw_task->tlen >= 16300 || small>=1024) {
                    chain_context->soft_sw_num++;
                    //造一个小的激励让fpga做，保持chain task数据和结果的完整性
                    sw_task_t* tmp_sw_task = create_sw_task(100, sw_context->query, 100, sw_context->target,
                                                            mat, opt->q, opt->e, opt->q2, opt->e2, 75, 
                                                            r->split_inv? opt->zdrop_inv : opt->zdrop,
                                                            opt->end_bonus, extra_flag|KSW_EZ_EXTZ_ONLY|KSW_EZ_RIGHT|KSW_EZ_REV_CIGAR, 
                                                            0);
                    add_sw_task(chain_task, tmp_sw_task);
                    //设置该sw的index值，用来进行结果替换
                    sw_task->sw_index = sw_index;
                    add_sw_task(soft_sw_task, sw_task);
                }
                else {
                    add_sw_task(chain_task, sw_task);
                }
                sw_index++;
                //mm_align_pair(km, opt, qe - qs, qseq, re - rs, tseq, mat, bw1, -1, opt->zdrop, extra_flag|KSW_EZ_APPROX_MAX, ez); // first pass: with approximate Z-drop
            }
			// test Z-drop and inversion Z-drop
			//if ((zdrop_code = mm_test_zdrop(km, opt, qseq, tseq, ez->n_cigar, ez->cigar, mat)) != 0)
			//	mm_align_pair(km, opt, qe - qs, qseq, re - rs, tseq, mat, bw1, -1, zdrop_code == 2? opt->zdrop_inv : opt->zdrop, extra_flag, ez); // second pass: lift approximate
			// update CIGAR
			/*if (ez->n_cigar > 0)
				mm_append_cigar(r, ez->n_cigar, ez->cigar);
			if (ez->zdropped) { // truncated by Z-drop; TODO: sometimes Z-drop kicks in because the next seed placement is wrong. This can be fixed in principle.
				for (j = i - 1; j >= 0; --j)
					if ((int32_t)a[as1 + j].x <= rs + ez->max_t){
						fprintf(stderr, "1.break\n");
                        break;
                    }
				dropped = 1;
				if (j < 0) j = 0;
				r->p->dp_score += ez->max;
				re1 = rs + (ez->max_t + 1);
				qe1 = qs + (ez->max_q + 1);
				if (cnt1 - (j + 1) >= opt->min_cnt) {
					mm_split_reg(r, r2, as1 + j + 1 - r->as, qlen, a);
					if (zdrop_code == 2) r2->split_inv = 1;
				}
                fprintf(stderr, "2.break\n");
				break;
			} else r->p->dp_score += ez->score;*/
			rs = re, qs = qe;
		}
	}

	if (!dropped && qe < qe0 && re < re0) { // right extension
		qseq = &qseq0[rev][qe];
		mm_idx_getseq(mi, rid, re, re0, tseq);
        //fprintf(stderr, "3.qlen=%d, tlen=%d, w=%d\n", qe0 - qe, re0 - re, bw);
        
        sw_context = (sw_context_t*)malloc(sizeof(sw_context_t));
        add_sw_context(chain_context, sw_context);
        //TODO 保存上下文
        sw_context->w = bw;
        sw_context->pos_flag = 2;
        sw_context->re = re;
        sw_context->qe = qe;
        sw_context->qe0 = qe0;

        sw_context->rs = rs;
        sw_context->qs = qs;
        sw_context->qs0 = qs0;
        
        sw_context->query = (uint8_t*)malloc((qe0 - qe) * sizeof(uint8_t));
        sw_context->target = (uint8_t*)malloc((re0 - re) * sizeof(uint8_t));
        memcpy(sw_context->query, qseq, (qe0 - qe) * sizeof(uint8_t));
        memcpy(sw_context->target, tseq, (re0 - re) * sizeof(uint8_t));
        sw_context->qlen = qe0 - qe;
        sw_context->tlen = re0 - re;
        sw_context->w = bw;
        sw_context->end_bonus = opt->end_bonus;
        sw_context->zdrop = opt->zdrop;
        sw_context->flag = extra_flag|KSW_EZ_EXTZ_ONLY;

        sw_task = create_sw_task(qe0 - qe, sw_context->query, re0 - re, sw_context->target,
                        mat, opt->q, opt->e, opt->q2, opt->e2, bw, 
                        opt->zdrop,
                        opt->end_bonus, extra_flag|KSW_EZ_EXTZ_ONLY, 
                        2);

        small = VSCMIN(sw_task->qlen, sw_task->tlen);
        small = VSCMIN(small, sw_task->w);
        if(sw_task->qlen >= 16300 || sw_task->tlen >= 16300 || small>=1024) {
            chain_context->soft_sw_num++;
            //造一个小的激励让fpga做，保持chain task数据和结果的完整性
            sw_task_t* tmp_sw_task = create_sw_task(100, sw_context->query, 100, sw_context->target,
                                                    mat, opt->q, opt->e, opt->q2, opt->e2, 75, 
                                                    r->split_inv? opt->zdrop_inv : opt->zdrop,
                                                    opt->end_bonus, extra_flag|KSW_EZ_EXTZ_ONLY|KSW_EZ_RIGHT|KSW_EZ_REV_CIGAR, 
                                                    0);
            add_sw_task(chain_task, tmp_sw_task);
            //设置该sw的index值，用来进行结果替换
            sw_task->sw_index = sw_index;
            add_sw_task(soft_sw_task, sw_task);
        }
        else {
            add_sw_task(chain_task, sw_task);
        }
        sw_index++;
		/*mm_align_pair(km, opt, qe0 - qe, qseq, re0 - re, tseq, mat, bw, opt->end_bonus, opt->zdrop, extra_flag|KSW_EZ_EXTZ_ONLY, ez);
		if (ez->n_cigar > 0) {
			mm_append_cigar(r, ez->n_cigar, ez->cigar);
			r->p->dp_score += ez->max;
		}
		re1 = re + (ez->reach_end? ez->mqe_t + 1 : ez->max_t + 1);
		qe1 = qe + (ez->reach_end? qe0 - qe : ez->max_q + 1);*/
	}
	//assert(qe1 <= qlen);

    chain_context->rid = rid;
    chain_context->re0 = re0;
    chain_context->rs0 = rs0;
    chain_context->tseq = (uint8_t*)malloc(tseq_len);
    memcpy(chain_context->tseq, tseq, tseq_len);
    chain_context->qseq0[0] = (uint8_t*)malloc(qlen * 2);
    chain_context->qseq0[1] = chain_context->qseq0[0] + qlen;
    memcpy(chain_context->qseq0[0], qseq0[0], qlen * 2);
    
	/*r->rs = rs1, r->re = re1;
	if (rev) r->qs = qlen - qe1, r->qe = qlen - qs1;
	else r->qs = qs1, r->qe = qe1;

	assert(re1 - rs1 <= re0 - rs0);
	if (r->p) {
		mm_idx_getseq(mi, rid, rs1, re1, tseq);
		mm_update_extra(r, &qseq0[r->rev][qs1], tseq, mat, opt->q, opt->e);
		if (rev && r->p->trans_strand)
			r->p->trans_strand ^= 3; // flip to the read strand
	}*/
    //先做软件的sw任务，保证硬件处理完成后软件已经有了结果
    if(soft_sw_task->sw_num > 0) {
        int index = 0;
        sw_result_t* result = create_result();
        for(index = 0; index < soft_sw_task->sw_num; index++) {
            ksw_extz_t* ez = (ksw_extz_t*)malloc(sizeof(ksw_extz_t));
            memset(ez, 0, sizeof(ksw_extz_t));
            sw_task_t* task = soft_sw_task->sw_tasks[index];
            
            ksw_extd2_sse(NULL, task->qlen, task->query, task->tlen, task->target, 5, task->mat, task->q, task->e, task->q2, task->e2, task->w, task->zdrop, task->end_bonus, task->flag, ez);

            add_result(result, ez, soft_sw_task->sw_tasks[index]->sw_index);
        }
        result->read_id = soft_sw_task->read_id;
        result->chain_id = soft_sw_task->chain_id;

        destroy_chain_sw_task(soft_sw_task);
        chain_context->soft_sw_result = result;     //将软件处理的sw结果挂在chain context上
    }
    if(chain_task->sw_num > 0) {

        //process_task(send_task, chain_task, read_index);  //将任务放到待发送队列
        send_to_fpga(&chain_task, 1, chain_task->data_size);

    }
    
	kfree(km, tseq);
}