#ifndef MMPRIV2_H
#define MMPRIV2_H

#include <assert.h>
#include "minimap.h"
#include "bseq.h"

#include "soft_sw.h"

#define MM_PARENT_UNSET   (-1)
#define MM_PARENT_TMP_PRI (-2)

#define MM_DBG_NO_KALLOC     0x1
#define MM_DBG_PRINT_QNAME   0x2
#define MM_DBG_PRINT_SEED    0x4
#define MM_DBG_PRINT_ALN_SEQ 0x8

#define MM_SEED_LONG_JOIN  (1ULL<<40)
#define MM_SEED_IGNORE     (1ULL<<41)
#define MM_SEED_TANDEM     (1ULL<<42)
#define MM_SEED_SELF       (1ULL<<43)

#define MM_SEED_SEG_SHIFT  48
#define MM_SEED_SEG_MASK   (0xffULL<<(MM_SEED_SEG_SHIFT))

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

#define mm_seq4_set(s, i, c) ((s)[(i)>>3] |= (uint32_t)(c) << (((i)&7)<<2))
#define mm_seq4_get(s, i)    ((s)[(i)>>3] >> (((i)&7)<<2) & 0xf)

#ifdef __cplusplus
extern "C" {
#endif

#ifndef KSTRING_T
#define KSTRING_T kstring_t
typedef struct __kstring_t {
	unsigned l, m;
	char *s;
} kstring_t;
#endif

typedef struct {
	int n_u, n_a;
	uint64_t *u;
	mm128_t *a;
} mm_seg_t;

double cputime(void);
double realtime(void);
long peakrss(void);

void radix_sort_128x(mm128_t *beg, mm128_t *end);
void radix_sort_64(uint64_t *beg, uint64_t *end);
uint32_t ks_ksmall_uint32_t(size_t n, uint32_t arr[], size_t kk);

void mm_sketch(void *km, const char *str, int len, int w, int k, uint32_t rid, int is_hpc, mm128_v *p);

void mm_write_sam_hdr(const mm_idx_t *mi, const char *rg, const char *ver, int argc, char *argv[]);
void mm_write_paf(kstring_t *s, const mm_idx_t *mi, const mm_bseq1_t *t, const mm_reg1_t *r, void *km, int opt_flag);
void mm_write_sam(kstring_t *s, const mm_idx_t *mi, const mm_bseq1_t *t, const mm_reg1_t *r, int n_regs, const mm_reg1_t *regs);
void mm_write_sam2(kstring_t *s, const mm_idx_t *mi, const mm_bseq1_t *t, int seg_idx, int reg_idx, int n_seg, const int *n_regs, const mm_reg1_t *const* regs, void *km, int opt_flag);

void mm_idxopt_init(mm_idxopt_t *opt);
const uint64_t *mm_idx_get(const mm_idx_t *mi, uint64_t minier, int *n);
int32_t mm_idx_cal_max_occ(const mm_idx_t *mi, float f);
mm128_t *mm_chain_dp(int max_dist_x, int max_dist_y, int bw, int max_skip, int min_cnt, int min_sc, int is_cdna, int n_segs, int64_t n, mm128_t *a, int *n_u_, uint64_t **_u, void *km);
void mm_align_skeleton(void *km, const mm_mapopt_t *opt, const mm_idx_t *mi, int qlen, const char *qstr, int *n_regs_, mm_reg1_t *regs, mm128_t *a, long read_index, context_t* context, user_params_t* params, int tid);
mm_reg1_t * mm_align_skeleton_ori(void *km, const mm_mapopt_t *opt, const mm_idx_t *mi, int qlen, const char *qstr, int *n_regs_, mm_reg1_t *regs, int n_a, mm128_t *a);
int mm_test_zdrop(void *km, const mm_mapopt_t *opt, const uint8_t *qseq, const uint8_t *tseq, uint32_t n_cigar, uint32_t *cigar, const int8_t *mat);
void mm_align_pair(void *km, const mm_mapopt_t *opt, int qlen, const uint8_t *qseq, int tlen, const uint8_t *tseq, const int8_t *mat, int w, int end_bonus, int zdrop, int flag, ksw_extz_t *ez);
void mm_append_cigar(mm_reg1_t *r, uint32_t n_cigar, uint32_t *cigar);
void mm_update_extra(mm_reg1_t *r, const uint8_t *qseq, const uint8_t *tseq, const int8_t *mat, int8_t q, int8_t e);
int mm_align1_inv(void *km, const mm_mapopt_t *opt, const mm_idx_t *mi, int qlen, uint8_t *qseq0[2], const mm_reg1_t *r1, const mm_reg1_t *r2, mm_reg1_t *r_inv, ksw_extz_t *ez);

mm_reg1_t *mm_gen_regs(void *km, uint32_t hash, int qlen, int n_u, uint64_t *u, mm128_t *a);
void mm_split_reg(mm_reg1_t *r, mm_reg1_t *r2, int n, int qlen, mm128_t *a);
void mm_sync_regs(void *km, int n_regs, mm_reg1_t *regs);
int mm_squeeze_a(void *km, int n_regs, mm_reg1_t *regs, mm128_t *a);
int mm_set_sam_pri(int n, mm_reg1_t *r);
void mm_set_parent(void *km, float mask_level, int n, mm_reg1_t *r, int sub_diff);
void mm_select_sub(void *km, float pri_ratio, int min_diff, int best_n, int *n_, mm_reg1_t *r);
void mm_select_sub_multi(void *km, float pri_ratio, float pri1, float pri2, int max_gap_ref, int min_diff, int best_n, int n_segs, const int *qlens, int *n_, mm_reg1_t *r);
void mm_filter_regs(void *km, const mm_mapopt_t *opt, int qlen, int *n_regs, mm_reg1_t *regs);
void mm_join_long(void *km, const mm_mapopt_t *opt, int qlen, int *n_regs, mm_reg1_t *regs, mm128_t *a);
void mm_hit_sort_by_dp(void *km, int *n_regs, mm_reg1_t *r);
void mm_set_mapq(void *km, int n_regs, mm_reg1_t *regs, int min_chain_sc, int match_sc, int rep_len, int is_sr);

void mm_est_err(const mm_idx_t *mi, int qlen, int n_regs, mm_reg1_t *regs, const mm128_t *a, int32_t n, const uint64_t *mini_pos);

mm_seg_t *mm_seg_gen(void *km, uint32_t hash, int n_segs, const int *qlens, int n_regs0, const mm_reg1_t *regs0, int *n_regs, mm_reg1_t **regs, const mm128_t *a);
void mm_seg_free(void *km, int n_segs, mm_seg_t *segs);
void mm_pair(void *km, int max_gap_ref, int dp_bonus, int sub_diff, int match_sc, const int *qlens, int *n_regs, mm_reg1_t **regs);

void mm_err_puts(const char *str);

#ifdef __cplusplus
}
#endif

#endif
