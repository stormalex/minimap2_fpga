#ifndef _CONTEXT_DATA_H_
#define _CONTEXT_DATA_H_

#include <stdint.h>


#include "minimap.h"
#include "user_to_drv.h"

#define  NT_LAST  (0)
#define  IS_LAST  (1)


//chaindp context
typedef struct {
    const mm_idx_t *mi;
    mm_reg1_t **regs;
    int *n_regs;
    mm_tbuf_t *b;
    mm128_t *mva;
    uint64_t *mini_pos;
    //uint64_t mvn,mvm;
    uint32_t hash;
    int n_segs;
    int qlen_sum;
    int gap_ref;
    int gap_qry;
    int rep_len;
    int n_mini_pos;
    int qlen[MM_MAX_SEG];
    char *qseqs[MM_MAX_SEG];
    char qname[MM_MAX_SEG];
} dp_context_t ;

typedef struct {
    //void *km;
    int32_t rev,rid,qs,qe,rs,re;
    int32_t bw1,extra_flag;
} zdrop_context_t ;

typedef struct {
    int32_t rs, re, qs, qe;
    int32_t rs0, rs1, qs1, re0, re1, qe1;
    int32_t rid,rev;
    int32_t llen, rlen, qlen;
    int32_t as1,cnt1;
    mm_reg1_t r2;
    zdrop_context_t *zdropctx;
} reg_context_t ;

typedef struct {
    const mm_idx_t *mi;
    void *newkm;
    mm_tbuf_t *b;
    mm128_t *a;
    mm_reg1_t *regs0, **regs;
    uint8_t *qseq0[2];//kfree after fpga comeback
    char **qseqs;
    int *n_regs, *qlens;
    int  n_regs0, rep_len, n_a, n_regs0_swdone, n_regs0_cap, dppos;
    int enable;
    reg_context_t *regctx;
} sw_context_t ;

/*
typedef struct {
    uint32_t swpos,swseq,enable;
} left_swkey_t ;
*/

extern thread_ctrl_t snd_ctrl, rcv_ctrl, x86_ctrl;
extern vsc_ring_buf_t *left_sw_pos[MAX_DP_THRD];

int add_req(int v);
int add_rsp(int v);
void set_reqrsp_zero();
void set_last_flag(int v);
int lastblk_comeback();

int add_req_sw(int v);
int add_rsp_sw(int v);
void set_reqrsp_zero_sw();
void set_last_flag_sw(int v);
int lastblk_comeback_sw();

void init_dpctx(void);
void init_swctx(void);
//void init_leftk(void);
void init_leftp(void);
void init_data_lock(void);
void deinit_data_lock(void);
void deinit_leftp(void);

void debug_context(void);
dp_to_ring_t * get_seed_data(void);
void chaindp_put_last(void);

//left_swkey_t *get_skdat_by_pos(uint32_t dpos, int tid);
dp_context_t *get_dpctx_by_pos(uint32_t dpos);
sw_context_t *get_swctx_by_pos(uint32_t dpos);
int get_free_swpos(void);
//int get_free_skpos(int tid);
void put_free_swpos(int pos);
void put_free_dppos(int pos);
//void put_free_skpos(int pos, int tid);

int sw_ring_submit(void *ptr, vsc_ring_buf_t *ring_buf);
sw_to_ring_t *sw_ring_init(void *km, const void *opt, uint16_t readnum, uint16_t tid, uint8_t swseq);


void save_dpctx(int dpos, const mm_idx_t *mi,int *n_regs,mm_reg1_t **regs,mm_tbuf_t *b,mm128_t *mva, /*uint64_t mvn,uint64_t mvm, uint64_t *mpos,*/ uint32_t hash, int n_segs, int qsum, int gapref, int gapqry, /*int replen, int nmpos,*/ const int *qlen, const char **seqs, const char *qname);

void save_swctx(int swpos, void *newkm, const mm_idx_t *mi, int *n_regs,mm_reg1_t **regs, mm_reg1_t *regs0, mm128_t *a, mm_tbuf_t *b, int *qlens, uint8_t **qseq0, char **qseqs, int n_regs0, int rep_len, int n_a, int dppos);

int save_seed_data(uint32_t ref, uint32_t qry, /*uint32_t seednum,*/ int flag, int occ, int qsum, uint16_t n_segs, /*mm128_t *seed,*/ void *km, const char *qname, const mm_idx_t *mi, mm128_v *mv);

uint32_t mm_align_skeleton_pre(void *km, const mm_mapopt_t *opt, const mm_idx_t *mi, uint8_t *qseq0[2], int qlen, const char *qstr, int *n_regs_, mm_reg1_t *regs, mm128_t *a, reg_context_t *regctx, sw_readhdr_t *read, uint32_t ctxpos, int regseq, int n_a);

#endif//_CONTEXT_DATA_H_
