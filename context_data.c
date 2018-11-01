#include <unistd.h>
#include <pthread.h>

#include "minimap.h"
#include "mmpriv.h"
#include "ksw2.h"

#include "kfifo.h"
#include "fpga_sim.h"
#include "context_data.h"
#include "user_to_drv.h"



thread_ctrl_t snd_ctrl, rcv_ctrl, x86_ctrl;

static dp_to_ring_t *seed_data = 0;
static pthread_spinlock_t seed_data_lock;

static vsc_ring_buf_t *dp_free_ctx  = 0;
static dp_context_t dp_ctx[DP_CTX_NUM];

static vsc_ring_buf_t *sw_free_ctx = 0;
static sw_context_t sw_ctx[SW_CTX_NUM];

//static vsc_ring_buf_t *left_sw[MAX_DP_THRD] = {0};
//static left_swkey_t lf_key[MAX_DP_THRD][SW_CTX_NUM]={{0,0,0}};

vsc_ring_buf_t *left_sw_pos[MAX_DP_THRD] = {0};


static int reqnum = 0;
static int rspnum = 0;
static int latblk = 0;//0 not last, 1=last,means submit by timeout

static int sw_reqnum = 0;
static int sw_rspnum = 0;
static int sw_latblk = 0;//0 not last, 1=last,means submit by timeout


void set_reqrsp_zero()
{
    __sync_lock_test_and_set(&reqnum, 0);
    __sync_lock_test_and_set(&rspnum, 0);
}

void set_last_flag(int v)
{
    __sync_lock_test_and_set(&latblk, v);
}

int lastblk_comeback()
{
    return __sync_add_and_fetch(&latblk, 0) == IS_LAST;
}


int add_req(int v)
{
    return __sync_add_and_fetch(&reqnum, v);
}

int add_rsp(int v)
{
    return __sync_add_and_fetch(&rspnum, v);
}


void set_reqrsp_zero_sw()
{
    __sync_lock_test_and_set(&sw_reqnum, 0);
    __sync_lock_test_and_set(&sw_rspnum,0);
}

void set_last_flag_sw(int v)
{
    __sync_lock_test_and_set(&sw_latblk, v);
}

int lastblk_comeback_sw()
{
    return __sync_add_and_fetch(&sw_latblk, 0) == IS_LAST;
}


int add_req_sw(int v)
{
    return __sync_add_and_fetch(&sw_reqnum, v);
}

int add_rsp_sw(int v)
{
    return __sync_add_and_fetch(&sw_rspnum, v);
}


void init_context_ring(vsc_ring_buf_t **ring, int qnum)
{
    *ring = ring_alloc(qnum);
    for (int k=0; k<qnum; k++)
    {
        put_to_ring_buf(*ring, (void*)k);
    }
}


void init_dpctx(void)
{
    init_context_ring(&dp_free_ctx, DP_CTX_NUM);
}

void init_swctx(void)
{
    init_context_ring(&sw_free_ctx, SW_CTX_NUM);
}
/*
void init_leftk(void)
{
    for (int k=0; k<MAX_DP_THRD; k++) 
        init_context_ring(&left_sw[k], SW_CTX_NUM);
}
*/
void init_leftp(void)
{
    for (int k=0; k<MAX_DP_THRD; k++)
        left_sw_pos[k] = ring_alloc(SW_LFT_NUM);
}

void deinit_leftp(void)
{
    for (int k=0; k<MAX_DP_THRD; k++)
        ring_free(left_sw_pos[k]);
}

dp_to_ring_t * get_seed_data(void)
{
    return seed_data;
}

int get_free_dppos(void)
{
    int *pos = (int*)0xdead;

    int get = 0;
    do {
        get = get_from_ring_buf(dp_free_ctx, (void**)&pos);
    } while(get < 0);

    return (int)pos;
}

int get_free_swpos(void)
{
    int *pos = (int*)0xdead;

    int get = 0;
    do {
        get = get_from_ring_buf(sw_free_ctx, (void**)&pos);
    } while(get < 0);

    return (int)pos;
}
/*
int get_free_skpos(int tid)
{
    int *pos = (int*)0xdead;

    int get = 0;
    do {
        get = get_from_ring_buf_st(left_sw[tid], (void**)&pos);
    } while(get < 0);

    return (int)pos;
}
*/
void put_free_dppos(int pos)
{
    put_to_ring_buf(dp_free_ctx, (void *)pos);
}

void put_free_swpos(int pos)
{
    put_to_ring_buf(sw_free_ctx, (void *)pos);
}
/*
void put_free_skpos(int pos, int tid)
{
    put_to_ring_buf_st(left_sw[tid], (void *)pos);
}
*/
dp_context_t *get_dpctx_by_pos(uint32_t dpos)
{
    return &dp_ctx[dpos];
}

sw_context_t *get_swctx_by_pos(uint32_t dpos)
{
    return &sw_ctx[dpos];
}

/*
left_swkey_t *get_skdat_by_pos(uint32_t dpos, int tid)
{
    return &lf_key[tid][dpos];
}
*/

void debug_context(void)
{
    for (int k=0; k<DP_CTX_NUM; k++) {
        dp_context_t *curr = &dp_ctx[k];
        if (curr->b && curr->b->km) fprintf(stderr,"ERR:...ctxpos %d km is %p\n",k,curr->b->km);
    }
}


void save_dpctx(int dpos, const mm_idx_t *mi,int *n_regs,mm_reg1_t **regs,mm_tbuf_t *b,mm128_t *mva, /*uint64_t mvn,uint64_t mvm, uint64_t *mpos,*/ uint32_t hash, int n_segs, int qsum, int gapref, int gapqry,/*int replen, int nmpos,*/ const int *qlen, const char **seqs, const char *qname)
{
    dp_context_t *curr = &dp_ctx[dpos];
    curr->mi = mi;
    curr->n_regs = n_regs;
    curr->regs = regs;
    curr->b = b;
    curr->mva = mva;
    //curr->mvn = mvn;
    //curr->mvm = mvm;
    //curr->mini_pos = mpos;
    curr->hash = hash;
    curr->n_segs = n_segs;
    curr->qlen_sum = qsum;
    curr->gap_ref = gapref;
    curr->gap_qry = gapqry;
    //curr->rep_len = replen;
    //curr->n_mini_pos = nmpos;

    memcpy(curr->qlen, qlen, MM_MAX_SEG*sizeof(int));
    memcpy(curr->qseqs, seqs, MM_MAX_SEG*sizeof(char *));
    memset(curr->qname, 0, MM_MAX_SEG);
    memcpy(curr->qname, qname, strnlen(qname,MM_MAX_SEG-1));
}

void chaindp_put_last(void)
{
    int stop = 0;

    do {
    //TODO: at the time, all worker thread are exit, lock is not necessary?
        pthread_spin_lock(&seed_data_lock);

        if (seed_data && seed_data->factnum) {//if is not necessary?

            //update memory, after then put to ring
            seed_data->lastblk = IS_LAST;
            add_req(seed_data->factnum);
            int idx = -1;
            do {
                idx = put_to_ring_buf(snd_ctrl.dp_ring_buf, seed_data);
            } while(idx < 0);

            stop = 1;
            seed_data = 0;//init again, or next trunk size will use this value and access error memory, cause segfault
        }//if
        else fprintf(stderr, "BUG: factnum or seed_data is null\n");

        pthread_spin_unlock(&seed_data_lock);
    }while(!stop);

}


sw_to_ring_t *sw_ring_init(void *km, const void *opt, uint16_t readnum, uint16_t tid, uint8_t swseq)
{
    //kmalloc sw_to_ring_t, no lock
    int memsize = sizeof(sw_to_ring_t);
    sw_to_ring_t *swring = (sw_to_ring_t *)kmalloc(km, memsize);
    //init
    swring->km  = km;
    swring->opt = opt;
    swring->tid = tid;
    swring->factnum = readnum;
    swring->lastblk = swseq;
    swring->total_size = 0;

    return swring;
}

int sw_ring_submit(void *ptr, vsc_ring_buf_t *ring_buf)
{
    int idx = -1;
    do {
        idx = put_to_ring_buf(ring_buf, ptr);
    } while(idx < 0);

    add_req_sw(1);
    return idx;
}

int get_seed_pos(void *km, const mm_idx_t *mi, uint32_t na)
{
    //pthread_spin_lock(&seed_data_lock);

    int last = -1;
    int newblock = 0;

    if (seed_data && 
      (DP_BATCHSIZE == seed_data->factnum)) {

        int idx = -1;
        do {
            idx = put_to_ring_buf(snd_ctrl.dp_ring_buf, seed_data);
        } while(idx < 0);

        newblock = 1;
        add_req(DP_BATCHSIZE);
        //fprintf(stderr, "put to ring, total %u, addr %p\n",add_req(0), seed_data);
    }

    if (!seed_data || (1 == newblock)) {
        seed_data = (dp_to_ring_t *)kmalloc(km, sizeof(dp_to_ring_t));
        seed_data->km  = km;
        seed_data->mi  = mi;
        seed_data->tid = 123;
        seed_data->total_size = 0;
        seed_data->factnum = 0;
        seed_data->lastblk = NT_LAST;

        struct timespec tv;
        clock_gettime(CLOCK_MONOTONIC, &tv);
        seed_data->magic = tv.tv_nsec;
    }

    seed_data->total_size += ALIGN_BYTE_N(64, na*sizeof(mm128_t));//64B align
    last = seed_data->factnum++;

    //pthread_spin_unlock(&seed_data_lock);

    return last;
}



int save_seed_data(uint32_t ref, uint32_t qry, /*uint32_t seednum,*/ int flag, int occ, int qsum, uint16_t n_segs,/*mm128_t *seed,*/ void *km, const char *qname, const mm_idx_t *mi, mm128_v *mv)
{
    int ctxpos = get_free_dppos();

    pthread_spin_lock(&seed_data_lock);

    int seedpos = get_seed_pos(km, mi, mv->n);

    chaindp_x86_t *dest = &seed_data->data[seedpos];
    dest->gap_ref = ref;
    dest->gap_qry = qry;
    dest->seednum = mv->n;//seednum;
    dest->mid_occ = occ;
    dest->optflag = flag;
    dest->qlen_sum= qsum;
    dest->n_segs  = n_segs;
    dest->pos  = ctxpos;
    dest->seed = (sim_mm128_t*)mv->a;//seed;
    //dest->seed_km = km;
    //dest->mv.n = mv->n;
    //dest->mv.m = mv->m;
    //dest->mv.a = mv->a;
    memset(dest->qname, 0, sizeof(dest->qname));
    memcpy(dest->qname, qname, strnlen(qname, QUERY_NAME_MAXLEN-1));

    pthread_spin_unlock(&seed_data_lock);

    return ctxpos;
}


void save_swctx(int swpos, void *newkm, const mm_idx_t *mi, int *n_regs,mm_reg1_t **regs, mm_reg1_t *regs0, mm128_t *a, mm_tbuf_t *b, int *qlens, uint8_t **qseq0, char **qseqs, int n_regs0, int rep_len, int n_a, int dppos)
{
    sw_context_t *curr = &sw_ctx[swpos];
    curr->newkm = newkm;
    curr->mi = mi;
    curr->a  = a;
    curr->b  = b;
    curr->regs0   = regs0;
    curr->regs    = regs;
    curr->n_regs  = n_regs;
    curr->qlens   = qlens;
    curr->n_regs0 = n_regs0;
    curr->rep_len = rep_len;
    curr->n_a     = n_a; 
    curr->qseqs   = qseqs;
    curr->qseq0[0] = qseq0[0];
    curr->qseq0[1] = qseq0[1];
    curr->n_regs0_swdone = 0;
    curr->dppos = dppos;
}


void init_data_lock(void)
{
    pthread_spin_init(&seed_data_lock, PTHREAD_PROCESS_PRIVATE);
}

void deinit_data_lock(void)
{
    pthread_spin_destroy(&seed_data_lock);
}


