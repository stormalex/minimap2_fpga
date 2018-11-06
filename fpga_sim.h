#ifndef _FPGA_SIM_H_
#define _FPGA_SIM_H_

#include <stdint.h>
#include "minimap.h"


#define QUERY_NAME_MAXLEN   (252)
#define DP_BATCHSIZE        (127)
#define DP_QUEUE_NUM        (1024)
#define SW_QUEUE_NUM        DP_QUEUE_NUM
#define DP_CTX_NUM          (DP_QUEUE_NUM * DP_BATCHSIZE)
#define SW_CTX_NUM          DP_CTX_NUM
#define SW_LFT_NUM          DP_CTX_NUM
//#define SW_SUBMIT_MEMSIZE   (2^28)//256MB

#define MAX_DP_THRD         (DP_QUEUE_NUM/DP_BATCHSIZE)

#define UNAVAILABLE_ADDR    (0xdeaddeaddeaddead)
//#define REV_CIGAR_FPGA      (0xA0B1C2D3)
#define FREE_CIGAR_X86      (0xDEAD)

//to debug use API,not thread
//#define DEBUG_API
//#define DUMP_FILES
#define RUN_ON_X86_DP
//#define RUN_ON_X86_SW
#define RUN_ON_X86_OT

#define FPGA_ON    1       //1:hardware    0:soft

typedef enum
{
    GLOBAL_INIT = 0,
    CHAINDP,
    KSW,
    COLLECT_SEED,
} DATA_TYPE_E;



typedef struct {uint64_t x,y;} sim_mm128_t; 


typedef struct
{
    sim_mm128_t *seed;//seed addr @x86
    void *seed_km;//alloc seed km before chaindp
    uint32_t gap_ref;
    uint32_t gap_qry;
    uint32_t seednum;
    int32_t  mid_occ,optflag,qlen_sum;
    uint32_t pos;//context pos@x86 queue
    uint16_t n_segs;
    //mm128_v mv;
    char qname[QUERY_NAME_MAXLEN];
} chaindp_x86_t;


typedef struct
{
    //uint64_t mvn,mvm;
    //uint32_t mvoffset;//from seed[0]
    uint32_t splitflag;
    uint32_t gap_ref;
    uint32_t gap_qry;
    uint32_t seednum;
    int32_t  mid_occ,optflag,qlen_sum;
    uint32_t pos;//context pos@x86 queue
    uint16_t n_segs;
    char qname[QUERY_NAME_MAXLEN];
    uint16_t pad2B;
    uint32_t pad32B[8];
    sim_mm128_t seed[0];//seed addr @x86
} chaindp_to_fpga_t;


typedef struct
{
    uint32_t offset;
    uint32_t size;
} chaindp_sndsubhdr_t;


typedef struct 
{
    uint32_t magic;
    uint32_t size;
    uint16_t tid;//also queue id
    uint16_t num;
    uint8_t  type;
    uint8_t  lat;
    uint16_t pad2B;
    uint32_t pad48B[12];
    chaindp_sndsubhdr_t subhdr[0];
} chaindp_sndhdr_t;

typedef struct {
    uint32_t offset_kb;
    uint32_t offset_nu;
    uint16_t kb;
    uint16_t nu;
    uint32_t ctxpos;
} chaindp_rcvsubhdr_t ;


typedef struct {
    uint32_t size;
    uint16_t tid;//also queue id
    uint16_t num;
    uint8_t  type;
    uint8_t  pad;
    chaindp_rcvsubhdr_t subhdr[0];
} chaindp_rcvhdr_t ; 

/////////sw

typedef struct 
{
    int16_t qlen;
    int16_t tlen;
    //int32_t qoff;
    //int32_t toff;
    int16_t flag;
    int16_t zdrop;
    int16_t bw;
    int8_t  end_bonus;
    int8_t  pad1B;
    int16_t qlen_align;
    int16_t tlen_align;
    //uint8_t swpos;//one of left,mid,right
} sw_to_fpga_t;


typedef struct
{
    uint32_t offset;//first sw segment, from head always
    uint32_t size;//total ez size of this reg
    uint32_t midnum;//middle cnt in x86
    uint32_t regpos;
    //uint8_t  hasleft;//0 not has, non 0 has
    //uint8_t  hasright;
} sw_reghdr_t;


typedef struct 
{
    uint64_t head;//node list's head
    uint32_t ctxpos;//curr read's ctx position
    uint16_t regnum;//regnum for curr read'sw
    uint16_t longsw;
    sw_reghdr_t reg;
} sw_readhdr_t;


typedef struct 
{
    void *km;
    uint32_t magic;
    uint32_t size;
    uint16_t tid;//thread id
    uint16_t num;//fact readnum
    uint8_t  type;
    uint8_t  lat;//last or not
    uint16_t freecigar;
    uint16_t swcount;
    uint8_t pad8B[6];
    sw_readhdr_t data[DP_BATCHSIZE];
} sw_sndhdr_t;



//for test only

typedef struct {
    int32_t p,v;
} pv_t ;

typedef struct {
    uint32_t n_a,rep_len,n_mini_pos;
    uint32_t moffset,aoffset;
    uint32_t subsize,ctxpos;
    uint32_t pad36B[9];
    uint8_t  pvf[0];
} dptest_rcvsubhdr_t ;


typedef struct {
    void *   km;
    uint32_t magic;//from dpring magic
    uint32_t size;
    uint16_t tid;
    uint16_t num;
    uint8_t  type;
    uint8_t  lat;
    uint16_t pad2B;
    uint32_t pad48B[12];
    dptest_rcvsubhdr_t subhdr[0];
} dptest_rcvhdr_t ; 

/*
typedef struct {
    void *   km;
    uint32_t magic;//from swring magic
    uint32_t size;
    uint16_t tid;//threadid
    uint16_t num;
    uint8_t  type;
    uint8_t  lat;
    uint16_t freecigar;
    uint8_t pad8B[8];
    sw_readhdr_t data[DP_BATCHSIZE];
} swtest_rcvhdr_t ;
*/

#define swtest_rcvhdr_t sw_sndhdr_t


static uint32_t merge_to_midnum(uint32_t midnum, int hasleft, int hasright)
{
    uint32_t ret = 0;
    if (hasleft)  ret |= midnum|0x80000000;//set bit 31 to left
    if (hasright) ret |= midnum|0x40000000;//set bit 30 to right
    return ret;
}

static int get_left(uint32_t midnum)
{
    return midnum&0x80000000;//get bit 31
}

static int get_right(uint32_t midnum)
{
    return midnum&0x40000000;//get bit 30
}

static uint32_t clear_midnum(uint32_t midnum)
{
    return midnum&0x3FFFFFFF;//clear bit 31-30
}


void dumpfile(int memsize, void *dest, char *fpname, char *drname);
//sim fpga calc
void calc_chaindp_batch(void *addr, void *km, const mm_idx_t *mi);
int64_t fpgasim_chain_dp(int max_dist_x, int max_dist_y, int bw, int max_skip, int min_cnt, int min_sc, int is_cdna, int n_segs, int64_t n, mm128_t *a, int *n_u_, uint64_t **_u, int32_t *n_v_, int32_t **_v, void *km);
mm128_t *collect_seed_hits(void *km, const mm_mapopt_t *opt, 
                           int max_occ, const mm_idx_t *mi, 
                           const char *qname, const mm128_v *mv, 
                           int qlen, int64_t *n_a, int *rep_len,
                           int *n_mini_pos, uint64_t **mini_pos);
mm128_t *mm_chain_dp_fpos(int min_cnt, int min_sc, int32_t n_u, int64_t n, int32_t *f, pv_t *pv, mm128_t *a, int *n_u_, uint64_t **_u, void *km);

void calc_sw(void *dest, void *km, const void *topt);

void debug_kmalloc(void *km, void *addr,int size, char *func, int line);
void debug_kfree(void *km, void *addr,int size, char *func, int line);



//#define MY_FPGA_API

#ifdef MY_FPGA_API
//FPGA real API
int fpga_init(int flag, int is_cdna, int maxskip, int minsc, int m, int longthres, int longdiff, int8_t mat0, int8_t mat1);
void fpga_finalize();

void* fpga_get_retbuf(int* size, uint8_t type);
void fpga_release_retbuf(void* addr);

void* fpga_get_writebuf(int size);
void fpga_writebuf_submit(void* addr, int size, uint8_t type);

#endif

#endif//_FPGA_SIM_H_
