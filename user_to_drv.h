#ifndef _USER_TO_DRIVER_H_
#define _USER_TO_DRIVER_H_

#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <pthread.h>

#include "netlink/list.h"

#include "kfifo.h"
#include "fpga_sim.h"

enum {
    SW_POS_LEFT = 0xab,
    SW_POS_MIDL,
    SW_POS_RIGH,
};

typedef struct
{
    uint16_t flag;
    uint16_t bw;
    uint16_t max_occ;
    uint16_t mid_occ;
    uint8_t  maxskip;
    uint8_t  mincnt;
    uint8_t  minscore;
    uint8_t  matchValue/*a*/,mismatchPenalty/*b*/,gapOpenPenalty/*q*/,gapExtendPenalty/*e*/,gapOpenPenalty2/*q2*/,gapExtendPenalty2/*e2*/;
    int8_t   end_bonus, zdrop;
} minimap2_global_t ;


typedef struct
{
    void *km;
    const mm_idx_t *mi;
    uint32_t magic;//for debug
    uint32_t total_size;//total seed size of factnum
    uint16_t tid;
    uint16_t factnum;
    uint8_t lastblk;//0=not,1=is last
    uint8_t type;
    chaindp_x86_t data[DP_BATCHSIZE];
} dp_to_ring_t;

typedef struct
{
    void *km;
    const void *opt;
    uint32_t magic;//for debug
    uint32_t total_size;//total memsize of factnum
    uint16_t tid;
    uint16_t factnum;//readnum
    uint8_t lastblk;//0=not,1=is last
    uint8_t type;//
    sw_readhdr_t data[DP_BATCHSIZE];
} sw_to_ring_t;


typedef struct
{
    void *qseq;
    uint32_t *S;//mi->S
    int64_t  offset;
    int32_t  qlen,tlen;
    int32_t  rid,n_seq,slen,st,en;
    int16_t  flag;
    int16_t  zdrop,zdrop_inv,bw;
    int16_t  end_bonus;
    int16_t  swpos;//one of left,mid,right
    struct nl_list_head entry;
} sw_node_t;



typedef struct
{
    int stop;//run or exit. 0=running
    vsc_ring_buf_t *dp_ring_buf;
    vsc_ring_buf_t *sw_ring_buf;
} thread_ctrl_t ;

typedef struct {
    int id;
    pthread_t tid;
    thread_ctrl_t* ring_buf;
}thread_param_t;

typedef void* (*CALLBACK)(void *);

thread_param_t *start_task_thread(int nthread, CALLBACK callback, void *param);
void* dp_task_sender(void *param);
void* sw_task_sender(void *param);
void* task_recver(void *param);
void* task_sender(void *param);
void* task_sender_sw(void *param);

//int dma_chaindp_prepare(uint16_t tid, uint8_t lat, uint16_t num, uint32_t size, void *dest, chaindp_x86_t *data);
//void dma_chaindp_submit();


#endif//_USER_TO_DRIVER_H_
