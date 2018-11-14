#ifndef __SOFT_SW_H__
#define __SOFT_SW_H__

#include "minimap.h"
#include "ksw2.h"

#define CONTEXT_NUM     1024    //保存上下文数量
#define CHAIN_TASK_NUM     1024    //保存chain任务的数量
#define CHAIN_RESULT_NUM     1024    //保存chain结果的数量

typedef struct {
    //原始数据
    int n_regs0;
    mm_reg1_t *regs0;       //开辟空间保存原始的regs，用完后需要释放
    int32_t n_a;
    mm128_t *a;             //开辟空间保存原始的a，用完后需要释放
    
    int i;  //表示第几条chain
   
    //align_regs context
    mm_tbuf_t *b;
    const mm_mapopt_t *opt;
    int qlen;
    const char* seq;
    long read_index;

    uint8_t* tseq;

    int rep_len;
    int *n_regs;
    mm_reg1_t **regs;

    int32_t rs;
    int32_t qs;
    int32_t qs0;
    int8_t mat[25];
    int32_t extra_flag;
    mm_reg1_t *r;
    mm_reg1_t *r2;
    int32_t as1;

    int32_t rev;
}context_t;

typedef struct {
    int32_t i;
}chain_context_t;

typedef struct {
    int qlen;
    uint8_t* query;
    int tlen;
    uint8_t* target;
    int8_t mat[25];
    int8_t q;
    int8_t e;
    int8_t q2;
    int8_t e2;
    int w;
    int zdrop;
    int end_bonus;
    int flag;
    char position;   //0-left 1-middle 2-right
}sw_task_t;

typedef struct {
    char pos_flag;
    int qlen;
    int tlen;
    int w;
    
    uint8_t* qseq;
    uint8_t* tseq;
    
    int32_t i;
    int32_t cnt1;
    
    int32_t re;
    int32_t qe;
    int32_t qe0;
}sw_context_t;

typedef struct {
    context_t* context;
    int sw_num;
    int sw_size;
    sw_task_t** sw_tasks;
    sw_context_t** sw_contexts;
    chain_context_t* chain_context;
}chain_sw_task_t;

typedef struct {
    context_t* context;
    int result_num;
    int result_size;
    ksw_extz_t** ezs;
    sw_context_t** sw_contexts;
    chain_context_t* chain_context;
}sw_result_t;


int send_task(chain_sw_task_t* chain_sw_tasks);
sw_result_t* get_result();

//创建一个sw任务
sw_task_t* create_sw_task(int qlen,
                        const uint8_t* query,
                        int tlen,
                        const uint8_t* target,
                        const int8_t* mat,
                        int q,
                        int e,
                        int q2,
                        int e2,
                        int w,
                        int zdrop,
                        int end_bonus,
                        int flag,
                        char position);


//销毁一个sw任务
void destroy_sw_task(sw_task_t* sw_task);

//创建一个chain的sw任务对象
chain_sw_task_t* create_chain_sw_task(context_t* context, chain_context_t* chain_context);

//向这个chain的sw任务对象中添加一个sw的任务
void add_sw_task(chain_sw_task_t* chain_sw_task, sw_task_t* sw_task, sw_context_t* sw_context);

//销毁一个chain的sw任务对象，包括它里面存储的所有sw的任务
void destroy_chain_sw_task(chain_sw_task_t* chain_sw_task);

//创建一个chain的结果
sw_result_t* create_result();

//向result中加入一个ez结果
void add_result(sw_result_t* result, ksw_extz_t* ez);

//销毁一个chain的所有result
void destroy_results(sw_result_t* result);

//sw处理的线程
void* sw_thread(void* arg);

//停止sw处理线程
void stop_sw_thread();
#endif //__SOFT_SW_H__
