#ifndef __SOFT_SW_H__
#define __SOFT_SW_H__

#include "minimap.h"
#include "ksw2.h"

#define CONTEXT_NUM     4096    //保存上下文数量
#define CHAIN_TASK_NUM     4096    //保存chain任务的数量
#define CHAIN_RESULT_NUM     4096    //保存chain结果的数量

typedef struct _sw_context_t sw_context_t;
typedef struct _chain_context_t chain_context_t;

typedef struct {
    //原始数据
    mm_reg1_t *regs0_ori;       //开辟空间保存原始的regs，用完后需要释放
    mm128_t *a_ori;             //开辟空间保存原始的a，用完后需要释放

    int n_regs0;
    mm_reg1_t *regs0;
    int32_t n_a;
    mm128_t *a;

    int i;  //表示第几条chain
   
    //align_regs context
    mm_tbuf_t *b;
    const mm_idx_t *mi;
    const mm_mapopt_t *opt;
    int qlen;
    char* seq;

    uint8_t* tseq;

    int rep_len;
    int *n_regs;        //存放结果
    mm_reg1_t **regs;   //存放结果

    int8_t mat[25];
    mm_reg1_t *r;

    int32_t rev;
    
    uint8_t *qseq0;     //结果处理完要释放

    long read_index;                    //read编号
    int chain_num;
    int chain_size;
    chain_context_t** chain_contexts;    //保存该read下每条chain的上下文
}context_t;

typedef struct _chain_context_t{
    context_t* context;

    int32_t i;
    int32_t rs;
    int32_t qs;
    int32_t qs0;
    
    int32_t rid;
    uint8_t* tseq;  //一个chain做完要释放
    uint8_t* qseq0[2]; //一个chain做完要释放
    
    int32_t re0;
    int32_t rs0;

    int sw_num;
    int sw_size;
    sw_context_t** sw_contexts;     //保存该chain下每条sw的上下文
}chain_context_t;

typedef struct {

    int qlen;
    const uint8_t* query;     //指向sw上下文的query
    int tlen;
    const uint8_t* target;    //指向sw上下文的target
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

typedef struct _sw_context_t{
    chain_context_t *chain_context;

    //保存做sw的现场
    int qlen;
    uint8_t* query;
    int tlen;
    uint8_t* target;
    int w;
    int zdrop;
    int end_bonus;
    int flag;
    int zdrop_flag;
    
    char pos_flag;  //0:left 1:mid 2:right
    
    int32_t i;
    int32_t cnt1;
    
    int32_t re;
    int32_t qe;
    int32_t qe0;
    
    int32_t rs;
    int32_t qs;
    int32_t qs0;
    
    int32_t as1;
}sw_context_t;

typedef struct {
    long read_id;
    int chain_id;

    int sw_num;
    int sw_size;
    sw_task_t** sw_tasks;
}chain_sw_task_t;

typedef struct {
    long read_id;
    int chain_id;

    int result_num;
    int result_size;
    ksw_extz_t** ezs;
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
chain_sw_task_t* create_chain_sw_task(long read_id, int chain_id);

//向这个chain的sw任务对象中添加一个sw的任务
void add_sw_task(chain_sw_task_t* chain_sw_task, sw_task_t* sw_task);

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

//创建一个read的上下文
context_t* create_context(long id);
//创建一个chain的上下文
chain_context_t* create_chain_context();
//向一个read上下文添加一个chain上下文
void add_chain_context(context_t* context, chain_context_t* chain_context);
//向一个chain上下文添加一个sw上下文
void add_sw_context(chain_context_t* chain_context, sw_context_t* sw_context);

//销毁一个chain及其所有sw的上下文
void destroy_chain_context(chain_context_t* chain_context);

int result_empty();

void init_result_array();
#endif //__SOFT_SW_H__
