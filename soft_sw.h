#ifndef __SOFT_SW_H__
#define __SOFT_SW_H__

#include "minimap.h"
#include "ksw2.h"

#define CONTEXT_NUM     1024    //保存上下文数量
#define CHAIN_TASK_NUM     1024    //保存chain任务的数量
#define CHAIN_RESULT_NUM     1024    //保存chain结果的数量

typedef struct {
    
}context_t;


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
    context_t* context;
    int sw_num;
    int sw_size;
    sw_task_t** sw_tasks;
}chain_sw_task_t;

typedef struct {
    int result_num;
    int result_size;
    char pos_flag[1024];
    ksw_extz_t** ezs;
}sw_result_t;

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
chain_sw_task_t* create_chain_sw_task(context_t* context);

//向这个chain的sw任务对象中添加一个sw的任务
void add_sw_task(chain_sw_task_t* chain_sw_task, sw_task_t* sw_task);

//销毁一个chain的sw任务对象，包括它里面存储的所有sw的任务
void destroy_chain_sw_task(chain_sw_task_t* chain_sw_task);

#endif //__SOFT_SW_H__
