#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <pthread.h>

#include "minimap.h"
#include "mmpriv.h"
#include "ksw2.h"
#include "kthread.h"
#include "kvec.h"
#include "kalloc.h"
#include "sdust.h"
#include "mmpriv.h"
#include "bseq.h"
#include "khash.h"

#include "soft_sw.h"

int sw_stop_flag = 1;

//保存上下文的数组
context_t* context_array[CONTEXT_NUM];
pthread_mutex_t context_mutex = PTHREAD_MUTEX_INITIALIZER;

void init_context()
{
    memset(context_array, 0, sizeof(context_array));
    return;
}

int insert_context(context_t* context)
{
    int i = 0;
    pthread_mutex_lock(&context_mutex);
    for(i = 0; i < CONTEXT_NUM; i++) {
        if(context_array[i] == NULL) {
            pthread_mutex_unlock(&context_mutex);
            return i;
        }
    }
    pthread_mutex_unlock(&context_mutex);
    return -1;
}

context_t* delete_context(int i)
{
    context_t* tmp;
    pthread_mutex_lock(&context_mutex);
    if(context_array[i] != NULL) {
        tmp = context_array[i];
        pthread_mutex_unlock(&context_mutex);
        return tmp;
    }
    pthread_mutex_unlock(&context_mutex);
    return NULL;
}

//保存任务的数组
chain_sw_task_t* chain_tasks_array[CHAIN_TASK_NUM];
pthread_mutex_t tasks_mutex = PTHREAD_MUTEX_INITIALIZER;
static int tasks_insert_index = 0;
static int tasks_delete_index = 0;

void init_task_array()
{
    memset(chain_tasks_array, 0, sizeof(chain_tasks_array));
    tasks_insert_index = 0;
    tasks_delete_index = 0;
    return;
}

void send_task(chain_sw_task_t* chain_sw_tasks)
{
    pthread_mutex_lock(&tasks_mutex);
    if((tasks_insert_index+1)%CHAIN_TASK_NUM == tasks_delete_index) {
        fprintf(stderr, "have no position to save task\n");
        pthread_mutex_unlock(&tasks_mutex);
        return;
    }
    tasks_insert_index++;
    tasks_insert_index = tasks_insert_index%CHAIN_TASK_NUM;
    if(chain_tasks_array[tasks_insert_index] != NULL) {
        fprintf(stderr, "NULL, index:%d, delete=%d\n", tasks_insert_index, tasks_delete_index);
        pthread_mutex_unlock(&tasks_mutex);
        return;
    }
    chain_tasks_array[tasks_insert_index] = chain_sw_tasks;
    fprintf(stderr, "insert to %d\n", tasks_insert_index);
    pthread_mutex_unlock(&tasks_mutex);
    return;
}

chain_sw_task_t* get_task()
{
    chain_sw_task_t* tmp;
    pthread_mutex_lock(&tasks_mutex);
    if(tasks_delete_index == tasks_insert_index) {
        fprintf(stderr, "have no task\n");
        pthread_mutex_unlock(&tasks_mutex);
        return NULL;
    }
    
    tmp = chain_tasks_array[tasks_delete_index];
    if(tmp == NULL) {
        fprintf(stderr, "NULL, index:%d, delete=%d\n", tasks_insert_index, tasks_delete_index);
        pthread_mutex_unlock(&tasks_mutex);
        return NULL;
    }
    fprintf(stderr, "get from %d\n", tasks_delete_index);
    chain_tasks_array[tasks_delete_index] = NULL;
    tasks_delete_index++;
    tasks_delete_index = tasks_delete_index%CHAIN_TASK_NUM;
    
    pthread_mutex_unlock(&tasks_mutex);
    
    return tmp;
}


//保存结果的数组
sw_result_t* results_array[CHAIN_RESULT_NUM];
pthread_mutex_t results_mutex = PTHREAD_MUTEX_INITIALIZER;
static int results_insert_index = 0;
static int results_delete_index = 0;

void init_result_array()
{
    memset(results_array, 0, sizeof(results_array));
    results_insert_index = 0;
    results_delete_index = 0;
    return;
}

void send_result(sw_result_t* results)
{
    pthread_mutex_lock(&results_mutex);
    if((results_insert_index+1)%CHAIN_RESULT_NUM == results_delete_index) {
        fprintf(stderr, "have no position to save result\n");
        pthread_mutex_unlock(&results_mutex);
        return;
    }
    results_insert_index++;
    results_insert_index = results_insert_index%CHAIN_RESULT_NUM;
    if(results_array[results_insert_index] != NULL) {
        fprintf(stderr, "NULL, index:%d, delete=%d\n", results_insert_index, results_delete_index);
        pthread_mutex_unlock(&results_mutex);
        return;
    }
    results_array[results_insert_index] = results;
    fprintf(stderr, "insert to %d\n", results_insert_index);
    pthread_mutex_unlock(&results_mutex);
    return;
}

sw_result_t* get_result()
{
    sw_result_t* tmp;
    pthread_mutex_lock(&results_mutex);
    if(results_delete_index == results_insert_index) {
        fprintf(stderr, "have no result\n");
        pthread_mutex_unlock(&results_mutex);
        return NULL;
    }
    
    tmp = results_array[results_delete_index];
    if(tmp == NULL) {
        fprintf(stderr, "NULL, index:%d, delete=%d\n", results_insert_index, results_delete_index);
        pthread_mutex_unlock(&results_mutex);
        return NULL;
    }
    results_array[results_delete_index] = NULL;
    fprintf(stderr, "get from %d\n", results_delete_index);
    results_delete_index++;
    results_delete_index = results_delete_index%CHAIN_RESULT_NUM;
    
    pthread_mutex_unlock(&results_mutex);
    return tmp;
}

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
                        char position)
{
    sw_task_t* new_task = (sw_task_t*)malloc(sizeof(sw_task_t));
    if(new_task != NULL) {
        new_task->qlen = qlen;
        new_task->query = (uint8_t *)malloc(qlen + 1);
        new_task->qlen = tlen;
        new_task->target = (uint8_t *)malloc(qlen + 1);
        memcpy(new_task->query, query, qlen * sizeof(uint8_t));
        memcpy(new_task->target, target, tlen * sizeof(uint8_t));
        new_task->query[qlen] = '\0';
        new_task->target[tlen] = '\0';
        memcpy(new_task->mat, mat, sizeof(new_task->mat));
        new_task->q = q;
        new_task->e = e;
        new_task->q2 = q2;
        new_task->e2 = e2;
        new_task->w = w;
        new_task->zdrop = zdrop;
        new_task->end_bonus = end_bonus;
        new_task->flag = flag;
        new_task->position = position;
    }
    return new_task;
}

void destroy_sw_task(sw_task_t* sw_task)
{
    if(sw_task != NULL) {
        if(sw_task->query != NULL)
            free(sw_task->query);
        if(sw_task->target != NULL)
            free(sw_task->target);
        free(sw_task);
    }
    return;
}

chain_sw_task_t* create_chain_sw_task(context_t* context)
{
    chain_sw_task_t* new_chain_sw_task = (chain_sw_task_t *)malloc(sizeof(chain_sw_task_t));
    if(new_chain_sw_task != NULL) {
        new_chain_sw_task->context = context;
        new_chain_sw_task->sw_num = 0;
        new_chain_sw_task->sw_size = 0;
        new_chain_sw_task->sw_tasks = NULL;
    }
    return new_chain_sw_task;
}

void destroy_chain_sw_task(chain_sw_task_t* chain_sw_task)
{
    int i = 0;
    if(chain_sw_task) {
        for(i = 0; i < chain_sw_task->sw_num; i++) {
            destroy_sw_task(chain_sw_task->sw_tasks[i]);
        }
        free(chain_sw_task->sw_tasks);
        free(chain_sw_task);
    }
    return;
}

void add_sw_task(chain_sw_task_t* chain_sw_task, sw_task_t* sw_task)
{
    if(chain_sw_task != NULL && sw_task != NULL) {
        if(chain_sw_task->sw_size == 0)
            chain_sw_task->sw_size = 32;
        else if(chain_sw_task->sw_size == chain_sw_task->sw_num)
            chain_sw_task->sw_size = chain_sw_task->sw_size*2;
        fprintf(stderr, "chain_sw_task->sw_size=%d\n", chain_sw_task->sw_size);
        chain_sw_task->sw_tasks = (sw_task_t**)realloc(chain_sw_task->sw_tasks, chain_sw_task->sw_size * sizeof(sw_task_t*));
        chain_sw_task->sw_num++;
        chain_sw_task->sw_tasks[chain_sw_task->sw_num] = sw_task;
    }
    return;
}

sw_result_t* create_result()
{
    sw_result_t* new_result = (sw_result_t*)malloc(sizeof(sw_result_t));
    if(new_result != NULL) {
        new_result->result_num = 0;
        new_result->result_size = 32;
        memset(new_result->pos_flag, 0xff, sizeof(new_result->pos_flag));
        new_result->ezs = (ksw_extz_t**)malloc(new_result->result_size * sizeof(ksw_extz_t*));
    }
    return new_result;
}

void add_result(sw_result_t* result, ksw_extz_t* ez)
{
    if(result != NULL && ez != NULL) {
        if(result->result_size == 0)
            result->result_size = 32;
        else if(result->result_size == result->result_num)
            result->result_size = result->result_size*2;
        fprintf(stderr, "result->result_size=%d\n", result->result_size);
        result->ezs = (ksw_extz_t**)realloc(result->ezs, result->result_size * sizeof(ksw_extz_t*));
        result->result_num++;
        result->ezs[result->result_num] = ez;
    }
    return;
}

void destroy_results(sw_result_t* result)
{
    int i = 0;
    if(result) {
        for(i = 0; i < result->result_num; i++) {
            free(result->ezs[i]->cigar);
            free(result->ezs[i]);
        }
        free(result->ezs);
        free(result);
    }
    return;
}

void sw_thread(void* arg)
{
    int i = 0;
    init_task_array();
    init_result_array();
    
    while(sw_stop_flag)
    {
        chain_sw_task_t* chain_tasks = get_task();
        if(chain_tasks == NULL)
            continue;
        sw_result_t* result = create_result();
        for(i = 0; i < chain_tasks->sw_num; i++) {
            ksw_extz_t* ez = (ksw_extz_t*)malloc(sizeof(ksw_extz_t));
            
            sw_task_t* task = chain_tasks->sw_tasks[i];
            
            ksw_extd2_sse(NULL, task->qlen, task->query, task->tlen, task->target, 5, task->mat, task->q, task->e, task->q2, task->e2, task->w, task->zdrop, task->end_bonus, task->flag, ez);
            result->pos_flag[i] = task->position;
            add_result(result, ez);
        }
        send_result(result);
    }
}