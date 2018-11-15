#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <pthread.h>
#include <unistd.h>

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

static int sw_stop_flag = 1;

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
static int tasks_head = 0;
static int tasks_tail = 0;

void init_task_array()
{
    memset(chain_tasks_array, 0, sizeof(chain_tasks_array));
    tasks_head = 0;
    tasks_tail = 0;
    return;
}

int task_is_full()
{
    return ((tasks_tail + 1) % CHAIN_TASK_NUM == tasks_head);
}
int task_is_empty()
{
    return (tasks_head == tasks_tail);
}

int send_task(chain_sw_task_t* chain_sw_tasks)
{
    //sleep(1);
    pthread_mutex_lock(&tasks_mutex);
    if (task_is_full()) {
        pthread_mutex_unlock(&tasks_mutex);
        //fprintf(stderr, "task ringbuf is full\n");
        return -1;
    }
    
    chain_tasks_array[tasks_tail] = chain_sw_tasks;
    //fprintf(stderr, "task insert to %d\n", tasks_tail);
    tasks_tail = (tasks_tail + 1) % CHAIN_TASK_NUM;
        
    pthread_mutex_unlock(&tasks_mutex);
    return 0;
}

chain_sw_task_t* get_task()
{
    chain_sw_task_t* tmp;
    //sleep(1);
    pthread_mutex_lock(&tasks_mutex);
    if (task_is_empty()) {
        pthread_mutex_unlock(&tasks_mutex);
        //fprintf(stderr, "task ringbuf is empty\n");
        return NULL;
    }
    tmp = chain_tasks_array[tasks_head];
    //fprintf(stderr, "task take from %d\n", tasks_head);
    tasks_head = (tasks_head + 1) % CHAIN_TASK_NUM;
    pthread_mutex_unlock(&tasks_mutex);
    return tmp;
}


//保存结果的数组
sw_result_t* results_array[CHAIN_RESULT_NUM];
pthread_mutex_t results_mutex = PTHREAD_MUTEX_INITIALIZER;
static int results_head = 0;
static int results_tail = 0;

void init_result_array()
{
    memset(results_array, 0, sizeof(results_array));
    results_head = 0;
    results_tail = 0;
    return;
}

int result_is_full()
{
    return ((results_tail + 1) % CHAIN_RESULT_NUM == results_head);
}
int result_is_empty()
{
    return (results_head == results_tail);
}

int send_result(sw_result_t* results)
{
    pthread_mutex_lock(&results_mutex);
    if (result_is_full()) {
        pthread_mutex_unlock(&results_mutex);
        //fprintf(stderr, "result ringbuf is full\n");
        return -1;
    }
    
    results_array[results_tail] = results;
    //fprintf(stderr, "result insert to %d\n", results_tail);
    results_tail = (results_tail + 1) % CHAIN_RESULT_NUM;
        
    pthread_mutex_unlock(&results_mutex);
    return 0;
}

sw_result_t* get_result()
{
    sw_result_t* tmp;
    //sleep(1);
    pthread_mutex_lock(&results_mutex);
    if (result_is_empty()) {
        pthread_mutex_unlock(&results_mutex);
        //fprintf(stderr, "result ringbuf is empty\n");
        return NULL;
    }
    tmp = results_array[results_head];
    results_head = (results_head + 1) % CHAIN_RESULT_NUM;
    
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
        new_task->query = query;
        new_task->tlen = tlen;
        new_task->target = target;
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
        free(sw_task);
    }
    return;
}

chain_sw_task_t* create_chain_sw_task(long read_id, int chain_id)
{
    chain_sw_task_t* new_chain_sw_task = (chain_sw_task_t *)malloc(sizeof(chain_sw_task_t));
    if(new_chain_sw_task != NULL) {
        new_chain_sw_task->read_id = read_id;
        new_chain_sw_task->chain_id = chain_id;
        new_chain_sw_task->sw_num = 0;
        new_chain_sw_task->sw_size = 32;
        new_chain_sw_task->sw_tasks = (sw_task_t**)malloc(new_chain_sw_task->sw_size * sizeof(sw_task_t*));
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
        else if(chain_sw_task->sw_size == chain_sw_task->sw_num) {
            chain_sw_task->sw_size = chain_sw_task->sw_size*2;
            //fprintf(stderr, "chain_sw_task(%p)->sw_size=%d\n", chain_sw_task, chain_sw_task->sw_size);
            chain_sw_task->sw_tasks = (sw_task_t**)realloc(chain_sw_task->sw_tasks, chain_sw_task->sw_size * sizeof(sw_task_t*));
        }
        chain_sw_task->sw_tasks[chain_sw_task->sw_num] = sw_task;
        chain_sw_task->sw_num++;
    }
    return;
}

sw_result_t* create_result()
{
    sw_result_t* new_result = (sw_result_t*)malloc(sizeof(sw_result_t));
    if(new_result != NULL) {
        new_result->result_num = 0;
        new_result->result_size = 32;
        new_result->ezs = (ksw_extz_t**)malloc(new_result->result_size * sizeof(ksw_extz_t*));
    }
    return new_result;
}

void add_result(sw_result_t* result, ksw_extz_t* ez)
{
    if(result != NULL && ez != NULL) {
        if(result->result_size == 0)
            result->result_size = 32;
        else if(result->result_size == result->result_num) {
            result->result_size = result->result_size*2;
            //fprintf(stderr, "result->result_size=%d\n", result->result_size);
            result->ezs = (ksw_extz_t**)realloc(result->ezs, result->result_size * sizeof(ksw_extz_t*));
        }
        result->ezs[result->result_num] = ez;
        result->result_num++;
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

context_t* create_context(long id)
{
    context_t* context = (context_t *)malloc(sizeof(context_t));
    if(context != NULL) {
        context->read_index = id;
        context->chain_num = 0;
        context->chain_size = 16;
        context->chain_contexts = (chain_context_t**)malloc(context->chain_size * sizeof(chain_context_t*));
    }
    return context;
}

void add_chain_context(context_t* context, chain_context_t* chain_context)
{
    if(context != NULL && chain_context != NULL) {
        if(context->chain_size == 0)
            context->chain_size = 32;
        else if(context->chain_size == context->chain_num) {
            context->chain_size = context->chain_size*2;
            context->chain_contexts = (chain_context_t**)realloc(context->chain_contexts, context->chain_size * sizeof(chain_context_t*));
        }
        context->chain_contexts[context->chain_num] = chain_context;
        context->chain_num++;
        chain_context->context = context;   //指向了这个read的上下文
    }
}

chain_context_t* create_chain_context()
{
    chain_context_t* chain_context = (chain_context_t *)malloc(sizeof(chain_context_t));
    if(chain_context != NULL) {
        chain_context->sw_num = 0;
        chain_context->sw_size = 16;
        chain_context->sw_contexts = (sw_context_t**)malloc(chain_context->sw_size * sizeof(sw_context_t*));
    }
    return chain_context;
}

void add_sw_context(chain_context_t* chain_context, sw_context_t* sw_context)
{
    if(chain_context != NULL && sw_context != NULL) {
        if(chain_context->sw_size == 0)
            chain_context->sw_size = 32;
        else if(chain_context->sw_size == chain_context->sw_num) {
            chain_context->sw_size = chain_context->sw_size*2;
            chain_context->sw_contexts = (sw_context_t**)realloc(chain_context->sw_contexts, chain_context->sw_size * sizeof(sw_context_t*));
        }
        chain_context->sw_contexts[chain_context->sw_num] = sw_context;
        chain_context->sw_num++;
        sw_context->chain_context = chain_context;   //指向了这个chain的上下文
    }
}

void stop_sw_thread()
{
    sw_stop_flag = 0;
}

void* sw_thread(void* arg)
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
        //fprintf(stderr, "sw task num=%d\n", chain_tasks->sw_num);
        for(i = 0; i < chain_tasks->sw_num; i++) {
            ksw_extz_t* ez = (ksw_extz_t*)malloc(sizeof(ksw_extz_t));
            memset(ez, 0, sizeof(ksw_extz_t));
            sw_task_t* task = chain_tasks->sw_tasks[i];
            
            ksw_extd2_sse(NULL, task->qlen, task->query, task->tlen, task->target, 5, task->mat, task->q, task->e, task->q2, task->e2, task->w, task->zdrop, task->end_bonus, task->flag, ez);

            add_result(result, ez);
        }
        result->read_id = chain_tasks->read_id;
        result->chain_id = chain_tasks->chain_id;

        destroy_chain_sw_task(chain_tasks);
        while(send_result(result));
    }
    return NULL;
}