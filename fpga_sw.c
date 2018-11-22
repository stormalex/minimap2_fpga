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

#include "fpga_sw.h"
#include "soft_sw.h"
#include "list.h"

static int send_thread_flag = 1;
static chain_sw_task_t* task_array[SEND_ARRAY_SIZE];

void stop_fpga_thread()
{
    send_thread_flag = 0;
}

void* send_task_thread(void* arg)
{
    int chain_num = 0;
    int data_size = 0;
    char* fpga_buf = NULL;
    char* p = NULL;
    fpga_task_t* head = NULL;
    fpga_task_id_t* chain_head = NULL;
    fpga_sw_task* sw_task = NULL;
    int wait_number = 0;
    int send_flag = 0;
    unsigned long long tag = 0;

    while(send_thread_flag)
    {
        chain_sw_task_t* chain_tasks = get_task();
        if(chain_tasks == NULL) {
            wait_number++;
            if(wait_number >= 1000) {
                wait_number = 0;
                send_flag = 1;
            }
        }
        else {
            wait_number = 0;
            task_array[chain_num++] = chain_tasks;
            data_size += chain_tasks->data_size;
            if(chain_num == SEND_ARRAY_SIZE) {
                send_flag = 1;
            }
        }

        if(send_flag == 1) {      //队列已满，或者发送条件满足
            send_flag = 0;

            if(data_size > 4 * 1024 * 1024) {
                fprintf(stderr, "data_size too large:%d\n", data_size);
                exit(0);
            }
            fpga_buf = (char*)fpga_get_writebuf(data_size, BUF_TYPE_SW);
            if(fpga_buf == NULL) {
                fprintf(stderr, "fpga_get_writebuf sw error\n");
                exit(0);
            }

            head = (fpga_task_t*)fpga_buf;        //指向buff首地址
            head->sw_num = 0;
            chain_head = (fpga_task_id_t*)((char*)fpga_buf + sizeof(fpga_task_t));
            sw_task = (fpga_sw_task*)((char*)fpga_buf + 4096);

            //设置头信息
            head->tag = tag++;
            head->chain_task_num = chain_num;

            int chain_index = 0;
            int sw_index = 0;
            for(chain_index = 0; chain_index < chain_num; chain_index++) {
                if((char*)chain_head >= (char*)(fpga_buf + 4096)) {
                    fprintf(stderr, "ERROR: too much chain, chain num=%d\n", chain_num);
                    exit(0);
                }

                //设置chain的信息
                chain_head->read_id = task_array[chain_index]->read_id;
                chain_head->chain_id = task_array[chain_index]->chain_id;
                chain_head->sw_num = task_array[chain_index]->sw_num;

                //设置sw任务的数据
                head->sw_num += chain_head->sw_num;     //统计一次调用sw任务数
                for(sw_index = 0; sw_index < chain_head->sw_num; sw_index++) {
                    sw_task->qlen = task_array[chain_index]->sw_tasks[sw_index]->qlen;
                    sw_task->tlen = task_array[chain_index]->sw_tasks[sw_index]->tlen;
                    sw_task->flag = task_array[chain_index]->sw_tasks[sw_index]->flag;
                    sw_task->zdrop = task_array[chain_index]->sw_tasks[sw_index]->zdrop;
                    sw_task->bw = task_array[chain_index]->sw_tasks[sw_index]->w;
                    sw_task->end_bonus = task_array[chain_index]->sw_tasks[sw_index]->end_bonus;
                    char *qseq = (char*)sw_task + sizeof(fpga_sw_task);    //设置qseq的地址
                    memcpy(qseq, task_array[chain_index]->sw_tasks[sw_index]->query, sw_task->qlen);
                    char *tseq = qseq + ADDR_ALIGN(sw_task->qlen, 16);
                    memcpy(tseq, task_array[chain_index]->sw_tasks[sw_index]->target, sw_task->tlen);

                    p = (char*)sw_task;
                    sw_task = (fpga_sw_task*)(p + (sizeof(fpga_sw_task) + ADDR_ALIGN(sw_task->qlen, 16) + ADDR_ALIGN(sw_task->tlen, 16)));    //指向下一个sw任务的地址
                }

                p = (char*)chain_head;
                chain_head = (fpga_task_id_t*)(p + sizeof(fpga_task_id_t));   //指向下一个chain的头
                destroy_chain_sw_task(task_array[chain_index]);     //销毁chain任务的数据
                task_array[chain_index] = NULL;
            }

            if((data_size + 4096) != ((char*)sw_task - (char*)fpga_buf)) {
                fprintf(stderr, "fpga_buf=%p, sw_task=%p, data_size=0x%x\n", fpga_buf, sw_task, data_size);
                exit(0);
            }
            fpga_writebuf_submit(fpga_buf, data_size + 4096, TYPE_SW);

            chain_num = 0;
            data_size = 0;
        }
    }
    return NULL;
}

void* recv_task_thread(void *arg)
{
    void* fpga_buf = NULL;
    int fpga_len;
    fpga_task_t* head = NULL;
    fpga_task_id_t* chain_head = NULL;
    int chain_index;
    int sw_index;
    char* p = NULL;
    ksw_extz_t* ez_addr = NULL;
    uint32_t *cigar_addr = NULL;

    while(send_thread_flag) {
        fpga_buf = fpga_get_retbuf(&fpga_len, RET_TYPE_SW);
        if(fpga_buf == NULL) {
            fprintf(stderr,"ERROR:fpga_len:%d\n",fpga_len);
            exit(1);
        }

        head = (fpga_task_t*)fpga_buf;
        int chain_num = head->chain_task_num;
        chain_head = (fpga_task_id_t*)((char*)fpga_buf + sizeof(fpga_task_t));
        ez_addr = (ksw_extz_t*)((char*)fpga_buf + 4096);
        cigar_addr = (uint32_t *)(ez_addr + head->sw_num);

        //处理每一个chain
        for(chain_index = 0; chain_index < chain_num; chain_index++) {
            sw_result_t* result = create_result();
            result->read_id = chain_head[chain_index].read_id;
            result->chain_id = chain_head[chain_index].chain_id;
            for(sw_index = 0; sw_index < chain_head[chain_index].sw_num; sw_index++) {

                //生成一个ez结果
                ksw_extz_t* ez = (ksw_extz_t*)malloc(sizeof(ksw_extz_t));
                memset(ez, 0, sizeof(ksw_extz_t));
                *ez = *ez_addr;
                ez->cigar = (uint32_t *)malloc(ez->n_cigar * sizeof(uint32_t));
                memcpy(ez->cigar, cigar_addr, ez->n_cigar * sizeof(uint32_t));
                ez->m_cigar = ez->n_cigar;

                ez_addr += 1;   //移动ez到下一个

                p = (char*)cigar_addr;
                cigar_addr = (uint32_t *)(p + ADDR_ALIGN((uint64_t)cigar_addr, 16));  //移动到下一个cigar开头

                add_result(result, ez); //将ez结果加入到chain结果中
            }
            
            while(send_result(result)); //发送给结果处理线程

            p = (char*)chain_head;
            chain_head = (fpga_task_id_t*)(p + sizeof(fpga_task_id_t));   //指向下一个chain的头
        }
        fpga_release_retbuf(fpga_buf);
    }
    return NULL;
}