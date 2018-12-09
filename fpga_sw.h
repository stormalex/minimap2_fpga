#ifndef __FPGA_SW_H__
#define __FPGA_SW_H__

#include "fpga.h"

#define DUMP_FILE   0   //0:fpga  1:soft

#define SEND_ARRAY_SIZE     128

#define ADDR_ALIGN(addr, align)   (((addr)+align-1)&~( align-1))

#define CHECK_ID    0x1234ABCD

typedef struct {
    unsigned long long tag;     //包的标记，可以用来标记唯一一个包
    int check_id;
    int chain_task_num;         //标记任务中有几个chain
    int sw_num;                 //标记任务中有几个sw任务
}fpga_task_t;

typedef struct {
    unsigned long read_id;      //标记该任务的read编号
    unsigned int chain_id;      //标记该任务的chain编号
    unsigned int sw_num;        //标记该chain下有几个sw任务
}fpga_task_id_t;

typedef struct {
    unsigned short qlen;
    unsigned short tlen;
    unsigned short flag;
    unsigned short zdrop;
    unsigned short bw;
    unsigned char end_bonus;
    unsigned char reserve1;
    unsigned char reserve2[4];
}fpga_sw_task;

typedef struct {
    void* data;
    unsigned int size;
}task_t;
typedef struct {
    void* data;
    unsigned int size;
}result_t;
void init_fpga_task_array();
void stop_fpga_send_thread();
void* send_task_thread(void* arg);
void stop_fpga_thread();
int send_fpga_task(task_t task);


void init_fpga_result_array();
void stop_fpga_recv_thread();
int send_fpga_result(result_t result);
int get_fpga_result(result_t* result);
void* recv_task_thread(void* arg);
#endif //__FPGA_SW_H__
