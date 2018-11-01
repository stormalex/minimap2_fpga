#ifndef __FPGA_H__
#define __FPGA_H__

#ifdef __cplusplus

extern "C"{

#endif

#define BLOCK       0
#define NOBLOCK     1

#define	FPGA_ERR_OK         (0)
#define FPGA_ERR_NOMEM      (-1)

#define TYPE_SW 0
#define TYPE_CD 1
#define TYPE_CS 2

typedef enum {
    BUF_TYPE_SW = 0,
    BUF_TYPE_CD = 1,
    BUF_TYPE_CS = 3,
}BUF_TYPE;

typedef enum {
    RET_TYPE_SW = 0,
    RET_TYPE_CD = 1,
    RET_TYPE_CS = 3,
}RET_TYPE;

int fpga_init(int flag);
void fpga_finalize();

void* fpga_get_retbuf(int* len, RET_TYPE type);         //获取一个结果缓冲区，数据可用
int fpga_release_retbuf(void* addr);               //释放一个结果缓冲区

void* fpga_get_writebuf(unsigned long size, BUF_TYPE type);                   //获得一个写缓冲区
int fpga_writebuf_submit(void* addr, unsigned int size, unsigned int type);    //向FPGA提交写缓冲区数据

int fpga_init_sw(void* parameters);

int fpga_send_sw(int id, int qlen, char* q, int tlen, char* t);


void fpga_test(void);

#ifdef __cplusplus

};

#endif

#endif //__FPGA_H__