#ifndef _VSC_KFIFO_H_
#define _VSC_KFIFO_H_

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <pthread.h>


typedef struct {
  pthread_spinlock_t lock;          /* protects concurrent modifications */
  uint64_t ridx;
  uint64_t widx;
  uint64_t size;
  uint32_t max_num;
  uint32_t num;
  uint64_t start_time;
  uint32_t total_num;
  uint64_t buff[0];
} vsc_ring_buf_t ;

#define likely(x)   __builtin_expect(!!(x), 1)
#define unlikely(x) __builtin_expect(!!(x), 0)


#define VSCMIN(x,y) ({ \
  typeof(x) _x = (x);     \
  typeof(y) _y = (y);     \
  (void) (&_x == &_y);            \
  _x < _y ? _x : _y; })

#define VSCMAX(x,y) ({ \
  typeof(x) _x = (x);     \
  typeof(y) _y = (y);     \
  (void) (&_x == &_y);            \
  _x > _y ? _x : _y; })


//#define RING_BUF_CELL_SIZE 1//sizeof(uint64_t)
static const uint64_t RING_BUF_CELL_SIZE = 1;

#define ALIGN_BYTE_N(n, v) ((v+n-1)&(~(n-1)))
#define is_power_of_2(x) ((x) != 0 && (((x) & ((x) - 1)) == 0))


static uint32_t get_power_2(uint32_t i){
  i|=i>>1;
  i|=i>>2;
  i|=i>>4;
  i|=i>>8;
  i|=i>>16;
  i+=1;
  i=i>>1;
  return i;
}

static uint32_t get_power_2_big(uint32_t i){
  i|=i>>1;
  i|=i>>2;
  i|=i>>4;
  i|=i>>8;
  i|=i>>16;
  i+=1;
  //i=i>>1;
  return i;
}


static void debug_ring(vsc_ring_buf_t *ring, int size)
{
  for (int k=0; k<size; k++)
    {
      fprintf(stderr,"buf[%d]=%u/0x%x\n",k,(int)ring->buff[k],(int)ring->buff[k]);
    }
}

static int put_to_ring_buf_st(vsc_ring_buf_t *ring, void *dataptr)
{
  int idx = -1;

  int len = VSCMIN(RING_BUF_CELL_SIZE, ring->size - ring->widx + ring->ridx);
  if (likely(len)) {
    idx = (ring->widx & (ring->size - 1));
    ring->buff[idx] = (uint64_t)dataptr;
    ring->widx += RING_BUF_CELL_SIZE;
  }

  return idx;
}

static int get_from_ring_buf_st(vsc_ring_buf_t *ring, void **dataptr)
{
  int idx = -1;

  int len = VSCMIN(RING_BUF_CELL_SIZE, ring->widx - ring->ridx);
  if (likely(len)) {
    idx = (ring->ridx & (ring->size - 1));
    *dataptr = (void *)ring->buff[idx];
    ring->ridx += RING_BUF_CELL_SIZE;
  }

  return idx;
}

static int put_to_ring_buf(vsc_ring_buf_t *ring, void *dataptr)
{
  int idx = -1;

  pthread_spin_lock(&ring->lock);

  int len = VSCMIN(RING_BUF_CELL_SIZE, ring->size - ring->widx + ring->ridx);
  if (likely(len)) {
    idx = (ring->widx & (ring->size - 1));
    ring->buff[idx] = (uint64_t)dataptr;
    ring->widx += RING_BUF_CELL_SIZE;
    ring->num++;
    ring->total_num++;
    if(ring->num > ring->max_num)
        ring->max_num = ring->num;
    if(ring->start_time == 0) {
        struct timespec t;
        clock_gettime(CLOCK_MONOTONIC, &t);
        ring->start_time = t.tv_sec * 1000000000 + t.tv_nsec;
    }
  }

  pthread_spin_unlock(&ring->lock);

  return idx;
}

static int get_from_ring_buf(vsc_ring_buf_t *ring, void **dataptr)
{
  int idx = -1;

  pthread_spin_lock(&ring->lock);

  int len = VSCMIN(RING_BUF_CELL_SIZE, ring->widx - ring->ridx);
  if (likely(len)) {
    idx = (ring->ridx & (ring->size - 1));
    *dataptr = (void *)ring->buff[idx];
    ring->ridx += RING_BUF_CELL_SIZE;
    ring->num--;
  }

  pthread_spin_unlock(&ring->lock);

  return idx;
}


static void init_ring_buf(vsc_ring_buf_t *ring, uint64_t blknum)
{
  //memset(ring->buff, 0, size);
  for (uint64_t k = 0; k < blknum; k++) ring->buff[k] = 0;
  ring->size = blknum;//<<3
  ring->widx = 0;
  ring->ridx = 0;
  ring->max_num = 0;
  ring->num = 0;
  ring->total_num = 0;
  ring->start_time = 0;

  pthread_spin_init(&ring->lock, PTHREAD_PROCESS_PRIVATE);
}


static vsc_ring_buf_t * ring_alloc(uint64_t blknum)
{
  if (!is_power_of_2(blknum)) blknum = get_power_2(blknum);
  uint64_t totalbyte = (blknum<<3)/*must have ()*/ + sizeof(vsc_ring_buf_t);

  vsc_ring_buf_t *ring = (vsc_ring_buf_t*)malloc(totalbyte);
  if (!ring) return 0;

  init_ring_buf(ring, blknum);
  return ring;
}

static void print_ring_info(vsc_ring_buf_t *ring)
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    unsigned long long time = t.tv_sec*1000000000 + t.tv_nsec;
    time = (time - ring->start_time)/1000000000;
    
    fprintf(stderr, "average queue depth=%d\n", ring->total_num/time);
}

static void ring_free(vsc_ring_buf_t *ring)
{
  if (ring) {
    pthread_spin_destroy(&ring->lock);
    free((void*)ring);
  }
}


#endif//_VSC_KFIFO_H_
