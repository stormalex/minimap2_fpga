#include "kfifo.h"

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <string.h>
#include <errno.h>

#include <time.h>
#include <unistd.h>


typedef struct _student_info
{
  uint64_t stu_id;
  uint32_t age;
  uint32_t score;
}student_info;

uint64_t alloc_num = 0;
uint64_t free_num = 0;
uint64_t failed_num = 0;

int thrdnum = 1, runnum = 10;

void print_student_info(const student_info *stu_info)
{
  return;
  printf("id:%lu\t",stu_info->stu_id);
  printf("age:%u\t",stu_info->age);
  printf("score:%u\n",stu_info->score);
}

student_info * get_student_info(time_t timer)
{
  student_info *stu_info = (student_info *)malloc(sizeof(student_info));
  if (!stu_info)
    {
      fprintf(stderr, "Failed to malloc memory.\n");
      return NULL;
    }
  int curr = __sync_fetch_and_add(&alloc_num, 1);
  srand(timer);
  stu_info->stu_id = curr;
  stu_info->age = rand() % 30;
  stu_info->score = rand() % 101;
  print_student_info(stu_info);
  return stu_info;
}

void * consumer_proc(void *arg)
{
  vsc_ring_buf_t *ring_buf = (vsc_ring_buf_t *)arg;
  student_info *stu_info; 
  while(1)
    {
      //sleep(2);
      stu_info = 0;
      //printf("------------------------------------------\n");
      //printf("get a student info from ring buffer.\n");
      get_from_ring_buf(ring_buf, (void **)&stu_info);
      if (stu_info) {
      print_student_info(stu_info);
      free(stu_info);
      free_num++;

      }
      //printf("------------------------------------------\n");
    }
  return (void *)ring_buf;
}

void * producer_proc(void *arg)
{
  int localnum = 0;
  time_t cur_time;
  vsc_ring_buf_t *ring_buf = (vsc_ring_buf_t *)arg;
  while(1)
    {
      time(&cur_time);
      srand(cur_time);
      int seed = rand() % 11111;
      //printf("\n\n");
      //printf("******************************************\n");
      student_info *stu_info = get_student_info(cur_time + seed);
      //printf("put a student info to ring buffer.\n");
      int succ = put_to_ring_buf(ring_buf, (void *)stu_info);
      if (succ < 0) {free(stu_info);__sync_fetch_and_add(&failed_num, 1);}

      //printf("******************************************\n");
      //printf("\n\n");
      if (++localnum == runnum) break;
      //sleep(1);
    }
  return (void *)ring_buf;
}

pthread_t consumer_thread(void *arg)
{
  int err;
  pthread_t tid;
  err = pthread_create(&tid, NULL, consumer_proc, arg);
  if (err != 0)
    {
      fprintf(stderr, "Failed to create consumer thread.errno:%u, reason:%s\n",
	      errno, strerror(errno));
      return -1;
    }
  return tid;
}
pthread_t producer_thread(void *arg)
{
  int err;
  pthread_t tid;
  for (int k = 0; k < thrdnum; k++)
    err = pthread_create(&tid, NULL, producer_proc, arg);
  if (err != 0)
    {
      fprintf(stderr, "Failed to create consumer thread.errno:%u, reason:%s\n",
	      errno, strerror(errno));
      return -1;
    }
  return tid;
}

void * stat_thread(void*)
{
  while(1)
    {
      sleep(3);
      printf("alloc_num %llu, free_num %llu, failed_num %llu\n",
	     alloc_num ,free_num , failed_num );
    }
}

int main(int argc, char *argv[])
{
  void * buffer = NULL;
  uint32_t size = 0;
  vsc_ring_buf_t *ring_buf = NULL;
  pthread_t consume_pid, produce_pid;

  ring_buf = ring_alloc(256);
  if (!ring_buf)
    {
      fprintf(stderr, "Failed to init ring buffer.\n");
      return -1;
    }
  thrdnum = atoi(argv[1]);
  runnum = atoi(argv[2]);

  printf("multi thread test..threadnum %u, runnum %u.....\n",thrdnum, runnum);
  produce_pid  = producer_thread((void*)ring_buf);
  consume_pid  = consumer_thread((void*)ring_buf);

  pthread_t stat;
  pthread_create(&stat, 0, stat_thread, 0);
  pthread_join(produce_pid, NULL);
  pthread_join(consume_pid, NULL);

  ring_free(ring_buf);

  return 0;
}
