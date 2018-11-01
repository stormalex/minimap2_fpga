#ifndef KTHREAD_H
#define KTHREAD_H

#ifdef __cplusplus
extern "C" {
#endif

void kt_for(int n_threads, void (*func)(void*,long,int), void *data, long n);
void kt_pipeline(int n_threads, void *(*func)(void*, int, void*), void *shared_data, int n_steps);

void kt_for_nowait(int n_threads, void (*func)(void*,long,int), void *data, long n, void (*last)(void), pthread_t **tidptr, void **wptr, void **kfptr);

void kt_for_nowait_clone(int n_threads, int (*func)(void*, int), void *data, pthread_t **tidptr, void**wptr, void **kfptr);


#ifdef __cplusplus
}
#endif

#endif
