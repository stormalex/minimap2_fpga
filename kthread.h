#ifndef KTHREAD_H
#define KTHREAD_H

#ifdef __cplusplus
extern "C" {
#endif

void kt_for(int n_threads, void (*func)(void*,long,int), void *data, long n);
void kt_for_map(int n_threads, void (*func)(void*,long,int,void*), void *data, long n, void* params, void (*last_func)(void*, int), void* (*extra_func)(void*));
void kt_pipeline(int n_threads, void *(*func)(void*, int, void*), void *shared_data, int n_steps);

#ifdef __cplusplus
}
#endif

#endif
