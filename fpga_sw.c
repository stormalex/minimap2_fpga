#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <pthread.h>
#include <unistd.h>
#include <errno.h>

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

void sw(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat,
				   int8_t q, int8_t e, int8_t q2, int8_t e2, int w, int zdrop, int end_bonus, int flag, ksw_extz_t *ez)
{
    int r, t, qe = q + e, n_col_, *off = 0, *off_end = 0, last_st, last_en, wl, wr, max_sc, min_sc, long_thres, long_diff;
	int with_cigar = !(flag&KSW_EZ_SCORE_ONLY), approx_max = !!(flag&KSW_EZ_APPROX_MAX);
	int32_t *H = 0, H0 = 0, last_H0_t = 0;
	uint8_t *qr, *sf, *mem2 = 0;
	
	
	int8_t q_, q2_, qe_, qe2_, sc_mch_, sc_mis_, m1_, sc_N_;
	int8_t *u, *v, *x, *y, *x2, *y2, *s, *p = 0;
//init ez
	ksw_reset_extz(ez);

	//改到上层调用
	if (m <= 1 || qlen <= 0 || tlen <= 0) return;
	//软件上层也需要做一次该判断和刷新
	if (q2 + e2 < q + e) t = q, q = q2, q2 = t, t = e, e = e2, e2 = t; // make sure q+e no larger than q2+e2
// set  the vector value of penalty or bonus
	//zero_ = 0;
	q_ = q;
	q2_ = q2;
	qe_ = q + e;
	qe2_ = q2 + e2;
	sc_mch_ = mat[0];
	sc_mis_ = mat[1];
	sc_N_ = -e2;
	m1_ = m- 1;

	if (w < 0) w = tlen > qlen? tlen : qlen;
	//bandwidth
	wl = wr = w;
	n_col_ = qlen < tlen? qlen : tlen;
	n_col_ = (n_col_ < w + 1? n_col_ : w + 1);

	//初始化SW的参数时，完成for循环步骤，求得max_sc min_sc并在初始化阶段发送到硬件
	for (t = 1, max_sc = mat[0], min_sc = mat[1]; t < m * m; ++t) {
		max_sc = max_sc > mat[t]? max_sc : mat[t];
		min_sc = min_sc < mat[t]? min_sc : mat[t];
	}
	//改到上层调用判断
	if (-min_sc > 2 * (q + e)) return; // otherwise, we won't see any mismatches
	//long_diff long_thres为软件计算好后，在init时传入硬件
	long_thres = e != e2? (q2 - q) / (e - e2) - 1 : 0;
	if (q2 + e2 + long_thres * e2 > q + e + long_thres * e)
		++long_thres;
	long_diff = long_thres * (e - e2) - (q2 - q) - e2;
	//fprintf(stderr,"alloc u\n");
	u = (int8_t*)calloc(tlen,sizeof(int8_t));
	v = (int8_t*)calloc(tlen,sizeof(int8_t));
	x = (int8_t*)calloc(tlen,sizeof(int8_t));
	y = (int8_t*)calloc(tlen,sizeof(int8_t));
	x2 = (int8_t*)calloc(tlen,sizeof(int8_t));
	y2 = (int8_t*)calloc(tlen,sizeof(int8_t));
	s = (int8_t*)calloc(tlen,sizeof(int8_t));
	sf = (uint8_t*)calloc(tlen,sizeof(int8_t));
	qr = (uint8_t*)calloc(qlen,sizeof(int8_t));
	//initialize u/v/x/x2/y/y2
	memset(u,  -q  - e,  tlen);
	memset(v,  -q  - e,  tlen);
	memset(x,  -q  - e,  tlen);
	memset(y,  -q  - e,  tlen);
	memset(x2, -q2 - e2, tlen);
	memset(y2, -q2 - e2, tlen);
	//approx_max = 0
	if (!approx_max) {
		H = (int32_t*)malloc( tlen * sizeof(int));
		for (t = 0; t < tlen; ++t) H[t] = KSW_NEG_INF;
	}

	if (with_cigar) {
		mem2 = (uint8_t*)malloc( ((qlen + tlen - 1) * n_col_ + 1));
		p = (int8_t*)((size_t)mem2);
		off = (int*)malloc( (qlen + tlen - 1) * sizeof(int) * 2);
		off_end = off + qlen + tlen - 1;
		p[0] = 0;
	}
	//traverse query seq
	for (t = 0; t < qlen; ++t) qr[t] = query[qlen - 1 - t];//qr存入query序列，并反转
	memcpy(sf, target, tlen); //sf 存入target序列
	//fprintf(stderr,"sc_mch_ is %d,sc_mis_ is %d\n",sc_mch_,sc_mis_);
	for (r = 0, last_st = last_en = -1; r < qlen + tlen - 1; ++r) {
		int st = 0, en = tlen - 1;
		int8_t x1, x21, v1;
		int32_t max_H = KSW_NEG_INF , max_t = -1, Hen1 = KSW_NEG_INF;
		//fprintf(stderr,"get query base and r = %d\n",r);
		uint8_t *qrr = qr + (qlen - 1 - r);
		//int8_t *u = (int8_t*)u, *v = (int8_t*)v, *x = (int8_t*)x, *x2 = (int8_t*)x2;
		//int8_t x1_, x21_, v1_;
		// find the boundaries
		//if r > qlen ,then st >0
		if (st < r - qlen + 1) st = r - qlen + 1;
		if (en > r) en = r;
		//wr = 751
		if (st < (r-wr+1)>>1) st = (r-wr+1)>>1; // take the ceil
		if (en > (r+wl)>>1) en = (r+wl)>>1; // take the floor
		if (st > en) {
			ez->zdropped = 1;
			break;
		}		
		//st0 = st, en0 = en;		
		//initialize x1 v1 
		if (st > 0) {
			//if (st - 1 >= last_st && st - 1 <= last_en) {
			if(st -1 > last_st || st -1 > last_en) fprintf(stderr,"error!\n");
			if (st - 1 == last_st) {
				//fprintf(stderr,"visit x/x2/u/v/y/y2 array\n");
				x1 = x[st - 1], x21 = x2[st - 1], v1 = v[st - 1]; // (r-1,s-1) calculated in the last round
			} else {
				x1 = -q - e, x21 = -q2 - e2;
				v1 = -q - e;
			}
		} else {
			x1 = -q - e, x21 = -q2 - e2;
			v1 = r == 0? -q - e : r < long_thres? -e : r == long_thres? long_diff : -e2;
		}
		if (en >= r) {
			((int8_t*)y)[r] = -q - e, ((int8_t*)y2)[r] = -q2 - e2;
			u[r] = r == 0? -q - e : r < long_thres? -e : r == long_thres? long_diff : -e2;
		}
		
		// loop fission: set scores first
		if (!(flag & KSW_EZ_GENERIC_SC)) {
			
			for (t = st; t <= en; t++) {
				uint8_t strq, strt,  mask,sqtmp,sttmp;
				int8_t tmp;
				strq = sf[t];
				strt = qrr[t];
				//mask = _mm_or_si128(_mm_cmpeq_epi8(strq, m1_), _mm_cmpeq_epi8(st, m1_));
				sqtmp = strq == m1_?1:0;
				sttmp = strt == m1_?1:0;
				mask =  sqtmp|sttmp;
				tmp = strq == strt?1:0;
				tmp = tmp ==1?sc_mch_:sc_mis_;
				tmp = mask == 1?sc_N_:tmp;
				s[t] = tmp;
				//if(t>16 && !((t-st)%16))fprintf(stderr,"target is %d,query is %d,s[t] is %d\n",sf[t],qrr[t],s[t]);
			}
		} 
		/*else {
			for (t = st; t <= en; ++t){
				((uint8_t*)s)[t] = mat[sf[t] * m + qrr[t]];
			}				
		}*/
		//x1_ = x1;//初始化值；
		//x21_ = x21;
		//v1_ = v1;
		//st_ = st , en_ = en ;
		assert(en - st + 1 <= n_col_);
		//max_H = H[en] = en > 0? H[en-1] + u[en] : H[en] + v[en]; // special casing the last element
		//max_t = en;
		//xt1 = x[r-1][t-1] 
		//vt1 = v[r-1][t-1]
		//x1_ = x[r-1][t]
		//v1_ = v[r-1][t]
		//a = x[r-1][t-1] + v[r-1][t-1]
		//b = y[r-1][t] + u[r-1][t]
		//a2 = x2[r-1][t-1] + v[r-1][t-1]
		//b2 = y2[r-1][t] + u[r-1][t]
		//fprintf(stderr,"\nr=%d  st=%d en=%d \npr: ",r,st,en);
		
		/*if (!with_cigar) { // score only
			if(!approx_max) Hen1 = en>0?H[en-1]:Hen1;
			for (t = st; t <= en; ++t) {	
				int8_t z, a, b, a2, b2, xt1, x2t1, vt1, ut, tmp;
				z = s[t];
				xt1 = x1;
				x1 = x[t];
				vt1 = v1;
				v1 = v[t];
				ut = u[t];
				b = y[t] + ut;
				x2t1 = x21;
				x21 = x2[t];
				b2 = y2[t] + ut;
				a = xt1 + vt1;
				a2 = x2t1 + vt1;								
				z = z>a?z:a;
				z = z>b?z:b;
				z = z>a2?z:a2;
				z = z>b2?z:b2;
				z = z<sc_mch_?z:sc_mch_;
				u[t] = z - vt1;
				v[t] = z - ut;
				tmp = z - q_;
				a = a - tmp;
				b = b - tmp;
				tmp = z - q2_;
				a2 = a2 - tmp;
				b2 = b2 - tmp;
				x[t] = a>0?a:0 - qe_;
				y[t] = b>0?b:0 - qe_;
				x2[t] = a2>0?a2:0 - qe2_;
				y2[t] = b2>0?b2:0 - qe2_;
				if(!approx_max){					
					if (r > 0) {	
						//Hen1 = en>0?H[en-1]:Hen1;
						if(t != en){
							H[t] += (int32_t)v[t];
							if (H[t] > max_H)
								max_H = H[t], max_t = t;
						}						
					}else H[0] = v[0] - qe, max_H = H[0], max_t = 0,Hen1 = H[0]; // special casing r==0
				}				
			}
		} else*/
		{
			int8_t *pr = p + r * n_col_ - st;
			off[r] = st, off_end[r] = en;
			if(!approx_max) Hen1 = en>0?H[en-1]:Hen1;
			//fprintf(stderr,"\nr=%d  st=%d en=%d\npr: ",r,st,en);
			for (t = st; t <= en; ++t) {				
				int8_t d, z, a, b, a2, b2, xt1, x2t1, vt1, ut, tmp;
				z = s[t];
				xt1 = x1;
				x1 = x[t];
				vt1 = v1;
				v1 = v[t];
				ut = u[t];
				b = y[t] + ut;
				x2t1 = x21;
				x21 = x2[t];
				b2 = y2[t] + ut;
				a = xt1 + vt1;
				a2 = x2t1 + vt1;		
				if(!(flag&KSW_EZ_RIGHT)){//right-extention and gap-filling
					d = a > z?1:0;
					z = z >a?z:a;
					d = b >z?2:d;
					z = z >b?z:b;
					d = a2 >z?3:d;
					z = z >a2?z:a2;
					d = b2 >z?4:d;
					z = z >b2?z:b2;
					z = z < sc_mch_?z:sc_mch_;
				}else{  //left extention
					d = z>a?0:1;
					z = z>a?z:a;
					d = z>b?d:2;
					z = z>b?z:b;
					d = z>a2?d:3;
					z = z>a2?z:a2;
					d = z>b2?d:4;
					z = z>b2?z:b2;
					z = z<sc_mch_?z:sc_mch_;
				}
				//block2
				u[t] = z - vt1;
				v[t] = z - ut;
				/*tmp = z - q_;
				a = a - tmp; 
				b = b - tmp;
				tmp = z - q2_;
				a2 = a2 - tmp;
				b2 = b2 - tmp;
				*/
				a = a + q_; a = a -z;
				b = b + q_; b = b -z;
				a2 = a2 + q2_;a2 = a2 -z;
				b2 = b2 + q2_;b2 = b2 -z;
				if(!(flag&KSW_EZ_RIGHT)){//right-extention and gap-filling
					d |= a>0?1<<3:0;
					d |= b>0?1<<4:0;
					d |= a2>0?1<<5:0;
					d |= b2>0?1<<6:0;					
				}else{ //left extention
					d |= a<0?0:1<<3;
					d |= b<0?0:1<<4;
					d |= a2<0?0:1<<5;
					d |= b2<0?0:1<<6;
				}
				x[t] = a>0?a - qe_:0 - qe_;					
				y[t] = b>0?b - qe_:0 - qe_;					
				x2[t] = a2>0?a2 - qe2_:0 - qe2_;					
				y2[t] = b2>0?b2 - qe2_:0 - qe2_;
				pr[t] = d;
				//fprintf(stderr,"%x ",pr[t]);
				if((x[t] -e) > 0 || x[t]<(0 - qe_)) 
				fprintf(stderr,"x over bounded and x is %d\n",x[t]);
				if((y[t] -e) > 0 || y[t]<(0 - qe_)) 
				fprintf(stderr,"y over bounded and y is %d\n",y[t]);
				if((x2[t] -e2) > 0 || x2[t]<(0 - qe2_)) 
				fprintf(stderr,"x2 over bounded and x2 is %d\n",x2[t]);
				if((y2[t] -e2) > 0 || y2[t]<(0 - qe2_)) 
				fprintf(stderr,"y2 over bounded and y2 is %d\n",y2[t]);
				if((v[t] - sc_mch_ -qe_) > 0 || v[t] < (0 - qe_)) 
				fprintf(stderr,"v over bounded and v is %d\n",v[t]);
				if((u[t] - sc_mch_ -qe_) > 0 || u[t] < (0 - qe_)) 
				fprintf(stderr,"u over bounded and u is %d\n",u[t]);
				if(!approx_max){					
					if (r > 0) {
						//Hen1 = H[en-1];//Hen1 = H[en-1]r-1,t-1
						if(t != en){						
							H[t] += (int32_t)v[t];
							if (H[t] > max_H)
								max_H = H[t], max_t = t;
						}						
					}else H[0] = v[0] - qe, max_H = H[0], max_t = 0; // special casing r==0
				}
			}			 
		}
		if(!approx_max){// update ez
			if(r > 0){
				H[en] = en > 0? Hen1 + u[en] : H[en] + v[en]; // special casing the last element
				if (H[en] >= max_H)
						max_H = H[en], max_t = en;
			}
			//fprintf(stderr,"r = %d,en = %d,u[en] = %d,H[en] = %d,Hen1 is %d,max_H is %d\n",r,en,u[en],H[en],Hen1,max_H);
			if (en == tlen - 1 && H[en] > ez->mte)
				ez->mte = H[en], ez->mte_q = r - en;
			if (r - st == qlen - 1 && H[st] > ez->mqe)
				ez->mqe = H[st], ez->mqe_t = st;
			if (ksw_apply_zdrop(ez, 1, max_H, r, max_t, zdrop, e2)) break;
			if (r == qlen + tlen - 2 && en == tlen - 1)
				ez->score = H[tlen - 1];
		}else{ // find approximate max; Z-drop might be inaccurate, too.
			if (r > 0) {
				if (last_H0_t >= st && last_H0_t <= en && last_H0_t + 1 >= st && last_H0_t + 1 <= en) {
					int32_t d0 = (int32_t)v[last_H0_t];
					int32_t d1 = (int32_t)u[last_H0_t + 1];
					if (d0 > d1) H0 += d0;
					else H0 += d1, ++last_H0_t;
				} else if (last_H0_t >= st && last_H0_t <= en) {
					H0 += v[last_H0_t];
				} else {
					//++last_H0_t, H0 += u[last_H0_t];
					if(st - last_H0_t > 1) {fprintf(stderr,"error!error!error!\n");}
					++last_H0_t,H0 += u[st];
				}
			} else H0 = v[0] - qe, last_H0_t = 0;
			if ((flag & KSW_EZ_APPROX_DROP) && ksw_apply_zdrop(ez, 1, H0, r, last_H0_t, zdrop, e2)) break;
			if (r == qlen + tlen - 2 && en == tlen - 1)
				ez->score = H0;
		}
		
		last_st = st, last_en = en;
	}
	free(u);
	free(v);
	free(x);
	free(x2);
	free(y);
	free(y2);
	free(s);
	free(qr);
	free(sf);

	if (!approx_max) free(H);
	if (with_cigar) { // backtrack
		int rev_cigar = !!(flag & KSW_EZ_REV_CIGAR);
		if (!ez->zdropped && !(flag&KSW_EZ_EXTZ_ONLY)) {
			ksw_backtrack(km, 1, rev_cigar, 0, (uint8_t*)p, off, off_end, n_col_, tlen-1, qlen-1, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
		} else if (!ez->zdropped && (flag&KSW_EZ_EXTZ_ONLY) && ez->mqe + end_bonus > ez->max) {
			ez->reach_end = 1;
			ksw_backtrack(km, 1, rev_cigar, 0, (uint8_t*)p, off, off_end, n_col_, ez->mqe_t, qlen-1, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
		} else if (ez->max_t >= 0 && ez->max_q >= 0) {
			ksw_backtrack(km, 1, rev_cigar, 0, (uint8_t*)p, off, off_end, n_col_, ez->max_t, ez->max_q, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
		}
        if (!rev_cigar) ez->revcigar = 0xA0B1C2D3;
		free(mem2); free(off);
		//fprintf(stderr,"free mem2,off\n");
	}
}

//保存结果的数组
static void* results_array[4096];
static pthread_mutex_t results_mutex = PTHREAD_MUTEX_INITIALIZER;
static int results_head = 0;
static int results_tail = 0;
static int result_num = 0;

static void init_queue()
{
    memset(results_array, 0, sizeof(results_array));
    results_head = 0;
    results_tail = 0;
    result_num = 0;
    return;
}

static int result_is_full()
{
    return ((results_tail + 1) % 4096 == results_head);
}
static int result_is_empty()
{
    return (results_head == results_tail);
}

static int send(void* results)
{
    pthread_mutex_lock(&results_mutex);
    if (result_is_full()) {
        pthread_mutex_unlock(&results_mutex);
        //fprintf(stderr, "result ringbuf is full\n");
        return -1;
    }
    
    results_array[results_tail] = results;
    result_num++;
    //fprintf(stderr, "result insert to %d\n", results_tail);
    results_tail = (results_tail + 1) % 4096;
        
    pthread_mutex_unlock(&results_mutex);
    return 0;
}

static void* get()
{
    void* tmp;
    //sleep(1);
    pthread_mutex_lock(&results_mutex);
    if (result_is_empty()) {
        pthread_mutex_unlock(&results_mutex);
        //fprintf(stderr, "result ringbuf is empty\n");
        return NULL;
    }
    tmp = results_array[results_head];
    results_array[results_head] = NULL;
    result_num--;
    results_head = (results_head + 1) % 4096;
    
    pthread_mutex_unlock(&results_mutex);
    
    return tmp;
}

void* sim_fpga(char* buf_4k, sw_task_t** sw_tasks, int num, int* size)
{
    int i = 0;
    fpga_task_t* head = NULL;
    ksw_extz_t* ez_addr;
    uint32_t *cigar_addr;
    char* p = NULL;
    ksw_extz_t ez;
    
    char* buf = (char*)malloc(4*10224*1024);
    
    memcpy(buf, buf_4k, 4096);
    head = (fpga_task_t*)buf;
    if(head->check_id != CHECK_ID) {
        fprintf(stderr, "head->check_id=%x, error buf\n", head->check_id);
        exit(0);
    }

    ez_addr = (ksw_extz_t*)((char*)buf + 4096);
    cigar_addr = (uint32_t *)((char*)ez_addr + head->sw_num * sizeof(ksw_extz_t));
    
    for(i = 0; i < num; i++) {
        sw_task_t* task = sw_tasks[i];
        memset(&ez, 0, sizeof(ez));
        sw(NULL, task->qlen, task->query, task->tlen, task->target, 5, task->mat, task->q, task->e, task->q2, task->e2, task->w, task->zdrop, task->end_bonus, task->flag, &ez);
        //ksw_extd2_sse(NULL, task->qlen, task->query, task->tlen, task->target, 5, task->mat, task->q, task->e, task->q2, task->e2, task->w, task->zdrop, task->end_bonus, task->flag, &ez);
        *ez_addr = ez;
        memcpy(cigar_addr, ez.cigar, ez.n_cigar * sizeof(uint32_t));
        free(ez.cigar);
        
        p = (char*)cigar_addr;
        cigar_addr = (uint32_t *)(p + ADDR_ALIGN(ez_addr->n_cigar*sizeof(uint32_t), 16));
        ez_addr += 1;
        
        free(task);
    }
    *size = ((char*)cigar_addr - buf);
    return buf;
}

static int send_thread_flag = 1;
//static chain_sw_task_t* task_array[SEND_ARRAY_SIZE];

void stop_fpga_thread()
{
    send_thread_flag = 0;
}
#define SW_TASK_NUM 20000

void* send_task_thread(void* arg)
{
    int tid = (int)arg;
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

    char* in_file = (char*)malloc(100);
    char* out_file = (char*)malloc(100);
    chain_sw_task_t** task_array = (chain_sw_task_t**)malloc(SEND_ARRAY_SIZE * sizeof(chain_sw_task_t*));
    
    init_queue();
    while(send_thread_flag)
    {
        chain_sw_task_t* chain_tasks = get_task();
        if(chain_tasks == NULL) {
            wait_number++;
            sleep(1);
            if(wait_number >= 3) {
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

            if(chain_num == 0) {
                continue;
            }

            if((data_size + 4096) > 4 * 1024 * 1024) {
                fprintf(stderr, "data_size too large:%d\n", data_size);
                exit(0);
            }
#if DUMP_FILE
            fpga_buf = (char*)malloc(4 * 1024 * 1024);
            if(fpga_buf == NULL) {
                fprintf(stderr, "malloc error\n");
                exit(0);
            }
            int all_sw_index = 0;
            sw_task_t** all_sw_tasks = (sw_task_t**)malloc(SW_TASK_NUM * sizeof(sw_task_t*));
            if(all_sw_tasks == NULL) {
                fprintf(stderr, "malloc sw tasks error\n");
                exit(0);
            }
#else
            fpga_buf = (char*)fpga_get_writebuf(data_size + 4096, BUF_TYPE_SW);
            if(fpga_buf == NULL) {
                fprintf(stderr, "fpga_get_writebuf sw error\n");
                exit(0);
            }
            
#endif
            head = (fpga_task_t*)fpga_buf;        //指向buff首地址
            head->check_id = CHECK_ID;
            head->sw_num = 0;
            chain_head = (fpga_task_id_t*)(fpga_buf + sizeof(fpga_task_t));
            sw_task = (fpga_sw_task*)(fpga_buf + 4096);
            //fprintf(stderr, "fpga_buf=%p\n", fpga_buf);
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
                //fprintf(stderr, "chain_head=%p\n", chain_head);
                chain_head->read_id = task_array[chain_index]->read_id;
                chain_head->chain_id = task_array[chain_index]->chain_id;
                chain_head->sw_num = task_array[chain_index]->sw_num;

                //设置sw任务的数据
                head->sw_num += chain_head->sw_num;     //统计一次调用sw任务数
                for(sw_index = 0; sw_index < chain_head->sw_num; sw_index++) {
#if DUMP_FILE
                    all_sw_tasks[all_sw_index++] = task_array[chain_index]->sw_tasks[sw_index];
#endif
                    sw_task->qlen = task_array[chain_index]->sw_tasks[sw_index]->qlen;
                    sw_task->tlen = task_array[chain_index]->sw_tasks[sw_index]->tlen;
                    sw_task->flag = task_array[chain_index]->sw_tasks[sw_index]->flag;
                    sw_task->zdrop = task_array[chain_index]->sw_tasks[sw_index]->zdrop;
                    sw_task->bw = task_array[chain_index]->sw_tasks[sw_index]->w;
                    sw_task->end_bonus = task_array[chain_index]->sw_tasks[sw_index]->end_bonus;
                    
                    uint8_t *qseq = (uint8_t*)sw_task + sizeof(fpga_sw_task);    //设置qseq的地址
                    memcpy(qseq, task_array[chain_index]->sw_tasks[sw_index]->query, sw_task->qlen);
                    uint8_t *tseq = qseq + ADDR_ALIGN(sw_task->qlen, 16);
                    memcpy(tseq, task_array[chain_index]->sw_tasks[sw_index]->target, sw_task->tlen);
                    
                    p = (char*)sw_task;
                    sw_task = (fpga_sw_task*)(p + sizeof(fpga_sw_task) + ADDR_ALIGN(sw_task->qlen, 16) + ADDR_ALIGN(sw_task->tlen, 16));    //指向下一个sw任务的地址
                }

                
                chain_head += 1;   //指向下一个chain的头
#if DUMP_FILE
                ;
#else
                destroy_chain_sw_task(task_array[chain_index]);     //销毁chain任务的数据
                task_array[chain_index] = NULL;
#endif
            }

            if((data_size + 4096) != ((char*)sw_task - fpga_buf)) {
                fprintf(stderr, "fpga_buf=%p, sw_task=%p, data_size=0x%x\n", fpga_buf, sw_task, data_size);
                exit(0);
            }
            else
                //fprintf(stderr, "2.fpga_buf=%p, sw_task=%p, data_size=0x%x\n", fpga_buf, sw_task, data_size);
            //fprintf(stderr, "send:tag:%lld, chain num:%d, sw num:%d\n", head->tag, head->chain_task_num, head->sw_num);
#if DUMP_FILE
            {
                fprintf(stderr, "write data\n");
                int ret = 0;
                int real_size = data_size + 4096;
                memset(in_file, 0, sizeof(in_file));
                memset(out_file, 0, sizeof(out_file));
                sprintf(in_file, "sw_%d_in.bin", tid);
                sprintf(out_file, "sw_%d_out.bin", tid);
                FILE* fp = fopen(in_file, "a");
                if(fp == NULL) {
                    fprintf(stderr, "fopen %d failed:%s\n", in_file, strerror(errno));
                    exit(0);
                }
                ret = fwrite(&real_size, 1, sizeof(real_size), fp);
                if(ret != sizeof(real_size)) {
                    fprintf(stderr, "fwrite failed:%s\n", strerror(errno));
                    exit(0);
                }
                ret = fwrite(fpga_buf, 1, real_size, fp);
                if(ret != real_size) {
                    fprintf(stderr, "fwrite failed:%s\n", strerror(errno));
                    exit(0);
                }
                fclose(fp);
                void* result_addr = sim_fpga(fpga_buf, all_sw_tasks, all_sw_index, &real_size);
                free(all_sw_tasks);
                fp = fopen(out_file, "a");
                if(fp == NULL) {
                    fprintf(stderr, "fopen %s failed:%s\n", out_file, strerror(errno));
                    exit(0);
                }
                ret = fwrite(&real_size, 1, sizeof(real_size), fp);
                if(ret != sizeof(real_size)) {
                    fprintf(stderr, "fwrite failed:%s\n", strerror(errno));
                    exit(0);
                }
                ret = fwrite(result_addr, 1, real_size, fp);
                if(ret != real_size) {
                    fprintf(stderr, "fwrite failed:%s\n", strerror(errno));
                    exit(0);
                }
                fclose(fp);
                
                free(fpga_buf);
                
                while(send(result_addr));
            }
#else
            fpga_writebuf_submit(fpga_buf, data_size + 4096, TYPE_SW);
            
#endif
            chain_num = 0;
            data_size = 0;
        }
    }
    free(task_array);
    free(in_file);
    free(out_file);
    return NULL;
}

void* recv_task_thread(void *arg)
{
    char* fpga_buf = NULL;
    int fpga_len;
    fpga_task_t* head = NULL;
    fpga_task_id_t* chain_head = NULL;
    int chain_index;
    int sw_index;
    char* p = NULL;
    ksw_extz_t* ez_addr = NULL;
    uint32_t *cigar_addr = NULL;

    while(send_thread_flag) {
#if DUMP_FILE
        fpga_buf = (char*)get();
        if(fpga_buf == NULL) {
            continue;
        }
#else
        fpga_buf = (char*)fpga_get_retbuf(&fpga_len, RET_TYPE_SW);
        if(fpga_buf == NULL || fpga_len == 0) {
            fprintf(stderr,"fpga_len:%d\n",fpga_len);
            return NULL;
        }
        if(fpga_len > 4 * 1024 * 1024) {
            fprintf(stderr, "fpga_len=%d\n", fpga_len);
            exit(0);
        }
#endif
        head = (fpga_task_t*)fpga_buf;
        if(head->check_id != CHECK_ID) {
            fpga_release_retbuf(fpga_buf);
            fprintf(stderr, "head->check_id=%x, error buf\n", head->check_id);
            continue;
        }
        //fprintf(stderr, "recv:tag:%lld, chain num:%d, sw num:%d\n", head->tag, head->chain_task_num, head->sw_num);
        int chain_num = head->chain_task_num;
        chain_head = (fpga_task_id_t*)((char*)fpga_buf + sizeof(fpga_task_t));
        ez_addr = (ksw_extz_t*)((char*)fpga_buf + 4096);
        cigar_addr = (uint32_t *)((char*)ez_addr + head->sw_num * sizeof(ksw_extz_t));

        //处理每一个chain
        for(chain_index = 0; chain_index < chain_num; chain_index++) {
            sw_result_t* result = create_result();
            result->read_id = chain_head[chain_index].read_id;
            result->chain_id = chain_head[chain_index].chain_id;

            for(sw_index = 0; sw_index < chain_head[chain_index].sw_num; sw_index++) {
                
                //生成一个ez结果
                ksw_extz_t* ez = (ksw_extz_t*)malloc(sizeof(ksw_extz_t));
                *ez = *ez_addr;
                ez->m_cigar = ez_addr->n_cigar + 1;
                
                
                if(ez->n_cigar == 0) {
                    //fprintf(stderr, "ez->n_cigar=0\n");
                }
                else {
                    ez->cigar = (uint32_t *)malloc(ez->m_cigar * sizeof(uint32_t));
                    memcpy(ez->cigar, cigar_addr, ez->n_cigar * sizeof(uint32_t));
                }
                
                ez_addr += 1;   //移动ez到下一个

                p = (char*)cigar_addr;
                cigar_addr = (uint32_t *)(p + ADDR_ALIGN(ez->n_cigar*sizeof(uint32_t), 16));  //移动到下一个cigar开头
                add_result(result, ez); //将ez结果加入到chain结果中
            }
            
            while(send_result(result)); //发送给结果处理线程

        }
#if DUMP_FILE
        free(fpga_buf);
#else
        fpga_release_retbuf(fpga_buf);
#endif
    }
    return NULL;
}

