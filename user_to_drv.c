#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <pthread.h>
#include "ksw2.h"
#include "kalloc.h"

#include "stat.h"
#include "kfifo.h"
#include "context_data.h"
#include "user_to_drv.h"

#include "fpga_sim.h"
#include "fpga.h"

#include "user_common.h"

static inline void seq_rev(uint32_t len, uint8_t *seq)
{
    uint32_t i;
    uint8_t t;
    for (i = 0; i < len>>1; ++i)
        t = seq[i], seq[i] = seq[len - 1 - i], seq[len - 1 - i] = t;
}

#define seq4_get(s, i)    ((s)[(i)>>3] >> (((i)&7)<<2) & 0xf)


static int idx_getseq(uint32_t n_seq, uint32_t slen, uint32_t rid, uint32_t st, uint32_t en, uint64_t offset, uint32_t *S, uint8_t *seq)
{
    uint64_t i, st1, en1;
    if (rid >= n_seq/*mi->n_seq*/ || st >= slen/*mi->seq[rid].len*/) return -1;
    if (en > slen/*mi->seq[rid].len*/) en = slen;/*mi->seq[rid].len;*/
    st1 = offset/*mi->seq[rid].offset*/ + st;
    en1 = offset/*mi->seq[rid].offset*/ + en;
    for (i = st1; i < en1; ++i)
        seq[i - st1] = seq4_get(S/*mi->S*/, i);
    return en - st;
}


pthread_t *start_task_thread(int nthread, CALLBACK callback, void *param)
{
    pthread_t *tid = (pthread_t *)malloc(nthread*sizeof(pthread_t));
    for (int k=0; k<nthread; k++) 
    {
        int err = pthread_create(&tid[k], 0, callback, param);
        if (0 != err)
        {
            fprintf(stderr, "Failed to start thread, errno:%u, reason:[%s]\n",errno,strerror(errno));
        }
    }
    return tid;
} 



int dma_chaindp_prepare(uint32_t magic, uint16_t tid, uint8_t lat, uint16_t num, uint32_t size, void *dest, chaindp_x86_t *data)
{
    chaindp_sndhdr_t *hdr = (chaindp_sndhdr_t*)dest;
    hdr->magic= magic;
    hdr->size = size;
    hdr->type = CHAINDP;
    hdr->tid  = tid;
    hdr->num  = num;
    hdr->lat  = lat;

    int fixed = sizeof(chaindp_sndhdr_t)+ALIGN_BYTE_N(64, hdr->num*sizeof(chaindp_sndsubhdr_t));

    for (int k=0; k<hdr->num; k++)
    {
          if (!k) hdr->subhdr[k].offset = fixed;
          else    hdr->subhdr[k].offset = hdr->subhdr[k-1].offset + hdr->subhdr[k-1].size;
          hdr->subhdr[k].size   = sizeof(chaindp_to_fpga_t)+ ALIGN_BYTE_N(64, data[k].seednum*sizeof(sim_mm128_t));

          chaindp_to_fpga_t *curr = (chaindp_to_fpga_t*)((char *)hdr+hdr->subhdr[k].offset);
          curr->splitflag = 0xFFFFFFFF;
          curr->gap_ref = data[k].gap_ref;
          curr->gap_qry = data[k].gap_qry;
          curr->seednum = data[k].seednum;
          curr->mid_occ = data[k].mid_occ;
          curr->optflag = data[k].optflag;
          curr->qlen_sum= data[k].qlen_sum;
          curr->n_segs = data[k].n_segs;
          curr->pos = data[k].pos;

          int len=strnlen(data[k].qname, QUERY_NAME_MAXLEN-1);
          memcpy(curr->qname, data[k].qname, len);
          curr->qname[len] = 0;

          for(int m=0; m<curr->seednum; m++)
            curr->seed[m] = data[k].seed[m];
    }

  return 0;
}


int dma_sw_prepare(uint32_t magic, uint16_t tid, uint8_t lat, uint16_t num, uint32_t size, void *dest, sw_readhdr_t *data, void *km)
{
    sw_sndhdr_t *hdr = (sw_sndhdr_t *)dest;
    hdr->magic=magic;
    hdr->size = size;//
    hdr->type = KSW;
    hdr->tid  = tid;
    hdr->num  = num;
    hdr->lat  = lat;
    hdr->freecigar = 0;
    hdr->km   = 0;

    memcpy(hdr->data, data, sizeof(hdr->data));

    sw_to_fpga_t *fpghdr = (sw_to_fpga_t *)(hdr->data + DP_BATCHSIZE);
  
    for (int k = 0; k < num; k++) {
        sw_readhdr_t *read = &hdr->data[k];
        sw_reghdr_t  *curr = &hdr->data[k].reg;
        curr->size = 0;

        if (0 != read->longsw) {
            sw_node_t *node = 0;
            struct nl_list_head *head = (struct nl_list_head *)hdr->data[k].head;
            nl_list_for_each_entry(node, head, entry) {
                fpghdr->qlen = node->qlen;
                fpghdr->tlen = node->tlen;
                fpghdr->qlen_align = ALIGN_BYTE_N(16, fpghdr->qlen);
                fpghdr->tlen_align = ALIGN_BYTE_N(16, fpghdr->tlen);
                //fpghdr->qoff = sizeof(sw_to_fpga_t);
                //fpghdr->toff = fpghdr->qoff + fpghdr->qlen;
                //fpghdr->swpos = node->swpos;
                fpghdr->flag  = node->flag;
                fpghdr->zdrop = node->zdrop;
                //fpghdr->zdrop_inv = node->zdrop_inv;
                fpghdr->bw        = node->bw;
                fpghdr->end_bonus = node->end_bonus;

                uint8_t *qseq = (uint8_t *)fpghdr + sizeof(sw_to_fpga_t);
                uint8_t *tseq = (uint8_t *)fpghdr + sizeof(sw_to_fpga_t) + fpghdr->qlen_align;

                idx_getseq(node->n_seq,node->slen,node->rid,node->st,node->en,node->offset,node->S,tseq);
                memcpy(qseq, node->qseq, fpghdr->qlen);

                if (SW_POS_LEFT == node->swpos) {
                    seq_rev(fpghdr->tlen, tseq);
                    seq_rev(fpghdr->qlen, qseq);
                } else if (SW_POS_MIDL == node->swpos) {
                } else if (SW_POS_RIGH == node->swpos) {
                } else {
                    fprintf(stderr, "#######node->swpos error %x\n", node->swpos);
                }

                curr->size += sizeof(sw_to_fpga_t)+fpghdr->qlen_align+fpghdr->tlen_align;

                //fpghdr move to next
                fpghdr = (sw_to_fpga_t *)((uint8_t *)fpghdr + sizeof(sw_to_fpga_t) + fpghdr->qlen_align + fpghdr->tlen_align);
            }//node loop list
            sw_node_t *next = 0;
            nl_list_for_each_entry_safe(node, next, head, entry) {
                kfree(km, node);
            }
            kfree(km, head);
        }
        if (0 == k) curr->offset = sizeof(sw_sndhdr_t);
        else curr->offset = hdr->data[k-1].reg.offset + hdr->data[k-1].reg.size;

    }//for reg

    return 0;
}//end func


int dma_sw2_prepare(uint32_t magic, uint16_t tid, uint8_t lat, uint16_t num, uint32_t size, void *dest, sw_readhdr_t *data, void *km)
{
    sw_sndhdr_t *hdr = (sw_sndhdr_t *)dest;
    hdr->magic= magic;
    hdr->size = size;
    hdr->type = KSW;
    hdr->tid  = tid;
    hdr->num  = num;
    hdr->lat  = lat;
    hdr->freecigar = 0;
    hdr->km   = km;

    sw_to_fpga_t *fpghdr = (sw_to_fpga_t *)(hdr->data + DP_BATCHSIZE);
  
    int p = 0;
    for (int k = 0; k < num; k++) {
        sw_readhdr_t *read = &data[k];

        if (0 == read->longsw) {
            sw_readhdr_t *dtrd = &hdr->data[p];
            *dtrd = *read;
            sw_reghdr_t  *curr = &dtrd->reg;
            curr->size = 0;

            sw_node_t *node = 0;
            struct nl_list_head *head = (struct nl_list_head *)dtrd->head;
            nl_list_for_each_entry(node, head, entry) {
                fpghdr->qlen = node->qlen;
                fpghdr->tlen = node->tlen;
                fpghdr->qlen_align = ALIGN_BYTE_N(16, fpghdr->qlen);
                fpghdr->tlen_align = ALIGN_BYTE_N(16, fpghdr->tlen);
                //fpghdr->qoff = sizeof(sw_to_fpga_t);
                //fpghdr->toff = fpghdr->qoff + fpghdr->qlen;
                //fpghdr->swpos = node->swpos;
                fpghdr->flag  = node->flag;
                fpghdr->zdrop = node->zdrop;
                //fpghdr->zdrop_inv = node->zdrop_inv;
                fpghdr->bw        = node->bw;
                fpghdr->end_bonus = node->end_bonus;

                uint8_t *qseq = (uint8_t *)fpghdr + sizeof(sw_to_fpga_t);
                uint8_t *tseq = (uint8_t *)fpghdr + sizeof(sw_to_fpga_t) + fpghdr->qlen_align;

                idx_getseq(node->n_seq,node->slen,node->rid,node->st,node->en,node->offset,node->S,tseq);
                memcpy(qseq, node->qseq, fpghdr->qlen);

                if (SW_POS_LEFT == node->swpos) {
                    seq_rev(fpghdr->tlen, tseq);
                    seq_rev(fpghdr->qlen, qseq);
                } else if (SW_POS_MIDL == node->swpos) {
                } else if (SW_POS_RIGH == node->swpos) {
                } else {
                    fprintf(stderr, "#######node->swpos error %x\n", node->swpos);
                }

                curr->size += sizeof(sw_to_fpga_t)+fpghdr->qlen_align+fpghdr->tlen_align;

                //fpghdr move to next
                fpghdr = (sw_to_fpga_t *)((uint8_t *)fpghdr + sizeof(sw_to_fpga_t) + fpghdr->qlen_align + fpghdr->tlen_align);
            }//node loop list
            sw_node_t *next = 0;
            nl_list_for_each_entry_safe(node, next, head, entry) {
                kfree(km, node);
            }
            kfree(km, head);

            if (0 == p) curr->offset = sizeof(sw_sndhdr_t);
            else curr->offset = hdr->data[p-1].reg.offset + hdr->data[p-1].reg.size;
            p++;
        }//end longsw

    }//for reg

    return 0;
}//end func


void proc_dp(vsc_ring_buf_t *ring_buf)
{
    dp_to_ring_t *snd_data = 0;
    int idx = get_from_ring_buf(ring_buf, (void **)&snd_data);
    if (idx >= 0)
    {
        //call driver API to send data
        uint32_t subalignmem = ALIGN_BYTE_N(64, snd_data->factnum*sizeof(chaindp_sndsubhdr_t));
        uint32_t memsize = subalignmem+snd_data->total_size+snd_data->factnum*sizeof(chaindp_to_fpga_t)+sizeof(chaindp_sndhdr_t);
#ifdef RUN_ON_X86_DP
        void *dest = kmalloc(snd_data->km,memsize);
#else
        void *dest = fpga_get_writebuf(2*memsize);
#endif
        set_dpstat(memsize);

        dma_chaindp_prepare(snd_data->magic, snd_data->tid, snd_data->lastblk, snd_data->factnum, memsize,dest, snd_data->data);
#if 0
        {
            static int seq = 0;
            char fpn[64]={0};
            char drn[64]={0};
            snprintf(fpn, sizeof(fpn)-1, "dpfpga%d-in",seq);
            snprintf(drn, sizeof(drn)-1, "dpdriv%d-in",seq);
            seq++;
            dumpfile(memsize, dest, fpn, drn);
        }
#endif
#ifdef RUN_ON_X86_DP
        //sim fpga proc
        calc_chaindp_batch(dest, snd_data->km, snd_data->mi);
        kfree(snd_data->km, dest);
#else
        fpga_writebuf_submit(dest, 2*memsize, CHAINDP);
#endif

        for (int k=0; k<snd_data->factnum; k++) {
            chaindp_x86_t *curr = &snd_data->data[k];
            //kfree(curr->seed_km, curr->seed);
            kfree(snd_data->km, curr->seed);//maybe bug, not snd_data->km
        }
        kfree(snd_data->km, snd_data);

    }//end if (snd_data)
    return;
}


void* dp_task_sender(void *param)
{
    thread_ctrl_t *ctrl = (thread_ctrl_t*)param;
    vsc_ring_buf_t *ring_buf = ctrl->dp_ring_buf;

    do {
        proc_dp(ring_buf);
    } while(!ctrl->stop);

    return param;
}


void test_memalign(void)
{
    assert(64 == sizeof(chaindp_sndhdr_t));
    assert(320 == sizeof(chaindp_to_fpga_t));
    assert(16 == sizeof(sw_to_fpga_t));
    assert(16 == sizeof(sw_reghdr_t));
    assert(32 == sizeof(sw_readhdr_t));
    assert(4096 == sizeof(swtest_rcvhdr_t));
    assert(32+DP_BATCHSIZE*32 == sizeof(sw_sndhdr_t));
}

void proc_sw(vsc_ring_buf_t *ring_buf)
{
    sw_to_ring_t *snd_data = 0;
    int idx = get_from_ring_buf(ring_buf, (void **)&snd_data);
    if (idx >= 0)
    {
      //fprintf(stderr, "sw ring rcver     magic %d,factnum %d,size %d,tid %x,last %x\n",snd_data->magic,snd_data->factnum,snd_data->total_size,snd_data->tid,snd_data->lastblk);
        int hdrsize = sizeof(sw_sndhdr_t);
        //if (snd_data->factnum > 0) 
        {
            for (int k = 0; k < snd_data->factnum; k++)
            {
                sw_readhdr_t *read = snd_data->data + k;
                //if (read->regnum != 0)
                {
                    sw_reghdr_t *cr = &read->reg;
                    hdrsize += (!!get_left(cr->midnum) + !!get_right(cr->midnum) + clear_midnum(cr->midnum))*sizeof(sw_to_fpga_t);
                }
            }
        }
        int memsize = hdrsize + snd_data->total_size;
        //sim alloc mem from driver API
        void *dest = kmalloc(snd_data->km, memsize);
        set_swstat(snd_data->factnum, memsize);
        dma_sw_prepare(snd_data->magic, snd_data->tid, snd_data->lastblk, snd_data->factnum, memsize, dest, snd_data->data, snd_data->km);
        //sim fpga
        calc_sw(dest, snd_data->km, snd_data->opt);

        kfree(snd_data->km, dest);

        //free x86 mem
        //fprintf(stderr, "free ring ptr %p,idx %d\n", snd_data,idx);
        kfree(snd_data->km, snd_data);

        snd_data = 0;
    }//end if (snd_data)

    test_memalign();
    return;
}

pthread_mutex_t file_mutex = PTHREAD_MUTEX_INITIALIZER;
int cxx = 0;

void proc_sw2(vsc_ring_buf_t *ring_buf)
{
    sw_to_ring_t *snd_data = 0;
    int idx = get_from_ring_buf(ring_buf, (void **)&snd_data);
    if (idx >= 0)
    {
        fprintf(stderr, "sw2 ring rcver    magic %d,factnum %d,size %d,tid %x,last %x,km %p\n",snd_data->magic,snd_data->factnum,snd_data->total_size,snd_data->tid,snd_data->lastblk,snd_data->km);
        //prepare for long sw
        struct timespec tv;
        sw_to_ring_t *swring = sw_ring_init(snd_data->km,snd_data->opt,0,snd_data->tid,NT_LAST);
        memset(swring->data, 0x0, sizeof(swring->data));
        clock_gettime(CLOCK_MONOTONIC, &tv);
        swring->magic = tv.tv_nsec;

        int swcount = 0;
        int hdrsize = sizeof(sw_sndhdr_t);
        //if (snd_data->factnum > 0) 
        {
            for (int k = 0; k < snd_data->factnum; k++)
            {
                sw_readhdr_t *read = snd_data->data + k;
                if (0 == read->longsw)
                {
                    sw_reghdr_t *cr = &read->reg;
                    swcount += (!!get_left(cr->midnum) + !!get_right(cr->midnum) + clear_midnum(cr->midnum));
                    hdrsize += (!!get_left(cr->midnum) + !!get_right(cr->midnum) + clear_midnum(cr->midnum))*sizeof(sw_to_fpga_t);
                }
                else 
                {
                    sw_readhdr_t *rdrd = &swring->data[swring->factnum++];
                    *rdrd = *read;
                    swring->total_size += rdrd->reg.size;
                    snd_data->total_size -= rdrd->reg.size;
                    fprintf(stderr, "^^^^^^^^^^^^^magic %d, factnum %d, mid %x,swpos %d,regpos %d\n", swring->magic,swring->factnum,rdrd->reg.midnum,rdrd->ctxpos,rdrd->reg.regpos);
                }

            }
        }

        int memsize = hdrsize + snd_data->total_size;
        sw_sndhdr_t *dest = NULL;
        if(memsize > 4096) {
#ifdef RUN_ON_X86_SW
            //sim alloc mem from driver API
            dest = (sw_sndhdr_t *)kmalloc(snd_data->km, memsize);
#else
            dest = fpga_get_writebuf(memsize, BUF_TYPE_SW);
            fprintf(stderr, "fpga_get_writebuf dest=%p\n", dest);
#endif
            set_swstat(snd_data->factnum, memsize);
            dma_sw2_prepare(snd_data->magic,snd_data->tid, snd_data->lastblk, snd_data->factnum, memsize, dest, snd_data->data, snd_data->km);
            //update factnum
            dest->num -= swring->factnum;
            dest->swcount = swcount;
        }
        if (swring->factnum > 0) sw_ring_submit(swring, x86_ctrl.sw_ring_buf);
        else kfree(snd_data->km, swring);
#ifdef DUMP_FILES
        {
            static int seq = 0;
            char fpn[64]={0};
            char drn[64]={0};
            snprintf(fpn, sizeof(fpn)-1, "swfpga%d-in",seq);
            snprintf(drn, sizeof(drn)-1, "swdriv%d-in",seq);
            seq++;
            dumpfile(memsize, dest, fpn, drn);
        }
#endif
        if(memsize > 4096) {
#ifdef RUN_ON_X86_SW
            //sim fpga
            calc_sw(dest, snd_data->km, snd_data->opt);
            kfree(snd_data->km, dest);
#else 
            /*pthread_mutex_lock(&file_mutex);
            FILE* fp = fopen("./sw_in.bin", "ab+");
            fwrite(&memsize, 4, 1, fp);
            fwrite(dest, 1, memsize, fp);
            fclose(fp);
            pthread_mutex_unlock(&file_mutex);*/
            fpga_writebuf_submit(dest, memsize, TYPE_SW);
#endif
        }
        //free x86 mem
        //fprintf(stderr, "free ring ptr %p,idx %d\n", snd_data,idx);
        kfree(snd_data->km, snd_data);

        snd_data = 0;
    }//end if (snd_data)

    test_memalign();
    return;
}


void* sw_task_sender(void *param)
{
    thread_ctrl_t *ctrl = (thread_ctrl_t*)param;
    vsc_ring_buf_t *ring_buf = ctrl->sw_ring_buf;

    do {
        proc_sw(ring_buf);
    } while(0);//(!ctrl->stop);

    return param;
}


void* task_sender(void *param)
{
    thread_ctrl_t *ctrl = (thread_ctrl_t*)param;
    vsc_ring_buf_t *dpring_buf = ctrl->dp_ring_buf;
    vsc_ring_buf_t *swring_buf = ctrl->sw_ring_buf;
    do {
        proc_dp(dpring_buf);
#ifndef DEBUG_API
        proc_sw2(swring_buf);
#endif
    } while(!ctrl->stop);
    fprintf(stderr, "task sender thread exit!\n");
    return param;
}

void* task_sender_sw(void *param)
{
    thread_ctrl_t *ctrl = (thread_ctrl_t*)param;
    vsc_ring_buf_t *swring_buf = ctrl->sw_ring_buf;
    do {
#ifndef DEBUG_API
        proc_sw(swring_buf);
#endif
    } while(!ctrl->stop);
    fprintf(stderr, "task sender sw thread exit!\n");
    return param;
}


void update_value(uint64_t *data)
{
  //update offset,size,cigar where from FPGA
  swtest_rcvhdr_t *hdr = (swtest_rcvhdr_t *)data;
/*
  fprintf(stderr, "km = %p\n", hdr->km);
  fprintf(stderr, "magic = %d\n", hdr->magic);
  fprintf(stderr, "size  = %d\n", hdr->size);
  fprintf(stderr, "tid   = %d\n", hdr->tid);
  fprintf(stderr, "num   = %d\n", hdr->num);
  fprintf(stderr, "type  = %d\n", hdr->type);
  fprintf(stderr, "lat   = %d\n", hdr->lat);
  fprintf(stderr, "free  = %d\n", hdr->freecigar);
  fprintf(stderr, "swcnt = %d\n", hdr->swcount);
*/
  for (int x=0; x<hdr->num; x++) {
    sw_readhdr_t *rdhdr = hdr->data + x;
/*
    fprintf(stderr, "*****************\n");
    fprintf(stderr, "ctxpos %d\n", rdhdr->ctxpos);
    //fprintf(stderr, "regnum %d\n", rdhdr->regnum);
    fprintf(stderr, "longsw %d\n", rdhdr->longsw); */
    rdhdr->reg.size = (!!get_left(rdhdr->reg.midnum)+!!get_right(rdhdr->reg.midnum)+clear_midnum(rdhdr->reg.midnum))*sizeof(ksw_extz_t);
    if (0 == x) rdhdr->reg.offset = sizeof(swtest_rcvhdr_t);
    else rdhdr->reg.offset = hdr->data[x-1].reg.offset+hdr->data[x-1].reg.size;
/*
    fprintf(stderr, "reg.offset %d\n", rdhdr->reg.offset);
    fprintf(stderr, "reg.size   %d\n", rdhdr->reg.size);
    fprintf(stderr, "reg.midnum %x\n", rdhdr->reg.midnum);
    fprintf(stderr, "reg.regpos %d\n", rdhdr->reg.regpos);
*/
  }

  ksw_extz_t *ez = (ksw_extz_t *)(hdr + 1);
  uint32_t *cigar_start = (uint32_t *)(ez + hdr->swcount);
  for (int x=0; x<hdr->swcount; x++) {
    ez->cigar = cigar_start;
    cigar_start += ALIGN_BYTE_N(16, ez->n_cigar*4)>>2;
    {
	//fprintf(stderr, "ez->ncigar %d,max %d,zdrop %d,maxq %d,maxt %d,mqe %d,mqet %d,mte %d,mteq %d,sc %d,reach %d, revcig %x\n",ez->n_cigar,ez->max,ez->zdropped,ez->max_q,ez->max_t,ez->mqe,ez->mqe_t,ez->mte,ez->mte_q,ez->score,ez->reach_end,ez->revcigar);
    }
    ez++;
  }

  //debug
  static uint64_t totalsw = 0;
  totalsw += hdr->swcount;
  fprintf(stderr, "***************total sw %llu\n", totalsw);
}


void* task_recver(void *param)
{
    thread_ctrl_t *ctrl = (thread_ctrl_t*)param;
    vsc_ring_buf_t *dpring_buf = ctrl->dp_ring_buf;
    vsc_ring_buf_t *swring_buf = ctrl->sw_ring_buf;
    vsc_ring_buf_t *curr_ring = 0;
    int *dprcv_data = 0, *swrcv_data = 0;

    while(!ctrl->stop)
    {
        int size = 0;
        int type = RET_TYPE_SW;
        void *ret = fpga_get_retbuf(&size, RET_TYPE_SW);
        if (0 == ret || 0 == size) {
            continue;
        }

        uint64_t *dest = (uint64_t *)kmalloc(0, size+8/*8B is for km ptr*/);
        if (0 == dest) {fpga_release_retbuf(ret); continue;}
        memcpy(dest, ret, size);
        fpga_release_retbuf(ret);

        switch(type) 
        {
        case RET_TYPE_CD:
        case RET_TYPE_CS:
            curr_ring = ctrl->dp_ring_buf;
            break;
        case RET_TYPE_SW:
            update_value(dest);
            curr_ring = ctrl->sw_ring_buf;
            break;
        default:
            break;
        }

        {
            int idx = -1;
            do {//loop dead until
                idx = put_to_ring_buf(curr_ring, (void *)dest);
            }while(idx < 0);
        }
    }

    return param;
}



