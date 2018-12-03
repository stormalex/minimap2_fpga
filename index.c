#include <stdlib.h>
#include <assert.h>
#if defined(WIN32) || defined(_WIN32)
#include <io.h> // for open(2)
#else
#include <unistd.h>
#endif
#include <fcntl.h>
#include <stdio.h>
#define __STDC_LIMIT_MACROS
#include "kthread.h"
#include "bseq.h"
#include "minimap.h"
#include "mmpriv.h"
#include "kvec.h"
#include "khash.h"

#define idx_hash(a) ((a)>>1)
#define idx_eq(a, b) ((a)>>1 == (b)>>1)
KHASH_INIT(idx, uint64_t, uint64_t, 1, idx_hash, idx_eq)
typedef khash_t(idx) idxhash_t;

KHASH_MAP_INIT_STR(str, uint32_t)

#define kroundup64(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, (x)|=(x)>>32, ++(x))

#define BUFFSIZE (16 * 1024 * 1024)
extern char txtbuff[BUFFSIZE];
extern int ntxt;
extern char * gentxt(uint8_t *dat, int datlen);
extern FILE *fidxcall;
extern FILE *fbreq, *fhreq, *fvreq, *fpreq;
extern FILE *fback, *fhack, *fvack, *fpack;

typedef struct _BBB {
    uint64_t h:36,hvalid:1, padding:27;
    uint64_t p:36,n_buckets:28;
} BBB;
int szb,szh,szp,szv;
int write64(uint8_t *buff, int n, int nall, FILE* fp)
{
    int newn = nall + n;
    char s[64];
    for (int i = 0; i < n; i++) {
        fprintf(fp, "%02x", buff[i]);
    }
    if (newn % 64 == 0) {
        fwrite("\n", 1, 1, fp);
    }
    return newn;
}

int fit4k(int nall, FILE* fp)
{
    uint8_t zero[64] = {0};
    if (nall % 64 != 0) {
        int n = (nall+63)/64*64 - nall;
        nall = write64((uint8_t *)zero, n, nall, fp);
    }
    if (nall % 4096 != 0) {
        int n = (nall+4095)/4096*4096 - nall;
        int m = n / 64;
        for (int i = 0; i < m; i++) {
            nall = write64((uint8_t *)zero, 64, nall, fp);
        }
    }
    fprintf(stderr, "--- fit to %d\n", nall);
    return nall;
}


typedef struct mm_idx_bucket_s {
	mm128_v a;   // (minimizer, position) array
	int32_t n;   // size of the _p_ array
	uint64_t *p; // position array for minimizers appearing >1 times
	void *h;     // hash table indexing _p_ and minimizers appearing once
    uint64_t hoff, poff;
} mm_idx_bucket_t;

mm_idx_t *mm_idx_init(int w, int k, int b, int flag)
{
	mm_idx_t *mi;
	if (k*2 < b) b = k * 2;
	if (w < 1) w = 1;
	mi = (mm_idx_t*)calloc(1, sizeof(mm_idx_t));
	mi->w = w, mi->k = k, mi->b = b, mi->flag = flag;
	mi->B = (mm_idx_bucket_t*)calloc(1<<b, sizeof(mm_idx_bucket_t));
	if (!(mm_dbg_flag & 1)) mi->km = km_init();
	return mi;
}

void mm_idx_destroy(mm_idx_t *mi)
{
	int i;
	if (mi == 0) return;
	if (mi->h) kh_destroy(str, (khash_t(str)*)mi->h);
	for (i = 0; i < 1<<mi->b; ++i) {
		free(mi->B[i].p);
		free(mi->B[i].a.a);
		kh_destroy(idx, (idxhash_t*)mi->B[i].h);
	}
	if (!mi->km) {
		for (i = 0; i < mi->n_seq; ++i)
			free(mi->seq[i].name);
		free(mi->seq);
	} else km_destroy(mi->km);
	free(mi->B); free(mi->S); free(mi);
}

uint64_t gret;
const uint64_t *mm_idx_get(const mm_idx_t *mi, uint64_t minier, int *n)
{
	int mask = (1<<mi->b) - 1;
	khint_t k;
    uint64_t breq,hreq,preq,vreq;
    uint64_t back,hack,pack,vack;
    breq=hreq=preq=vreq=back=hack=pack=vack=-1;
	mm_idx_bucket_t *b = &mi->B[minier&mask];
    breq = minier&mask;
    fprintf(fbreq, "%09lx\n", breq);
	idxhash_t *h = (idxhash_t*)b->h;
	*n = 0;
	if (h == 0) {
        uint64_t B[2] = {0};
        B[1] = b->hoff << 28 | b->poff >> 8;
        B[0] =(((b->poff&0xFF)<<56) | 0<<24);
        fprintf(fidxcall, "breq %lx\nhoff %lx\npoff %lx\nnh %lx\n", breq, b->hoff, b->poff, 0UL);
        fprintf(fback, "%016lx\n", B[1]);
        fprintf(fback, "%016lx\n", B[0]);
        return 0;
    } else {
        uint64_t n_buckets = kh_end(h);
        uint64_t B[2] = {0};
        B[1] = b->hoff << 28 | b->poff >> 8;
        B[0] =(((b->poff&0xFF)<<56) | n_buckets<<24);
        fprintf(fidxcall, "breq %lx\nhoff %lx\npoff %lx\nnh %lx\n", breq, b->hoff, b->poff, n_buckets);
        fprintf(fback, "%016lx\n", B[1]);
        fprintf(fback, "%016lx\n", B[0]);
    }
	k = kh_get(idx, h, minier>>mi->b<<1);
    uint64_t key = minier>>mi->b<<1;
    if (h->n_buckets) {
        khint_t k1, i, last, mask, step = 0;
        mask = h->n_buckets - 1;
        k1 = ((key)>>1); i = k1 & mask;
        last = i;
        int getk = 0;
        uint64_t h2[2] = {0};
        hreq = b->hoff + i;
        h2[1] = (uint64_t)h->flags[i>>4] << 32 | (uint32_t)(kh_key(h, i) >> 16);
        h2[0] = (kh_key(h, i) & 0xffff) << 48;
        fprintf(fhreq, "%09lx\n", hreq);
        fprintf(fhack, "%016lx\n", h2[1]);
        fprintf(fhack, "%016lx\n", h2[0]);
        while (!((h->flags[i>>4]>>((i&0xfU)<<1))&2) && (((h->flags[i>>4]>>((i&0xfU)<<1))&1) || !((h->keys[i])>>1 == (key)>>1))) {
            i = (i + (++step)) & mask;
            if (i == last) {
                k = h->n_buckets;
                getk = 1;
                break;
            }
            uint64_t h2[2] = {0};
            hreq = b->hoff + i;
            h2[1] = (uint64_t)h->flags[i>>4] << 32 | (uint32_t)(kh_key(h, i) >> 16);
            h2[0] = (kh_key(h, i) & 0xffff) << 48;
            fprintf(fhreq, "%09lx\n", hreq);
            fprintf(fhack, "%016lx\n", h2[1]);
            fprintf(fhack, "%016lx\n", h2[0]);
        }
        if (getk == 0) {
            k = ((h->flags[i>>4]>>((i&0xfU)<<1))&3)? h->n_buckets : i;
        }
    } else
        k = 0;
    hreq = minier>>mi->b<<1;
    hack = k;
    fprintf(fidxcall, "hreq %lx\nhack %lx\n", hreq, hack);
	if (k == kh_end(h)) return 0;
	if (kh_key(h, k)&1) { // special casing when there is only one k-mer
		*n = 1;
        uint64_t *v = &kh_val(h, k);
        vreq = b->hoff + k;
        vack = *v;
        fprintf(fidxcall, "vreq %lx\nvack %lx\n", vreq, vack);
        fprintf(fvreq, "%09lx\n", vreq);
        fprintf(fvack, "%016lx\n", vack);
        return v;
	} else {
		*n = (uint32_t)kh_val(h, k);
        vreq = b->hoff + k;
        vack = *n;
        fprintf(fidxcall, "vreq %lx\nvack %lx\n", vreq, vack);
        fprintf(fvreq, "%09lx\n", vreq);
        fprintf(fvack, "%016lx\n", vack);
		uint64_t *p = &b->p[kh_val(h, k)>>32];
        uint64_t poff = b->poff + (kh_val(h, k)>>32);
        preq = poff;
        fprintf(fidxcall, "preq %lx\n", preq);
        fprintf(fpreq, "%016lx\n", poff << 28 | (kh_val(h, k) & 0xffffffff));
        for (int i = 0; i < b->n; i++) {
            pack = p[i];
            fprintf(fidxcall, "pack %lx\n", pack);
            fprintf(fpack, "%016lx\n", pack);
        }
        return p;
	}
}

void mm_idx_stat(const mm_idx_t *mi)
{
	int i, n = 0, n1 = 0;
	uint64_t sum = 0, len = 0;
	fprintf(stderr, "[M::%s] kmer size: %d; skip: %d; is_hpc: %d; #seq: %d\n", __func__, mi->k, mi->w, mi->flag&MM_I_HPC, mi->n_seq);
	for (i = 0; i < mi->n_seq; ++i)
		len += mi->seq[i].len;
	for (i = 0; i < 1<<mi->b; ++i)
		if (mi->B[i].h) n += kh_size((idxhash_t*)mi->B[i].h);
	for (i = 0; i < 1<<mi->b; ++i) {
		idxhash_t *h = (idxhash_t*)mi->B[i].h;
		khint_t k;
		if (h == 0) continue;
		for (k = 0; k < kh_end(h); ++k)
			if (kh_exist(h, k)) {
				sum += kh_key(h, k)&1? 1 : (uint32_t)kh_val(h, k);
				if (kh_key(h, k)&1) ++n1;
			}
	}
	fprintf(stderr, "[M::%s::%.3f*%.2f] distinct minimizers: %d (%.2f%% are singletons); average occurrences: %.3lf; average spacing: %.3lf\n",
			__func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0), n, 100.0*n1/n, (double)sum / n, (double)len / sum);
}

int mm_idx_index_name(mm_idx_t *mi)
{
	khash_t(str) *h;
	uint32_t i;
	int has_dup = 0, absent;
	if (mi->h) return 0;
	h = kh_init(str);
	for (i = 0; i < mi->n_seq; ++i) {
		khint_t k;
		k = kh_put(str, h, mi->seq[i].name, &absent);
		if (absent) kh_val(h, k) = i;
		else has_dup = 1;
	}
	mi->h = h;
	if (has_dup && mm_verbose >= 2)
		fprintf(stderr, "[WARNING] some database sequences have identical sequence names\n");
	return has_dup;
}

int mm_idx_name2id(const mm_idx_t *mi, const char *name)
{
	khash_t(str) *h = (khash_t(str)*)mi->h;
	khint_t k;
	if (h == 0) return -2;
	k = kh_get(str, h, name);
	return k == kh_end(h)? -1 : kh_val(h, k);
}

int mm_idx_getseq(const mm_idx_t *mi, uint32_t rid, uint32_t st, uint32_t en, uint8_t *seq)
{
	uint64_t i, st1, en1;
	if (rid >= mi->n_seq || st >= mi->seq[rid].len) return -1;
	if (en > mi->seq[rid].len) en = mi->seq[rid].len;
	st1 = mi->seq[rid].offset + st;
	en1 = mi->seq[rid].offset + en;
	for (i = st1; i < en1; ++i)
		seq[i - st1] = mm_seq4_get(mi->S, i);
	return en - st;
}

int32_t mm_idx_cal_max_occ(const mm_idx_t *mi, float f)
{
	int i;
	size_t n = 0;
	uint32_t thres;
	khint_t *a, k;
	if (f <= 0.) return INT32_MAX;
	for (i = 0; i < 1<<mi->b; ++i)
		if (mi->B[i].h) n += kh_size((idxhash_t*)mi->B[i].h);
	a = (uint32_t*)malloc(n * 4);
	for (i = n = 0; i < 1<<mi->b; ++i) {
		idxhash_t *h = (idxhash_t*)mi->B[i].h;
		if (h == 0) continue;
		for (k = 0; k < kh_end(h); ++k) {
			if (!kh_exist(h, k)) continue;
			a[n++] = kh_key(h, k)&1? 1 : (uint32_t)kh_val(h, k);
		}
	}
	thres = ks_ksmall_uint32_t(n, a, (uint32_t)((1. - f) * n)) + 1;
	free(a);
	return thres;
}

int cmp(rname_rid_t *r0, rname_rid_t *r1)
{
    return strcmp(r0->rname, r1->rname);
}
/*********************************
 * Sort and generate hash tables *
 *********************************/

static void worker_post(void *g, long i, int tid)
{
	int j, start_a, start_p, n, n_keys;
	idxhash_t *h;
	mm_idx_t *mi = (mm_idx_t*)g;
	mm_idx_bucket_t *b = &mi->B[i];
	if (b->a.n == 0) return;

	// sort by minimizer
	radix_sort_128x(b->a.a, b->a.a + b->a.n);

	// count and preallocate
	for (j = 1, n = 1, n_keys = 0, b->n = 0; j <= b->a.n; ++j) {
		if (j == b->a.n || b->a.a[j].x>>8 != b->a.a[j-1].x>>8) {
			++n_keys;
			if (n > 1) b->n += n;
			n = 1;
		} else ++n;
	}
	h = kh_init(idx);
	kh_resize(idx, h, n_keys);
	b->p = (uint64_t*)calloc(b->n, 8);

	//modified by LQX

	// create the hash table
	/*
	n=1 val ----63...43|42...22| 21   |20...0
				 refid |refpos |strand|rankid
	n>1 p[k]----63...43|42...22| 21   |20...0
				 refid |refpos |strand|rankid
	*/
	for (j = 1, n = 1, start_a = start_p = 0; j <= b->a.n; ++j) {
		if (j == b->a.n || b->a.a[j].x>>8 != b->a.a[j-1].x>>8) {
			khint_t itr;
			int absent;
			mm128_t *p = &b->a.a[j-1];
			itr = kh_put(idx, h, p->x>>8>>mi->b<<1, &absent);
			assert(absent && j - start_a == n);
			if (n == 1) {
				kh_key(h, itr) |= 1;
				//uint64_t val = p->y;
				uint64_t refid = (p->y)>>32;
				uint64_t refpos_strand = (p->y)&0xFFFFFFFF;
				uint32_t rankid = mi->rever_rid[refid];
				uint64_t val = ((refid&0x1FFFFF)<<43)|((refpos_strand&0x3FFFFF)<<21)|((rankid&0x1FFFFF));
				kh_val(h, itr) = val;
				//fprintf(stderr,"finished work_post_n_1 \n");
			} else {
				int k;
				for (k = 0; k < n; ++k)
					b->p[start_p + k] = b->a.a[start_a + k].y;
				radix_sort_64(&b->p[start_p], &b->p[start_p + n]); // sort by position; needed as in-place radix_sort_128x() is not stable
				for (k = 0; k < n; ++k){
					uint64_t refid = (b->p[start_p+k])>>32;
					uint64_t refpos_strand = (b->p[start_p+k])&0xFFFFFFFF;
					uint32_t rankid = mi->rever_rid[refid];
					uint64_t pk = ((refid&0x1FFFFF)<<43)|((refpos_strand&0x3FFFFF)<<21)|((rankid&0x1FFFFF));
					b->p[start_p+k] = pk;
				}
				kh_val(h, itr) = (uint64_t)start_p<<32 | n;
				start_p += n;
				//fprintf(stderr,"finished work_post\n");
				//modified by LQX

			}
			start_a = j, n = 1;
		} else ++n;
	}
	b->h = h;
	assert(b->n == start_p);

	// deallocate and clear b->a
	kfree(0, b->a.a);
	b->a.n = b->a.m = 0, b->a.a = 0;
	//fprintf(stderr,"finished work_post\n");
}

static void mm_idx_post(mm_idx_t *mi, int n_threads)
{
	kt_for(n_threads, worker_post, mi, 1<<mi->b);
}

/******************
 * Generate index *
 ******************/

#include <string.h>
#include <zlib.h>
#include "bseq.h"

typedef struct {
	int mini_batch_size;
	uint64_t batch_size, sum_len;
	mm_bseq_file_t *fp;
	mm_idx_t *mi;
} pipeline_t;

typedef struct {
    int n_seq;
	mm_bseq1_t *seq;
	mm128_v a;
} step_t;

static void mm_idx_add(mm_idx_t *mi, int n, const mm128_t *a)
{
	int i, mask = (1<<mi->b) - 1;
	for (i = 0; i < n; ++i) {
		mm128_v *p = &mi->B[a[i].x>>8&mask].a;
		kv_push(mm128_t, 0, *p, a[i]);
	}
}

static void *worker_pipeline(void *shared, int step, void *in)
{
	int i;
    pipeline_t *p = (pipeline_t*)shared;
    if (step == 0) { // step 0: read sequences
        step_t *s;
		if (p->sum_len > p->batch_size) return 0;
        s = (step_t*)calloc(1, sizeof(step_t));
		s->seq = mm_bseq_read(p->fp, p->mini_batch_size, 0, &s->n_seq); // read a mini-batch
		if (s->seq) {
			uint32_t old_m, m;
			assert((uint64_t)p->mi->n_seq + s->n_seq <= UINT32_MAX); // to prevent integer overflow
			// make room for p->mi->seq
			old_m = p->mi->n_seq, m = p->mi->n_seq + s->n_seq;
			kroundup32(m); kroundup32(old_m);
			if (old_m != m)
				p->mi->seq = (mm_idx_seq_t*)krealloc(p->mi->km, p->mi->seq, m * sizeof(mm_idx_seq_t));
			// make room for p->mi->S
			if (!(p->mi->flag & MM_I_NO_SEQ)) {
				uint64_t sum_len, old_max_len, max_len;
				for (i = 0, sum_len = 0; i < s->n_seq; ++i) sum_len += s->seq[i].l_seq;
				old_max_len = (p->sum_len + 7) / 8;
				max_len = (p->sum_len + sum_len + 7) / 8;
				kroundup64(old_max_len); kroundup64(max_len);
				if (old_max_len != max_len) {
					p->mi->S = (uint32_t*)realloc(p->mi->S, max_len * 4);
					memset(&p->mi->S[old_max_len], 0, 4 * (max_len - old_max_len));
				}
			}
			// populate p->mi->seq
			for (i = 0; i < s->n_seq; ++i) {
				mm_idx_seq_t *seq = &p->mi->seq[p->mi->n_seq];
				uint32_t j;
				if (!(p->mi->flag & MM_I_NO_NAME)) {
					seq->name = (char*)kmalloc(p->mi->km, strlen(s->seq[i].name) + 1);
					strcpy(seq->name, s->seq[i].name);
				} else seq->name = 0;
				seq->len = s->seq[i].l_seq;
				seq->offset = p->sum_len;
				// copy the sequence
				if (!(p->mi->flag & MM_I_NO_SEQ)) {
					for (j = 0; j < seq->len; ++j) { // TODO: this is not the fastest way, but let's first see if speed matters here
						uint64_t o = p->sum_len + j;
						int c = seq_nt4_table[(uint8_t)s->seq[i].seq[j]];
						mm_seq4_set(p->mi->S, o, c);
					}
				}
				// update p->sum_len and p->mi->n_seq
				p->sum_len += seq->len;
				s->seq[i].rid = p->mi->n_seq++;
			}
			return s;
		} else free(s);
    } else if (step == 1) { // step 1: compute sketch
        step_t *s = (step_t*)in;
		for (i = 0; i < s->n_seq; ++i) {
			mm_bseq1_t *t = &s->seq[i];
			if (t->l_seq > 0)
				mm_sketch(0, t->seq, t->l_seq, p->mi->w, p->mi->k, t->rid, p->mi->flag&MM_I_HPC, &s->a);
			else if (mm_verbose >= 2)
				fprintf(stderr, "[WARNING] the length database sequence '%s' is 0\n", t->name);
			free(t->seq); free(t->name);
		}
		free(s->seq); s->seq = 0;
		return s;
    } else if (step == 2) { // dispatch sketch to buckets
        step_t *s = (step_t*)in;
		mm_idx_add(p->mi, s->a.n, s->a.a);
		kfree(0, s->a.a); free(s);
	}
    return 0;
}

mm_idx_t *mm_idx_gen(mm_bseq_file_t *fp, int w, int k, int b, int flag, int mini_batch_size, int n_threads, uint64_t batch_size)
{
	pipeline_t pl;
	if (fp == 0 || mm_bseq_eof(fp)) return 0;
	memset(&pl, 0, sizeof(pipeline_t));
	pl.mini_batch_size = mini_batch_size < batch_size? mini_batch_size : batch_size;
	pl.batch_size = batch_size;
	pl.fp = fp;
	pl.mi = mm_idx_init(w, k, b, flag);

	kt_pipeline(n_threads < 3? n_threads : 3, worker_pipeline, &pl, 3);
	if (mm_verbose >= 3)
		fprintf(stderr, "[M::%s::%.3f*%.2f] collected minimizers\n", __func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0));
	#if 1
	{
		mm_idx_t *mi = pl.mi;
		mi->rname_rid = (rname_rid_t *)calloc(mi->n_seq, sizeof(rname_rid_t));
		mi->rever_rid = (int *)calloc(mi->n_seq, sizeof(int));

		//fprintf(stderr, "---------------------\n");
		for (int k=0; k<mi->n_seq; k++)
			{
				const mm_idx_seq_t *s = &mi->seq[k];
				strncpy(mi->rname_rid[k].rname, s->name, 255);
				mi->rname_rid[k].rname[strlen(s->name)] = 0;
				mi->rname_rid[k].rid = k;
			}

		qsort(mi->rname_rid, mi->n_seq, sizeof(rname_rid_t), cmp);
		/*
		fprintf(stderr, "---------------------\n");
		for (int k=0; k<mi->n_seq; k++)
			{
				fprintf(stderr, "%s\n", mi->rname_rid[k].rname);
			}
		*/
		for (int k=0; k<mi->n_seq; k++)
			{
				mi->rever_rid[mi->rname_rid[k].rid] = k;
			}
		/*
		for (int k=0; k<mi->n_seq; k++)
			{
				fprintf(stderr, "table,%d:%d\n", k, mi->rever_rid[k]);
			}
		*/
	}
	#endif
	mm_idx_post(pl.mi, n_threads);
	if (mm_verbose >= 3)
		fprintf(stderr, "[M::%s::%.3f*%.2f] sorted minimizers\n", __func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0));

    uint8_t buff[1024];
    mm_idx_t *mi = pl.mi;
    int i, j;
    FILE *fidxb, *fidxh, *fidxp, *fidxv;
    FILE *fidxbt, *fidxht, *fidxpt, *fidxvt;
    if ((fidxb = fopen("idxb.txt", "w")) == NULL) {
        perror("open data file failed!");
        exit(1);
    }
    if ((fidxh = fopen("idxh.txt", "w")) == NULL) {
        perror("open data file failed!");
        exit(1);
    }
    if ((fidxp = fopen("idxp.txt", "w")) == NULL) {
        perror("open data file failed!");
        exit(1);
    }
    if ((fidxv = fopen("idxv.txt", "w")) == NULL) {
        perror("open data file failed!");
        exit(1);
    }
    if ((fidxbt = fopen("idxbt.txt", "w")) == NULL) {
        perror("open data file failed!");
        exit(1);
    }

    //gen idx data
    uint64_t allh = 0;
    uint64_t allp = 0;
	uint32_t loop_num = 0;
    for(i = 0;i < (1<<mi->b);++i){
        kh_idx_t* h = (idxhash_t*)mi->B[i].h;
        uint64_t B[2] = {0};
        if (mi->B[i].h) {
            uint64_t values;
            //Get B Array
            int n_buckets;
            n_buckets = kh_end((idxhash_t*)mi->B[i].h);
            B[1] = allh << 28 | allp >> 8;
            B[0] =(((allp&0xFF)<<56) | n_buckets<<24);
            fprintf(fidxbt, "%lx %lx %x\n", allh, allp, n_buckets);
            szb = write64((uint8_t *)B, 16, szb, fidxb);

            //set B
            mi->B[i].hoff = allh;
            mi->B[i].poff = allp;

            allh += n_buckets;
            allp += mi->B[i].n;
            //将 n_buckets h p 合并到128bit中，实际是拆分两个64bit 顺序分别为h-36bit|p-28bit|p-8bit|n_buckets-32bit|padding
            for(int j = 0;j < n_buckets;++j){
                uint32_t flags;
                uint64_t keys;
                flags = h->flags[j>>4];
                keys = h->keys[j];
                uint64_t a[2]={0};

                //Get H Array
                a[1] |= flags;
                a[1] = (a[1] << 32) |((keys&0xffffffffffffUL)>>16);
                a[0] |= keys&0xffff;
                a[0] = a[0]<<48;
                szh = write64((uint8_t *)a, 16, szh, fidxh);

                //Get Values Array
                values = h->vals[j];//value是后期可能要优化为48bit    21-bit|22-bit|padding
                szv = write64((uint8_t *)&values, 8, szv, fidxv);
				loop_num++;
            }

            //Get p Array
            int n = mi->B[i].n;
            for(int j = 0;j < n;++j){
                values = mi->B[i].p[j];
                szp = write64((uint8_t *)&values, 8, szp, fidxp);
            }
        } else {
            fprintf(fidxbt, "%lx %lx %x\n", 0L, 0L, 0);
            szb = write64((uint8_t *)B, 16, szb, fidxb);
        }
    }

    //align to 4KB
    fit4k(szb, fidxb);
    fit4k(szh, fidxh);
    fit4k(szp, fidxp);
    fit4k(szv, fidxv);

    fclose(fidxb);
    fclose(fidxh);
    fclose(fidxp);
    fclose(fidxv);

    fclose(fidxbt);

	return pl.mi;
}

mm_idx_t *mm_idx_build(const char *fn, int w, int k, int flag, int n_threads) // a simpler interface; deprecated
{
	mm_bseq_file_t *fp;
	mm_idx_t *mi;
	fp = mm_bseq_open(fn);
	if (fp == 0) return 0;
	mi = mm_idx_gen(fp, w, k, 14, flag, 1<<18, n_threads, UINT64_MAX);
	mm_bseq_close(fp);
	return mi;
}

mm_idx_t *mm_idx_str(int w, int k, int is_hpc, int bucket_bits, int n, const char **seq, const char **name)
{
	uint64_t sum_len = 0;
	mm128_v a = {0,0,0};
	mm_idx_t *mi;
	int i, flag = 0;
	if (n <= 0) return 0;
	for (i = 0; i < n; ++i) // get the total length
		sum_len += strlen(seq[i]);
	if (is_hpc) flag |= MM_I_HPC;
	if (name == 0) flag |= MM_I_NO_NAME;
	if (bucket_bits < 0) bucket_bits = 14;
	mi = mm_idx_init(w, k, bucket_bits, flag);
	mi->n_seq = n;
	mi->seq = (mm_idx_seq_t*)kcalloc(mi->km, n, sizeof(mm_idx_seq_t)); // ->seq is allocated from km
	mi->S = (uint32_t*)calloc((sum_len + 7) / 8, 4);
	for (i = 0, sum_len = 0; i < n; ++i) {
		const char *s = seq[i];
		mm_idx_seq_t *p = &mi->seq[i];
		uint32_t j;
		if (name && name[i]) {
			p->name = (char*)kmalloc(mi->km, strlen(name[i]) + 1);
			strcpy(p->name, name[i]);
		}
		p->offset = sum_len;
		p->len = strlen(s);
		for (j = 0; j < p->len; ++j) {
			int c = seq_nt4_table[(uint8_t)s[j]];
			uint64_t o = sum_len + j;
			mm_seq4_set(mi->S, o, c);
		}
		sum_len += p->len;
		if (p->len > 0) {
			a.n = 0;
			mm_sketch(0, s, p->len, w, k, i, is_hpc, &a);
			mm_idx_add(mi, a.n, a.a);
		}
	}
	free(a.a);
	mm_idx_post(mi, 1);
	return mi;
}

/*************
 * index I/O *
 *************/

void mm_idx_dump(FILE *fp, const mm_idx_t *mi)
{
	uint64_t sum_len = 0;
	uint32_t x[5];
	int i;

	x[0] = mi->w, x[1] = mi->k, x[2] = mi->b, x[3] = mi->n_seq, x[4] = mi->flag;
	fwrite(MM_IDX_MAGIC, 1, 4, fp);
	fwrite(x, 4, 5, fp);
	for (i = 0; i < mi->n_seq; ++i) {
		uint8_t l;
		l = strlen(mi->seq[i].name);
		fwrite(&l, 1, 1, fp);
		fwrite(mi->seq[i].name, 1, l, fp);
		fwrite(&mi->seq[i].len, 4, 1, fp);
		sum_len += mi->seq[i].len;
	}
	for (i = 0; i < 1<<mi->b; ++i) {
		mm_idx_bucket_t *b = &mi->B[i];
		khint_t k;
		idxhash_t *h = (idxhash_t*)b->h;
		uint32_t size = h? h->size : 0;
		fwrite(&b->n, 4, 1, fp);
		fwrite(b->p, 8, b->n, fp);
		fwrite(&size, 4, 1, fp);
		if (size == 0) continue;
		for (k = 0; k < kh_end(h); ++k) {
			uint64_t x[2];
			if (!kh_exist(h, k)) continue;
			x[0] = kh_key(h, k), x[1] = kh_val(h, k);
			fwrite(x, 8, 2, fp);
		}
	}
	if (!(mi->flag & MM_I_NO_SEQ))
		fwrite(mi->S, 4, (sum_len + 7) / 8, fp);
	fflush(fp);
}

mm_idx_t *mm_idx_load(FILE *fp)
{
	int i;
	char magic[4];
	uint32_t x[5];
	uint64_t sum_len = 0;
	mm_idx_t *mi;

	if (fread(magic, 1, 4, fp) != 4) return 0;
	if (strncmp(magic, MM_IDX_MAGIC, 4) != 0) return 0;
	if (fread(x, 4, 5, fp) != 5) return 0;
	mi = mm_idx_init(x[0], x[1], x[2], x[4]);
	mi->n_seq = x[3];
	mi->seq = (mm_idx_seq_t*)kcalloc(mi->km, mi->n_seq, sizeof(mm_idx_seq_t));
	for (i = 0; i < mi->n_seq; ++i) {
		uint8_t l;
		mm_idx_seq_t *s = &mi->seq[i];
		fread(&l, 1, 1, fp);
		s->name = (char*)kmalloc(mi->km, l + 1);
		fread(s->name, 1, l, fp);
		s->name[l] = 0;
		fread(&s->len, 4, 1, fp);
		s->offset = sum_len;
		sum_len += s->len;
	}
	for (i = 0; i < 1<<mi->b; ++i) {
		mm_idx_bucket_t *b = &mi->B[i];
		uint32_t j, size;
		khint_t k;
		idxhash_t *h;
		fread(&b->n, 4, 1, fp);
		b->p = (uint64_t*)malloc(b->n * 8);
		fread(b->p, 8, b->n, fp);
		fread(&size, 4, 1, fp);
		if (size == 0) continue;
		b->h = h = kh_init(idx);
		kh_resize(idx, h, size);
		for (j = 0; j < size; ++j) {
			uint64_t x[2];
			int absent;
			fread(x, 8, 2, fp);
			k = kh_put(idx, h, x[0], &absent);
			assert(absent);
			kh_val(h, k) = x[1];
		}
	}
	if (!(mi->flag & MM_I_NO_SEQ)) {
		mi->S = (uint32_t*)malloc((sum_len + 7) / 8 * 4);
		fread(mi->S, 4, (sum_len + 7) / 8, fp);
	}
	return mi;
}

int64_t mm_idx_is_idx(const char *fn)
{
	int fd, is_idx = 0;
	off_t ret, off_end;
	char magic[4];

	if (strcmp(fn, "-") == 0) return 0; // read from pipe; not an index
	fd = open(fn, O_RDONLY);
	if (fd < 0) return -1; // error
	if ((off_end = lseek(fd, 0, SEEK_END)) >= 4) {
		lseek(fd, 0, SEEK_SET);
		ret = read(fd, magic, 4);
		if (ret == 4 && strncmp(magic, MM_IDX_MAGIC, 4) == 0)
			is_idx = 1;
	}
	close(fd);
	return is_idx? off_end : 0;
}

mm_idx_reader_t *mm_idx_reader_open(const char *fn, const mm_idxopt_t *opt, const char *fn_out)
{
	int64_t is_idx;
	mm_idx_reader_t *r;
	is_idx = mm_idx_is_idx(fn);
	if (is_idx < 0) return 0; // failed to open the index
	r = (mm_idx_reader_t*)calloc(1, sizeof(mm_idx_reader_t));
	r->is_idx = is_idx;
	if (opt) r->opt = *opt;
	else mm_idxopt_init(&r->opt);
	if (r->is_idx) {
		r->fp.idx = fopen(fn, "rb");
		r->idx_size = is_idx;
	} else r->fp.seq = mm_bseq_open(fn);
	if (fn_out) r->fp_out = fopen(fn_out, "wb");
	return r;
}

void mm_idx_reader_close(mm_idx_reader_t *r)
{
	if (r->is_idx) fclose(r->fp.idx);
	else mm_bseq_close(r->fp.seq);
	if (r->fp_out) fclose(r->fp_out);
	free(r);
}

mm_idx_t *mm_idx_reader_read(mm_idx_reader_t *r, int n_threads)
{
	mm_idx_t *mi;
	if (r->is_idx) {
		mi = mm_idx_load(r->fp.idx);
		if (mi && mm_verbose >= 2 && (mi->k != r->opt.k || mi->w != r->opt.w || (mi->flag&MM_I_HPC) != (r->opt.flag&MM_I_HPC)))
			fprintf(stderr, "[WARNING]\033[1;31m Indexing parameters (-k, -w or -H) overridden by parameters used in the prebuilt index.\033[0m\n");
	} else
		mi = mm_idx_gen(r->fp.seq, r->opt.w, r->opt.k, r->opt.bucket_bits, r->opt.flag, r->opt.mini_batch_size, n_threads, r->opt.batch_size);
	if (mi) {
		if (r->fp_out) mm_idx_dump(r->fp_out, mi);
		++r->n_parts;
	}
	return mi;
}

int mm_idx_reader_eof(const mm_idx_reader_t *r) // TODO: in extremely rare cases, mm_bseq_eof() might not work
{
	return r->is_idx? (feof(r->fp.idx) || ftell(r->fp.idx) == r->idx_size) : mm_bseq_eof(r->fp.seq);
}
