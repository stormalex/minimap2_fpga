#include <stdint.h>
#include <string.h>
#include <stdio.h>

#include "minimap.h"
#include "mmpriv.h"
#include "kalloc.h"
#include "ksw2.h"

#include "kfifo.h"
#include "context_data.h"
#include "fpga_sim.h"



static const char LogTable256[256] = {
#define LT(n) n, n, n, n, n, n, n, n, n, n, n, n, n, n, n, n
	-1, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
	LT(4), LT(5), LT(5), LT(6), LT(6), LT(6), LT(6),
	LT(7), LT(7), LT(7), LT(7), LT(7), LT(7), LT(7), LT(7)
};

static inline int ilog2_32(uint32_t v)
{
	register uint32_t t, tt;
	if ((tt = v>>16)) return (t = tt>>8) ? 24 + LogTable256[t] : 16 + LogTable256[tt];
	return (t = v>>8) ? 8 + LogTable256[t] : LogTable256[v];
}


void dumpfile(int memsize, void *dest, char *fpname, char *drname)
{
    char txtname[64] = {0};
    char hexname[64] = {0};
    unsigned char r64[64] ={0};
    unsigned char d64[129]={0};
    int left = memsize;
    
    uint8_t *d = (uint8_t *)dest;

    snprintf(txtname, sizeof(txtname)-1, "%s.txt", fpname);
    snprintf(hexname, sizeof(hexname)-1, "%s.bin", drname);

    FILE *fp = fopen(txtname, "w");
    for (int k=0; k<memsize;) {
        if (left < 64) {
            for (int x=0; x<left; x++) r64[x]=d[k+left-1-x];
            k += left;
            for (int x=0; x<left; x++)
                snprintf(d64+2*x, 3, "%02x", r64[x]);
            
            fwrite(d64, 2*left, 1, fp);
            fwrite("\r\n",2,1,fp);
        } else {
            for (int x=0; x<64; x++) r64[x]=d[k+63-x];
            left -= 64;
            k+= 64;
            for (int x=0; x<64; x++)
                snprintf(d64+2*x, 3, "%02x", r64[x]);
            
            fwrite(d64, 2*64, 1, fp);
            fwrite("\r\n",2,1,fp);
        }
    }
    //for (int x=0; x<64; x++) fprintf(stderr, "%02x", d[x]);            
    fclose(fp);
    
    FILE *fd = fopen(hexname, "wb");
    //fwrite(&memsize, 4, 1, fd);
    for (int x=0; x<memsize; x++) fwrite(d+x, 1, 1, fd);
    fclose(fd);
}


int64_t fpgasim_chain_dp(int max_dist_x, int max_dist_y, int bw, int max_skip, int min_cnt, int min_sc, int is_cdna, int n_segs, int64_t n, mm128_t *a, int *n_u_, uint64_t **_u, int32_t *n_v_, int32_t **_v, void *km)
{ // TODO: make sure this works when n has more than 32 bits
	int32_t k, *f, *p, *t, *v, n_u, n_v;
	int64_t i, j, st = 0;
	uint64_t *u, /* *u2,*/ sum_qspan = 0;
	float avg_qspan;
	//mm128_t *b, *w;

	if (_u) *_u = 0, *n_u_ = 0;
	f = (int32_t*)kmalloc(km, n * 4);
	p = (int32_t*)kmalloc(km, n * 4);
	t = (int32_t*)kmalloc(km, n * 4);
	v = (int32_t*)kmalloc(km, n * 4);
	memset(t, 0, n * 4);

	for (i = 0; i < n; ++i) sum_qspan += a[i].y>>32&0xff;
	avg_qspan = (float)sum_qspan / n;

	// fill the score and backtrack arrays
	for (i = 0; i < n; ++i) {
		uint64_t ri = a[i].x;
		int64_t max_j = -1;
		int32_t qi = (int32_t)a[i].y, q_span = a[i].y>>32&0xff; // NB: only 8 bits of span is used!!!
		int32_t max_f = q_span, n_skip = 0, min_d;
		int32_t sidi = (a[i].y & MM_SEED_SEG_MASK) >> MM_SEED_SEG_SHIFT;
		while (st < i && ri - a[st].x > max_dist_x) ++st;
		for (j = i - 1; j >= st; --j) {
			int64_t dr = ri - a[j].x;
			int32_t dq = qi - (int32_t)a[j].y, dd, sc, log_dd;
			int32_t sidj = (a[j].y & MM_SEED_SEG_MASK) >> MM_SEED_SEG_SHIFT;
			if ((sidi == sidj && dr == 0) || dq <= 0) continue; // don't skip if an anchor is used by multiple segments; see below
			if ((sidi == sidj && dq > max_dist_y) || dq > max_dist_x) continue;
			dd = dr > dq? dr - dq : dq - dr;
			if (sidi == sidj && dd > bw) continue;
			if (n_segs > 1 && !is_cdna && sidi == sidj && dr > max_dist_y) continue;
			min_d = dq < dr? dq : dr;
			sc = min_d > q_span? q_span : dq < dr? dq : dr;
			log_dd = dd? ilog2_32(dd) : 0;
			if (is_cdna || sidi != sidj) {
				int c_log, c_lin;
				c_lin = (int)(dd * .01 * avg_qspan);
				c_log = log_dd;
				if (sidi != sidj && dr == 0) ++sc; // possibly due to overlapping paired ends; give a minor bonus
				else if (dr > dq || sidi != sidj) sc -= c_lin < c_log? c_lin : c_log;
				else sc -= c_lin + (c_log>>1);
			} else sc -= (int)(dd * .01 * avg_qspan) + (log_dd>>1);
			sc += f[j];
			if (sc > max_f) {
				max_f = sc, max_j = j;
				if (n_skip > 0) --n_skip;
			} else if (t[j] == i) {
				if (++n_skip > max_skip)
					break;
			}
			if (p[j] >= 0) t[p[j]] = i;
		}
		f[i] = max_f, p[i] = max_j;
		v[i] = max_j >= 0 && v[max_j] > max_f? v[max_j] : max_f; // v[] keeps the peak score up to i; f[] is the score ending at i, not always the peak
	}

	// find the ending positions of chains
	memset(t, 0, n * 4);
	for (i = 0; i < n; ++i)
		if (p[i] >= 0) t[p[i]] = 1;
	for (i = n_u = 0; i < n; ++i)
		if (t[i] == 0 && v[i] >= min_sc)
			++n_u;
	if (n_u == 0) {
		kfree(km, a); kfree(km, f); kfree(km, p); kfree(km, t); kfree(km, v);
		return 0;
	}
	u = (uint64_t*)kmalloc(km, n_u * 8);
	for (i = n_u = 0; i < n; ++i) {
		if (t[i] == 0 && v[i] >= min_sc) {
			j = i;
			while (j >= 0 && f[j] < v[j]) j = p[j]; // find the peak that maximizes f[]
			if (j < 0) j = i; // TODO: this should really be assert(j>=0)
			u[n_u++] = (uint64_t)f[j] << 32 | j;
		}
	}
	radix_sort_64(u, u + n_u);
	for (i = 0; i < n_u>>1; ++i) { // reverse, s.t. the highest scoring chain is the first
		uint64_t t = u[i];
		u[i] = u[n_u - i - 1], u[n_u - i - 1] = t;
	}

	// backtrack
	memset(t, 0, n * 4);
	for (i = n_v = k = 0; i < n_u; ++i) { // starting from the highest score
		int32_t n_v0 = n_v, k0 = k;
		j = (int32_t)u[i];
		do {
			v[n_v++] = j;
			t[j] = 1;
			j = p[j];
		} while (j >= 0 && t[j] == 0);
		if (j < 0) {
			if (n_v - n_v0 >= min_cnt) u[k++] = u[i]>>32<<32 | (n_v - n_v0);
		} else if ((int32_t)(u[i]>>32) - f[j] >= min_sc) {
			if (n_v - n_v0 >= min_cnt) u[k++] = ((u[i]>>32) - f[j]) << 32 | (n_v - n_v0);
		}
		if (k0 == k) n_v = n_v0; // no new chain added, reset
	}

	*n_v_ = n_v, *_v = v;
	*n_u_ = n_u = k, *_u = u; // NB: note that u[] may not be sorted by score here

	// free temporary arrays
	kfree(km, f); kfree(km, p); kfree(km, t);

	return 1;
}

void debug_kmalloc(void *km, void *addr,int size, char *func, int line)
{
  size = 0;
  //fprintf(stderr, "call kmalloc in [%s:%d], size %u->%p:%p\n", func,line,size,km,addr);
}

void debug_kfree(void *km, void *addr,int size, char *func, int line)
{
  //fprintf(stderr, "call kfee in [%s:%d], size %d->%p:%p\n", func,line,size,km,addr);
}

int32_t simsc_chain_dp(int max_dist_x, int max_dist_y, int n_segs, int64_t n, mm128_t *a, void *km, int32_t *f, pv_t *pv)
{ // TODO: make sure this works when n has more than 32 bits
	int32_t k;
	int32_t n_u;
	int64_t i, j, st = 0;
	uint64_t sum_qspan = 0;
	float avg_qspan;

	int bw = 500;
	int max_skip = 25;
	int min_cnt = 3;
	int min_sc = 40;
	int is_cdna = 0;

    int32_t *t = (int32_t*)kmalloc(km, n * 4);
	memset(t, 0, n * 4);

	for (i = 0; i < n; ++i) sum_qspan += a[i].y>>32&0xff;
	avg_qspan = (float)sum_qspan / n;

	// fill the score and backtrack arrays
	for (i = 0; i < n; ++i) {
		uint64_t ri = a[i].x;
		int64_t max_j = -1;
		int32_t qi = (int32_t)a[i].y, q_span = a[i].y>>32&0xff; // NB: only 8 bits of span is used!!!
		int32_t max_f = q_span, n_skip = 0, min_d;
		int32_t sidi = (a[i].y & MM_SEED_SEG_MASK) >> MM_SEED_SEG_SHIFT;
		while (st < i && ri - a[st].x > max_dist_x) ++st;
		for (j = i - 1; j >= st; --j) {
			int64_t dr = ri - a[j].x;
			int32_t dq = qi - (int32_t)a[j].y, dd, sc, log_dd;
			int32_t sidj = (a[j].y & MM_SEED_SEG_MASK) >> MM_SEED_SEG_SHIFT;
			if ((sidi == sidj && dr == 0) || dq <= 0) continue; // don't skip if an anchor is used by multiple segments; see below
			if ((sidi == sidj && dq > max_dist_y) || dq > max_dist_x) continue;
			dd = dr > dq? dr - dq : dq - dr;
			if (sidi == sidj && dd > bw) continue;
			if (n_segs > 1 && !is_cdna && sidi == sidj && dr > max_dist_y) continue;
			min_d = dq < dr? dq : dr;
			sc = min_d > q_span? q_span : dq < dr? dq : dr;
			log_dd = dd? ilog2_32(dd) : 0;
			if (is_cdna || sidi != sidj) {
				int c_log, c_lin;
				c_lin = (int)(dd * .01 * avg_qspan);
				c_log = log_dd;
				if (sidi != sidj && dr == 0) ++sc; // possibly due to overlapping paired ends; give a minor bonus
				else if (dr > dq || sidi != sidj) sc -= c_lin < c_log? c_lin : c_log;
				else sc -= c_lin + (c_log>>1);
			} else sc -= (int)(dd * .01 * avg_qspan) + (log_dd>>1);
			sc += f[j];
			if (sc > max_f) {
				max_f = sc, max_j = j;
				if (n_skip > 0) --n_skip;
			} else if (t[j] == i) {
				if (++n_skip > max_skip)
					break;
			}
			if ((pv[j].p<<1) >= 0) t[pv[j].p<<1>>1] = i;
		}
		f[i] = max_f, pv[i].p = max_j&0x7fffffff;
		pv[i].v = max_j >= 0 && pv[max_j].v > max_f? pv[max_j].v : max_f; // v[] keeps the peak score up to i; f[] is the score ending at i, not always the peak
        if (max_j >= 0) pv[max_j].p |= 0x80000000;
	}
	// find the ending positions of chains
	memset(t, 0, n * 4);
	//for (i = 0; i < n; ++i)
	//	if (pv->p[i] >= 0) t[pv->p[i]] = 1;
	for (i = n_u = 0; i < n; ++i)
		if ((pv[i].p>>31) == 0 && pv[i].v >= min_sc)
			++n_u;
    pv[n].p = n_u, pv[n].v = n_u;

    kfree(km, t);

	return n_u;
}

mm128_t * simall_chain_dp(int max_dist_x, int max_dist_y, int n_segs, int64_t n, mm128_t *a, int32_t *outk, int *n_u_, uint64_t **_u, void *km,void *sdkm)
{ // TODO: make sure this works when n has more than 32 bits
	int32_t k, *f, *p, *t, *v, n_u, n_v;
	int64_t i, j, st = 0;
	uint64_t *u, *u2, sum_qspan = 0;
	float avg_qspan;
	mm128_t *b, *w;


	int bw = 500;
	int max_skip = 25;
	int min_cnt = 3;
	int min_sc = 40;
	int is_cdna = 0;


	if (_u) *_u = 0, *n_u_ = 0;
	f = (int32_t*)kmalloc(km, n * 4);
	p = (int32_t*)kmalloc(km, n * 4);
	t = (int32_t*)kmalloc(km, n * 4);
	v = (int32_t*)kmalloc(km, n * 4);
	memset(t, 0, n * 4);

	for (i = 0; i < n; ++i) sum_qspan += a[i].y>>32&0xff;
	avg_qspan = (float)sum_qspan / n;

	// fill the score and backtrack arrays
	for (i = 0; i < n; ++i) {
		uint64_t ri = a[i].x;
		int64_t max_j = -1;
		int32_t qi = (int32_t)a[i].y, q_span = a[i].y>>32&0xff; // NB: only 8 bits of span is used!!!
		int32_t max_f = q_span, n_skip = 0, min_d;
		int32_t sidi = (a[i].y & MM_SEED_SEG_MASK) >> MM_SEED_SEG_SHIFT;
		while (st < i && ri - a[st].x > max_dist_x) ++st;
		for (j = i - 1; j >= st; --j) {
			int64_t dr = ri - a[j].x;
			int32_t dq = qi - (int32_t)a[j].y, dd, sc, log_dd;
			int32_t sidj = (a[j].y & MM_SEED_SEG_MASK) >> MM_SEED_SEG_SHIFT;
			if ((sidi == sidj && dr == 0) || dq <= 0) continue; // don't skip if an anchor is used by multiple segments; see below
			if ((sidi == sidj && dq > max_dist_y) || dq > max_dist_x) continue;
			dd = dr > dq? dr - dq : dq - dr;
			if (sidi == sidj && dd > bw) continue;
			if (n_segs > 1 && !is_cdna && sidi == sidj && dr > max_dist_y) continue;
			min_d = dq < dr? dq : dr;
			sc = min_d > q_span? q_span : dq < dr? dq : dr;
			log_dd = dd? ilog2_32(dd) : 0;
			if (is_cdna || sidi != sidj) {
				int c_log, c_lin;
				c_lin = (int)(dd * .01 * avg_qspan);
				c_log = log_dd;
				if (sidi != sidj && dr == 0) ++sc; // possibly due to overlapping paired ends; give a minor bonus
				else if (dr > dq || sidi != sidj) sc -= c_lin < c_log? c_lin : c_log;
				else sc -= c_lin + (c_log>>1);
			} else sc -= (int)(dd * .01 * avg_qspan) + (log_dd>>1);
			sc += f[j];
			if (sc > max_f) {
				max_f = sc, max_j = j;
				if (n_skip > 0) --n_skip;
			} else if (t[j] == i) {
				if (++n_skip > max_skip)
					break;
			}
			if (p[j] >= 0) t[p[j]] = i;
		}
		f[i] = max_f, p[i] = max_j;
		v[i] = max_j >= 0 && v[max_j] > max_f? v[max_j] : max_f; // v[] keeps the peak score up to i; f[] is the score ending at i, not always the peak
	}

	// find the ending positions of chains
	memset(t, 0, n * 4);
	for (i = 0; i < n; ++i)
		if (p[i] >= 0) t[p[i]] = 1;
	for (i = n_u = 0; i < n; ++i)
		if (t[i] == 0 && v[i] >= min_sc)
			++n_u;
	if (n_u == 0) {
	  /*kfree(sdkm, a);*/ kfree(km, f); kfree(km, p); kfree(km, t); kfree(km, v);
	  return 0;
	}
	u = (uint64_t*)kmalloc(km, n_u * 8);
	for (i = n_u = 0; i < n; ++i) {
		if (t[i] == 0 && v[i] >= min_sc) {
			j = i;
			while (j >= 0 && f[j] < v[j]) j = p[j]; // find the peak that maximizes f[]
			if (j < 0) j = i; // TODO: this should really be assert(j>=0)
			u[n_u++] = (uint64_t)f[j] << 32 | j;
		}
	}
	radix_sort_64(u, u + n_u);
	for (i = 0; i < n_u>>1; ++i) { // reverse, s.t. the highest scoring chain is the first
		uint64_t t = u[i];
		u[i] = u[n_u - i - 1], u[n_u - i - 1] = t;
	}

	// backtrack
	memset(t, 0, n * 4);
	for (i = n_v = k = 0; i < n_u; ++i) { // starting from the highest score
		int32_t n_v0 = n_v, k0 = k;
		j = (int32_t)u[i];
		do {
			v[n_v++] = j;
			t[j] = 1;
			j = p[j];
		} while (j >= 0 && t[j] == 0);
		if (j < 0) {
			if (n_v - n_v0 >= min_cnt) u[k++] = u[i]>>32<<32 | (n_v - n_v0);
		} else if ((int32_t)(u[i]>>32) - f[j] >= min_sc) {
			if (n_v - n_v0 >= min_cnt) u[k++] = ((u[i]>>32) - f[j]) << 32 | (n_v - n_v0);
		}
		if (k0 == k) n_v = n_v0; // no new chain added, reset
	}
	*n_u_ = n_u = k, *_u = u; // NB: note that u[] may not be sorted by score here

	// free temporary arrays
	kfree(km, f); kfree(km, p); kfree(km, t);

	// write the result to b[]
	b = (mm128_t*)kmalloc(km, n_v * sizeof(mm128_t));
	for (i = 0, k = 0; i < n_u; ++i) {
		int32_t k0 = k, ni = (int32_t)u[i];
		for (j = 0; j < ni; ++j)
			b[k] = a[v[k0 + (ni - j - 1)]], ++k;
	}
	kfree(km, v);

	// sort u[] and a[] by a[].x, such that adjacent chains may be joined (required by mm_join_long)
	w = (mm128_t*)kmalloc(km, n_u * sizeof(mm128_t));
	for (i = k = 0; i < n_u; ++i) {
		w[i].x = b[k].x, w[i].y = (uint64_t)k<<32|i;
		k += (int32_t)u[i];
	}
	radix_sort_128x(w, w + n_u);
	u2 = (uint64_t*)kmalloc(km, n_u * 8);
	for (i = k = 0; i < n_u; ++i) {
		int32_t j = (int32_t)w[i].y, n = (int32_t)u[j];

		u2[i] = u[j];
		memcpy(&a[k], &b[w[i].y>>32], n * sizeof(mm128_t));
		k += n;
	}
	memcpy(u, u2, n_u * 8);
	memcpy(b, a, k * sizeof(mm128_t)); // write _a_ to _b_ and deallocate _a_ because _a_ is oversized, sometimes a lot
	/*kfree(sdkm, a);*/ kfree(km, w); kfree(km, u2);

	*outk = k;

	return b;
}


void calc_chaindp_batch(void *addr, void *km, const mm_idx_t *mi)
{
  int oldsize = sizeof(dptest_rcvhdr_t);
  chaindp_sndhdr_t *hdr = (chaindp_sndhdr_t *)addr;

  dptest_rcvhdr_t *ret = (dptest_rcvhdr_t*)kmalloc(km, oldsize);
  ret->km  = km;
  ret->size= oldsize;
  ret->type= hdr->type;
  ret->tid = hdr->tid;
  ret->num = hdr->num;
  ret->lat = hdr->lat;
  ret->magic= hdr->magic;

  uint32_t preoff = 0;
  for (int x=0; x<hdr->num; x++) {
    int n_u,rep_len,n_mini_pos;
    int64_t n_a;
    uint64_t *mini_pos;

    chaindp_to_fpga_t *data = (chaindp_to_fpga_t*)((char *)hdr + hdr->subhdr[x].offset);
    int qlen_sum = data->qlen_sum;
    mm_mapopt_t opt;
    opt.flag = data->optflag;

    mm128_v mv;
    mv.n = data->seednum;
    mv.a = (mm128_t*)data->seed;

    mm128_t *a = collect_seed_hits(km, &opt, data->mid_occ, mi, data->qname, &mv, qlen_sum, &n_a, &rep_len, &n_mini_pos, &mini_pos);

    uint32_t subsize = sizeof(dptest_rcvsubhdr_t) + ((n_a<<5)-(n_a<<2)+8)/* 3*int/pvf+mm128_t/a+8/n_u */ + (n_mini_pos<<3)/*mini_pos is uint64*/;
    oldsize += subsize;
    ret = (dptest_rcvhdr_t *)krealloc(km, ret, oldsize);

    dptest_rcvsubhdr_t *sub = (dptest_rcvsubhdr_t *)((uint8_t *)ret->subhdr + preoff);
    preoff += subsize;

    sub->n_a = n_a;
    sub->rep_len = rep_len;
    sub->n_mini_pos = n_mini_pos;
    sub->aoffset = (n_a<<3)+(n_a<<2)+8;
    sub->moffset = sub->aoffset+(n_a<<4);
    sub->subsize = subsize;
    sub->ctxpos  = data->pos;

    pv_t *pv = (pv_t*)sub->pvf;
    int32_t *f = (int32_t *)(pv + n_a + 1);

    n_u = simsc_chain_dp(data->gap_ref,data->gap_qry,data->n_segs,n_a,a,km,f,pv);

    if (n_u != 0) {
      mm128_t *da = (mm128_t *)(sub->pvf + sub->aoffset);
      memcpy(da, a, n_a<<4);

      uint64_t *dm = (uint64_t *)(sub->pvf + sub->moffset);
      memcpy(dm, mini_pos, n_mini_pos<<3);
    }

    ret->size += oldsize;

    //free a  b u...
    kfree(km, a);
    kfree(km, mini_pos);
  }

  //put to rcv ring
  int idx = -1;
  do {
    idx = put_to_ring_buf(rcv_ctrl.dp_ring_buf, ret);
  } while(idx < 0);

}


void calc_sw(void *dest, void *km, const void *topt)
{
  //int k=0, a=2,b=4,q=4,e=2,q2=24,e2=1;
  const mm_mapopt_t *opt = (const mm_mapopt_t *)topt;

  int8_t mat[25] = {0};
  ksw_gen_simple_mat(5, mat, opt->a, opt->b);

  sw_sndhdr_t *hdr = (sw_sndhdr_t *)dest;

  int hdrsize = sizeof(swtest_rcvhdr_t);
  //if (hdr->num > 0) 
  {
      for (int k=0; k<hdr->num; k++)
      {
          sw_readhdr_t *read = hdr->data + k;
          //if (read->regnum != 0)
          {
              sw_reghdr_t *cr = &read->reg;
              hdrsize += (!!get_left(cr->midnum) + !!get_right(cr->midnum) + clear_midnum(cr->midnum))*sizeof(ksw_extz_t);
          }
      }
  }
  
  swtest_rcvhdr_t *ret = (swtest_rcvhdr_t *)kmalloc(km, hdrsize);
  ret->km   = km;
  ret->magic= hdr->magic;
  ret->size = 67890;
  ret->type = hdr->type;
  ret->tid  = hdr->tid;
  ret->num  = hdr->num;
  ret->lat  = hdr->lat;
  ret->freecigar = FREE_CIGAR_X86;

  for (int k=0; k<hdr->num; k++) {
    sw_readhdr_t *rd = hdr->data + k;
    sw_reghdr_t  *cr = &rd->reg;
    sw_to_fpga_t *dt = (sw_to_fpga_t *)((uint8_t *)hdr + cr->offset);

    sw_readhdr_t *dr = ret->data + k;
    sw_reghdr_t  *rr = &dr->reg;
    rr->midnum       = cr->midnum;
    rr->regpos       = cr->regpos;
    //rr->hasleft      = cr->hasleft;
    //rr->hasright     = cr->hasright;

    //dr->regnum       = rd->regnum;
    dr->ctxpos       = rd->ctxpos;
    dr->longsw       = rd->longsw;

    //int total = cr->hasleft + cr->midnum + cr->hasright;//total segment of curr reg
    rr->size = /*rd->regnum * */ (!!get_left(rr->midnum) + !!get_right(rr->midnum) + clear_midnum(rr->midnum)) * sizeof(ksw_extz_t);
    if (0 == k) rr->offset = sizeof(swtest_rcvhdr_t);
    else rr->offset = ret->data[k-1].reg.offset + ret->data[k-1].reg.size;

    //if (0 == rd->regnum) continue;//should be continue after calc size&offset, or will bug when last regnum=0

    ksw_extz_t *ez = (ksw_extz_t *)((uint8_t *)ret + rr->offset);
    memset(ez, 0, (!!get_left(rr->midnum) + !!get_right(rr->midnum) + clear_midnum(rr->midnum))*sizeof(ksw_extz_t));

    const uint8_t *qseq = 0;
    const uint8_t *tseq = 0;

    if (get_left(cr->midnum)) {
        qseq = (uint8_t *)dt + sizeof(sw_to_fpga_t);
        tseq = (uint8_t *)dt + sizeof(sw_to_fpga_t) + dt->qlen_align;

        //ksw_extd2_sse(km, dt->qlen, qseq, dt->tlen, tseq, 5, mat, q, e, q2, e2, dt->bw, dt->zdrop, dt->end_bonus, dt->flag, ez);
        //fprintf(stderr, "left before align, qlen %d, tlen %d\n", dt->qlen,dt->tlen);
        mm_align_pair(km, opt, dt->qlen, qseq, dt->tlen, tseq, mat, dt->bw, dt->end_bonus, dt->zdrop, dt->flag, ez);
        //fprintf(stderr, "left,qlen %d,tlen %d,bw %d,endbos %d,zdrop %d,flag %x,ez max %d,zdrp %d,maxq %d,maxt %d,mqe %d,mqet %d,score %d,ncig %d,reaend %d\n",dt->qlen,dt->tlen,dt->bw,dt->end_bonus,dt->zdrop,dt->flag,ez->max,ez->zdropped,ez->max_q,ez->max_t,ez->mqe,ez->mqe_t,ez->score,ez->n_cigar,ez->reach_end);
        //int x=0;fprintf(stderr, "\nsimleft,ez %p,cigar %p:",ez,ez->cigar);for(;x<ez->n_cigar;x++)fprintf(stderr, " %02x,",ez->cigar[x]);fprintf(stderr, "\n");
        dt = (sw_to_fpga_t *)((uint8_t *)dt + sizeof(sw_to_fpga_t) + dt->qlen_align + dt->tlen_align);
        ez++;
    }
    int crmidnum = clear_midnum(cr->midnum);
    for (int x = 0; x < crmidnum; x++) {
        qseq = (uint8_t *)dt + sizeof(sw_to_fpga_t);
        tseq = (uint8_t *)dt + sizeof(sw_to_fpga_t) + dt->qlen_align;
        //ksw_extd2_sse(km, dt->qlen, qseq, dt->tlen, tseq, 5, mat, q, e, q2, e2, dt->bw, dt->zdrop, dt->end_bonus, dt->flag|KSW_EZ_APPROX_MAX, ez);
        mm_align_pair(km, opt, dt->qlen, qseq, dt->tlen, tseq, mat, dt->bw, dt->end_bonus, dt->zdrop, dt->flag, ez);
        //fprintf(stderr, "mid,qlen %d,tlen %d,bw %d,endbos %d,zdrop %d,flag %x,ez max %d,zdrp %d,maxq %d,maxt %d,mqe %x,mqet %d,score %x,ncig %d,reaend %d\n",dt->qlen,dt->tlen,dt->bw,dt->end_bonus,dt->zdrop,dt->flag|KSW_EZ_APPROX_MAX,ez->max,ez->zdropped,ez->max_q,ez->max_t,ez->mqe,ez->mqe_t,ez->score,ez->n_cigar,ez->reach_end);
#if 0
        int zdrop_code;
        if ((zdrop_code = mm_test_zdrop(km, opt, qseq, tseq, ez->n_cigar, ez->cigar, mat)) != 0) {
            //ksw_extd2_sse(km, dt->qlen, qseq, dt->tlen, tseq, 5, mat, q, e, q2, e2, dt->bw, zdrop_code == 2? dt->zdrop_inv : dt->zdrop, dt->end_bonus, dt->flag, fstez);
            mm_align_pair(km, opt, dt->qlen, qseq, dt->tlen, tseq, mat, dt->bw, dt->end_bonus, zdrop_code == 2? dt->zdrop_inv : dt->zdrop, dt->flag, ez);
            //fprintf(stderr, "mid zdrop,qlen %d,tlen %d,bw %d,endbos %d,zcode %d, zinv %d,zdrop %d,factzdrop %d,flag %x,ez max %d,zdrp %d,maxq %d,maxt %d,mqe %d,mqet %d,score %d,ncig %d,reaend %d\n",dt->qlen,dt->tlen,dt->bw,dt->end_bonus,zdrop_code,dt->zdrop_inv,dt->zdrop,zdrop_code == 2?dt->zdrop_inv : dt->zdrop,dt->flag,ez->max,ez->zdropped,ez->max_q,ez->max_t,ez->mqe,ez->mqe_t,ez->score,ez->n_cigar,ez->reach_end);
            
        }
        ez->m_cigar = zdrop_code;//save code to mcigar, as mcigar is free
#endif
        dt = (sw_to_fpga_t *)((uint8_t *)dt + sizeof(sw_to_fpga_t) + dt->qlen_align + dt->tlen_align);
        ez++;
    } 
    if (get_right(cr->midnum)) {
        qseq = (uint8_t *)dt + sizeof(sw_to_fpga_t);
        tseq = (uint8_t *)dt + sizeof(sw_to_fpga_t) + dt->qlen_align;
        //ksw_extd2_sse(km, dt->qlen, qseq, dt->tlen, tseq, 5, mat, q, e, q2, e2, dt->bw, dt->zdrop, dt->end_bonus, dt->flag, ez);
        mm_align_pair(km, opt, dt->qlen, qseq, dt->tlen, tseq, mat, dt->bw, dt->end_bonus, dt->zdrop, dt->flag, ez);
        //fprintf(stderr, "rig,qlen %d,tlen %d,bw %d,endbos %d,zdrop %d,flag %x,ez max %d,zdrp %d,maxq %d,maxt %d,mqe %d,mqet %d,score %d,ncig %d,reaend %d\n",dt->qlen,dt->tlen,dt->bw,dt->end_bonus,dt->zdrop,dt->flag,ez->max,ez->zdropped,ez->max_q,ez->max_t,ez->mqe,ez->mqe_t,ez->score,ez->n_cigar,ez->reach_end);

        dt = (sw_to_fpga_t *)((uint8_t *)dt + sizeof(sw_to_fpga_t) + dt->qlen_align + dt->tlen_align);
        ez++;
    }
  }//for regnum

#ifdef DUMP_FILES
  {
      //for debug, dump output memory
      int eznum = 0;
      int memsize = hdrsize;
      for (int k=0; k<hdr->num; k++) {
          sw_readhdr_t *dr = ret->data + k;
          sw_reghdr_t  *rr = &dr->reg;
          ksw_extz_t *ez = (ksw_extz_t *)((uint8_t *)ret + rr->offset);
          if (get_left(rr->midnum)) {
              memsize += ALIGN_BYTE_N(128, ez->n_cigar);
              ez++;
              eznum++;
          }
          int crmidnum = clear_midnum(rr->midnum);
          for (int x = 0; x < crmidnum; x++) {
              memsize += ALIGN_BYTE_N(128, ez->n_cigar);
              ez++;
              eznum++;
          }
          if (get_right(rr->midnum)) {
              memsize += ALIGN_BYTE_N(128, ez->n_cigar);
              ez++;
              eznum++;
          }
      }

      uint8_t *filemem = (uint8_t *)kmalloc(km, memsize);
      uint8_t *orig = filemem;
      memcpy(filemem, ret, sizeof(swtest_rcvhdr_t));
      filemem += sizeof(swtest_rcvhdr_t);
      memcpy(filemem, (uint8_t *)ret + sizeof(swtest_rcvhdr_t), eznum*sizeof(ksw_extz_t));
      filemem += eznum*sizeof(ksw_extz_t);

      for (int k=0; k<hdr->num; k++) {
          sw_readhdr_t *dr = ret->data + k;
          sw_reghdr_t  *rr = &dr->reg;
          ksw_extz_t *ez = (ksw_extz_t *)((uint8_t *)ret + rr->offset);

          if (dr->longsw) continue;
          if (get_left(rr->midnum)) {
              memcpy(filemem, ez->cigar, ez->n_cigar);
              filemem += ALIGN_BYTE_N(128, ez->n_cigar);
              ez++;
          }
          int crmidnum = clear_midnum(rr->midnum);
          for (int x = 0; x < crmidnum; x++) {
              memcpy(filemem, ez->cigar, ez->n_cigar);
              filemem += ALIGN_BYTE_N(128, ez->n_cigar);
              ez++;
          }
          if (get_right(rr->midnum)) {
              memcpy(filemem, ez->cigar, ez->n_cigar);
              filemem += ALIGN_BYTE_N(128, ez->n_cigar);
              ez++;
          }
      }//end for

      static int seq = 0;
      char fpn[64]={0};
      char drn[64]={0};
      snprintf(fpn, sizeof(fpn)-1, "swfpga%d-out",seq);
      snprintf(drn, sizeof(drn)-1, "swdriv%d-out",seq);
      seq++;
      dumpfile(memsize, orig, fpn, drn);
      kfree(km, orig);
  }
#endif
  //put to rcv ring
  int idx = -1;
  do {
    idx = put_to_ring_buf(rcv_ctrl.sw_ring_buf, ret);
  } while(idx < 0);

  //fprintf(stderr, "send sw back....ret %p,last=%d,tid=%d,num=%d\n",ret, ret->lat,ret->tid,ret->num);
}



#ifdef MY_FPGA_API

int fpga_init(int flag, int is_cdna, int maxskip, int minsc, int m, int longthres, int longdiff, int8_t mat0, int8_t mat1)
{
}

void fpga_finalize()
{
}

void* fpga_get_retbuf(int* len, uint8_t type)
{
}

void fpga_release_retbuf(void* addr)
{
}

void* fpga_get_writebuf(int size)
{
}

void fpga_writebuf_submit(void* addr, int size, uint8_t type)
{
}



#endif//MY_FPGA_API

