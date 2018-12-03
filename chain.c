#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include "minimap.h"
#include "mmpriv.h"
#include "kalloc.h"

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

struct new_seed
{
    mm128_t seed;
    int32_t p;
    int32_t f;
};

extern pthread_mutex_t counter_mutex;

typedef struct __attribute__((__packed__)) _FPGAHDR {
    uint32_t magic;
    uint32_t size;
    uint16_t tid;
    uint16_t num;
    uint8_t type;
    uint8_t lat;
    uint8_t padding[50];
} FPGAHDR;
typedef struct __attribute__((__packed__)) _DPHDR {
    uint32_t gap_ref;
    uint32_t gap_qry;
    uint32_t seednum;
    uint32_t qlensum;
    uint32_t ctxpos;
    uint32_t bid;
    uint16_t n_segs;
    uint8_t b;
    uint8_t padding[37];
} DPHDR;
typedef struct __attribute__((__packed__)) _ODPHDR {
    uint32_t err_flag;
    uint32_t ctxpos;
    uint32_t subsize;
    uint32_t n_a;
    uint32_t n_minipos;
    uint32_t rep_len;
    uint8_t padding[40];
} ODPHDR;

#define BUFFSIZE (16 * 1024 * 1024)
extern char txtbuff[BUFFSIZE];
extern int ntxt;
extern char * gentxt(uint8_t *dat, int datlen);
extern FILE *fcshi, *fcsho, *fdpi, *fdpo;
extern FILE *fbcshi, *fbcsho, *fbdpi, *fbdpo;

static uint8_t gbuff[BUFFSIZE];
static int gbufflen = 0;
static uint8_t gobuff[BUFFSIZE];
static int gobufflen = 0;
static int maxdp = 8;
static int maxcsh = 8;
static int gn = 0;
static uint32_t gmagic = 0;

mm128_t *mm_chain_dp(int max_dist_x, int max_dist_y, int bw, int max_skip, int min_cnt, int min_sc, int is_cdna, int n_segs, int64_t n, mm128_t *a, int *n_u_, uint64_t **_u, void *km)
{ // TODO: make sure this works when n has more than 32 bits
	int32_t k, *f, *p, *t, *v, n_u, n_v;
	int64_t i, j, st = 0;
	uint64_t *u, *u2, sum_qspan = 0;
	float avg_qspan;
	mm128_t *b, *w;
    struct new_seed* fpga_a = NULL;
    int32_t *fpga_id = NULL;

    // sim input
    int orglen;
    if (gn == 0) {
        memset(gbuff, 0, sizeof(FPGAHDR));
        memset(gobuff, 0, sizeof(FPGAHDR));
        FPGAHDR *a = (FPGAHDR *)gbuff;
        a->magic = gmagic;
        a->type = 1;
        a->tid = 0x66;
        a->lat = 0x77;
        gbufflen = sizeof(FPGAHDR);
        FPGAHDR *o = (FPGAHDR *)gobuff;
        o->magic = gmagic;
        o->type = 1;
        o->tid = 0x66;
        o->lat = 0x77;
        gobufflen = sizeof(FPGAHDR);
    }
    DPHDR *input = (DPHDR *)(gbuff + gbufflen);
    memset(input, 0, sizeof(DPHDR));
    input->gap_ref = max_dist_x;
    input->gap_qry = max_dist_y;
    input->seednum = n;
    input->ctxpos = 0x88;
    input->n_segs = n_segs;
    gbufflen += sizeof(DPHDR);
    memcpy(gbuff + gbufflen, a, n * 16);
    gbufflen += n * 16;
    if (gbufflen % 64 != 0) {
        gbufflen = (gbufflen >> 6 << 6) + 64;
    }
    assert(gbufflen < BUFFSIZE);
    FILE *fp;
    if ((fp = fopen("DP_IN.dat", "a")) != NULL) {
        fprintf(fp, "%x\n", max_dist_x);
        fprintf(fp, "%x\n", max_dist_y);
        fprintf(fp, "%lx\n", n);
        fprintf(fp, "%x\n", n_segs);
        fclose(fp);
    }
    if ((fp = fopen("DP_SEED.dat", "a")) != NULL) {
        fprintf(fp, "\n\n%x\n", gn);
        for (i = 0; i < n; i++) {
            fprintf(fp, "%016lx%016lx\n", a[i].x, a[i].y);
        }
        fclose(fp);
    }

	if (_u) *_u = 0, *n_u_ = 0;
	f = (int32_t*)kmalloc(km, n * 4);
	p = (int32_t*)kmalloc(km, n * 4);
	t = (int32_t*)kmalloc(km, n * 4);
	v = (int32_t*)kmalloc(km, n * 4);
    fpga_id = (int32_t*)kmalloc(km, n * sizeof(int32_t));
    fpga_a = (struct new_seed*)kmalloc(km, n * sizeof(struct new_seed));
	memset(t, 0, n * 4);
    //memset(fpga_id, 0xff, n * sizeof(int32_t));
    for(i = 0; i < n; i++)
        fpga_id[i] = -1;


	for (i = 0; i < n; ++i) sum_qspan += a[i].y>>32&0xff;
	avg_qspan = (float)sum_qspan / n;

    uint32_t new_i = 0;
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

        //
        if(p[i] >= 0) {                 //
            if(fpga_id[p[i]] != -1) {   //seed已经在fpga_a中了
                fpga_a[new_i].p = fpga_id[p[i]]<<2;
            }
            else {
                fpga_a[new_i].seed = a[p[i]];
                fpga_a[new_i].f = f[p[i]];
                fpga_id[p[i]] = new_i;

                fpga_a[new_i].p = (-1)<<2;
                fpga_a[new_i].p |= (v[p[i]] >= min_sc);
                fpga_a[new_i].p |= ((f[p[i]] < v[p[i]]) << 1);

                new_i++;
            }
        }
        if((v[i] >= min_sc) || (p[i] >= 0)) {
            fpga_a[new_i].seed = a[i];
            fpga_a[new_i].f = f[i];

            fpga_id[i] = new_i;

            if(p[i]>=0) fpga_a[new_i].p = fpga_id[p[i]]<<2;
            else fpga_a[new_i].p = (-1)<<2;

            fpga_a[new_i].p |= (v[i] >= min_sc);
            fpga_a[new_i].p |= (f[i] < v[i]) << 1;

            new_i++;
        }
    }

    // sim output
    orglen = gobufflen;
    ODPHDR *output = (ODPHDR *)(gobuff + gobufflen);
    memset(output, 0, sizeof(ODPHDR));
    output->ctxpos = 0x88;
    output->n_a = new_i;
    output->n_minipos = 0;
    gobufflen += sizeof(ODPHDR);
    memcpy(gobuff + gobufflen, fpga_a, new_i * 24);
    gobufflen += new_i * 24;
    if (gobufflen % 64 != 0) {
        gobufflen = (gobufflen >> 6 << 6) + 64;
    }
    output->subsize = gobufflen - orglen;
    assert(gobufflen < BUFFSIZE);

    gn++;
    if (gn == maxdp) {
        assert(gbufflen % 64 == 0);
        assert(gobufflen % 64 == 0);

        //set nread,size
        FPGAHDR *a = (FPGAHDR *)gbuff;
        a->size = gbufflen;
        a->num = gn;
        FPGAHDR *o = (FPGAHDR *)gobuff;
        o->size = gobufflen;
        o->num = gn;

        //write sim
        gentxt(gbuff, gbufflen);
        fwrite(txtbuff, 1, ntxt, fdpi);
        gentxt(gobuff, gobufflen);
        fwrite(txtbuff, 1, ntxt, fdpo);

        //write
        fwrite(gbuff, 1, gbufflen, fbdpi);
        fwrite(gobuff, 1, gobufflen, fbdpo);

        gn = 0;
        gmagic++;
    }

    //lvjingbang
    uint64_t counter = 0;

	// find the ending positions of chains
	memset(t, 0, n * 4);
	for (i = 0; i < new_i; ++i)
		//if (fpga_a[i].p>>2 >= 0) t[fpga_a[i].p>>2] = 1;
        if (fpga_a[i].p >= 0) t[fpga_a[i].p>>2] = 1;
	for (i = n_u = 0; i < new_i; ++i)
		//if (t[i] == 0 && v[i] >= min_sc)
        if (((fpga_a[i].p & 0x01) == 1) && (t[i] == 0))
			++n_u;
	if (n_u == 0) {
		kfree(km, a); kfree(km, f); kfree(km, p); kfree(km, t); kfree(km, v);
        kfree(km, fpga_a);
        kfree(km, fpga_id);
		return 0;
	}
    kfree(km, f);
    kfree(km, p);

	u = (uint64_t*)kmalloc(km, n_u * 8);
	for (i = n_u = 0; i < new_i; ++i) {
		//if (t[i] == 0 && v[i] >= min_sc) {
        if (((fpga_a[i].p & 0x01) == 1) && (t[i] == 0)) {
			j = i;
			//while (j >= 0 && f[j] < v[j]) j = p[j]; // find the peak that maximizes f[]
            while (j >= 0 && (fpga_a[j].p & 2)) j = fpga_a[j].p>>2;
			if (j < 0) j = i; // TODO: this should really be assert(j>=0)
			u[n_u++] = (uint64_t)fpga_a[j].f << 32 | j;
		}
	}
	radix_sort_64(u, u + n_u);  //TODO 修改为从大到小
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
			//j = p[j];
            j = fpga_a[j].p>>2;
		} while (j >= 0 && t[j] == 0);
		if (j < 0) {
			if (n_v - n_v0 >= min_cnt) u[k++] = u[i]>>32<<32 | (n_v - n_v0);
		} else if ((int32_t)(u[i]>>32) - fpga_a[j].f >= min_sc) {
			if (n_v - n_v0 >= min_cnt) u[k++] = ((u[i]>>32) - fpga_a[j].f) << 32 | (n_v - n_v0);
		}
		if (k0 == k) n_v = n_v0; // no new chain added, reset
	}
	*n_u_ = n_u = k, *_u = u; // NB: note that u[] may not be sorted by score here

	// free temporary arrays
	//kfree(km, f); kfree(km, p);
    kfree(km, t);

	// write the result to b[]
	b = (mm128_t*)kmalloc(km, n_v * sizeof(mm128_t));
	for (i = 0, k = 0; i < n_u; ++i) {
		int32_t k0 = k, ni = (int32_t)u[i];
		for (j = 0; j < ni; ++j)
			b[k] = fpga_a[v[k0 + (ni - j - 1)]].seed, ++k;
	}
	kfree(km, v);

    kfree(km, fpga_a);
    kfree(km, fpga_id);

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
	kfree(km, a); kfree(km, w); kfree(km, u2);



	return b;
}
