#include "stat.h"
#include "fpga_sim.h"




#include <stdio.h>
#include <stdint.h>

static uint64_t min_factnum = 0;
static uint64_t max_factnum = 0;
static uint64_t avg_factnum = 0;

static uint64_t min_memsize = 0;
static uint64_t max_memsize = 0;
static uint64_t avg_memsize = 0;

static uint64_t callnums = 0;
static uint64_t nullnums = 0;


void set_swstat(int factnum, int memsize)
{
    callnums++;
    if (0 == factnum) {nullnums++;return;}

    if (0 == min_memsize) min_memsize = memsize;
    else if (memsize < min_memsize) min_memsize = memsize;
    if (0 == max_memsize) max_memsize = memsize;
    else if (memsize > max_memsize) max_memsize = memsize;
    avg_memsize += memsize;

    if (0 == min_factnum) min_factnum = factnum;
    else if (factnum < min_factnum) min_factnum = factnum;
    if (0 == max_factnum) max_factnum = factnum;
    else if (factnum > max_factnum) max_factnum = factnum;
    avg_factnum += factnum;
}


void out_swstat(void)
{
    if (0 == callnums) callnums = 1;
    fprintf(stderr, "SW STAT: total callnum %lu, nullnum %lu, max tasknum %lu, min %lu, avg %lu; max mem %lu, min %lu, avg %lu\n", callnums, nullnums, max_factnum, min_factnum, avg_factnum/callnums, max_memsize, min_memsize, avg_memsize/callnums);
}

static uint64_t dp_callnums = 0;
static uint64_t dp_min_memsize = 0;
static uint64_t dp_max_memsize = 0;
static uint64_t dp_avg_memsize = 0;

void set_dpstat(int memsize)
{
    dp_callnums++;
    if (0 == dp_min_memsize) dp_min_memsize = memsize;
    else if (memsize < dp_min_memsize) dp_min_memsize = memsize;
    if (0 == dp_max_memsize) dp_max_memsize = memsize;
    else if (memsize > dp_max_memsize) dp_max_memsize = memsize;
    dp_avg_memsize += memsize;
}


void out_dpstat(void)
{
    if(dp_callnums != 0)
        fprintf(stderr, "DP STAT: batch %u, max mem %lu, min %lu, avg %lu\n", DP_BATCHSIZE, dp_max_memsize, dp_min_memsize, dp_avg_memsize/dp_callnums);
}
