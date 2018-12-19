#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include "bseq.h"
#include "minimap.h"
#include "mmpriv.h"
#ifdef HAVE_GETOPT
#include <getopt.h>
#else
#include "getopt.h"
#endif

#include <pthread.h>
#include "soft_sw.h"
#include "fpga_sw.h"

#define MM_VERSION "2.10-r761"

#ifdef __linux__
#include <sys/resource.h>
#include <sys/time.h>
void liftrlimit()
{
	struct rlimit r;
	getrlimit(RLIMIT_AS, &r);
	r.rlim_cur = r.rlim_max;
	setrlimit(RLIMIT_AS, &r);
}
#else
void liftrlimit() {}
#endif

static struct option long_options[] = {
	{ "bucket-bits",    required_argument, 0, 0 },
	{ "mb-size",        required_argument, 0, 'K' },
	{ "seed",           required_argument, 0, 0 },
	{ "no-kalloc",      no_argument,       0, 0 },
	{ "print-qname",    no_argument,       0, 0 },
	{ "no-self",        no_argument,       0, 'D' },
	{ "print-seeds",    no_argument,       0, 0 },
	{ "max-chain-skip", required_argument, 0, 0 },
	{ "min-dp-len",     required_argument, 0, 0 },
	{ "print-aln-seq",  no_argument,       0, 0 },
	{ "splice",         no_argument,       0, 0 },
	{ "cost-non-gt-ag", required_argument, 0, 'C' },
	{ "no-long-join",   no_argument,       0, 0 },
	{ "sr",             no_argument,       0, 0 },
	{ "frag",           required_argument, 0, 0 },
	{ "secondary",      required_argument, 0, 0 },
	{ "cs",             optional_argument, 0, 0 },
	{ "end-bonus",      required_argument, 0, 0 },
	{ "no-pairing",     no_argument,       0, 0 },
	{ "splice-flank",   required_argument, 0, 0 },
	{ "idx-no-seq",     no_argument,       0, 0 },
	{ "end-seed-pen",   required_argument, 0, 0 },   // 21
	{ "for-only",       no_argument,       0, 0 },   // 22
	{ "rev-only",       no_argument,       0, 0 },   // 23
	{ "heap-sort",      required_argument, 0, 0 },   // 24
	{ "all-chain",      no_argument,       0, 'P' },
	{ "dual",           required_argument, 0, 0 },   // 26
	{ "max-clip-ratio", required_argument, 0, 0 },   // 27
	{ "min-occ-floor",  required_argument, 0, 0 },   // 28
	{ "MD",             no_argument,       0, 0 },   // 29
	{ "help",           no_argument,       0, 'h' },
	{ "max-intron-len", required_argument, 0, 'G' },
	{ "version",        no_argument,       0, 'V' },
	{ "min-count",      required_argument, 0, 'n' },
	{ "min-chain-score",required_argument, 0, 'm' },
	{ "mask-level",     required_argument, 0, 'M' },
	{ "min-dp-score",   required_argument, 0, 's' },
	{ "sam",            no_argument,       0, 'a' },
	{ 0, 0, 0, 0}
};

static inline int64_t mm_parse_num(const char *str)
{
	double x;
	char *p;
	x = strtod(optarg, &p);
	if (*p == 'G' || *p == 'g') x *= 1e9;
	else if (*p == 'M' || *p == 'm') x *= 1e6;
	else if (*p == 'K' || *p == 'k') x *= 1e3;
	return (int64_t)(x + .499);
}

static inline void yes_or_no(mm_mapopt_t *opt, int flag, int long_idx, const char *arg, int yes_to_set)
{
	if (yes_to_set) {
		if (strcmp(arg, "yes") == 0 || strcmp(arg, "y") == 0) opt->flag |= flag;
		else if (strcmp(arg, "no") == 0 || strcmp(arg, "n") == 0) opt->flag &= ~flag;
		else fprintf(stderr, "[WARNING]\033[1;31m option '--%s' only accepts 'yes' or 'no'.\033[0m\n", long_options[long_idx].name);
	} else {
		if (strcmp(arg, "yes") == 0 || strcmp(arg, "y") == 0) opt->flag &= ~flag;
		else if (strcmp(arg, "no") == 0 || strcmp(arg, "n") == 0) opt->flag |= flag;
		else fprintf(stderr, "[WARNING]\033[1;31m option '--%s' only accepts 'yes' or 'no'.\033[0m\n", long_options[long_idx].name);
	}
}

void init_task_array();
void init_result_array();
void stop_sw_thread();

double chaindp_time[100];
double step0_first = 0;
double step0[100];
double step1[100];
double step2[100];

double process_task_time[100];
double get_mem_time[100];
double package_time[100];
double send_time[100];
double chaindp_sw_time[100];
double sw_time[100];
double sw_task_time[100];
double sw_soft_time[100];
double create_time;
double mem_init_time;

double chaindp_1[100];
double chaindp_2[100];
double chaindp_3[100];
double chaindp_4[100];

double multi_thread_time[100];

double fun_1_time[100];
double fun_2_time[100];
double fun_3_time[100];

int main(int argc, char *argv[])
{
	const char *opt_str = "2aSDw:k:K:t:r:f:Vv:g:G:I:d:XT:s:x:Hcp:M:n:z:A:B:O:E:m:N:Qu:R:hF:LC:y";
	mm_mapopt_t opt;
	mm_idxopt_t ipt;
	int i, c, n_threads = 3, long_idx;
	char *fnw = 0, *rg = 0, *s;
	FILE *fp_help = stderr;
	mm_idx_reader_t *idx_rdr;
	mm_idx_t *mi;

    create_time = 0;
    memset(chaindp_time, 0, sizeof(chaindp_time));
    memset(step0, 0, sizeof(step0));
    memset(step1, 0, sizeof(step1));
    memset(step2, 0, sizeof(step2));
    memset(process_task_time, 0, sizeof(process_task_time));
    memset(get_mem_time, 0, sizeof(get_mem_time));
    memset(package_time, 0, sizeof(package_time));
    memset(send_time, 0, sizeof(send_time));
    memset(chaindp_sw_time, 0, sizeof(chaindp_sw_time));
    memset(sw_time, 0, sizeof(sw_time));
    memset(sw_task_time, 0, sizeof(sw_task_time));
    memset(sw_soft_time, 0, sizeof(sw_soft_time));
    memset(chaindp_1, 0, sizeof(chaindp_1));
    memset(chaindp_2, 0, sizeof(chaindp_2));
    memset(chaindp_3, 0, sizeof(chaindp_3));
    memset(chaindp_4, 0, sizeof(chaindp_4));
    memset(multi_thread_time, 0, sizeof(multi_thread_time));
    memset(fun_1_time, 0, sizeof(fun_1_time));
    memset(fun_2_time, 0, sizeof(fun_2_time));
    memset(fun_3_time, 0, sizeof(fun_3_time));
    
	mm_verbose = 3;
	liftrlimit();
	mm_realtime0 = realtime();
	mm_set_opt(0, &ipt, &opt);

	while ((c = getopt_long(argc, argv, opt_str, long_options, &long_idx)) >= 0) // apply option -x/preset first
		if (c == 'x') {
			if (mm_set_opt(optarg, &ipt, &opt) < 0) {
				fprintf(stderr, "[ERROR] unknown preset '%s'\n", optarg);
				return 1;
			}
			break;
		}
	optind = 0; // for musl getopt, optind=0 has the same effect as optreset=1; older libc doesn't have optreset

	while ((c = getopt_long(argc, argv, opt_str, long_options, &long_idx)) >= 0) {
		if (c == 'w') ipt.w = atoi(optarg);
		else if (c == 'k') ipt.k = atoi(optarg);
		else if (c == 'H') ipt.flag |= MM_I_HPC;
		else if (c == 'd') fnw = optarg; // the above are indexing related options, except -I
		else if (c == 'r') opt.bw = (int)mm_parse_num(optarg);
		else if (c == 't') n_threads = atoi(optarg);
		else if (c == 'v') mm_verbose = atoi(optarg);
		else if (c == 'g') opt.max_gap = (int)mm_parse_num(optarg);
		else if (c == 'G') mm_mapopt_max_intron_len(&opt, (int)mm_parse_num(optarg));
		else if (c == 'F') opt.max_frag_len = (int)mm_parse_num(optarg);
		else if (c == 'N') opt.best_n = atoi(optarg);
		else if (c == 'p') opt.pri_ratio = atof(optarg);
		else if (c == 'M') opt.mask_level = atof(optarg);
		else if (c == 'c') opt.flag |= MM_F_OUT_CG | MM_F_CIGAR;
		else if (c == 'D') opt.flag |= MM_F_NO_DIAG;
		else if (c == 'P') opt.flag |= MM_F_ALL_CHAINS;
		else if (c == 'X') opt.flag |= MM_F_ALL_CHAINS | MM_F_NO_DIAG | MM_F_NO_DUAL | MM_F_NO_LJOIN; // -D -P --no-long-join --dual=no
		else if (c == 'a') opt.flag |= MM_F_OUT_SAM | MM_F_CIGAR;
		else if (c == 'Q') opt.flag |= MM_F_NO_QUAL;
		else if (c == 'Y') opt.flag |= MM_F_SOFTCLIP;
		else if (c == 'L') opt.flag |= MM_F_LONG_CIGAR;
		else if (c == 'y') opt.flag |= MM_F_COPY_COMMENT;
		else if (c == 'T') opt.sdust_thres = atoi(optarg);
		else if (c == 'n') opt.min_cnt = atoi(optarg);
		else if (c == 'm') opt.min_chain_score = atoi(optarg);
		else if (c == 'A') opt.a = atoi(optarg);
		else if (c == 'B') opt.b = atoi(optarg);
		else if (c == 's') opt.min_dp_max = atoi(optarg);
		else if (c == 'C') opt.noncan = atoi(optarg);
		else if (c == 'I') ipt.batch_size = mm_parse_num(optarg);
		else if (c == 'K') opt.mini_batch_size = (int)mm_parse_num(optarg);
		else if (c == 'R') rg = optarg;
		else if (c == 'h') fp_help = stdout;
		else if (c == '2') opt.flag |= MM_F_2_IO_THREADS;
		else if (c == 0 && long_idx == 0) ipt.bucket_bits = atoi(optarg); // --bucket-bits
		else if (c == 0 && long_idx == 2) opt.seed = atoi(optarg); // --seed
		else if (c == 0 && long_idx == 3) mm_dbg_flag |= MM_DBG_NO_KALLOC; // --no-kalloc
		else if (c == 0 && long_idx == 4) mm_dbg_flag |= MM_DBG_PRINT_QNAME; // --print-qname
		else if (c == 0 && long_idx == 6) mm_dbg_flag |= MM_DBG_PRINT_QNAME | MM_DBG_PRINT_SEED, n_threads = 1; // --print-seed
		else if (c == 0 && long_idx == 7) opt.max_chain_skip = atoi(optarg); // --max-chain-skip
		else if (c == 0 && long_idx == 8) opt.min_ksw_len = atoi(optarg); // --min-dp-len
		else if (c == 0 && long_idx == 9) mm_dbg_flag |= MM_DBG_PRINT_QNAME | MM_DBG_PRINT_ALN_SEQ, n_threads = 1; // --print-aln-seq
		else if (c == 0 && long_idx ==10) opt.flag |= MM_F_SPLICE; // --splice
		else if (c == 0 && long_idx ==12) opt.flag |= MM_F_NO_LJOIN; // --no-long-join
		else if (c == 0 && long_idx ==13) opt.flag |= MM_F_SR; // --sr
		else if (c == 0 && long_idx ==17) opt.end_bonus = atoi(optarg); // --end-bonus
		else if (c == 0 && long_idx ==18) opt.flag |= MM_F_INDEPEND_SEG; // --no-pairing
		else if (c == 0 && long_idx ==20) ipt.flag |= MM_I_NO_SEQ; // --idx-no-seq
		else if (c == 0 && long_idx ==21) opt.anchor_ext_shift = atoi(optarg); // --end-seed-pen
		else if (c == 0 && long_idx ==22) opt.flag |= MM_F_FOR_ONLY; // --for-only
		else if (c == 0 && long_idx ==23) opt.flag |= MM_F_REV_ONLY; // --rev-only
		else if (c == 0 && long_idx ==27) opt.max_clip_ratio = atof(optarg); // --max-clip-ratio
		else if (c == 0 && long_idx ==28) opt.min_mid_occ = atoi(optarg); // --min-occ-floor
		else if (c == 0 && long_idx ==29) opt.flag |= MM_F_OUT_MD; // --MD
		else if (c == 0 && long_idx == 14) { // --frag
			yes_or_no(&opt, MM_F_FRAG_MODE, long_idx, optarg, 1);
		} else if (c == 0 && long_idx == 15) { // --secondary
			yes_or_no(&opt, MM_F_NO_PRINT_2ND, long_idx, optarg, 0);
		} else if (c == 0 && long_idx == 16) { // --cs
			opt.flag |= MM_F_OUT_CS | MM_F_CIGAR;
			if (optarg == 0 || strcmp(optarg, "short") == 0) {
				opt.flag &= ~MM_F_OUT_CS_LONG;
			} else if (strcmp(optarg, "long") == 0) {
				opt.flag |= MM_F_OUT_CS_LONG;
			} else if (strcmp(optarg, "none") == 0) {
				opt.flag &= ~MM_F_OUT_CS;
			} else if (mm_verbose >= 2) {
				fprintf(stderr, "[WARNING]\033[1;31m --cs only takes 'short' or 'long'. Invalid values are assumed to be 'short'.\033[0m\n");
			}
		} else if (c == 0 && long_idx == 19) { // --splice-flank
			yes_or_no(&opt, MM_F_SPLICE_FLANK, long_idx, optarg, 1);
		} else if (c == 0 && long_idx == 24) { // --heap-sort
			yes_or_no(&opt, MM_F_HEAP_SORT, long_idx, optarg, 1);
		} else if (c == 0 && long_idx == 26) { // --dual
			yes_or_no(&opt, MM_F_NO_DUAL, long_idx, optarg, 0);
		} else if (c == 'S') {
			opt.flag |= MM_F_OUT_CS | MM_F_CIGAR | MM_F_OUT_CS_LONG;
			if (mm_verbose >= 2)
				fprintf(stderr, "[WARNING]\033[1;31m option -S is deprecated and may be removed in future. Please use --cs=long instead.\033[0m\n");
		} else if (c == 'V') {
			puts(MM_VERSION);
			return 0;
		} else if (c == 'f') {
			double x;
			char *p;
			x = strtod(optarg, &p);
			if (x < 1.0) opt.mid_occ_frac = x, opt.mid_occ = 0;
			else opt.mid_occ = (int)(x + .499);
			if (*p == ',') opt.max_occ = (int)(strtod(p+1, &p) + .499);
		} else if (c == 'u') {
			if (*optarg == 'b') opt.flag |= MM_F_SPLICE_FOR|MM_F_SPLICE_REV; // both strands
			else if (*optarg == 'f') opt.flag |= MM_F_SPLICE_FOR, opt.flag &= ~MM_F_SPLICE_REV; // match GT-AG
			else if (*optarg == 'r') opt.flag |= MM_F_SPLICE_REV, opt.flag &= ~MM_F_SPLICE_FOR; // match CT-AC (reverse complement of GT-AG)
			else if (*optarg == 'n') opt.flag &= ~(MM_F_SPLICE_FOR|MM_F_SPLICE_REV); // don't try to match the GT-AG signal
			else {
				fprintf(stderr, "[ERROR]\033[1;31m unrecognized cDNA direction\033[0m\n");
				return 1;
			}
		} else if (c == 'z') {
			opt.zdrop = opt.zdrop_inv = strtol(optarg, &s, 10);
			if (*s == ',') opt.zdrop_inv = strtol(s + 1, &s, 10);
		} else if (c == 'O') {
			opt.q = opt.q2 = strtol(optarg, &s, 10);
			if (*s == ',') opt.q2 = strtol(s + 1, &s, 10);
		} else if (c == 'E') {
			opt.e = opt.e2 = strtol(optarg, &s, 10);
			if (*s == ',') opt.e2 = strtol(s + 1, &s, 10);
		}
	}
	if ((opt.flag & MM_F_SPLICE) && (opt.flag & MM_F_FRAG_MODE)) {
		fprintf(stderr, "[ERROR]\033[1;31m --splice and --frag should not be specified at the same time.\033[0m\n");
		return 1;
	}
	if (!fnw && !(opt.flag&MM_F_CIGAR))
		ipt.flag |= MM_I_NO_SEQ;
	if (mm_check_opt(&ipt, &opt) < 0)
		return 1;

	if (argc == optind || fp_help == stdout) {
		fprintf(fp_help, "Usage: minimap2 [options] <target.fa>|<target.idx> [query.fa] [...]\n");
		fprintf(fp_help, "Options:\n");
		fprintf(fp_help, "  Indexing:\n");
		fprintf(fp_help, "    -H           use homopolymer-compressed k-mer\n");
		fprintf(fp_help, "    -k INT       k-mer size (no larger than 28) [%d]\n", ipt.k);
		fprintf(fp_help, "    -w INT       minizer window size [%d]\n", ipt.w);
		fprintf(fp_help, "    -I NUM       split index for every ~NUM input bases [4G]\n");
		fprintf(fp_help, "    -d FILE      dump index to FILE []\n");
		fprintf(fp_help, "  Mapping:\n");
		fprintf(fp_help, "    -f FLOAT     filter out top FLOAT fraction of repetitive minimizers [%g]\n", opt.mid_occ_frac);
		fprintf(fp_help, "    -g NUM       stop chain enlongation if there are no minimizers in INT-bp [%d]\n", opt.max_gap);
		fprintf(fp_help, "    -G NUM       max intron length (effective with -xsplice; changing -r) [200k]\n");
		fprintf(fp_help, "    -F NUM       max fragment length (effective with -xsr or in the fragment mode) [800]\n");
		fprintf(fp_help, "    -r NUM       bandwidth used in chaining and DP-based alignment [%d]\n", opt.bw);
		fprintf(fp_help, "    -n INT       minimal number of minimizers on a chain [%d]\n", opt.min_cnt);
		fprintf(fp_help, "    -m INT       minimal chaining score (matching bases minus log gap penalty) [%d]\n", opt.min_chain_score);
//		fprintf(fp_help, "    -T INT       SDUST threshold; 0 to disable SDUST [%d]\n", opt.sdust_thres); // TODO: this option is never used; might be buggy
		fprintf(fp_help, "    -X           skip self and dual mappings (for the all-vs-all mode)\n");
		fprintf(fp_help, "    -p FLOAT     min secondary-to-primary score ratio [%g]\n", opt.pri_ratio);
		fprintf(fp_help, "    -N INT       retain at most INT secondary alignments [%d]\n", opt.best_n);
		fprintf(fp_help, "  Alignment:\n");
		fprintf(fp_help, "    -A INT       matching score [%d]\n", opt.a);
		fprintf(fp_help, "    -B INT       mismatch penalty [%d]\n", opt.b);
		fprintf(fp_help, "    -O INT[,INT] gap open penalty [%d,%d]\n", opt.q, opt.q2);
		fprintf(fp_help, "    -E INT[,INT] gap extension penalty; a k-long gap costs min{O1+k*E1,O2+k*E2} [%d,%d]\n", opt.e, opt.e2);
		fprintf(fp_help, "    -z INT[,INT] Z-drop score and inversion Z-drop score [%d,%d]\n", opt.zdrop, opt.zdrop_inv);
		fprintf(fp_help, "    -s INT       minimal peak DP alignment score [%d]\n", opt.min_dp_max);
		fprintf(fp_help, "    -u CHAR      how to find GT-AG. f:transcript strand, b:both strands, n:don't match GT-AG [n]\n");
		fprintf(fp_help, "  Input/Output:\n");
		fprintf(fp_help, "    -a           output in the SAM format (PAF by default)\n");
		fprintf(fp_help, "    -Q           don't output base quality in SAM\n");
		fprintf(fp_help, "    -L           write CIGAR with >65535 ops at the CG tag\n");
		fprintf(fp_help, "    -R STR       SAM read group line in a format like '@RG\\tID:foo\\tSM:bar' []\n");
		fprintf(fp_help, "    -c           output CIGAR in PAF\n");
		fprintf(fp_help, "    --cs[=STR]   output the cs tag; STR is 'short' (if absent) or 'long' [none]\n");
		fprintf(fp_help, "    --MD         output the MD tag\n");
		fprintf(fp_help, "    -Y           use soft clipping for supplementary alignments\n");
		fprintf(fp_help, "    -t INT       number of threads [%d]\n", n_threads);
		fprintf(fp_help, "    -K NUM       minibatch size for mapping [500M]\n");
//		fprintf(fp_help, "    -v INT       verbose level [%d]\n", mm_verbose);
		fprintf(fp_help, "    --version    show version number\n");
		fprintf(fp_help, "  Preset:\n");
		fprintf(fp_help, "    -x STR       preset (always applied before other options) []\n");
		fprintf(fp_help, "                 map-pb: -Hk19 (PacBio vs reference mapping)\n");
		fprintf(fp_help, "                 map-ont: -k15 (Oxford Nanopore vs reference mapping)\n");
		fprintf(fp_help, "                 asm5: -k19 -w19 -A1 -B19 -O39,81 -E3,1 -s200 -z200 (asm to ref mapping; break at 5%% div.)\n");
		fprintf(fp_help, "                 asm10: -k19 -w19 -A1 -B9 -O16,41 -E2,1 -s200 -z200 (asm to ref mapping; break at 10%% div.)\n");
		fprintf(fp_help, "                 ava-pb: -Hk19 -Xw5 -m100 -g10000 --max-chain-skip 25 (PacBio read overlap)\n");
		fprintf(fp_help, "                 ava-ont: -k15 -Xw5 -m100 -g10000 -r2000 --max-chain-skip 25 (ONT read overlap)\n");
		fprintf(fp_help, "                 splice: long-read spliced alignment (see minimap2.1 for details)\n");
		fprintf(fp_help, "                 sr: short single-end reads without splicing (see minimap2.1 for details)\n");
		fprintf(fp_help, "\nSee `man ./minimap2.1' for detailed description of command-line options.\n");
		return fp_help == stdout? 0 : 1;
	}

	if ((opt.flag & MM_F_SR) && argc - optind > 3) {
		fprintf(stderr, "[ERROR] incorrect input: in the sr mode, please specify no more than two query files.\n");
		return 1;
	}
	idx_rdr = mm_idx_reader_open(argv[optind], &ipt, fnw);
	if (idx_rdr == 0) {
		fprintf(stderr, "[ERROR] failed to open file '%s'\n", argv[optind]);
		return 1;
	}
	if (!idx_rdr->is_idx && fnw == 0 && argc - optind < 2) {
		fprintf(stderr, "[ERROR] missing input: please specify a query file to map or option -d to keep the index\n");
		mm_idx_reader_close(idx_rdr);
		return 1;
	}
	if (opt.best_n == 0 && (opt.flag&MM_F_CIGAR) && mm_verbose >= 2)
		fprintf(stderr, "[WARNING]\033[1;31m `-N 0' reduces alignment accuracy. Please use --secondary=no to suppress secondary alignments.\033[0m\n");



    fprintf(stderr, "50 align with 16 is %d\n", ADDR_ALIGN(50, 16));
    fprintf(stderr, "64 align with 16 is %d\n", ADDR_ALIGN(64, 16));
    
    fprintf(stderr, "sizeof(fpga_task_t)=%ld\n", sizeof(fpga_task_t));
    fprintf(stderr, "sizeof(fpga_task_id_t)=%ld\n", sizeof(fpga_task_id_t));
    fprintf(stderr, "sizeof(fpga_sw_task)=%ld\n", sizeof(fpga_sw_task));
    fprintf(stderr, "sizeof(ksw_extz_t)=%ld\n", sizeof(ksw_extz_t));
    
    assert(sizeof(fpga_sw_task) == 16);
    
    init_task_array();
    init_result_array();
    init_fpga_task_array();
    init_fpga_result_array();

#if !DUMP_FILE
    int ret = fpga_init(BLOCK);
    if(ret) {
        printf("fpga_init failed\n");
        return -1;
    }
#endif
    /*pthread_t send_tid[10];
    pthread_t recv_tid[10];
    int thread_i = 0;
    int tid[10];
    for(thread_i = 0; thread_i<1; thread_i++) {
        tid[thread_i] = thread_i;
        pthread_create(&send_tid[thread_i], NULL, send_task_thread, (void*)tid[thread_i]);
        pthread_create(&recv_tid[thread_i], NULL, recv_task_thread, tid[thread_i]);
    }*/
    
    /*sleep(3);
    uint64_t version;
    int count = 10;
    while(count--) {
        version = fpga_read_reg(0);
        fprintf(stderr, "version = 0x%016lx\n", version);
    }
    exit(1);*/
    
    fpga_set_block();
    pthread_t send_tid, recv_tid;
    pthread_create(&send_tid, NULL, send_task_thread, NULL);
    pthread_create(&recv_tid, NULL, recv_task_thread, NULL);

    while ((mi = mm_idx_reader_read(idx_rdr, n_threads)) != 0) {
		if ((opt.flag & MM_F_CIGAR) && (mi->flag & MM_I_NO_SEQ)) {
			fprintf(stderr, "[ERROR] the prebuilt index doesn't contain sequences.\n");
			mm_idx_destroy(mi);
			mm_idx_reader_close(idx_rdr);
			return 1;
		}
		if ((opt.flag & MM_F_OUT_SAM) && idx_rdr->n_parts == 1) {
			if (mm_idx_reader_eof(idx_rdr)) {
				mm_write_sam_hdr(mi, rg, MM_VERSION, argc, argv);
			} else {
				mm_write_sam_hdr(0, rg, MM_VERSION, argc, argv);
				if (mm_verbose >= 2)
					fprintf(stderr, "[WARNING]\033[1;31m For a multi-part index, no @SQ lines will be outputted.\033[0m\n");
			}
		}
		if (mm_verbose >= 3)
			fprintf(stderr, "[M::%s::%.3f*%.2f] loaded/built the index for %d target sequence(s)\n",
					__func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0), mi->n_seq);
		if (argc != optind + 1) mm_mapopt_update(&opt, mi);
		if (mm_verbose >= 3) mm_idx_stat(mi);
		if (!(opt.flag & MM_F_FRAG_MODE)) {
			for (i = optind + 1; i < argc; ++i)
				mm_map_file(mi, argv[i], &opt, n_threads);
		} else {
			mm_map_file_frag(mi, argc - (optind + 1), (const char**)&argv[optind + 1], &opt, n_threads);
		}
		mm_idx_destroy(mi);
	}
	mm_idx_reader_close(idx_rdr);

	if (fflush(stdout) == EOF) {
		fprintf(stderr, "[ERROR] failed to write the results\n");
		exit(EXIT_FAILURE);
	}
    
    stop_fpga_recv_thread();
    stop_fpga_send_thread();
    fpga_exit_block();
    pthread_join(recv_tid, NULL);
    pthread_join(send_tid, NULL);
#if !DUMP_FILE
    fpga_exit_block();
#endif
    /*for(thread_i = 0; thread_i<1; thread_i++) {
        pthread_join(send_tid[thread_i], NULL);
        pthread_join(recv_tid[thread_i], NULL);
    }*/
#if !DUMP_FILE
    //fpga_set_block();
    fpga_finalize();
#endif


	if (mm_verbose >= 3) {
		fprintf(stderr, "[M::%s] Version: %s\n", __func__, MM_VERSION);
		fprintf(stderr, "[M::%s] CMD:", __func__);
		for (i = 0; i < argc; ++i)
			fprintf(stderr, " %s", argv[i]);
		fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__, realtime() - mm_realtime0, cputime());
	}
    fprintf(stderr, "first step 0 time:%.3f msec\n", step0_first);

    double chain_sw_total = 0;
    for(i = 0; i < sizeof(chaindp_sw_time)/sizeof(chaindp_sw_time[0]); i++) {
        if(chaindp_sw_time[i] == 0)
            break;
        chain_sw_total += chaindp_sw_time[i];
    }
    fprintf(stderr, "chain sw total time:      %.3f msec i=%d\n", chain_sw_total/i, i);
    
    double chain_dp_total = 0;
    for(i = 0; i < sizeof(chaindp_time)/sizeof(chaindp_time[0]); i++) {
        if(chaindp_time[i] == 0)
            break;
        chain_dp_total += chaindp_time[i];
    }
    fprintf(stderr, "chain dp time:      %.3f msec i=%d chain_dp_total=%.3f\n", chain_dp_total/i, i, chain_dp_total);
    
    //chaindp 1
    double chaindp_1_total = 0;
    for(i = 0; i < sizeof(chaindp_1)/sizeof(chaindp_1[0]); i++) {
        if(chaindp_1[i] == 0)
            break;
        chaindp_1_total += chaindp_1[i];
    }
    fprintf(stderr, "chain dp 1 time:      %.3f msec i=%d\n", chaindp_1_total/i, i);
    
    //chaindp 2
    double chaindp_2_total = 0;
    for(i = 0; i < sizeof(chaindp_2)/sizeof(chaindp_2[0]); i++) {
        if(chaindp_2[i] == 0)
            break;
        chaindp_2_total += chaindp_2[i];
    }
    fprintf(stderr, "chain dp 2 time:      %.3f msec i=%d\n", chaindp_2_total/i, i);
    
    //chaindp 3
    double chaindp_3_total = 0;
    for(i = 0; i < sizeof(chaindp_3)/sizeof(chaindp_3[0]); i++) {
        if(chaindp_3[i] == 0)
            break;
        chaindp_3_total += chaindp_3[i];
    }
    fprintf(stderr, "chain dp 3 time:      %.3f msec i=%d\n", chaindp_3_total/i, i);
    
    //chaindp 4
    double chaindp_4_total = 0;
    for(i = 0; i < sizeof(chaindp_4)/sizeof(chaindp_4[0]); i++) {
        if(chaindp_4[i] == 0)
            break;
        chaindp_4_total += chaindp_4[i];
    }
    fprintf(stderr, "chain dp 4 time:      %.3f msec i=%d\n", chaindp_4_total/i, i);
    
    double sw_total = 0;
    for(i = 0; i < sizeof(sw_time)/sizeof(sw_time[0]); i++) {
        if(sw_time[i] == 0)
            break;
        sw_total += sw_time[i];
    }
    fprintf(stderr, "sw time:      %.3f msec i=%d\n", sw_total/i, i);
    
    double process_total = 0;
    for(i = 0; i < sizeof(process_task_time)/sizeof(process_task_time[0]); i++) {
        if(process_task_time[i] == 0)
            break;
        process_total += process_task_time[i];
    }
    fprintf(stderr, "process task time:      %.3f msec i=%d\n", process_total/i, i);
    
    double get_mem_total = 0;
    for(i = 0; i < sizeof(get_mem_time)/sizeof(get_mem_time[0]); i++) {
        if(get_mem_time[i] == 0)
            break;
        get_mem_total += get_mem_time[i];
    }
    fprintf(stderr, "get mem time:      %.3f msec i=%d\n", get_mem_total/i, i);
    
    double package_total = 0;
    for(i = 0; i < sizeof(package_time)/sizeof(package_time[0]); i++) {
        if(package_time[i] == 0)
            break;
        package_total += package_time[i];
    }
    fprintf(stderr, "package time:      %.3f msec i=%d\n", package_total/i, i);
    
    double send_total = 0;
    for(i = 0; i < sizeof(send_time)/sizeof(send_time[0]); i++) {
        if(send_time[i] == 0)
            break;
        send_total += send_time[i];
    }
    fprintf(stderr, "submit time:      %.3f msec i=%d\n", send_total/i, i);
    
    double sw_task_total = 0;
    for(i = 0; i < sizeof(sw_task_time)/sizeof(sw_task_time[0]); i++) {
        if(sw_task_time[i] == 0)
            break;
        sw_task_total += sw_task_time[i];
    }
    fprintf(stderr, "sw task time:      %.3f msec i=%d\n", sw_task_total/i, i);
    
    double sw_soft_total = 0;
    for(i = 0; i < sizeof(sw_soft_time)/sizeof(sw_soft_time[0]); i++) {
        if(sw_soft_time[i] == 0)
            break;
        sw_soft_total += sw_soft_time[i];
    }
    fprintf(stderr, "sw soft time:      %.3f msec i=%d\n", sw_soft_total/i, i);
    
    fprintf(stderr, "mem init time:           %.3f msec\n", mem_init_time);
    fprintf(stderr, "create thread time:      %.3f msec\n", create_time);
    
    for(i = 0; i < sizeof(multi_thread_time)/sizeof(multi_thread_time[0]); i++) {
            if(multi_thread_time[i] == 0)
                    break;
            fprintf(stderr, "multi_thread_time time[%d]:      %.3f msec\n", i, multi_thread_time[i]);
    }
    fprintf(stderr, "\n");
    
    
    fprintf(stderr, "\n\n\n");
    for(i = 0; i < sizeof(step0)/sizeof(step0[0]); i++) {
            if(step0[i] == 0)
                    break;
            fprintf(stderr, "step0 time[%d]:      %.3f msec\n", i, step0[i]);
    }
    fprintf(stderr, "\n");
    double step1_total = 0;
    for(i = 0; i < sizeof(step1)/sizeof(step1[0]); i++) {
            if(step1[i] == 0)
                    break;
            step1_total += step1[i];
            fprintf(stderr, "step1 time[%d]:      %.3f msec\n", i, step1[i]);
    }
    fprintf(stderr, "step1 total time:      %.3f msec\n", step1_total);
    fprintf(stderr, "\n");
    for(i = 0; i < sizeof(step2)/sizeof(step2[0]); i++) {
            if(step2[i] == 0)
                    break;
            fprintf(stderr, "step2 time[%d]:      %.3f msec\n", i, step2[i]);
    }
    
    fprintf(stderr, "\n");
    double fun1_total = 0;
    for(i = 0; i < sizeof(fun_1_time)/sizeof(fun_1_time[0]); i++) {
            if(fun_1_time[i] == 0)
                    break;
            fun1_total += fun_1_time[i];
    }
    fprintf(stderr, "fun_1_time time:      %.3f msec\n", fun1_total/i);
    fprintf(stderr, "\n");
    
    double fun2_total = 0;
    for(i = 0; i < sizeof(fun_2_time)/sizeof(fun_2_time[0]); i++) {
            if(fun_2_time[i] == 0)
                    break;
            fun2_total += fun_2_time[i];
    }
    fprintf(stderr, "fun_2_time time:      %.3f msec\n", fun2_total/i);
    fprintf(stderr, "\n");
    
    double fun3_total = 0;
    for(i = 0; i < sizeof(fun_3_time)/sizeof(fun_3_time[0]); i++) {
            if(fun_3_time[i] == 0)
                    break;
            fun3_total += fun_3_time[i];
    }
    fprintf(stderr, "fun_3_time time:      %.3f msec\n", fun3_total/i);
    fprintf(stderr, "\n");
    
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    
    fprintf(stderr, "ru_utime=%.3f\n", (usage.ru_utime.tv_sec*1000 + usage.ru_utime.tv_usec*1e-3));
    fprintf(stderr, "ru_stime=%.3f\n", (usage.ru_stime.tv_sec*1000 + usage.ru_stime.tv_usec*1e-3));
    fprintf(stderr, "ru_maxrss=%ld\n", usage.ru_maxrss);
    fprintf(stderr, "ru_ixrss=%ld\n", usage.ru_ixrss);
    fprintf(stderr, "ru_idrss=%ld\n", usage.ru_idrss);
    fprintf(stderr, "ru_isrss=%ld\n", usage.ru_isrss);
    fprintf(stderr, "ru_minflt=%ld\n", usage.ru_minflt);
    fprintf(stderr, "ru_majflt=%ld\n", usage.ru_majflt);
    fprintf(stderr, "ru_nswap=%ld\n", usage.ru_nswap);
    fprintf(stderr, "ru_inblock=%ld\n", usage.ru_inblock);
    fprintf(stderr, "ru_oublock=%ld\n", usage.ru_oublock);
    fprintf(stderr, "ru_msgsnd=%ld\n", usage.ru_msgsnd);
    fprintf(stderr, "ru_msgrcv=%ld\n", usage.ru_msgrcv);
    fprintf(stderr, "ru_nsignals=%ld\n", usage.ru_nsignals);
    fprintf(stderr, "ru_nvcsw=%ld\n", usage.ru_nvcsw);
    fprintf(stderr, "ru_nivcsw=%ld\n", usage.ru_nivcsw);
    
    
	return 0;
}
