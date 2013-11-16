/**
 * @file ccd.c
 *
 * @author andy.cj.wu@gmail.com
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2010 Washington State University. All rights reserved.
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 */
#include "config.h"

#include <cassert>
#include <climits>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>

#include <time.h>
#include <sys/time.h>

#define USE_SET 1
#if USE_SET
#include <set>
#include <vector>
#include <utility>
using std::set;
using std::vector;
using std::make_pair;
using std::pair;
#endif
#include <iostream>
using std::cout;
using std::endl;

#include <unistd.h>

#include "bucket.hpp"
#include "combinations.h"
#include "constants.h"
#include "csequence.h"
#include "stree.hpp"

#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#endif


using namespace pgraph;

static void parse_command_line(int argc, char **argv,
        char *sequence_file, char *config_file, size_t *n_sequences);
static inline size_t zpower(size_t base, size_t n);

static inline double wtime() {
    struct timespec ts;
#ifdef __MACH__ // OS X does not have clock_gettime, use clock_get_time
    clock_serv_t cclock;
    mach_timespec_t mts;
    host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
    clock_get_time(cclock, &mts);
    mach_port_deallocate(mach_task_self(), cclock);
    ts.tv_sec = mts.tv_sec;
    ts.tv_nsec = mts.tv_nsec;
#else
    int retval = clock_gettime(CLOCK_REALTIME, &ts);
    assert(0 == retval);
#endif
    return double(ts.tv_sec) + double(ts.tv_nsec) / 1000000000.0;
}


int main(int argc, char *argv[])
{
    char sequence_file[FILENAME_MAX];/* path to sequence (fasta) file */
    char config_file[FILENAME_MAX]; /* path to config file */
    size_t n_sequences = 0;         /* number of sequences (cmd line param) */
    sequences_t *sequences = NULL;  /* all sequences parsed from fasta file */
    param_t param;                  /* config file parameters */
    time_t t1 = 0;                  /* start timer */
    time_t t2 = 0;                  /* stop timer */
    suffix_buckets_t *suffix_buckets = NULL;
//#ifdef DEBUG
    size_t i = 0;                   /* for loop index */
//#endif
    ssize_t d = 0;                  /* signed for loop index */
    size_t n_triangular = 0;        /* number of possible pairs */
#if USE_SET
    set<pair<unsigned int,unsigned int> > pairs;
#else
    int *dup = NULL;                /* track duplicate pairs */
#endif

    setbuf(stdout, (char *)0);

    /* reading in pamameters */
    parse_command_line(argc, argv, sequence_file, config_file, &n_sequences);

    /* read in configurations */
    get_params(config_file, &param);
    
    printf("----------------------------------------------\n");
    printf("%-20s: %.78s\n", "fasta seq", sequence_file);
    printf("%-20s: %.78s\n", "config file", config_file);
    printf("%-20s: %d\n", "slide size", param.window_size);
    printf("%-20s: %d\n", "exactMatch len", param.exact_match_len);
    printf("%-20s: %d\n", "AlignOverLongerSeq", param.AOL);
    printf("%-20s: %d\n", "MatchSimilarity", param.SIM);
    printf("%-20s: %d\n", "OptimalScoreOverSelfScore", param.OS);
    printf("%-20s: %zu\n", "nseqs", n_sequences);
    printf("----------------------------------------------\n");

    (void) time(&t1);
    sequences = pg_load_fasta(sequence_file, DOLLAR);
    (void) time(&t2);
    printf("<%zu> amino acids (<%zu> chars) are loaded in <%lld> secs\n",
            sequences->size, sequences->n_chars, (long long)(t2-t1));
    printf("max_seq_size %zu\n", sequences->max_seq_size);
    unsigned long ntasks = binomial_coefficient(n_sequences, 2);
    printf("ntasks %lu\n", ntasks);

    (void) time(&t1);
    suffix_buckets = create_suffix_buckets_old(sequences, param);
    (void) time(&t2);
    printf("Bucketing finished in <%lld> secs\n", (long long)(t2-t1));
    
    #ifdef DEBUG
    {
        for (i=0; i<100; ++i) {
            printf("n_buckets[%d] - %d\n", i, buckets[i].size);
        }
        exit(0);
    }
    #endif

    n_triangular = sequences->size * (sequences->size + 1U) / 2U;
#if USE_SET
#else
    dup = new int[n_triangular];
#pragma omp parallel for
    for (d = 0; d < (ssize_t)n_triangular; ++d) {
        dup[d] = MAYBE;
    }
#endif
    
    /* suffix tree construction & processing */
    printf("n_buckets=%zu\n", suffix_buckets->buckets_size);
    printf("n_suffixes=%zu\n", suffix_buckets->suffixes_size);
    (void) time(&t1);
    size_t count = 0;
    cout << "tree" << ",";
    stats_t::print_header(cout, "fanout");
    cout << ",";
    stats_t::print_header(cout, "depth");
    cout << ",";
    stats_t::print_header(cout, "sequence_length");
    cout << ",";
    stats_t::print_header(cout, "suffix_length");
    cout << ",";
    cout << "time_build";
    cout << ",";
    cout << "time_process";
    cout << endl;
    double total_time_build = 0;
    double total_time_process = 0;
    for (d = 0; d < (long)suffix_buckets->buckets_size; ++d) {
        if (NULL != suffix_buckets->buckets[d].suffixes) {
            double btimer = wtime();
            stree_t *tree = build_tree(
                    sequences, &(suffix_buckets->buckets[d]), param);
            btimer = wtime() - btimer;
            double ptimer = wtime();
#if USE_SET
            generate_pairs(tree, sequences, pairs, param);
#else
            generate_pairs(tree, sequences, dup, param);
#endif
            ptimer = wtime() - ptimer;
            cout << d << ",";
            stats_t::print(cout, tree->fanout);
            cout << ",";
            stats_t::print(cout, tree->depth);
            cout << ",";
            stats_t::print(cout, tree->sequence_length);
            cout << ",";
            stats_t::print(cout, tree->suffix_length);
            cout << ",";
            cout << btimer;
            cout << ",";
            cout << ptimer;
            cout << endl;
            free_tree(tree);
            count += 1;
            total_time_build += btimer;
            total_time_process += ptimer;
        }
    }
    (void) time(&t2);
    printf("%zu non-empty trees constructed and processed in <%lld> secs\n",
            count, (long long)(t2-t1));
    cout << "total time build  " << total_time_build << endl;
    cout << "total time process" << total_time_process << endl;
#if USE_SET
    printf("%zu/%lu pairs generated\n", pairs.size(), ntasks);
    {
        /* generate statistics on how much was saved by filtering */
        set<pair<unsigned int,unsigned int> >::iterator it;
        vector<pair<unsigned int,unsigned int> > pairs_vec(pairs.begin(), pairs.end());
        unsigned long work_total = 0;
        unsigned long histo_width = 10000;
        unsigned long histo_size = sequences->max_seq_size * sequences->max_seq_size / histo_width;
        unsigned long *histo = new unsigned long[histo_size];
        (void)memset(histo, 0, sizeof(unsigned long)*histo_size);
        (void) time(&t1);
#pragma omp parallel for
        for (d = 0; d < (ssize_t)pairs_vec.size(); ++d) {
            size_t work = 0;
            unsigned long histo_index = 0;
            size_t s1Len = 0;
            size_t s2Len = 0;

            s1Len = sequences->seq[pairs_vec[d].first].size;
            s2Len = sequences->seq[pairs_vec[d].second].size;
            work = s1Len * s2Len;
            histo_index = work / histo_width;
            histo[histo_index] += 1;
#pragma omp atomic
            work_total += work;
        }
        (void) time(&t2);
        printf("pairs analyzed in <%lld> secs\n", (long long)(t2-t1));
        printf("work_total=%lu\n", work_total);
#if PRINT_HISTO
        printf("histogram\n");
        printf("%lu", histo[0]);
        for (i=1; i<histo_size; ++i) {
            printf(",%lu", histo[i]);
        }
        printf("\n");
#endif
        delete [] histo;
    }
    {
        /* generate statistics on how much was saved by filtering */
        unsigned long work_total = 0;
        unsigned long histo_width = 10000;
        unsigned long histo_size = sequences->max_seq_size * sequences->max_seq_size / histo_width;
        unsigned long *histo = new unsigned long[histo_size];
        (void)memset(histo, 0, sizeof(unsigned long)*histo_size);
        (void) time(&t1);
#pragma omp parallel for
        for (d = 0; d < (ssize_t)ntasks; ++d) {
            unsigned long seq_id[2];
            size_t work = 0;
            unsigned long histo_index = 0;
            size_t s1Len = 0;
            size_t s2Len = 0;

            k_combination2(d, seq_id);
            s1Len = sequences->seq[seq_id[0]].size;
            s2Len = sequences->seq[seq_id[1]].size;
            work = s1Len * s2Len;
            histo_index = work / histo_width;
            histo[histo_index] += 1;
#pragma omp atomic
            work_total += work;
        }
        (void) time(&t2);
        printf("filter analyzed in <%lld> secs\n", (long long)(t2-t1));
        printf("work_total=%lu\n", work_total);
#if PRINT_HISTO
        printf("histogram\n");
        printf("%lu", histo[0]);
        for (i=1; i<histo_size; ++i) {
            printf(",%lu", histo[i]);
        }
        printf("\n");
        delete [] histo;
#endif
    }
#else
    /* generate statistics on how much was saved by filtering */
    unsigned long work_no = 0;
    unsigned long work_yes = 0;
    unsigned long skipped = 0;
    unsigned long histo_width = 10000;
    unsigned long histo_size = sequences->max_seq_size * sequences->max_seq_size / histo_width;
    unsigned long *histo = new unsigned long[histo_size];
    (void)memset(histo, 0, sizeof(unsigned long)*histo_size);
    int cutOff = param.AOL * param.SIM;
    printf("      ntasks=%lu\n", ntasks);
    printf("n_triangular=%zu\n", n_triangular);
    (void) time(&t1);
#pragma omp parallel for
    for (d = 0; d < (ssize_t)ntasks; ++d) {
        unsigned long seq_id[2];
        size_t index = 0;
        size_t work = 0;
        unsigned long histo_index = 0;
        size_t s1Len = 0;
        size_t s2Len = 0;

        k_combination2(d, seq_id);
        s1Len = sequences->seq[seq_id[0]].size;
        s2Len = sequences->seq[seq_id[1]].size;
        index = (n_sequences*seq_id[0]) + seq_id[1] - (seq_id[0]*(seq_id[0]+1)/2);
        work = s1Len * s2Len;
        histo_index = work / histo_width;
        if ((s1Len <= s2Len && (100 * s1Len < cutOff * s2Len))
                || (s2Len < s1Len && (100 * s2Len < cutOff * s1Len))) {
            dup[index] = NO;
        }
        switch (dup[index]) {
            case NO:
            case MAYBE:
#pragma omp atomic
                work_no += work;
#pragma omp atomic
                skipped += 1;
#pragma omp atomic
                histo[histo_index] += 1;
                break;
            case YES:
#pragma omp atomic
                work_yes += work;
                break;
            default:
                assert(0);
        }
    }
    (void) time(&t2);
    printf("filter analyzed in <%lld> secs\n", (long long)(t2-t1));
    printf("skipped %zu alignments out of %zu (%f%%)\n",
            skipped, ntasks, 100.0*double(skipped)/double(ntasks));
    printf(" work_no=%zu\n", work_no);
    printf("work_yes=%zu\n", work_yes);
    printf("work saved is %f%%\n", 100.0*double(work_no)/double(work_no+work_yes));
#if PRINT_HISTO
    printf("histogram\n");
    printf("%lu", histo[0]);
    for (i=1; i<histo_size; ++i) {
        printf(",%lu", histo[i]);
    }
    printf("\n");
#endif

    /* now use a simple length-based cutoff filter */
    printf("using a simple length-based cutoff filter\n");
    ntasks = binomial_coefficient(n_sequences, 2);
    work_no = 0;
    work_yes = 0;
    skipped = 0;
    (void)memset(histo, 0, sizeof(unsigned long)*histo_size);
    (void) time(&t1);
    for (d = 0; d < (ssize_t)ntasks; ++d) {
        unsigned long seq_id[2];
        size_t index = 0;
        size_t work = 0;
        unsigned long histo_index = 0;
        size_t s1Len = 0;
        size_t s2Len = 0;

        k_combination2(d, seq_id);
        s1Len = sequences->seq[seq_id[0]].size;
        s2Len = sequences->seq[seq_id[1]].size;
        index = (n_sequences*seq_id[0]) + seq_id[1] - (seq_id[0]*(seq_id[0]+1)/2);
        work = s1Len * s2Len;
        histo_index = work / histo_width;
        if ((s1Len <= s2Len && (100 * s1Len < cutOff * s2Len))
                || (s2Len < s1Len && (100 * s2Len < cutOff * s1Len))) {
            work_no += work;
            skipped += 1;
            histo[histo_index] += 1;
        }
        else {
            work_yes += work;
        }
    }
    (void) time(&t2);
    printf("filter analyzed in <%lld> secs\n", (long long)(t2-t1));
    printf("skipped %zu alignments out of %zu (%f%%)\n",
            skipped, ntasks, 100.0*double(skipped)/double(ntasks));
    printf(" work_no=%zu\n", work_no);
    printf("work_yes=%zu\n", work_yes);
    printf("work saved is %f%%\n", 100.0*double(work_no)/double(work_no+work_yes));
#if PRINT_HISTO
    printf("histogram\n");
    printf("%lu", histo[0]);
    for (i=1; i<histo_size; ++i) {
        printf(",%lu", histo[i]);
    }
    printf("\n");
#endif

    delete [] dup;
#endif
    pg_free_sequences(sequences);
    free_suffix_buckets(suffix_buckets);

    return EXIT_SUCCESS; 
}


static void parse_command_line(int argc, char **argv,
        char *sequence_file, char *config_file, size_t *n_sequences)
{
	int option;
	int cnt = 0;
	
	while (-1 != (option = getopt(argc, argv, "f:n:c:"))) {
		switch (option) {
		    case '?':
                printf("Usage : %s -f {fasta} -n {#seqs} -c {cfg file}\n",
                        argv[0]);
			    exit(-1);
		    case 'f':
			    strcpy(sequence_file, optarg);
			    cnt++;
			    break;
            case 'c':
                strcpy(config_file, optarg);
                cnt++;
                break;
		    case 'n':
                *n_sequences = atol(optarg);
			    cnt++;
			    break;
		}
	}

    if (cnt != 3) {
        printf("Usage : %s -f {fasta} -n {#seqs} -c {cfg file}\n", argv[0]);
        exit(0);
    }
}


static inline size_t zpower(size_t base, size_t n)
{
    size_t p = 1;

    for (/**/; n > 0; --n) {
        assert(p < SIZE_MAX / base);
        p *= base;
    }

    return p;
}

