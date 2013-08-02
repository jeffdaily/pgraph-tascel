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

#include <unistd.h>

#include "bucket.hpp"
#include "constants.h"
#include "csequence.h"
#include "stree.hpp"

using namespace pgraph;

static void parse_command_line(int argc, char **argv,
        char *sequence_file, char *config_file, size_t *n_sequences);
static inline size_t zpower(size_t base, size_t n);


int main(int argc, char *argv[])
{
    char sequence_file[FILENAME_MAX];/* path to sequence (fasta) file */
    char config_file[FILENAME_MAX]; /* path to config file */
    size_t n_sequences = 0;         /* number of sequences (cmd line param) */
    sequences_t *sequences = NULL;  /* all sequences parsed from fasta file */
    size_t maxSeqLen = 0;           /* longest sequence length parsed */
    bucket_t *buckets = NULL;       /* all buckets */
    size_t n_buckets;               /* number of buckets */
    suffix_t *suffixes = NULL;      /* all suffixes */
    size_t n_suffixes = 0;          /* number of suffixes */
    param_t param;                  /* config file parameters */
    time_t t1 = 0;                  /* start timer */
    time_t t2 = 0;                  /* stop timer */
    suffix_buckets_t *suffix_buckets = NULL;
    size_t i = 0;                   /* for loop index */
    size_t n_triangular = 0;        /* number of possible pairs */
    int *dup = NULL;                /* track duplicate pairs */

    setbuf(stdout, (char *)0);

    /* reading in pamameters */
    parse_command_line(argc, argv, sequence_file, config_file, &n_sequences);

    /* read in configurations */
    get_params(config_file, &param);
    
    printf("----------------------------------------------\n");
    printf("%-15s: %.78s\n", "fasta seq", sequence_file);
    printf("%-15s: %.78s\n", "config file", config_file);
    printf("%-15s: %d\n", "slide size", param.window_size);
    printf("%-15s: %d\n", "exactMatch len", param.exact_match_len);
    printf("%-15s: %zu\n", "nseqs", n_sequences);
    printf("----------------------------------------------\n");

    (void) time(&t1);
    sequences = pg_load_fasta(sequence_file, DOLLAR);
    (void) time(&t2);
    printf("<%zu> amino acids are loaded in <%lld> secs\n",
            sequences->size, (long long)(t2-t1));

    (void) time(&t1);
    suffix_buckets = create_suffix_buckets(sequences, param);
    (void) time(&t2);
    printf("Bucketing finished in <%lld> secs\n", (long long)(t2-t1));
    
#if 0
    (void) time(&t1);
    count_sort_buckets(buckets, n_buckets);
    //qsort(buckets, n_buckets, sizeof *buckets, bktCmp);
    (void) time(&t2);
    printf("Bucketing sorted in <%lld> secs\n", (long long)(t2-t1));
#endif

    #ifdef DEBUG
    {
        int i;
        for (i=0; i<100; ++i) {
            printf("n_buckets[%d] - %d\n", i, buckets[i].size);
        }
        exit(0);
    }
    #endif

    n_triangular = sequences->size * (sequences->size + 1U) / 2U;
    dup = new int[n_triangular];
    for (i = 0; i < n_triangular; ++i) {
        dup[i] = 2;
    }
    
    /* suffix tree construction & processing */
    (void) time(&t1);
    size_t count = 0;
    for (i = 0; i < suffix_buckets->buckets_size; ++i) {
        if (NULL != suffix_buckets->buckets[i].suffixes) {
            stree_t *tree = build_tree(
                    sequences, &(suffix_buckets->buckets[i]), param);
            generate_pairs(tree, sequences, dup, param);
            free_tree(tree);
            ++count;
        }
    }
    (void) time(&t2);
    free(dup);
    printf("%zu non-empty trees constructed and processed in <%lld> secs\n",
            count, (long long)(t2-t1));
    
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

