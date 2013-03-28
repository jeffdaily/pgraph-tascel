/**
 * @file ccd.c
 *
 * @author andy.cj.wu@gmail.com
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2010 Washington State University. All rights reserved.
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include "bucket.h"
#include "cfg.h"
#include "elib.h"
#include "lib.h"
#include "loadseq.h"
#include "pairs.h"
#include "search.h"
#include "stree.h"
#include "uFind.h"


static void parse_command_line(int argc, char **argv,
        char *sequence_file, char *config_file, size_t *n_sequences);


int main(int argc, char *argv[])
{
    char sequence_file[FILENAME_MAX];/* path to sequence (fasta) file */
    char config_file[FILENAME_MAX]; /* path to config file */
    size_t n_sequences = 0;         /* number of sequences (cmd line param) */
    sequence_t *sequences = NULL;   /* all sequences parsed from fasta file */
    size_t maxSeqLen = 0;           /* longest sequence length parsed */
    size_t n_acids = 0;             /* number of amino acids parsed */
    bucket_t *buckets = NULL;       /* all buckets */
    size_t n_buckets;               /* number of buckets */
    suffix_t *suffixes = NULL;      /* all suffixes */
    size_t n_suffixes = 0;          /* number of suffixes */
    ufind_t *union_set = NULL;      /* TODO */
    int singletons = 2;             /* TODO */
    param_t param;                  /* config file parameters */
    time_t t1 = 0;                  /* start timer */
    time_t t2 = 0;                  /* stop timer */

    setbuf(stdout, (char *)0);

    /* reading in pamameters */
    parse_command_line(argc, argv, sequence_file, config_file, &n_sequences);

    /* read in configurations */
    param.window_size = get_config_val(config_file, "SlideWindowSize");
    param.exact_match_len = get_config_val(config_file, "ExactMatchLen");
    param.AOL = get_config_val(config_file, "AlignOverLongerSeq");
    param.SIM = get_config_val(config_file, "MatchSimilarity");
    param.OS = get_config_val(config_file, "OptimalScoreOverSelfScore");
    
    /* initialize dynamic programming blosum table */
    init_map(SIGMA);

    printf("----------------------------------------------\n");
    printf("%-15s: %.32s\n", "fasta seq", sequence_file);
    printf("%-15s: %.32s\n", "config file", config_file);
    printf("%-15s: %d\n", "slide size", param.window_size);
    printf("%-15s: %d\n", "exactMatch len", param.exact_match_len);
    printf("----------------------------------------------\n");

    sequences = emalloc(n_sequences*sizeof(sequence_t));

    /* load sequences */
    (void) time(&t1);
    load_all_sequences(sequence_file, n_sequences, sequences,
            &n_acids, &maxSeqLen);
    (void) time(&t2);
    printf("<%zu> amino acids are loaded in <%lld> secs\n",
            n_acids, (long long)(t2-t1));

    /* allocate memory for *buckets, *suffixes */
    n_buckets = zpower(SIGMA, param.window_size);
    buckets = emalloc(n_buckets*sizeof(bucket_t));
    /* n_acids includes '$' */
    n_suffixes = n_acids - n_sequences*param.window_size;
    suffixes = emalloc(n_suffixes*sizeof(suffix_t));

    /* build the buckets */
    (void) time(&t1);
    build_buckets(
            sequences, n_sequences,
            buckets,   n_buckets,
            suffixes,  n_suffixes,
            param.window_size); 
    (void) time(&t2);
    printf("Bucketing finished in <%lld> secs\n", (long long)(t2-t1));
    
    (void) time(&t1);
    count_sort_buckets(buckets, n_buckets);
    //qsort(buckets, n_buckets, sizeof *buckets, bktCmp);
    (void) time(&t2);
    printf("Bucketing sorted in <%lld> secs\n", (long long)(t2-t1));

    #ifdef DEBUG
    {
        int i;
        for (i=0; i<100; ++i) {
            printf("n_buckets[%d] - %d\n", i, buckets[i].size);
        }
        exit(0);
    }
    #endif

    union_set = init_union(n_sequences);

    /* suffix tree construction & processing */
    (void) time(&t1);
    build_forest(
            buckets,   n_buckets,
            sequences, n_sequences,
            maxSeqLen, union_set, &param);
    (void) time(&t2);
    printf("Tree constructed and processed in <%lld> secs\n",
            (long long)(t2-t1));
    print_pairs();
    
    disp_all_clusters(union_set, n_sequences, &singletons, ".", sequences); 

    /* free mem */
    free(suffixes);
    free(buckets);
    free_sequences(sequences, n_sequences);
    free(sequences);
    free_union(union_set);

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
