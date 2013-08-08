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
#include <sstream>

#include <unistd.h>

#include <mpi.h>

#include "constants.h"
#include "tascelx.hpp"
#include "mpix.hpp"
#include "SequenceDatabase.hpp"
#include "SequenceDatabaseArmci.hpp"
#include "SuffixBuckets.hpp"
#include "SuffixTree.hpp"

using namespace pgraph;
using namespace std;

static void parse_command_line(int argc, char **argv,
        char *sequence_file, char *config_file,
        size_t &n_sequences, size_t &budget);

static int comm_rank = 0;              /* communicator rank */
static int comm_size = 0;              /* communicator size */

int main(int argc, char *argv[])
{
    MPI_Comm comm = MPI_COMM_NULL;  /* communicator */
    int status = 0;                 /* return value from MPI routines */
    size_t budget = 0;              /* memory budget */
    char sequence_file[FILENAME_MAX];/* path to sequence (fasta) file */
    char config_file[FILENAME_MAX]; /* path to config file */
    size_t n_sequences = 0;         /* number of sequences (cmd line param) */
    SequenceDatabase *sequences = NULL; /* all sequences parsed from fasta file */
    size_t maxSeqLen = 0;           /* longest sequence length parsed */
    time_t t1 = 0;                  /* start timer */
    time_t t2 = 0;                  /* stop timer */
    SuffixBuckets *suffix_buckets = NULL;
    size_t i = 0;                   /* for loop index */
    size_t n_triangular = 0;        /* number of possible pairs */
    int *dup = NULL;                /* track duplicate pairs */

    /* initialize MPI */
    MPI_Init(&argc, &argv);
    status = MPI_Comm_dup(MPI_COMM_WORLD, &comm);
    assert(MPI_SUCCESS == status);
    status = MPI_Comm_rank(comm, &comm_rank);
    assert(MPI_SUCCESS == status);
    status = MPI_Comm_size(comm, &comm_size);
    assert(MPI_SUCCESS == status);

    setbuf(stdout, (char *)0);

    /* reading in pamameters */
    parse_command_line(argc, argv, sequence_file, config_file,
            n_sequences, budget);

    /* read in configurations */
    Parameters parameters(config_file, comm);
    
    if (0 == comm_rank) {
        printf("----------------------------------------------\n");
        printf("%-20s: %.78s\n", "fasta seq", sequence_file);
        printf("%-20s: %.78s\n", "config file", config_file);
        printf("%-20s: %d\n", "slide size", parameters.window_size);
        printf("%-20s: %d\n", "exactMatch len", parameters.exact_match_len);
        printf("%-20s: %d\n", "AlignOverLongerSeq", parameters.AOL);
        printf("%-20s: %d\n", "MatchSimilarity", parameters.SIM);
        printf("%-20s: %d\n", "OptimalScoreOverSelfScore", parameters.OS);
        printf("%-20s: %zu\n", "nseqs", n_sequences);
        printf("----------------------------------------------\n");
    }

    (void) time(&t1);
    sequences = new SequenceDatabaseArmci(sequence_file, budget, comm, 1, DOLLAR);
    (void) time(&t2);
    if (0 == comm_rank) {
        printf("<%zu> amino acids (<%zu> chars) are loaded in <%lld> secs\n",
                sequences->get_global_count(),
                sequences->get_global_size(),
                (long long)(t2-t1));
    }
    assert(sequences->get_global_count() == n_sequences);

    (void) time(&t1);
    suffix_buckets = new SuffixBuckets(sequences, parameters, comm);
    (void) time(&t2);
    if (0 == comm_rank) {
        printf("Bucketing finished in <%lld> secs\n", (long long)(t2-t1));
    }
    
    n_triangular = sequences->get_global_count() *
            (sequences->get_global_count() + 1U) / 2U;
    dup = new int[n_triangular];
    for (i = 0; i < n_triangular; ++i) {
        dup[i] = 2;
    }
    
    /* suffix tree construction & processing */
    (void) time(&t1);
    size_t count = 0;
    /* TODO */
    for (i = 0; i < suffix_buckets->buckets_size; ++i) {
        if (NULL != suffix_buckets->buckets[i].suffixes) {
            SuffixTree *tree = new SuffixTree(
                    sequences, &(suffix_buckets->buckets[i]), parameters);
            tree->generate_pairs(dup);
            delete tree;
            ++count;
        }
    }
    (void) time(&t2);
    mpix_reduce(count, MPI_SUM);
    delete [] dup;
    if (0 == comm_rank) {
        printf("%zu non-empty trees constructed and processed in <%lld> secs\n",
                count, (long long)(t2-t1));
    }
    
    delete sequences;
#if 0
    free_suffix_buckets(suffix_buckets);
#endif

    return EXIT_SUCCESS; 
}


const char *USAGE = "Usage : %s -f {fasta} -n {#seqs} -c {cfg file} -m {memory budget}\n";


static void parse_command_line(int argc, char **argv,
        char *sequence_file, char *config_file, size_t &n_sequences,
        size_t &budget)
{
	int option;
	int cnt = 0;
    char budget_multiplier = 'b';
	
	while (-1 != (option = getopt(argc, argv, "f:n:c:m:"))) {
		switch (option) {
		    case '?':
                if (0 == comm_rank) {
                    printf(USAGE, argv[0]);
                }
			    exit(EXIT_FAILURE);
		    case 'f':
			    strcpy(sequence_file, optarg);
			    cnt++;
			    break;
            case 'c':
                strcpy(config_file, optarg);
                cnt++;
                break;
		    case 'n':
                n_sequences = atol(optarg);
			    cnt++;
			    break;
            case 'm':
                istringstream iss(optarg);
                iss >> budget >> budget_multiplier;
                cnt++;
                if (budget <= 0) {
                    if (0 == comm_rank) {
                        cerr << "memory budget must be positive real number" << endl;
                    }
                    exit(EXIT_FAILURE);
                }
                if (budget_multiplier == 'b' || budget_multiplier == 'B') {
                    budget *= 1; /* byte */
                }
                else if (budget_multiplier == 'k' || budget_multiplier == 'K') {
                    budget *= 1024; /* kilobyte */
                }
                else if (budget_multiplier == 'm' || budget_multiplier == 'M') {
                    budget *= 1048576; /* megabyte */
                }
                else if (budget_multiplier == 'g' || budget_multiplier == 'G') {
                    budget *= 1073741824; /* gigabyte */
                }
                break;
        }
    }

    if (cnt != 4) {
        if (0 == comm_rank) {
            printf(USAGE, argv[0]);
        }
        exit(EXIT_FAILURE);
    }
}

