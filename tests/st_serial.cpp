/**
 * @author jeff.daily@pnnl.gov
 *
 * Alignment of an input dataset using work stealing. Each MPI task reads the
 * input file.
 */
#include "config.h"

#include <mpi.h>
#include <tascel.h>
#include <tascel/UniformTaskCollection.h>
#include <tascel/UniformTaskCollIter.h>
#include <tascel/UniformTaskCollectionSplit.h>
#if THREADED
#include "tascelx.hpp"
#endif

#include <algorithm>
#include <cassert>
#include <cctype>
#include <cfloat>
#include <cmath>
#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "alignment.hpp"
#include "AlignStats.hpp"
#include "Bootstrap.hpp"
#include "combinations.h"
#include "EdgeResult.hpp"
#include "mpix.hpp"
#include "mpix-types.hpp"
#include "PairCheck.hpp"
#include "PairCheckGlobal.hpp"
#include "PairCheckGlobalServer.hpp"
#include "PairCheckLocal.hpp"
#include "PairCheckSemiLocal.hpp"
#include "PairCheckSmp.hpp"
#include "Parameters.hpp"
#include "Sequence.hpp"
#include "SequenceDatabase.hpp"
#include "SequenceDatabaseReplicated.hpp"
#include "Stats.hpp"
#include "SuffixBuckets.hpp"
#include "SuffixBucketsTascel.hpp"
#include "SuffixTree.hpp"
#include "TreeStats.hpp"

#define ENABLE_ARMCI 0
#if HAVE_ARMCI && ENABLE_ARMCI
#include "SequenceDatabaseArmci.hpp"
#include "SuffixBucketsArmci.hpp"
#endif

#define printf(...) fprintf(stdout, __VA_ARGS__); fflush(stdout);

using namespace std;
using namespace tascel;
using namespace pgraph;


int main(int argc, char **argv)
{
    double time_main;
    int rank;
    int nprocs;
    vector<string> all_argv;
    TreeStats *stats_tree = NULL;
    SequenceDatabase *sequences = NULL;
    Parameters *parameters = NULL;
    char delimiter = '\0';

    pgraph::initialize(argc, argv);
    rank = mpix::comm_rank(pgraph::comm);
    nprocs = mpix::comm_size(pgraph::comm);
    time_main = MPI_Wtime();

    assert(nprocs == 1);

    /* initialize tascel */
    TascelConfig::initialize(NUM_WORKERS_DEFAULT, pgraph::comm);

    /* initialize global data */
    stats_tree = new TreeStats[NUM_WORKERS];
    parameters = new Parameters;

    /* MPI standard does not guarantee all procs receive argc and argv */
    all_argv = mpix::bcast(argc, argv, pgraph::comm);

    /* print the command line arguments */
    if (0 == rank) {
        ostringstream header;
        header.fill('-');
        header << left << setw(79) << "--- Command Line ";
        cout << header.str() << endl;
        for (size_t j=0; j<all_argv.size(); ++j) {
            cout << "argv[" << j << "]='" << all_argv[j] << "'" << endl;
        }
    }

    /* sanity check that we got the correct number of arguments */
    if (all_argv.size() <= 1 || all_argv.size() >= 4) {
        if (0 == rank) {
            if (all_argv.size() <= 1) {
                cout << "missing input file" << endl;
            }
            else if (all_argv.size() >= 4) {
                cout << "too many arguments" << endl;
            }
            cout << "usage: " << all_argv[0]
                << " sequence_file <config_file>" << endl;
        }
        TascelConfig::finalize();
        pgraph::finalize();
        return 1;
    }
    else if (all_argv.size() >= 3) {
        parameters->parse(all_argv[2].c_str(), pgraph::comm);
    }
    else if (all_argv.size() >= 2) {
        /* do nothing */
    }

    /* print parameters */
    if (0 == rank) {
        ostringstream header;
        header.fill('-');
        header << left << setw(79) << "-- Parameters ";
        cout << header.str() << endl;
        //cout << *parameters << endl;
        cout << "ExactMatchLength: " << parameters->exact_match_length << endl;
        cout << " SlideWindowSize: " << parameters->window_size << endl;
    }

    delimiter = parameters->alphabet_dollar;

    double time_seq = MPI_Wtime();
    sequences = new SequenceDatabaseReplicated(all_argv[1],
            parameters->memory_sequences, pgraph::comm, delimiter);
    time_seq = MPI_Wtime() - time_seq;
    if (0 == rank) {
        cout << "time sequence db open " << time_seq << " seconds" << endl;
    }

    /* create a single bucket containing all suffixes */
    double time_bucket = MPI_Wtime();
    Bucket *bucket = new Bucket;
    //size_t n_suffixes = sequences->char_size() - sequences->size();
    size_t n_suffixes = sequences->char_size();
    vector<Suffix> suffixes(n_suffixes);
    size_t suffix_index = 0;
    if (0 == rank) {
        cout << "n_suffixes = " << n_suffixes << endl;
    }
    for (size_t i=0; i<sequences->size(); ++i) {
        size_t stop_index = 0;
        const char *sequence_data = NULL;
        size_t sequence_length = 0;
        Sequence *sequence = sequences->get_sequence(i);

        sequence->get_sequence(sequence_data, sequence_length);
        /* stop_index stops before the assumed DOLLAR terminal character */
        stop_index = sequence_length - 2;
        for (size_t j = 0; j <= stop_index; ++j) {
            suffixes[suffix_index].sid = i;
            suffixes[suffix_index].pid = j;
            suffixes[suffix_index].bid = 0;
            suffixes[suffix_index].k = 0;
            suffixes[suffix_index].next = &suffixes[suffix_index+1];
            suffix_index++;
        }
        delete sequence;
    }
    if (0 == rank) {
        cout << "suffix_index = " << suffix_index << endl;
    }
    suffixes[suffix_index-1].next = NULL;
    bucket->suffixes = &suffixes[0];
    bucket->size = suffix_index;
    bucket->bid = 0;
    bucket->k = 0;
    bucket->owner = 0;
    time_bucket = MPI_Wtime() - time_bucket;
    if (0 == rank) {
        cout << "time bucket = " << time_bucket << " seconds" << endl;
    }

    double time_tree = MPI_Wtime();
    SuffixTree *tree = new SuffixTree(sequences, bucket, *parameters, 0);
    time_tree = MPI_Wtime() - time_tree;
#if PRINT
    tree->print();
#endif
    if (0 == rank) {
        cout << "time tree = " << time_tree << " seconds" << endl;
    }
    cout << "tree size " << tree->get_size() << endl;
    cout << "tree size internal " << tree->get_size_internal() << endl;

    double time_pairs = MPI_Wtime();
    vector<pair<size_t,size_t> > pairs;
    tree->generate_pairs(pairs);
    time_pairs = MPI_Wtime() - time_pairs;
#if 0
    double time_dup = MPI_Wtime();
    set<pair<size_t,size_t> > pairs_set(pairs.begin(), pairs.end());
    time_dup = MPI_Wtime() - time_dup;
    if (0 == rank) {
        cout << "time pairs = " << time_pairs << " seconds" << endl;
        cout << "time dup = " << time_dup << " seconds" << endl;
        for (set<pair<size_t,size_t> >::iterator pit=pairs_set.begin();
                pit!=pairs_set.end(); ++pit) {
            cout << "pair " << pit->first << "," << pit->second << endl;
        }
        cout << "pairs_set.size() = " << pairs_set.size() << endl;
        cout << "pairs.size() = " << pairs.size() << endl;
    }
#else
    if (0 == rank) {
        cout << "time pairs = " << time_pairs << " seconds" << endl;
        cout << "pairs generated = " << pairs.size() << endl;
        set<pair<size_t,size_t> > set_pairs(pairs.begin(),pairs.end());
        cout << "pairs unique    = " << set_pairs.size() << endl;
#if PRINT
        for (set<pair<size_t,size_t> >::iterator it=set_pairs.begin();
                it!=set_pairs.end(); ++it) {
            cout << it->first << "," << it->second << endl;
        }
#endif
    }
#endif

    delete [] stats_tree;
    delete bucket;
    delete parameters;
    delete sequences;
    delete tree;

    time_main = MPI_Wtime() - time_main;
    if (0 == rank) {
        cout << "time_main " << time_main << " seconds" << endl;
    }
    TascelConfig::finalize();
    pgraph::finalize();

    return 0;
}

