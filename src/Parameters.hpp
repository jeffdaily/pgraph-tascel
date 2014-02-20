/**
 * @file Parameters.hpp
 *
 * @author andy.cj.wu@gmail.com
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2010 Washington State University. All rights reserved.
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 */
#ifndef _PGRAPH_PARAMETERS_H_
#define _PGRAPH_PARAMETERS_H_

#include <iostream>
#include <string>
#include <vector>

#include <mpi.h>

using std::ostream;
using std::string;
using std::vector;

namespace pgraph {

/**
 * Parameters for various parts of the pgraph software.
 */
class Parameters
{
public:
    /* Keys */
    static const string KEY_ALIGN_OVER_LONGER_SEQUENCE;
    static const string KEY_MATCH_SIMILARITY;
    static const string KEY_OPTIMAL_SCORE_OVER_SELF_SCORE;
    static const string KEY_EXACT_MATCH_LENGTH;
    static const string KEY_SLIDE_WINDOW_SIZE;
    static const string KEY_OPEN;
    static const string KEY_GAP;
    static const string KEY_MEMORY_WORKER;
    static const string KEY_MEMORY_SEQUENCES;
    static const string KEY_SKIP_PREFIXES;
    static const string KEY_OUTPUT_ALL;
    static const string KEY_OUTPUT_TO_DISK;
    static const string KEY_DISTRIBUTE_SEQUENCES;
    static const string KEY_USE_LENGTH_FILTER;
    static const string KEY_USE_ITERATOR;
    static const string KEY_PRINT_STATS;

    /* Defaults */
    static const int DEF_ALIGN_OVER_LONGER_SEQUENCE;
    static const int DEF_MATCH_SIMILARITY;
    static const int DEF_OPTIMAL_SCORE_OVER_SELF_SCORE;
    static const int DEF_EXACT_MATCH_LENGTH;
    static const int DEF_SLIDE_WINDOW_SIZE;
    static const int DEF_OPEN;
    static const int DEF_GAP;
    static const size_t DEF_MEMORY_WORKER;
    static const size_t DEF_MEMORY_SEQUENCES;
    static const vector<string> DEF_SKIP_PREFIXES;
    static const bool DEF_OUTPUT_ALL;
    static const bool DEF_OUTPUT_TO_DISK;
    static const bool DEF_DISTRIBUTE_SEQUENCES;
    static const bool DEF_USE_LENGTH_FILTER;
    static const bool DEF_USE_ITERATOR;
    static const bool DEF_PRINT_STATS;

    /**
     * Constructs empty (default) parameters.
     */
    Parameters();

    /**
     * Constructs and parses the given parameters file.
     *
     * @param[in] parameters_file the parameters file
     * @param[in] comm the communicator
     */
    Parameters(const char *parameters_file, MPI_Comm comm);

    ~Parameters() {
        skip_prefixes.clear();
    }

    /**
     * Parses the given parameters file.
     *
     * @param[in] parameters_file the parameters file
     */
    void parse(const char *parameters_file, MPI_Comm comm);

    int AOL;            /**< AlignOverLongerSequence */
    int SIM;            /**< MatchSimilarity */
    int OS;             /**< OptimalScoreOverSelfScore */
    int exact_match_length;/**< exact match length cutoff */
    int window_size;    /**< slide window size */
    int open;           /**< open penalty for affine gap alignment */
    int gap;            /**< gap extension penalty for affine gap alignment */
    size_t memory_worker; /**< memory budget per worker task pool */
    size_t memory_sequences; /**< memory budget for sequence database */
    vector<string> skip_prefixes; /**< buckets to skip */
    bool output_all;    /**< whether to output all results instead of only edges */
    bool output_to_disk;/**< whether to write results to disk (no==debugging) */
    bool distribute_sequences;/**< whether to allow sequence DB to distribute */
    bool use_length_filter; /**< whether to skip obviously bad alignments */
    bool use_iterator;  /**< whether to use tascel iterator */
    bool print_stats;   /**< whether to print detailed timings */
};

ostream& operator<< (ostream &os, const Parameters &p);

}; /* namespace pgraph */

#endif /* _PGRAPH_PARAMETERS_H_ */
