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

#include <sys/types.h>
#include <regex.h>

#include <ostream>
#include <string>
#include <vector>

#include <mpi.h>

using ::std::ostream;
using ::std::string;
using ::std::vector;

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
    static const string KEY_USE_COUNTER;
    static const string KEY_USE_TREE;
    static const string KEY_USE_TREE_DYNAMIC;
    static const string KEY_USE_TREE_HYBRID;
    static const string KEY_PRINT_STATS;
    static const string KEY_ALPHABET;
    static const string KEY_ALPHABET_BEGIN;
    static const string KEY_ALPHABET_DOLLAR;
    static const string KEY_DUPLICATES_LOCAL;
    static const string KEY_DUPLICATES_SEMILOCAL;
    static const string KEY_DUPLICATES_SMP;
    static const string KEY_DUPLICATES_GLOBAL;
    static const string KEY_BUCKET_CUTOFF;
    static const string KEY_SKIP_TREE;
    static const string KEY_PERFORM_ALIGNMENTS;

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
    static const bool DEF_USE_COUNTER;
    static const bool DEF_USE_TREE;
    static const bool DEF_USE_TREE_DYNAMIC;
    static const bool DEF_USE_TREE_HYBRID;
    static const bool DEF_PRINT_STATS;
    static const string DEF_ALPHABET;
    static const char DEF_ALPHABET_BEGIN;
    static const char DEF_ALPHABET_DOLLAR;
    static const bool DEF_DUPLICATES_LOCAL;
    static const bool DEF_DUPLICATES_SEMILOCAL;
    static const bool DEF_DUPLICATES_SMP;
    static const bool DEF_DUPLICATES_GLOBAL;
    static const size_t DEF_BUCKET_CUTOFF;
    static const bool DEF_SKIP_TREE;
    static const bool DEF_PERFORM_ALIGNMENTS;

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
        for (size_t i=0; i<re.size(); ++i) {
            regfree(re[i]);
            delete re[i];
            re[i] = NULL;
        }
        re.clear();
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
    vector<regex_t*> re; /**< regex(s) of buckets to skip */
    bool output_all;    /**< whether to output all results instead of only edges */
    bool output_to_disk;/**< whether to write results to disk (no==debugging) */
    bool distribute_sequences;/**< whether to allow sequence DB to distribute */
    bool use_length_filter; /**< whether to skip obviously bad alignments */
    bool use_iterator;  /**< whether to use tascel iterator */
    bool use_counter;   /**< whether to use ARMCI task counter */
    bool use_tree;      /**< whether to use tree filter as separate phase  */
    bool use_tree_dynamic; /**< whether to use tree filter as-needed */
    bool use_tree_hybrid;  /**< whether to allow tree tasks to be stolen */
    bool print_stats;   /**< whether to print detailed timings */
    string alphabet;    /**< the valid characters expected in input file */
    char alphabet_begin;/**< character representing left terminal */
    char alphabet_dollar;/**< character representing right terminal */
    bool dup_local;     /**< whether to check for duplicate pairs */
    bool dup_semilocal; /**< whether to check for duplicate pairs across invocations */
    bool dup_smp;       /**< whether to check for duplicate pairs across invocations, thread safe*/
    bool dup_global;    /**< whether to check for duplicate pairs distributed */
    size_t bucket_cutoff; /**< how many stddev above bucket size to discard */
    bool skip_tree;     /**< don't use tree if cutoff == exact match length */
    bool perform_alignments; /**< when debugging, sometimes useful to not align */
};

ostream& operator<< (ostream &os, const Parameters &p);

}; /* namespace pgraph */

#endif /* _PGRAPH_PARAMETERS_H_ */
