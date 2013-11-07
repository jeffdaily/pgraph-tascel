/**
 * @file stree.hpp
 *
 * @author andy.cj.wu@gmail.com
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2010 Washington State University. All rights reserved.
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 */
#ifndef _PGRAPH_STREE_H_
#define _PGRAPH_STREE_H_

#include "bucket.hpp"

#include <iostream>
#include <set>
#include <string>
#include <vector>
#include <utility>

using std::pair;
using std::ostream;
using std::set;
using std::string;
using std::vector;

namespace pgraph {

/**
 * suffix tree node
 */
typedef struct _stnode_t {
    int fanout;
    int depth; /**< depth since the root, not including the initial size k */
    size_t rLeaf; /**< right most leaf index */
    suffix_t **lset;  /**< subtree's nodes branched according to left
                              characters */
    struct _stnode_t **nset; /**< set of immediate child nodes */
} stnode_t;

typedef struct _stats_t_ {
    vector<double> values;
    double sum;
    double min;
    double max;
    double mean;
    double mode;
    double median;
    double stddev;

    static void print_header(ostream &os, const string &prefix) {
        os << prefix << "_sum" << ",";
        os << prefix << "_min" << ",";
        os << prefix << "_max" << ",";
        os << prefix << "_mean" << ",";
        os << prefix << "_mode" << ",";
        os << prefix << "_median" << ",";
        os << prefix << "_stddev";
    }
    static void print(ostream &os, const _stats_t_ &stats) {
        os << stats.sum << ",";
        os << stats.min << ",";
        os << stats.max << ",";
        os << stats.mean << ",";
        os << stats.mode << ",";
        os << stats.median << ",";
        os << stats.stddev;
    }
} stats_t;


/**
 * suffix tree containing suffix tree nodes and other data
 */
typedef struct {
    stnode_t *nodes;        /**< tree nodes */
    size_t size;            /**< number of nodes */
    size_t size_internal;   /**< number of internal nodes */
    size_t size_leaf;       /**< number of leaf nodes */
    suffix_t **lset_array;  /**< memory for all node's lsets (SIGMA*nnodes) */
    stnode_t **nset_array;  /**< memory for all node's nsets (SIGMA*nnodes) */
    stats_t fanout; 
    stats_t depth; 
    stats_t sequence_length; 
    stats_t suffix_length; 
} stree_t;


/**
 * Builds a tree for the given bucket (list of suffixes).
 *
 * Each bucket will be branched into SIGMA different sub-buckets in a DFS way,
 * the total space will be (SIGMA*4*depth) bytes
 *
 * @note NOTE: 'suffixes' attribute of bucket parameter is modified after this
 * function, which also means the previous bucket does not exist at all
 *
 * @param[in] sequences all fasta sequences
 * @param[in] bucket linked list of suffixes for this bucket
 * @param[in] param alignment parameters
 * @return the suffix tree
 */
stree_t* build_tree(
        sequences_t *sequences, bucket_t *bucket, param_t param);


/**
 * TODO
 *
 * @param[in] tree TODO
 */
void free_tree(stree_t *tree);


/**
 * Generate promising pairs for alignment.
 *
 * @param[in] tree
 * @param[in] sequences
 * @param[in] dup
 * @param[in] param
 */
void generate_pairs(
        stree_t *tree, sequences_t *sequences, int *dup, param_t param);


/**
 * Generate promising pairs for alignment.
 *
 * @param[in] tree
 * @param[in] sequences
 * @param[in] dup
 * @param[in] param
 */
void generate_pairs(
        stree_t *tree, sequences_t *sequences, set<pair<unsigned int,unsigned int> > &dup, param_t param);


}; /* namespace pgraph */

#endif /* _PGRAPH_STREE_H_ */
