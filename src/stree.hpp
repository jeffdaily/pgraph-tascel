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

namespace pgraph {

/**
 * suffix tree node
 */
typedef struct {
    int depth; /**< depth since the root, not including the initial size k */
    size_t rLeaf; /**< right most leaf index */
    suffix_t **lset;  /**< subtree's nodes branched according to left
                              characters */
} stnode_t;


/**
 * suffix tree containing suffix tree nodes and other data
 */
typedef struct {
    stnode_t *nodes;        /**< tree nodes */
    size_t size;            /**< number of nodes */
    suffix_t **lset_array;  /**< memory for all node's lsets (SIGMA*nnodes) */
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


}; /* namespace pgraph */

#endif /* _PGRAPH_STREE_H_ */
