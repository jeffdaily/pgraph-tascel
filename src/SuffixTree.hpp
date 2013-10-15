/**
 * @file SuffixTree.hpp
 *
 * @author andy.cj.wu@gmail.com
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2010 Washington State University. All rights reserved.
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 */
#ifndef _PGRAPH_SUFFIXTREE_H_
#define _PGRAPH_SUFFIXTREE_H_

#include <cstddef>
#include <set>
#include <utility>

using std::pair;
using std::set;
using std::size_t;

namespace pgraph {

#include "Parameters.hpp"

class Bucket;
class SequenceDatabase;
class Suffix;

/**
 * suffix tree node
 */
class SuffixTreeNode
{
    public:
        int depth;      /**< depth since the root, not including initial size k */
        size_t rLeaf;   /**< right most leaf index */
        Suffix **lset;  /**< subtree's nodes branched according to left */
};


/**
 * suffix tree containing suffix tree nodes and other data
 */
class SuffixTree
{
    public:
        /**
         * Builds a tree for the given bucket (list of suffixes).
         *
         * Each bucket will be branched into SIGMA different sub-buckets in a
         * DFS way, the total space will be (SIGMA*4*depth) bytes
         *
         * @note NOTE: 'suffixes' attribute of bucket parameter is modified
         * after this function, which also means the previous bucket does not
         * exist at all
         *
         * @param[in] sequences all fasta sequences
         * @param[in] bucket linked list of suffixes for this bucket
         * @param[in] param alignment parameters
         * @return the suffix tree
         */
        SuffixTree(SequenceDatabase *sequences,
                   Bucket *bucket,
                   const Parameters &param);

        /** Destructor. */
        ~SuffixTree();

        /**
         * Generate promising pairs for alignment.
         *
         * @param[out] pairs
         */
        void generate_pairs(set<pair<size_t,size_t> > &pairs);

    private:
        void create();

        SequenceDatabase *sequences;
        Bucket *bucket;
        Parameters param;
        SuffixTreeNode *nodes;  /**< tree nodes */
        size_t size;            /**< number of nodes */
        Suffix **lset_array;    /**< memory for all node's lsets (SIGMA*nnodes) */
};

}; /* namespace pgraph */

#endif /* _PGRAPH_SUFFIXTREE_H_ */
