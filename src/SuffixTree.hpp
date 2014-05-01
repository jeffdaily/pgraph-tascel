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
#include <vector>
#include <utility>

#include "Parameters.hpp"
#include "SequenceDatabase.hpp"
#include "Stats.hpp"
#include "Suffix.hpp"
#include "SuffixBuckets.hpp"

using ::std::pair;
using ::std::set;
using ::std::vector;
using ::std::size_t;

namespace pgraph {

/**
 * suffix tree node
 */
class SuffixTreeNode
{
    public:
        size_t depth;   /**< depth since the root, not including initial size k */
        size_t rLeaf;   /**< right most leaf index */
        Suffix **lset;  /**< subtree's nodes branched according to left */
};


typedef void(*SuffixTreePairCallback)(pair<size_t,size_t>);


/**
 * suffix tree containing suffix tree nodes and other data
 */
class SuffixTree
{
    public:
        struct VectorCallback {
            vector<pair<size_t,size_t> > &pairs;
            VectorCallback(vector<pair<size_t,size_t> > &pairs)
                : pairs(pairs) {}
            bool operator()(const pair<size_t,size_t> &p) {
                pairs.push_back(p);
                return false;
            }
        };

        struct SetCallback {
            set<pair<size_t,size_t> > &pairs;
            SetCallback(set<pair<size_t,size_t> > &pairs)
                : pairs(pairs) {}
            bool operator()(const pair<size_t,size_t> &p) {
                pairs.insert(p);
                return false;
            }
        };

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
         * @param[in] k kmer size, if different than param.window_size
         * @return the suffix tree
         */
        SuffixTree(SequenceDatabase *sequences,
                   Bucket *bucket,
                   const Parameters &param,
                   int k=-1);

        /** Destructor. */
        ~SuffixTree();

        /**
         * Generate promising pairs for alignment.
         *
         * @param[out] pairs
         */
        bool generate_pairs(set<pair<size_t,size_t> > &pairs) {
            return generate_pairs(SetCallback(pairs));
        }

        /**
         * Generate promising pairs for alignment.
         *
         * @param[out] pairs
         */
        bool generate_pairs(vector<pair<size_t,size_t> > &pairs) {
            return generate_pairs(VectorCallback(pairs));
        }

        /**
         * Generate promising pairs for alignment.
         *
         * @param[out] pairs
         */
        template <class Callback>
        bool generate_pairs(Callback callback);

        size_t get_size() const { return size; }

        size_t get_size_internal() const { return size_internal; }

        const Stats& get_fanout() const { return fanout_stats; }

        const Stats& get_depth() const { return depth_stats; }

        const Stats& get_suffix_length() const { return suffix_length_stats; }

        static bool length_filter(size_t s1Len, size_t s2Len, size_t cutOff);

    private:
        size_t build_tree_recursive(Suffix *suffixes, int depth);
        size_t build_tree_recursive2(Suffix *suffixes, int depth);
        void compute_lset(Suffix *suffixes, Suffix **lset);
        int next_diff_pos(Suffix *suffixes, int depth);
        int next_diff_pos2(Suffix *suffixes, int depth, int size);
        bool is_candidate(Suffix *p, Suffix *q);
        size_t* count_sort();
        void merge_lsets(size_t sIndex, size_t eIndex);

        template <class Callback>
        bool proc_leaf(Suffix **lset, Callback callback);

        template <class Callback>
        bool proc_node(size_t sIndex, size_t eIndex, Callback callback);

        SequenceDatabase *sequences;
        //Bucket *bucket;
        const Parameters &param;/**< user parameters */
        const size_t SIGMA;     /**< alphabet size */
        const char DOLLAR;      /**< terminal character */
        const char BEGIN;       /**< leftmost terminal character */
        const string alphabet;  /**< the alphabet */
        vector<size_t> alphabet_table;/** lookup table for alphabet index */
        int window_size;        /**< sliding window size */
        SuffixTreeNode *nodes;  /**< tree nodes */
        size_t size;            /**< number of nodes */
        size_t size_internal;   /**< number of internal nodes */
        Suffix **lset_array;    /**< memory for all node's lsets (SIGMA*nnodes) */
        Stats fanout_stats;     /**< average fanout */
        Stats depth_stats;      /**< node depth stats */
        Stats suffix_length_stats;/**< suffix length stats */
        Suffix **tails;

        static const size_t npos;
};


template <class Callback>
bool SuffixTree::generate_pairs(Callback callback)
{
    size_t *srtIndex = count_sort();
    size_t EM = param.exact_match_length;

    assert(param.exact_match_length >= 1);

    /* srtIndex maintain an order of NON-increasing depth of nodes[] */
    for (size_t i = 0; i < size; i++) {
        size_t sIndex = srtIndex[i];
        SuffixTreeNode *node = &nodes[sIndex];

#ifdef DEBUG
        printf("stNode->depth=%d, node->rLeaf=%ld, sIndex=%ld\n", node->depth, node->rLeaf, sIndex);
#endif

        if (node->depth >= EM - 1) {
            if (node->rLeaf == sIndex) { /* leaf node */
                if (proc_leaf(node->lset, callback)) return true;
            }
            else {                       /* internal node */
                size_t eIndex = node->rLeaf;
                if (proc_node(sIndex, eIndex, callback)) return true;
                merge_lsets(sIndex, eIndex);
            }
        }
        else {
            /* stnodes are sorted, so later part
             * will not satisfy EM cutoff */
            break;
        }
    }

    /* free */
    delete [] srtIndex;

    return false;
}


template <class Callback>
bool SuffixTree::proc_leaf(Suffix **lset, Callback callback)
{
    size_t i;
    size_t j;
    Suffix *p = NULL;
    Suffix *q = NULL;
    int cutOff;

    cutOff = param.AOL * param.SIM;

    for (i = 0; i < SIGMA; i++) {
        if (lset[i]) {
            if (i == alphabet_table[BEGIN]) { /* inter cross */
                for (p = lset[i]; p != NULL; p = p->next) {
                    for (q = p->next; q != NULL; q = q->next) {
                        if (is_candidate(p, q)) {
                            bool retval = false;
                            if (p->sid > q->sid) {
                                retval = callback(make_pair(q->sid, p->sid));
                            }
                            else {
                                retval = callback(make_pair(p->sid, q->sid));
                            }
                            if (retval) return true;
                        }
                    }
                }
            }

            /* intra cross */
            for (j = i + 1; j < SIGMA; j++) {
                if (lset[j]) {
                    for (p = lset[i]; p != NULL; p = p->next) {
                        for (q = lset[j]; q != NULL; q = q->next) {
                            if (is_candidate(p, q)) {
                                bool retval = false;
                                if (p->sid > q->sid) {
                                    retval = callback(make_pair(q->sid, p->sid));
                                }
                                else {
                                    retval = callback(make_pair(p->sid, q->sid));
                                }
                                if (retval) return true;
                            }
                        }
                    }
                }
            }
        }
    }

    return false;
}


template <class Callback>
bool SuffixTree::proc_node(size_t sIndex, size_t eIndex, Callback callback)
{
    size_t m;
    size_t n;
    size_t s;
    size_t t;
    Suffix *p = NULL;
    Suffix *q = NULL;

    /* pairs generation loop for internal node */
    for (m = sIndex + 1; m < eIndex; m = nodes[m].rLeaf+1) {
        for (n = nodes[m].rLeaf+1; n <= eIndex; n = nodes[n].rLeaf+1) {
            for (s = 0; s < SIGMA; s++) {
                if (nodes[m].lset[s]) {
                    for (t = 0; t < SIGMA; t++) {
                        if (nodes[n].lset[t]) {
                            if (s != t || s == alphabet_table[BEGIN]) {
                                for (p = nodes[m].lset[s]; p != NULL; p = p->next) {
                                    for (q = nodes[n].lset[t]; q != NULL; q = q->next) {
                                        if (is_candidate(p, q)) {
                                            bool retval = false;
                                            if (p->sid > q->sid) {
                                                retval = callback(make_pair(q->sid, p->sid));
                                            }
                                            else {
                                                retval = callback(make_pair(p->sid, q->sid));
                                            }
                                            if (retval) return true;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    return false;
}

}; /* namespace pgraph */

#endif /* _PGRAPH_SUFFIXTREE_H_ */
