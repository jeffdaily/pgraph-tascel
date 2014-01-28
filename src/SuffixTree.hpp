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

using std::pair;
using std::set;
using std::vector;
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


typedef void(*SuffixTreePairCallback)(pair<size_t,size_t>);


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

        /**
         * Generate promising pairs for alignment.
         *
         * @param[out] pairs
         */
        void generate_pairs(vector<pair<size_t,size_t> > &pairs);

        /**
         * Generate promising pairs for alignment.
         *
         * @param[out] pairs
         */
        template <class Callback>
        bool generate_pairs_cb(Callback callback);

        double get_size_internal() const { return size_internal; }

        double get_avg_fanout() const { return fanout; }

        double get_avg_depth() const { return avgdepth; }

        double get_deepest() const { return deepest; }

        double get_suffix_avg_length() const { return suffix_avg_length; }

        double get_suffix_max_length() const { return suffix_max_length; }

    private:
        void create();

        SequenceDatabase *sequences;
        Bucket *bucket;
        Parameters param;
        SuffixTreeNode *nodes;  /**< tree nodes */
        size_t size;            /**< number of nodes */
        Suffix **lset_array;    /**< memory for all node's lsets (SIGMA*nnodes) */
        double fanout;          /**< average fanout */
        double size_internal;   /**< number of internal nodes */
        double avgdepth;        /**< average node depth */
        double deepest;         /**< deepest node depth */
        double suffix_avg_length;   /**< average suffix length */
        double suffix_max_length;   /**< longest suffix */
};

/**
 * counting sort based on the depth of the tree nodes.
 *
 * @param stNodes -
 * @param srtIndex -
 * @param nStNodes -
 * @param maxSeqLen -
 * -----------------------------------------------------------*/
static inline void
count_sort(SuffixTreeNode *stNodes, int *srtIndex, size_t nStNodes, size_t maxSeqLen)
{
    size_t i;
    int *cnt = NULL;
    int depth;

    cnt = new int[maxSeqLen];
    if (NULL == cnt) {
        perror("count_sort: malloc cnt");
        exit(EXIT_FAILURE);
    }

    /* init counters */
    for (i = 0; i < maxSeqLen; i++) {
        cnt[i] = 0;
    }

    /* fill in counters */
    for (i = 0; i < nStNodes; i++) {
        depth = stNodes[i].depth; /* depth starts from 0 */
        assert(depth < maxSeqLen);
        cnt[depth]++;
    }

    /* suffix sum to sort in a reverse way */
    for (i = maxSeqLen - 2; i > 0; i--) {
        cnt[i] += cnt[i + 1];
    }
    cnt[0] += cnt[1];

    /* store the sorted index into srtIndex */
    for (i = 0; i < nStNodes; i++) {
        depth = stNodes[i].depth;
        srtIndex[cnt[depth] - 1] = i;
        cnt[depth]--;
    }

    delete [] cnt;
}

static inline int
is_candidate(SequenceDatabase *seqs, size_t nSeqs,
             Suffix *p, Suffix *q, Parameters param)
{
    size_t s1Len = 0;
    size_t s2Len = 0;
    int f1 = 0;
    int f2 = 0;
    int cutOff = param.AOL * param.SIM;
    char result = FALSE;

    f1 = p->sid;
    f2 = q->sid;
    assert(f1 < nSeqs);
    assert(f2 < nSeqs);

    if (f1 == f2) {
        result = FALSE;
    }
    else {
        if (f1 > f2) {
            int swap = f1;
            f1 = f2;
            f2 = swap;
        }
        s1Len = (*seqs)[f1].get_sequence_length() - 1;
        s2Len = (*seqs)[f2].get_sequence_length() - 1;
        if (s1Len <= s2Len) {
            if (100 * s1Len < cutOff * s2Len) {
                result = FALSE;
            }
            else {
                result = TRUE;
            }
        }
        else {
            if (100 * s2Len < cutOff * s1Len) {
                result = FALSE;
            }
            else {
                result = TRUE;
            }
        }
    }

    return result;
}


template <class Callback>
bool SuffixTree::generate_pairs_cb(Callback callback)
{
    SuffixTreeNode *stNodes = NULL;
    int *srtIndex = NULL;
    size_t nStNodes = 0;
    int nSeqs = 0;
    int maxSeqLen = 0;
    size_t i = 0;
    int j = 0;
    SuffixTreeNode *stnode = NULL;
    int sIndex;
    int eIndex;
    size_t m;
    size_t n;
    size_t s;
    size_t t;
    Suffix *p = NULL;
    Suffix *q = NULL;
    int EM;
    int cutOff; /* cut off value of filter 1 */

    srtIndex = new int[this->size];
    count_sort(this->nodes, srtIndex, this->size, sequences->get_max_length());
    stNodes = this->nodes;
    nStNodes = this->size;
    nSeqs = sequences->get_global_count();
    maxSeqLen = sequences->get_max_length();

    /* only two rows are allocated */
    assert(NROW == 2);

    EM = param.exact_match_len;
    cutOff = param.AOL * param.SIM;

    /* srtIndex maintain an order of NON-increasing depth of stNodes[] */
    for (i = 0; i < nStNodes; i++) {
        sIndex = srtIndex[i];
        stnode = &stNodes[sIndex];

#ifdef DEBUG
        printf("stNode->depth=%d, stnode->rLeaf=%ld, sIndex=%ld\n", stnode->depth, stnode->rLeaf, sIndex);
#endif

        bool retval = false;
        if (stnode->depth >= EM - 1) {
            if (stnode->rLeaf == sIndex) { /* leaf node */
                Suffix **lset = stnode->lset;
                SequenceDatabase *seqs = sequences;
                size_t i;
                size_t j;
                Suffix *p = NULL;
                Suffix *q = NULL;
                int cutOff;

                cutOff = param.AOL * param.SIM;

                for (i = 0; i < SIGMA; i++) {
                    if (lset[i]) {
                        if (i == BEGIN - 'A') { /* inter cross */
                            for (p = lset[i]; p != NULL; p = p->next) {
                                for (q = p->next; q != NULL; q = q->next) {
                                    if (TRUE == is_candidate(seqs, nSeqs, p, q, param)) {
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
                                        if (TRUE == is_candidate(seqs, nSeqs, p, q, param)) {
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
            else {                       /* internal node */
                eIndex = stnode->rLeaf;

                /* pairs generation loop for internal node */
                for (m = sIndex + 1; m < eIndex; m = stNodes[m].rLeaf+1) {
                    for (n = stNodes[m].rLeaf+1; n <= eIndex; n = stNodes[n].rLeaf+1) {
                        for (s = 0; s < SIGMA; s++) {
                            if (stNodes[m].lset[s]) {
                                for (t = 0; t < SIGMA; t++) {
                                    if (stNodes[n].lset[t]) {
                                        if (s != t || s == (BEGIN - 'A')) {
                                            for (p = stNodes[m].lset[s]; p != NULL; p = p->next) {
                                                for (q = stNodes[n].lset[t]; q != NULL; q = q->next) {
                                                    if (TRUE == is_candidate(
                                                                sequences, nSeqs, p, q, param)) {
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

                /* merge the lsets of subtree */
                for (m = 0; m < SIGMA; m++) {
                    p = NULL;
                    for (j = sIndex + 1; j <= eIndex; j = stNodes[j].rLeaf + 1) {
                        if ((q = stNodes[j].lset[m])) {

                            /* empty the subtree's ptrs array */
                            stNodes[j].lset[m] = NULL;
                            if (p == NULL) {
                                p = q;
                                stNodes[sIndex].lset[m] = q;
                            }
                            else {
                                p->next = q;
                            }

                            /* walk to the end */
                            while (p->next) {
                                p = p->next;
                            }
                        }
                    }
                }
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


}; /* namespace pgraph */

#endif /* _PGRAPH_SUFFIXTREE_H_ */
