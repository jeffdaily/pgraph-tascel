/**
 * @file SuffixArray.hpp
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2010 Washington State University. All rights reserved.
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 */
#ifndef _PGRAPH_SUFFIXARRAY_H_
#define _PGRAPH_SUFFIXARRAY_H_

#include <climits>
#include <cstddef>
#include <map>
#include <set>
#include <stack>
#include <vector>
#include <utility>

#include "Parameters.hpp"
#include "Pair.hpp"
#include "SequenceDatabase.hpp"
#include "Stats.hpp"
#include "Suffix.hpp"
#include "SuffixBuckets.hpp"

using ::std::map;
using ::std::pair;
using ::std::set;
using ::std::stack;
using ::std::vector;
using ::std::size_t;

namespace pgraph {

struct quad {
    int lcp;
    int lb;
    int rb;
    vector<quad> children;

    quad()
        : lcp(0), lb(0), rb(INT_MAX), children() {}
    quad(int lcp, int lb, int rb)
        : lcp(lcp), lb(lb), rb(rb), children() {}
    quad(int lcp, int lb, int rb, vector<quad> children)
        : lcp(lcp), lb(lb), rb(rb), children(children) {}

    bool empty() { return rb == INT_MAX; }
};

typedef void(*SuffixArrayPairCallback)(Pair);


/**
 * suffix tree containing suffix tree nodes and other data
 */
class SuffixArray
{
    public:
        struct VectorCallback {
            vector<Pair> &pairs;
            VectorCallback(vector<Pair> &pairs)
                : pairs(pairs) {}
            bool operator()(const Pair &p) {
                pairs.push_back(p);
                return false;
            }
        };

        struct SetCallback {
            set<Pair> &pairs;
            SetCallback(set<Pair> &pairs)
                : pairs(pairs) {}
            bool operator()(const Pair &p) {
                pairs.insert(p);
                return false;
            }
        };

        /**
         * Builds a suffix array for the given bucket (list of suffixes).
         *
         * @param[in] sequences all fasta sequences
         * @param[in] bucket linked list of suffixes for this bucket
         * @param[in] param alignment parameters
         * @param[in] k kmer size, if different than param.window_size
         * @return the suffix array
         */
        SuffixArray(SequenceDatabase *sequences,
                   Bucket *bucket,
                   const Parameters &param,
                   int k=-1);

        /** Destructor. */
        ~SuffixArray();

        /**
         * Generate promising pairs for alignment.
         *
         * @param[out] pairs
         */
        bool generate_pairs(set<Pair> &pairs) {
            return generate_pairs(SetCallback(pairs));
        }

        /**
         * Generate promising pairs for alignment.
         *
         * @param[out] pairs
         */
        bool generate_pairs(vector<Pair> &pairs) {
            return generate_pairs(VectorCallback(pairs));
        }

        /**
         * Generate promising pairs for alignment.
         *
         * @param[out] pairs
         */
        template <class Callback>
        bool generate_pairs(Callback callback);

        static bool length_filter(size_t s1Len, size_t s2Len, size_t cutOff);

        void print();

    private:
        bool is_candidate(Suffix *p, Suffix *q);

        template <class Callback>
        bool pair_check(
                int &count_generated,
                Callback callback,
                const int &i,
                const int &j);

        template <class Callback>
        bool process(
                int &count,
                int &count_generated,
                Callback callback,
                const quad &q,
                const int &cutoff);

        Sequence& get_sequence(size_t i) {
            if (sequences_cache.count(i)) {
                return *sequences_cache[i];
            }
            else {
                return *(sequences_cache[i] = sequences->get_sequence(i));
            }
        };

        SequenceDatabase *sequences;
        //Bucket *bucket;
        const Parameters &param;/**< user parameters */
        const size_t SIGMA;     /**< alphabet size */
        const char DOLLAR;      /**< terminal character */
        const char BEGIN;       /**< leftmost terminal character */
        const string alphabet;  /**< the alphabet */
        vector<size_t> alphabet_table;/** lookup table for alphabet index */
        int window_size;        /**< sliding window size */
        unsigned char *T;       /**< the concatenated suffix input */
        int *SA;                /**< the suffix array */
        int *LCP;               /**< the lcp array */
        unsigned char *BWT;     /**< the burrows wheeler array */
        int *SID;               /**< the sequence ID array */
        int *POS;               /**< the position array */
        vector<int> END;        /**< end index for each sequence */
        int n;                  /**< number of characters in the input */
        int sid;                /**< number of sequences in input */
        char sentinal;          /**< sentinal character */
        map<size_t,Sequence*> sequences_cache; /**< local sequence cache */

        static const size_t npos;
};


template <class Callback>
bool SuffixArray::generate_pairs(Callback callback)
{
    int count = 0;
    int count_generated = 0;
    int cutoff = param.exact_match_length;

    /* When counting k-mer buckets, the GSA we create will put all
     * sentinals either at the beginning or end of the SA. We don't want
     * to count all of the terminals, nor do we want to process them in
     * our bottom-up traversal. */
    /* do the sentinals appear at the beginning or end of SA? */
    fprintf(stderr, "n=%d sid=%d\n", n, sid);
    int bup_start = 1;
    int bup_stop = n;
    if (T[SA[0]] == sentinal) {
        fprintf(stderr, "sentinals at beginning\n");
        bup_start = sid+1;
        bup_stop = n;
    }
    else if (T[SA[n-1]] == sentinal) {
        fprintf(stderr, "sentinals at end\n");
        bup_start = 1;
        bup_stop = n-sid;
    }
    else {
        fprintf(stderr, "T[SA[0]]=%c T[SA[n-1]]=%c sentinal=%c\n",
                T[SA[0]], T[SA[n-1]], sentinal);
        fprintf(stderr, "sentinals not found at beginning or end of SA\n");
    }
    printf("bup_start=%d\n", bup_start);
    printf("bup_stop =%d\n", bup_stop);

    stack<quad> the_stack;
    quad last_interval;
    the_stack.push(quad());
    for (int i = bup_start; i <= bup_stop; ++i) {
        int lb = i - 1;
        while (LCP[i] < the_stack.top().lcp) {
            the_stack.top().rb = i - 1;
            last_interval = the_stack.top();
            the_stack.pop();
            process(count, count_generated, callback, last_interval, cutoff);
            lb = last_interval.lb;
            if (LCP[i] <= the_stack.top().lcp) {
                last_interval.children.clear();
                the_stack.top().children.push_back(last_interval);
                last_interval = quad();
            }
        }
        if (LCP[i] > the_stack.top().lcp) {
            if (!last_interval.empty()) {
                last_interval.children.clear();
                the_stack.push(quad(LCP[i],lb,INT_MAX,vector<quad>(1, last_interval)));
                last_interval = quad();
            }
            else {
                the_stack.push(quad(LCP[i],lb,INT_MAX));
            }
        }
    }
    the_stack.top().rb = bup_stop - 1;
    process(count, count_generated, callback, the_stack.top(), cutoff);

    return false;
}


template <class Callback>
bool SuffixArray::pair_check(
        int &count_generated,
        Callback callback,
        const int &i,
        const int &j)
{
    bool retval = false;
    const int &sidi = SID[SA[i]];
    const int &sidj = SID[SA[j]];
    if (BWT[i] != BWT[j] || BWT[i] == sentinal) {
        if (sidi != sidj) {
            ++count_generated;
            if (sidi < sidj) {
                retval = callback(make_pair(sidi,sidj));
            }
            else {
                retval = callback(make_pair(sidj,sidi));
            }
        }
    }
    return retval;
}


template <class Callback>
bool SuffixArray::process(
        int &count,
        int &count_generated,
        Callback callback,
        const quad &q,
        const int &cutoff)
{
    bool retval = false;
    const int n_children = q.children.size();
    int child_index = 0;

    ++count;

    if (q.lcp < cutoff) return false;

    if (n_children) {
        for (int i=q.lb; i<=q.rb; ++i) {
            int j = i+1;
            if (child_index < n_children) {
                if (i >= q.children[child_index].lb) {
                    j = q.children[child_index].rb+1;
                    if (i >= q.children[child_index].rb) {
                        ++child_index;
                    }
                }
            }
            for (/*nope*/; j<=q.rb; ++j) {
                retval = pair_check(count_generated, callback, i, j);
                if (retval) return true;
            }
        }
    }
    else {
        for (int i=q.lb; i<=q.rb; ++i) {
            for (int j=i+1; j<=q.rb; ++j) {
                retval = pair_check(count_generated, callback, i, j);
                if (retval) return true;
            }
        }
    }

    return false;
}

}; /* namespace pgraph */

#endif /* _PGRAPH_SUFFIXARRAY_H_ */
