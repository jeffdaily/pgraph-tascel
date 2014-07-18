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

#include <algorithm>
#include <climits>
#include <cstddef>
#include <map>
#include <set>
#include <stack>
#include <string>
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
using ::std::string;
using ::std::vector;
using ::std::size_t;

#ifndef SIZE_MAX
#define SIZE_MAX (size_t(-1))
#endif

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

        int get_size() const { return n; }

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
        Bucket *bucket;         /**< the bucket */
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
        vector<int> SIDmap;     /**< end index for each sequence */
        int n;                  /**< number of characters in the input */
        int sid;                /**< number of sequences in input */
        char sentinal;          /**< sentinal character */
        map<size_t,Sequence*> sequences_cache; /**< local sequence cache */

        static const size_t npos;
};


/* the suffix array we create will represent more than just the subtree we are
 * after, so we must find the "bup_start/bup_stop" indexes corresponding to the
 * subtree we want to generate pairs from. use a binary search with this
 * comparator */
struct SuffixArrayComparator {
    Bucket *bucket;
    unsigned char *T;
    int *SA;
    int bup_start;

    SuffixArrayComparator(Bucket *bucket, unsigned char *T, int *SA, int bup_start)
        : bucket(bucket)
        , T(T)
        , SA(SA)
        , bup_start(bup_start)
    {}

    bool operator() (const int &sa_val, const string &val) {
#if DEBUG
        cout << &sa_val - &SA[bup_start] << endl;
        cout << string((const char *)&T[sa_val], bucket->k) << " < " << val;
        cout << " = " << (string((const char*)&T[sa_val], bucket->k) < val) << endl;
#endif
        return string((const char*)&T[sa_val], bucket->k) < val;
    }
    bool operator() (const string &val, const int &sa_val) {
#if DEBUG
        cout << &sa_val - &SA[bup_start] << endl;
        cout << val << " < " << string((const char *)&T[sa_val], bucket->k);
        cout << " = " << (val < string((const char*)&T[sa_val], bucket->k)) << endl;
#endif
        return val < string((const char*)&T[sa_val], bucket->k);
    }
};


static size_t powz(size_t base, size_t n)
{
    size_t p = 1;

    for(/*empty*/; n > 0; --n) {
        assert(p < SIZE_MAX/base);
        p *= base;
    }

    return p;
}


template <class Callback>
bool SuffixArray::generate_pairs(Callback callback)
{
    int count = 0;
    int count_generated = 0;
    int cutoff = param.exact_match_length;

    string kmer;
    {
        size_t bid = bucket->bid;
        for (int i = bucket->k-1; i >= 0; --i) {
            size_t tmp = powz(param.alphabet.size(),i);
            size_t quo = bid / tmp;
            size_t rem = bid % tmp;
            assert(quo <= param.alphabet.size());
            kmer += param.alphabet[quo];
            bid = rem;
        }
    }
#if DEBUG
    cout << "kmer=" << kmer << endl;
#endif

    /* When counting k-mer buckets, the GSA we create will put all
     * sentinals either at the beginning or end of the SA. We don't want
     * to count all of the terminals, nor do we want to process them in
     * our bottom-up traversal. */
    /* do the sentinals appear at the beginning or end of SA? */
#if DEBUG
    fprintf(stderr, "n=%d sid=%d\n", n, sid);
#endif
    int bup_start = 0;
    int bup_stop = n;
    if (T[SA[0]] == sentinal) {
#if DEBUG
        fprintf(stderr, "sentinals at beginning\n");
#endif
        while (T[SA[bup_start]] == sentinal) {
            ++bup_start;
        }
    }
    else if (T[SA[n-1]] == sentinal) {
#if DEBUG
        fprintf(stderr, "sentinals at end\n");
#endif
        while (T[SA[bup_stop-1]] == sentinal) {
            --bup_stop;
        }
    }
#if DEBUG
    else {
        fprintf(stderr, "T[SA[0]]=%c T[SA[n-1]]=%c sentinal=%c\n",
                T[SA[0]], T[SA[n-1]], sentinal);
        fprintf(stderr, "sentinals not found at beginning or end of SA\n");
    }
#endif
#if DEBUG
    printf("bup_start=%d SA[%d]=%d SA[%d]=%d SA[%d]=%d SA[%d]=%d SA[%d]=%d\n",
            bup_start,
            bup_start-2, SA[bup_start-2],
            bup_start-1, SA[bup_start-1],
            bup_start-0, SA[bup_start-0],
            bup_start+1, SA[bup_start+1],
            bup_start+2, SA[bup_start+2]);
    printf("T[%d]=%c T[%d]=%c T[%d]=%c T[%d]=%c T[%d]=%c\n",
            SA[bup_start-2], T[SA[bup_start-2]],
            SA[bup_start-1], T[SA[bup_start-1]],
            SA[bup_start-0], T[SA[bup_start-0]],
            SA[bup_start+1], T[SA[bup_start+1]],
            SA[bup_start+2], T[SA[bup_start+2]]);
    printf("bup_stop=%d SA[%d]=%d SA[%d]=%d SA[%d]=%d SA[%d]=%d SA[%d]=%d\n",
            bup_stop,
            bup_stop-2, SA[bup_stop-2],
            bup_stop-1, SA[bup_stop-1],
            bup_stop-0, SA[bup_stop-0],
            bup_stop+1, SA[bup_stop+1],
            bup_stop+2, SA[bup_stop+2]);
    printf("T[%d]=%c T[%d]=%c T[%d]=%c T[%d]=%c T[%d]=%c\n",
            SA[bup_stop-2], T[SA[bup_stop-2]],
            SA[bup_stop-1], T[SA[bup_stop-1]],
            SA[bup_stop-0], T[SA[bup_stop-0]],
            SA[bup_stop+1], T[SA[bup_stop+1]],
            SA[bup_stop+2], T[SA[bup_stop+2]]);
#endif

    int new_bup_start = bup_start + ::std::lower_bound(
            &SA[bup_start], &SA[bup_stop], kmer,
            SuffixArrayComparator(bucket, T, SA, bup_start)) - &SA[bup_start];
#if DEBUG
    printf("new_bup_start=%d SA[%d]=%d SA[%d]=%d SA[%d]=%d SA[%d]=%d SA[%d]=%d\n",
            new_bup_start,
            new_bup_start-2, SA[new_bup_start-2],
            new_bup_start-1, SA[new_bup_start-1],
            new_bup_start-0, SA[new_bup_start-0],
            new_bup_start+1, SA[new_bup_start+1],
            new_bup_start+2, SA[new_bup_start+2]);
    printf("T[%d]=%c T[%d]=%c T[%d]=%c T[%d]=%c T[%d]=%c\n",
            SA[new_bup_start-2], T[SA[new_bup_start-2]],
            SA[new_bup_start-1], T[SA[new_bup_start-1]],
            SA[new_bup_start-0], T[SA[new_bup_start-0]],
            SA[new_bup_start+1], T[SA[new_bup_start+1]],
            SA[new_bup_start+2], T[SA[new_bup_start+2]]);
#endif

    int new_bup_stop = bup_start + ::std::upper_bound(
            &SA[bup_start], &SA[bup_stop+1], kmer,
            SuffixArrayComparator(bucket, T, SA, bup_start)) - &SA[bup_start];
#if DEBUG
    printf("new_bup_stop=%d SA[%d]=%d SA[%d]=%d SA[%d]=%d SA[%d]=%d SA[%d]=%d\n",
            new_bup_stop,
            new_bup_stop-2, SA[new_bup_stop-2],
            new_bup_stop-1, SA[new_bup_stop-1],
            new_bup_stop-0, SA[new_bup_stop-0],
            new_bup_stop+1, SA[new_bup_stop+1],
            new_bup_stop+2, SA[new_bup_stop+2]);
    printf("T[%d]=%c T[%d]=%c T[%d]=%c T[%d]=%c T[%d]=%c\n",
            SA[new_bup_stop-2], T[SA[new_bup_stop-2]],
            SA[new_bup_stop-1], T[SA[new_bup_stop-1]],
            SA[new_bup_stop-0], T[SA[new_bup_stop-0]],
            SA[new_bup_stop+1], T[SA[new_bup_stop+1]],
            SA[new_bup_stop+2], T[SA[new_bup_stop+2]]);
#endif

    stack<quad> the_stack;
    quad last_interval;
    the_stack.push(quad());
    for (int i = new_bup_start; i <= new_bup_stop; ++i) {
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
            if (SIDmap[sidi] < SIDmap[sidj]) {
                retval = callback(make_pair(SIDmap[sidi],SIDmap[sidj]));
            }
            else {
                retval = callback(make_pair(SIDmap[sidj],SIDmap[sidi]));
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
