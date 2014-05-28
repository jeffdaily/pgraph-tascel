/**
 * @file SuffixTree.cpp
 *
 * @author andy.cj.wu@gmail.com
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2010 Washington State University. All rights reserved.
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 */
#include "config.h"

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <set>
#include <utility>
#include <vector>

using ::std::cout;
using ::std::endl;
using ::std::make_pair;
using ::std::numeric_limits;
using ::std::pair;
using ::std::set;
using ::std::size_t;
using ::std::vector;

#include "Parameters.hpp"
#include "Sequence.hpp"
#include "SequenceDatabase.hpp"
#include "SuffixBuckets.hpp"
#include "SuffixTree.hpp"

#define ERROR -1
#define DOL_END -2
#ifndef SIZE_MAX
#define SIZE_MAX (size_t(-1))
#endif

namespace pgraph {

const size_t SuffixTree::npos(-1);

SuffixTree::SuffixTree(
        SequenceDatabase *sequences, Bucket *bucket, const Parameters &param, int k)
    :   sequences(sequences)
    //,   bucket(bucket)
    ,   param(param)
    ,   SIGMA(param.alphabet.size())
    ,   DOLLAR(param.alphabet_dollar)
    ,   BEGIN(param.alphabet_begin)
    ,   alphabet(param.alphabet)
    ,   alphabet_table(numeric_limits<unsigned char>::max(), npos)
    ,   window_size(k > 0 ? k : param.window_size)
    ,   nodes(NULL)
    ,   size(0U)
    ,   size_internal(0U)
    ,   lset_array(NULL)
    ,   fanout_stats()
    ,   depth_stats()
    ,   suffix_length_stats()
    ,   tails(NULL)
    ,   sequences_cache()
{
    size_t n_nodes = 2 * bucket->size;
    size_t i = 0;

#if DEBUG
    cout << "SuffixTree"
        << " bid=" << bucket->bid
        << " size=" << bucket->size
        << " bucket->k=" << bucket->k
        << " k=" << k
        << " param.window_size=" << param.window_size
        << " window_size=" << window_size
        << endl;
#endif

    for (i=0; i<SIGMA; ++i) {
        alphabet_table[(unsigned char)(param.alphabet[i])] = i;
    }

    /* allocate lset pointer memory for tree nodes */
    lset_array = new Suffix*[SIGMA * (n_nodes+1)];
    if (NULL == lset_array) {
        perror("build_tree: malloc lset array");
        exit(EXIT_FAILURE);
    }
    /* initialize lset pointer memory */
    for (i = 0; i < SIGMA * n_nodes; ++i) {
        lset_array[i] = NULL;
    }

    /* allocate tree nodes */
    nodes = new SuffixTreeNode[n_nodes];
    if (NULL == nodes) {
        perror("build_tree: malloc tree nodes");
        exit(EXIT_FAILURE);
    }

    /* initialize tree nodes */
    for (i = 0; i < n_nodes; ++i) {
        nodes[i].depth = 0;
        nodes[i].rLeaf = 0;
        nodes[i].lset = &(lset_array[i*SIGMA]);
    }
    tails = &(lset_array[SIGMA * n_nodes]);

    /* gather suffix statistics */
    Suffix *p = NULL;
    Suffix *q = NULL;
    for (p = bucket->suffixes; p != NULL; p = q) {
        q = p->next;
        Sequence &seq = get_sequence(p->sid);
        /* depth also includes the '$' */
        double len = (seq.size() - 1) - p->pid;
        suffix_length_stats.push_back(len);
    }

    /* the following bulk fetch of sequences caused significant slow
     * down */
#if 0
    /* make sure all required Sequence instances are cached */
    {
        set<size_t> sids;
        Suffix *p = NULL;
        for (p = bucket->suffixes; p != NULL; p = p->next) {
            sids.insert(p->sid);
        }
        sequences_cache = sequences->get_sequences(sids);
    }
#endif

    build_tree_recursive(bucket->suffixes, window_size - 1);
}


SuffixTree::~SuffixTree()
{
    delete [] nodes;
    delete [] lset_array;
    for (map<size_t,Sequence*>::iterator it=sequences_cache.begin();
            it!=sequences_cache.end(); ++it) {
        delete it->second;
    }
}



/**
 * TODO
 */
size_t SuffixTree::build_tree_recursive(Suffix *suffixes, int depth)
{
    int diffPos;

    diffPos = next_diff_pos(suffixes, depth);
    if (diffPos == ERROR) { /* can never happen */
        fprintf(stderr, "wrong when exploring diff chars!");
        exit(EXIT_FAILURE);
    }
    else if (diffPos == DOL_END) { /* leaf node */
        Sequence &seq = get_sequence(suffixes->sid);

        /* depth also includes the '$' */
        nodes[size].depth = (seq.size() - 1) - suffixes->pid;
        nodes[size].rLeaf = size; /* point to itself */

        compute_lset(suffixes, nodes[size].lset);

        return size++;
    }
    else {  /* internal node */
        size_t i;
        size_t j; /* store the index of the internal node */
        Suffix **heads = NULL;
        Suffix *p = NULL;
        Suffix *q = NULL;
        size_t rLeaf = SIZE_MAX;

        j = size; /* store st_index in the stack */
        size++;
        nodes[j].depth = diffPos - 1 ;

        depth_stats.push_back(nodes[j].depth);

        heads = nodes[j].lset;
        for (i = 0; i < SIGMA; i++) {
            heads[i] = NULL;
        }

        /* partition suffixes into SIGMA sub-buckets */
        for (p = suffixes; p != NULL; p = q) {
            char hChar = get_sequence(p->sid)[p->pid + diffPos];
            size_t hIndex = alphabet_table[(unsigned)hChar];
            if (hIndex == npos) {
                fprintf(stderr, "invalid h index caused by '%c'\n", hChar);
            }
            assert(hIndex != npos);
            assert(hIndex < SIGMA);

            q = p->next;
            p->next = heads[hIndex];
            heads[hIndex] = p;
        }

        /* compute fanout */
        size_internal++;
        for (i = 0; i < SIGMA; i++) {
            size_t fanout = 0;
            if (heads[i]) {
                ++fanout;
            }
            fanout_stats.push_back(fanout);
        }

        /* recursively construct the tree in DFS way */
        for (i = 0; i < SIGMA; i++) {
            if (heads[i]) {
                /* branching with '$' */
                if (i == (alphabet_table[DOLLAR])) {
                    nodes[size].depth =
                        (get_sequence(heads[i]->sid).size() - 1 - heads[i]->pid);
                    nodes[size].rLeaf = size;

                    compute_lset(heads[i], nodes[size].lset);
                    rLeaf = size++;
                }
                else {
                    rLeaf = build_tree_recursive(heads[i], diffPos);
                }

                /* put it back into NULL */
                heads[i] = NULL;
            }
        }

        /* store the right most leaf in the internal node */
        nodes[j].rLeaf = rLeaf;
        return nodes[j].rLeaf;
    }
}


/**
 * TODO
 */
size_t SuffixTree::build_tree_recursive2(Suffix *suffixes, int depth)
{
    int diffPos;

    diffPos = next_diff_pos2(suffixes, depth, size);
    if (diffPos == ERROR) { /* can never happen */
        fprintf(stderr, "wrong when exploring diff chars!");
        exit(EXIT_FAILURE);
    }
    else if (diffPos == DOL_END) { /* leaf node */
        suffixes = nodes[size].lset[0];
        nodes[size].lset[0] = NULL;
        Sequence &seq = get_sequence(suffixes->sid);

        /* depth also includes the '$' */
        nodes[size].depth = (seq.size() - 1) - suffixes->pid;
        nodes[size].rLeaf = size; /* point to itself */

        compute_lset(suffixes, nodes[size].lset);

        return size++;
    }
    else {  /* internal node */
        size_t i;
        size_t j; /* store the index of the internal node */
        Suffix **heads = NULL;
        Suffix *p = NULL;
        Suffix *q = NULL;
        size_t rLeaf = SIZE_MAX;

        j = size; /* store st_index in the stack */
        size++;
        nodes[j].depth = diffPos - 1 ;

        depth_stats.push_back(nodes[j].depth);

        heads = nodes[j].lset;

        /* compute fanout */
        size_internal++;
        for (i = 0; i < SIGMA; i++) {
            size_t fanout = 0;
            if (heads[i]) {
                ++fanout;
            }
            fanout_stats.push_back(fanout);
        }

        /* recursively construct the tree in DFS way */
        for (i = 0; i < SIGMA; i++) {
            if (heads[i]) {
                /* branching with '$' */
                if (i == (alphabet_table[DOLLAR])) {
                    nodes[size].depth =
                        (get_sequence(heads[i]->sid).size() - 1 - heads[i]->pid);
                    nodes[size].rLeaf = size;

                    compute_lset(heads[i], nodes[size].lset);
                    rLeaf = size++;
                }
                else {
                    rLeaf = build_tree_recursive2(heads[i], diffPos);
                }

                /* put it back into NULL */
                heads[i] = NULL;
            }
        }

        /* store the right most leaf in the internal node */
        nodes[j].rLeaf = rLeaf;
        return nodes[j].rLeaf;
    }
}


/**
 * compute lset according to the left character of pid to
 * ensure the left maximality. Mark the special case of
 * whole sequence. Use BEGIN='O' as the special character.
 *
 * @param suffixes -
 * @param seqs - all seqs info
 * @param lset -
 */
void SuffixTree::compute_lset(Suffix *suffixes, Suffix **lset)
{
    Suffix *p = NULL;
    Suffix *q = NULL;
    int lIndex;

    for (p = suffixes; p != NULL; p = q) {
        q = p->next;
        if (p->pid == 0) {
            lIndex = alphabet_table[BEGIN];
        }
        else {
            lIndex = alphabet_table[get_sequence(p->sid)[p->pid-1]];
        }
        p->next = lset[lIndex];
        lset[lIndex] = p;
    }
}


/**
 * find the next different position of suffixes.
 * NOTE: if all suffixes end at '$', then report it as
 *       isLeaf = YES
 *
 * return -2: all dollar ending
 *        -1: error, cannot happen
 *         i: next diff position from starting point
 *
 * @param[in] suffixes - bucket in linked list way
 * @param[in] depth - depth since root, which has 0 depth
 * ---------------------------------------------------*/
int SuffixTree::next_diff_pos(Suffix *suffixes, int depth)
{
    int i;
    Suffix *p = NULL;
    char pCh;
    char cCh;

    assert(suffixes != NULL);
    i = depth; /* position checked until this point */
    p = suffixes;

    while (1) {
        i++; /* step forward one more char for comparison */
        if ((i + p->pid) >= get_sequence(p->sid).size()) {
            cout << "i=" << i
                << " p->pid=" << p->pid
                << " i+p->pid=" << i+p->pid
                << " <? p->sid.size=" << get_sequence(p->sid).size()
                << endl;
            cout << get_sequence(p->sid);
            cout << string(p->pid, 'X') << endl;;
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
        assert((i + p->pid) < get_sequence(p->sid).size());
        pCh = get_sequence(p->sid)[p->pid + i];

        for (p = p->next; p != NULL; p = p->next) {
            assert((i + p->pid) < get_sequence(p->sid).size());
            cCh = get_sequence(p->sid)[p->pid + i];
            if (cCh != pCh) {
                return i;
            }
        }

        /* reset suffixes to next round of character comparison */
        p = suffixes;

        /* all suffixes ending in '$' */
        if (pCh == DOLLAR) {
            return DOL_END;
        }
    }
}


/**
 * find the next different position of suffixes.
 * NOTE: if all suffixes end at '$', then report it as
 *       isLeaf = YES
 *
 * return -2: all dollar ending
 *        -1: error, cannot happen
 *         i: next diff position from starting point
 *
 * @param[in] suffixes - bucket in linked list way
 * @param[in] depth - depth since root, which has 0 depth
 * ---------------------------------------------------*/
int SuffixTree::next_diff_pos2(Suffix *suffixes, int depth, int size)
{
    int i;
    size_t j;
    char pCh;
    char cCh;
    size_t hIndex;
    Suffix **heads = nodes[size].lset; /* temp space for heads */
    //Suffix **tails = new Suffix*[SIGMA]; /* temp space for tails */
    Suffix *p = NULL;
    Suffix *q = NULL;
    int result = -1;

    assert(suffixes != NULL);
    i = depth; /* position checked until this point */

    for (j = 0; j < SIGMA; j++) {
        heads[j] = NULL;
        tails[j] = NULL;
    }

    /* partition suffixes into SIGMA sub-buckets */
    p = suffixes;
    while (result == -1) {
        i++; /* step forward one more char for comparison */
        pCh = get_sequence(p->sid)[p->pid + i];
        hIndex = alphabet_table[(unsigned)pCh];
        if (hIndex == npos) {
            fprintf(stderr, "invalid h index caused by '%c'\n", pCh);
        }
        assert(hIndex != npos);
        assert(hIndex < SIGMA);
        q = p->next;
        p->next = heads[hIndex];
        heads[hIndex] = p;
        tails[hIndex] = p;

        for (p = q; p != NULL; p = q) {
            cCh = get_sequence(p->sid)[p->pid + i];
            hIndex = alphabet_table[(unsigned)cCh];
            if (hIndex == npos) {
                fprintf(stderr, "invalid h index caused by '%c'\n", cCh);
            }
            assert(hIndex != npos);
            assert(hIndex < SIGMA);
            q = p->next;
            p->next = heads[hIndex];
            heads[hIndex] = p;
            if (tails[hIndex] == NULL) {
                tails[hIndex] = p;
            }
            if (cCh != pCh) {
                result = i;
            }
        }

        if (result == -1) {
            /* went through all suffixes, all matched */
            /* stitch them back together */
            p = NULL;
            for (j = 0; j < SIGMA; ++j) {
                if (heads[j]) {
                    tails[j]->next = p;
                    p = heads[j];
                }
                heads[j] = NULL;
                tails[j] = NULL;
            }
            if (pCh == DOLLAR) {
                result = DOL_END;
                heads[0] = p;
            }
        }
    }

    //delete [] tails;

    return result;
}


bool SuffixTree::is_candidate(Suffix *p, Suffix *q)
{
    size_t s1Len = 0;
    size_t s2Len = 0;
    size_t f1 = 0;
    size_t f2 = 0;
    int cutOff = param.AOL * param.SIM / 100;
    bool result = false;

    f1 = p->sid;
    f2 = q->sid;
    assert(f1 < sequences->size());
    assert(f2 < sequences->size());

    if (f1 == f2) {
        result = false;
    }
    else {
        if (f1 > f2) {
            size_t swap = f1;
            f1 = f2;
            f2 = swap;
        }
        s1Len = get_sequence(f1).size() - 1;
        s2Len = get_sequence(f2).size() - 1;
        result = length_filter(s1Len, s2Len, cutOff);
    }

    return result;
}


/**
 * Returns true if the strings s1 and s2 should be considered for
 * alignment.
 */
bool SuffixTree::length_filter(size_t s1Len, size_t s2Len, size_t cutOff)
{
    bool result = true;

    if (s1Len <= s2Len) {
        if (100 * s1Len < cutOff * s2Len) {
            result = false;
        }
    }
    else {
        if (100 * s2Len < cutOff * s1Len) {
            result = false;
        }
    }

    return result;
}


/**
 * counting sort based on the depth of the tree nodes.
 * @return a sorted index of the SuffixTree nodes
 */
size_t* SuffixTree::count_sort()
{
    size_t i;
    size_t *cnt = NULL;
    size_t depth;
    size_t maxSeqLen = sequences->longest();
    size_t *srtIndex = new size_t[size];

    cnt = new size_t[maxSeqLen];
    if (NULL == cnt) {
        perror("count_sort: malloc cnt");
        exit(EXIT_FAILURE);
    }

    /* init counters */
    for (i = 0; i < maxSeqLen; i++) {
        cnt[i] = 0;
    }

    /* fill in counters */
    for (i = 0; i < size; i++) {
        depth = nodes[i].depth; /* depth starts from 0 */
        assert(((unsigned)depth) < maxSeqLen);
        cnt[depth]++;
    }

    /* suffix sum to sort in a reverse way */
    for (i = maxSeqLen - 2; i > 0; i--) {
        cnt[i] += cnt[i + 1];
    }
    cnt[0] += cnt[1]; /* because var i was unsigned */

    /* store the sorted index into srtIndex */
    for (i = 0; i < size; i++) {
        depth = nodes[i].depth;
        srtIndex[cnt[depth] - 1] = i;
        assert(cnt[depth] >= 1);
        cnt[depth]--;
    }

    delete [] cnt;

    return srtIndex;
}


void SuffixTree::merge_lsets(size_t sIndex, size_t eIndex)
{
    size_t m = 0;
    size_t j = 0;
    Suffix *p = NULL;
    Suffix *q = NULL;

    /* merge the lsets of subtree */
    for (m = 0; m < SIGMA; m++) {
        p = NULL;
        for (j = sIndex + 1; j <= eIndex; j = nodes[j].rLeaf + 1) {
            if ((q = nodes[j].lset[m])) {

                /* empty the subtree's ptrs array */
                nodes[j].lset[m] = NULL;
                if (p == NULL) {
                    p = q;
                    nodes[sIndex].lset[m] = q;
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

}; /* namespace pgraph */

