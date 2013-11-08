/**
 * @file stree.cpp
 *
 * @author andy.cj.wu@gmail.com
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2010 Washington State University. All rights reserved.
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 */
#include "config.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include <algorithm>
#include <map>
#include <set>
#include <utility>
#include <vector>

#include "bucket.hpp"
#include "constants.h"
#include "alignment.hpp"
#include "stree.hpp"

#define PRINT_EDGES 0

using std::make_pair;
using std::map;
using std::pair;
using std::set;

#define STATIC static inline

namespace pgraph {

STATIC void calc_stats(stats_t &stats)
{
    map<double,int> counts;
    int count = 0;

    stats.sum = 0.0;
    stats.min = 0.0;
    stats.max = 0.0;
    stats.mean = 0.0;
    stats.mode = 0.0;
    stats.median = 0.0;
    stats.stddev = 0.0;

    if (stats.values.empty()) {
        return;
    }

    // sort for median
    std::sort(stats.values.begin(), stats.values.end());

    // calc median
    if (stats.values.size() % 2U == 0) {
        // evenly sized
        stats.median =  stats.values[((stats.values.size())/2U)];
        stats.median += stats.values[((stats.values.size())/2U)-1U];
        stats.median =  stats.median / 2.0;
    }
    else {
        // oddly sized
        stats.median = stats.values[((stats.values.size()+1U)/2U)-1U];
    }

    // calc mean and mode together
    for (vector<double>::const_iterator it=stats.values.begin();
            it!=stats.values.end(); ++it) {
        stats.sum += *it;
        stats.min = *it < stats.min ? *it : stats.min;
        stats.max = *it > stats.max ? *it : stats.max;
        if (counts.count(*it) == 0) {
            counts[*it] = 1;
        }
        else {
            counts[*it] += 1;
        }
        if (counts[*it] > count) {
            count = counts[*it];
            stats.mode = *it;
        }
    }

    stats.mean = stats.sum / stats.values.size();

    // calc stddev
    for (vector<double>::const_iterator it=stats.values.begin();
            it!=stats.values.end(); ++it) {
        stats.stddev += std::pow(*it-stats.mean,2.0);
    }
    stats.stddev = std::pow(stats.stddev / stats.values.size(), 0.5);
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
STATIC void
compute_lset(suffix_t *suffixes, sequence_t *seqs, suffix_t **lset)
{
    suffix_t *p = NULL;
    suffix_t *q = NULL;
    int lIndex;

    for (p = suffixes; p != NULL; p = q) {
        q = p->next;
        if (p->pid == 0) {
            lIndex = BEGIN - 'A';
            p->next = lset[lIndex];
            lset[lIndex] = p;
        }
        else {
            lIndex = seqs[p->sid].str[p->pid - 1] - 'A';
            p->next = lset[lIndex];
            lset[lIndex] = p;
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
 * @param[in] seqs - all fasta seqs
 * @param[in] suffixes - bucket in linked list way
 * @param[in] depth - depth since root, which has 0 depth
 * @param[in] window_size - slide window size for bucketing
 * ---------------------------------------------------*/
STATIC int
nextDiffPos(sequence_t *seqs, suffix_t *suffixes, int depth, int window_size)
{
    int i;
    suffix_t *p = NULL;
    char pCh;
    char cCh;

    assert(suffixes != NULL);
    i = depth; /* position checked until this point */
    p = suffixes;

    while (1) {
        i++; /* step forward one more char for comparison */
        assert((i + p->pid) < seqs[p->sid].size);
        pCh = seqs[p->sid].str[p->pid + i];

        for (p = p->next; p != NULL; p = p->next) {
            assert((i + p->pid) < seqs[p->sid].size);
            cCh = seqs[p->sid].str[p->pid + i];
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
 * TODO
 */
STATIC size_t
build_tree_recursive(
    sequence_t *sequences, size_t n_sequences,
    suffix_t *suffixes, int depth, int window_size,
    stnode_t *st_nodes, size_t *st_index)
{
    int diffPos;
    size_t i;
    size_t j; /* store the index of the internal node */
    suffix_t **heads = NULL;
    suffix_t *p = NULL;
    suffix_t *q = NULL;
    int hIndex;
    sequence_t *seq = NULL;
    size_t rLeaf = SIZE_MAX;

    diffPos = nextDiffPos(sequences, suffixes, depth, window_size);
    if (diffPos == ERROR) { /* can never happen */
        fprintf(stderr, "wrong when exploring diff chars!");
        exit(EXIT_FAILURE);
    }
    else if (diffPos == DOL_END) { /* leaf node */
        seq = &sequences[suffixes[0].sid];

        /* depth also includes the '$' */
        st_nodes[*st_index].fanout = 0;
        st_nodes[*st_index].depth = (seq->size - 1) - suffixes[0].pid;
        st_nodes[*st_index].rLeaf = *st_index; /* point to itself */

        compute_lset(suffixes, sequences, st_nodes[*st_index].lset);

        return (*st_index)++;
    }
    else {  /* internal node */

        j = *st_index; /* store st_index in the stack */
        (*st_index)++;
        st_nodes[j].fanout = 0;
        st_nodes[j].depth = diffPos - 1 ;

        heads = st_nodes[j].lset;
        for (i = 0; i < SIGMA; i++) {
            heads[i] = NULL;
        }

#ifdef DEBUG
        printf("============================\n");
        print_suffixes(suffixes);
        printf("============================\n");
#endif

        /* partition suffixes into SIGMA sub-buckets */
        for (p = suffixes; p != NULL; p = q) {
            q = p->next;
            hIndex = sequences[p->sid].str[p->pid + diffPos] - 'A';
            assert(hIndex >= 0 && (unsigned)hIndex < SIGMA);

            p->next = heads[hIndex];
            heads[hIndex] = p;
        }

#ifdef DEBUG
        for (i = 0; i < SIGMA; i++) {
            if (heads[i]) {
                print_suffixes(heads[i]);
            }
        }
#endif
        for (i = 0; i < SIGMA; i++) {
            if (heads[i]) {
                st_nodes[j].fanout += 1;
            }
        }

        /* recursively construct the tree in DFS way */
        for (i = 0; i < SIGMA; i++) {
            if (heads[i]) {
                /* branching with '$' */
                if (i == (DOLLAR - 'A')) {
                    st_nodes[*st_index].fanout = 0;
                    st_nodes[*st_index].depth =
                        (sequences[heads[i]->sid].size - 1 - heads[i]->pid);
                    st_nodes[*st_index].rLeaf = *st_index;

                    compute_lset(heads[i], sequences, st_nodes[*st_index].lset);
                    rLeaf = (*st_index)++;
                }
                else {
                    rLeaf = build_tree_recursive(
                            sequences, n_sequences, heads[i], diffPos,
                            window_size, st_nodes, st_index);
                }

                /* put it back into NULL */
                heads[i] = NULL;
            }
        }

        /* store the right most leaf in the internal node */
        st_nodes[j].rLeaf = rLeaf;
        return st_nodes[j].rLeaf;
    }
}


stree_t* build_tree(
        sequences_t *sequences, bucket_t *bucket, param_t param)
{
    stree_t *tree = NULL;
    size_t n_nodes = 0;
    size_t i = 0;

    /* allocate tree */
    tree = new stree_t;
    if (NULL == tree) {
        perror("build_tree: malloc tree");
        exit(EXIT_FAILURE);
    }

    /* allocate tree nodes */
    n_nodes = 2 * bucket->size;
    tree->nodes = new stnode_t[n_nodes];
    if (NULL == tree->nodes) {
        perror("build_tree: malloc tree nodes");
        exit(EXIT_FAILURE);
    }

    /* allocate statistical stuff */
    tree->fanout.values.reserve(n_nodes);
    tree->depth.values.reserve(n_nodes);
    tree->sequence_length.values.reserve(bucket->size);
    tree->suffix_length.values.reserve(bucket->size);

    tree->size = 0;

    /* allocate lset pointer memory for tree nodes */
    tree->lset_array = new suffix_t*[SIGMA * n_nodes];
    if (NULL == tree->lset_array) {
        perror("build_tree: malloc lset array");
        exit(EXIT_FAILURE);
    }
    /* initialize lset pointer memory */
    for (i = 0; i < SIGMA * n_nodes; ++i) {
        tree->lset_array[i] = NULL;
    }

    /* initialize tree nodes */
    for (i = 0; i < n_nodes; ++i) {
        tree->nodes[i].depth = 0;
        tree->nodes[i].rLeaf = 0;
        tree->nodes[i].lset = &(tree->lset_array[i*SIGMA]);
    }

    /* get suffix and sequence statistics */
    {
        suffix_t *p = NULL;
        suffix_t *q = NULL;
        for (p = bucket->suffixes; p != NULL; p = q) {
            size_t len = sequences->seq[p->sid].size;
            tree->sequence_length.values.push_back(len);
            tree->suffix_length.values.push_back(len - p->pid);
            q = p->next;
        }
    }

    (void) build_tree_recursive(
            sequences->seq, sequences->size,
            bucket->suffixes, param.window_size - 1, param.window_size,
            tree->nodes, &(tree->size));

    /* get tree node statistics */
    for (i = 0; i < tree->size; ++i) {
        tree->fanout.values.push_back(tree->nodes[i].fanout);
        tree->depth.values.push_back(tree->nodes[i].fanout);
    }

    calc_stats(tree->fanout);
    calc_stats(tree->depth);
    calc_stats(tree->sequence_length);
    calc_stats(tree->suffix_length);

    return tree;
}


void free_tree(stree_t *tree)
{
    assert(NULL != tree);
    delete [] tree->nodes;
    delete [] tree->lset_array;
    delete tree;
}


STATIC int
is_candidate(sequence_t *seqs, size_t nSeqs,
             suffix_t *p, suffix_t *q, param_t param, int *dup)
{
    size_t s1Len = 0;
    size_t s2Len = 0;
    int f1 = 0;
    int f2 = 0;
    int cutOff = param.AOL * param.SIM;
    cell_t result;
    size_t index = 0;

    f1 = p->sid;
    f2 = q->sid;
    assert(size_t(f1) < nSeqs);
    assert(size_t(f2) < nSeqs);

    if (f1 > f2) {
        int swap = f1;
        f1 = f2;
        f2 = swap;
    }
    index = (nSeqs*f1) + f2 - (f1*(f1+1)/2);

    if (f1 == f2) {
        dup[index] = FALSE;
        return FALSE;
    }

    s1Len = seqs[f1].size - 1;
    s2Len = seqs[f2].size - 1;

    if (dup[index] == MAYBE) {
        if (s1Len <= s2Len) {
            if (100 * s1Len < cutOff * s2Len) {
                dup[index] = FALSE;
            }
            else {
                dup[index] = TRUE;
            }
        }
        else {
            if (100 * s2Len < cutOff * s1Len) {
                dup[index] = FALSE;
            }
            else {
                dup[index] = TRUE;
            }
        }
        /* check if it is an edge */
        if (dup[index] == MAYBE) {
            int ignore1;
            size_t ignore2;
            assert(0);
            dup[index] = is_edge_blosum(result,
                    seqs[f1].str, seqs[f1].size,
                    seqs[f2].str, seqs[f2].size,
                    param.AOL, param.SIM, param.OS,
                    ignore1, ignore2);
        }
        
        return dup[index];
    }
    else {
        return FALSE;
    }
}


STATIC int
is_candidate(sequence_t *seqs, size_t nSeqs,
             suffix_t *p, suffix_t *q, param_t param, 
             set<pair<unsigned int,unsigned int> > &dup)
{
    size_t s1Len = 0;
    size_t s2Len = 0;
    int f1 = 0;
    int f2 = 0;
    int cutOff = param.AOL * param.SIM;
    cell_t result;
    size_t index = 0;

    f1 = p->sid;
    f2 = q->sid;
    assert(size_t(f1) < nSeqs);
    assert(size_t(f2) < nSeqs);

    if (f1 > f2) {
        int swap = f1;
        f1 = f2;
        f2 = swap;
    }

    if (f1 == f2) {
        return FALSE;
    }

    s1Len = seqs[f1].size - 1;
    s2Len = seqs[f2].size - 1;

    if (s1Len <= s2Len) {
        if (100 * s1Len < cutOff * s2Len) {
            return FALSE;
        }
        else {
            return TRUE;
        }
    }
    else {
        if (100 * s2Len < cutOff * s1Len) {
            return FALSE;
        }
        else {
            return TRUE;
        }
    }

    assert(0);
    return FALSE;
}


/**
 * counting sort based on the depth of the tree nodes.
 *
 * @param stNodes -
 * @param srtIndex -
 * @param nStNodes -
 * @param maxSeqLen -
 * -----------------------------------------------------------*/
STATIC void
count_sort(stnode_t *stNodes, int *srtIndex, size_t nStNodes, size_t maxSeqLen)
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
        assert(depth < int(maxSeqLen));
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


STATIC void
procLeaf(suffix_t **lset, sequence_t *seqs, int nSeqs, param_t param, int *dup)
{
    size_t i;
    size_t j;
    suffix_t *p = NULL;
    suffix_t *q = NULL;
    int cutOff;

    cutOff = param.AOL * param.SIM;

    for (i = 0; i < SIGMA; i++) {
        if (lset[i]) {
            if (i == BEGIN - 'A') { /* inter cross */
                for (p = lset[i]; p != NULL; p = p->next) {
                    for (q = p->next; q != NULL; q = q->next) {
                        if (TRUE == is_candidate(seqs, nSeqs, p, q, param, dup)) {
#if PRINT_EDGES
                            if (p->sid > q->sid) {
                                printf("edge\t%zu\t%zu\n", q->sid, p->sid);
                            }
                            else {
                                printf("edge\t%zu\t%zu\n", p->sid, q->sid);
                            }
#endif
                        }
                    }
                }
            }

            /* intra cross */
            for (j = i + 1; j < SIGMA; j++) {
                if (lset[j]) {
                    for (p = lset[i]; p != NULL; p = p->next) {
                        for (q = lset[j]; q != NULL; q = q->next) {
                            if (TRUE == is_candidate(seqs, nSeqs, p, q, param, dup)) {
#if PRINT_EDGES
                                if (p->sid > q->sid) {
                                    printf("edge\t%zu\t%zu\n", q->sid, p->sid);
                                }
                                else {
                                    printf("edge\t%zu\t%zu\n", p->sid, q->sid);
                                }
#endif
                            }
                        }
                    }
                }
            }
        }
    }
}


STATIC void
procLeaf(suffix_t **lset, sequence_t *seqs, int nSeqs, param_t param, set<pair<unsigned int,unsigned int> > &dup)
{
    size_t i;
    size_t j;
    suffix_t *p = NULL;
    suffix_t *q = NULL;
    int cutOff;

    cutOff = param.AOL * param.SIM;

    for (i = 0; i < SIGMA; i++) {
        if (lset[i]) {
            if (i == BEGIN - 'A') { /* inter cross */
                for (p = lset[i]; p != NULL; p = p->next) {
                    for (q = p->next; q != NULL; q = q->next) {
                        if (TRUE == is_candidate(seqs, nSeqs, p, q, param, dup)) {
                            if (p->sid > q->sid) {
                                dup.insert(make_pair(q->sid,p->sid));
                            }
                            else {
                                dup.insert(make_pair(p->sid,q->sid));
                            }
                        }
                    }
                }
            }

            /* intra cross */
            for (j = i + 1; j < SIGMA; j++) {
                if (lset[j]) {
                    for (p = lset[i]; p != NULL; p = p->next) {
                        for (q = lset[j]; q != NULL; q = q->next) {
                            if (TRUE == is_candidate(seqs, nSeqs, p, q, param, dup)) {
                                if (p->sid > q->sid) {
                                    dup.insert(make_pair(q->sid,p->sid));
                                }
                                else {
                                    dup.insert(make_pair(p->sid,q->sid));
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}



STATIC void
procNode(stnode_t *stNodes, int sIndex, int eIndex, suffix_t **lset, sequence_t *seqs, int nSeqs, param_t param, set<pair<unsigned int,unsigned int> > &dup)
{
    size_t m;
    size_t n;
    size_t s;
    size_t t;
    suffix_t *p = NULL;
    suffix_t *q = NULL;

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
                                                    seqs, nSeqs, p, q,
                                                    param, dup)) {
                                            if (p->sid > q->sid) {
                                                dup.insert(make_pair(q->sid,p->sid));
                                            }
                                            else {
                                                dup.insert(make_pair(p->sid,q->sid));
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
    }
}


STATIC void
mergeLsets(stnode_t *stNodes, int sIndex, int eIndex)
{
    size_t m = 0;
    int j = 0;
    suffix_t *p = NULL;
    suffix_t *q = NULL;

    /* merge the lsets of subtree */
    for (m = 0; m < SIGMA; m++) {
        p = NULL;
        for (j = sIndex + 1; j <= eIndex; j++) {
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


void generate_pairs(
        stree_t *tree, sequences_t *sequences, int *dup, param_t param)
{
    stnode_t *stNodes = NULL;
    int *srtIndex = NULL;
    size_t nStNodes = 0;
    sequence_t *seqs = NULL;
    int nSeqs = 0;
    int maxSeqLen = 0;
    size_t i = 0;
    int j = 0;
    stnode_t *stnode = NULL;
    int sIndex;
    int eIndex;
    size_t m;
    size_t n;
    size_t s;
    size_t t;
    suffix_t *p = NULL;
    suffix_t *q = NULL;
    int EM;
    int cutOff; /* cut off value of filter 1 */

    srtIndex = new int[tree->size];
    count_sort(tree->nodes, srtIndex, tree->size, sequences->max_seq_size);
    stNodes = tree->nodes;
    nStNodes = tree->size;
    seqs = sequences->seq;
    nSeqs = sequences->size;
    maxSeqLen = sequences->max_seq_size;

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

        if (stnode->depth >= EM - 1) {
            if (stnode->rLeaf == sIndex) { /* leaf node */
                procLeaf(stnode->lset, seqs, nSeqs, param, dup);
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
                                                                seqs, nSeqs, p, q,
                                                                param, dup)) {
#if PRINT_EDGES
                                                        if (p->sid > q->sid) {
                                                            printf("edge\t%zu\t%zu\n",
                                                                    q->sid, p->sid);
                                                        }
                                                        else {
                                                            printf("edge\t%zu\t%zu\n",
                                                                    p->sid, q->sid);
                                                        }
#endif
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
                    for (j = sIndex + 1; j <= eIndex; j++) {
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
}


void generate_pairs(
        stree_t *tree, sequences_t *sequences, set<pair<unsigned int,unsigned int> > &dup, param_t param)
{
    stnode_t *stNodes = NULL;
    int *srtIndex = NULL;
    size_t nStNodes = 0;
    sequence_t *seqs = NULL;
    int nSeqs = 0;
    int maxSeqLen = 0;
    size_t i = 0;
    int j = 0;
    stnode_t *stnode = NULL;
    int sIndex;
    int eIndex;
    size_t m;
    size_t n;
    size_t s;
    size_t t;
    suffix_t *p = NULL;
    suffix_t *q = NULL;
    int EM;
    int cutOff; /* cut off value of filter 1 */

    srtIndex = new int[tree->size];
    count_sort(tree->nodes, srtIndex, tree->size, sequences->max_seq_size);
    stNodes = tree->nodes;
    nStNodes = tree->size;
    seqs = sequences->seq;
    nSeqs = sequences->size;
    maxSeqLen = sequences->max_seq_size;

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

        if (stnode->depth >= EM - 1) {
            if (stnode->rLeaf == sIndex) { /* leaf node */
                procLeaf(stnode->lset, seqs, nSeqs, param, dup);
            }
            else {                       /* internal node */
                eIndex = stnode->rLeaf;
                /* pairs generation loop for internal node */
                procNode(stNodes, sIndex, eIndex, stnode->lset, seqs, nSeqs, param, dup);
                /* merge the lsets of subtree */
                mergeLsets(stNodes, sIndex, eIndex);
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
}



}; /* namespace pgraph */

