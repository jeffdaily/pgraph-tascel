/**
 * @file stree.h
 *
 * @author andy.cj.wu@gmail.com
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2010 Washington State University. All rights reserved.
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 */
#ifndef STREE_H_
#define STREE_H_

#include "bucket.h"
#include "dynamic.h"
#include "loadseq.h"

enum {ERROR = -1, DOL_END = -2};

/**
 * suffix tree node
 */
typedef struct {
    int depth; /**< depth since the root, not including the initial size k */
    int rLeaf; /**< right most leaf index */
    suffix_t *lset[SIGMA];  /**< subtree's nodes branched according to left
                              characters */
} stnode_t;


/**
 * Builds a tree for the given bucket (list of suffixes).
 *
 * Each bucket will be branched into SIGMA different sub-buckets in a DFS way,
 * the total space will be (26*4*depth) bytes
 *
 * NOTE: 'suffixes' parameter is modified after this function, which also means
 * the previous bucket does not exist at all
 *
 * @param[in] sequences all fasta sequences
 * @param[in] n_sequences number of fasta sequences
 * @param[in,out] suffixes suffixes list, modified thereafter
 * @param[in] depth characters compared since root, which has 0 depth
 * @param[in] window_size - slide window size
 * @param[in] st_nodes stnode_t array in DFS order
 * @param[in] st_index @todo TODO
 * @return @todo TODO
 */
int build_tree(sequence_t *sequences, size_t n_sequences, suffix_t *suffixes,
               int depth, int window_size, stnode_t *st_nodes, int *st_index);



/**
 * compute lset according to the left character of pid to
 * ensure the left maximality. Mark the special case of
 * whole sequence. Use BEGIN='O' as the special character.
 *
 * @param suffixes -
 * @param seqs - all seqs info
 * @param lset -
 */
void compute_lset(suffix_t *suffixes, sequence_t *sequences, suffix_t **lset);


/**
 * Constructs the tree, sorts it, and generates pairs.
 *
 * @param[in] sequences -
 * @param[in] n_sequences -
 * @param[in] suffixes -
 * @param[in] n_suffixes -
 * @param[in] max_seq_len -
 * @param[in] union_set -
 * @param[in] param -
 * @param[in] ind -
 */
void
process_bucket(sequence_t *sequences, size_t n_sequences,
               suffix_t *suffixes, size_t n_suffixes,
               size_t max_seq_len, param_t *param, int *dup);


/**
 * Builds a suffix tree for each bucket.
 *
 * @param buckets bucket
 * @param n_buckets -
 * @param seqs global fasta seqs
 * @param nseqs  #seqs in fasta file
 * @param window_size slide window size
 * @param NN #chars in all seqs
 * @param nStNodes #stnodes in constructed forest
 */
void build_forest(bucket_t *buckets, size_t n_buckets, sequence_t *sequences,
                  size_t n_sequences, size_t max_seq_len, param_t *param);


int nextDiffPos(sequence_t *seqs, suffix_t *lset, int depth, int k);


int lsetSize(suffix_t *lset);
void freeLset(suffix_t *lset);
void freeStNodes(stnode_t *stNodes, int size);


void print_stnodes(stnode_t *stNodes, int stIndex, int blSize, int ind);
void printLset(suffix_t **lset);
void printCnt(int *cnt, int size);
void printStNode(stnode_t *stNodes, int i);
void printStNodes(stnode_t *stNodes, int offset);

#endif /* end of stree.h */
