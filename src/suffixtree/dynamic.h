/**
 * @file dynamic.h
 *
 * @author andy.cj.wu@gmail.com
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2010 Washington State University. All rights reserved.
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 *
 * @todo TODO remove SIGMA and NROW; dynamically assign alphabet
 */
#ifndef DYNAMIC_H_
#define DYNAMIC_H_

#ifdef __cplusplus
extern "C" {
#endif 

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <limits.h>

#define SIGMA 26U /**< size of alphabet */
#define NROW 2U /**< number of rows in dynamic programming table, always 2 */

enum {OPEN = -10, GAP = -1};
enum {NO = 0, YES = 1};
enum {FALSE = 0, TRUE = 1};

/**
 * cell_t for dynamic alignment result
 */
typedef struct {
    int score;  /**< alignment score */
    int ndig;   /**< #(matches) */
    int alen;   /**< alignment length */
} cell_t;

/**
 * parameters for is_edge test, packed in struct
 */
typedef struct {
    int AOL;            /**< AlignOverLongerSeq */
    int SIM;            /**< MatchSimilarity */
    int OS;             /**< OptimalScoreOverSelfScore */
    int exact_match_len;/**< exact match length cutoff */
    int window_size;    /**< slide window size */
} param_t;


/**
 * Initialize the protein index mapping used internally by the scoring matrix.
 *
 * @param[in] nAA the number of characters in the alphabet
 */
void init_map(int nAA);


/**
 * Calculates the score if the given sequence were aligned with itself.
 *
 * @param[in] s the sequence as a C string
 * @param[in] ns the length of s
 * @return the self score
 */
int self_score(const char *s, size_t ns);


/**
 * Allocates a 2-dimensional cell_t array for sequence alignment.
 *
 * @param[in] nrow number of rows (always 2)
 * @param[in] ncol length of longest sequence to align
 * @return the table
 */
cell_t **alloc_tbl(int nrow, int ncol);


/**
 * Allocates a 2-dimensional int array (for sequence alignment).
 *
 * @param nrow number of rows
 * @param ncol number of columns e.g. length of longest sequence to align
 * @return the table
 */
int **alloc_int(int nrow, int ncol);


/**
 * Frees the 2-dimensional cell_t table. Symmetric to alloc_tbl().
 *
 * @param[in] tbl the cell_t table
 * @param[in] the number of rows (always 2)
 */
void free_tbl(cell_t **tbl, int nrow);


/**
 * Frees the 2-dimensional int table. Symmetric to alloc_int().
 *
 * @param[in] tbl the int table
 * @param[in] the number of rows (always 2)
 */
void free_int(int **tbl, int nrow);


/**
 * Implementation of affine gap pairwise sequence alignment.
 *
 * It is a space efficient version: only two rows are required; also mem for
 * all dynamic tables are allocated ONLY ONCE outside of this function call
 * using alloc_tbl() and alloc_int() and passed as tbl, del, and ins arguments.
 *
 * @param[in] s1 sequence s1
 * @param[in] s1Len sequence length of <s1>, strlen(s1)
 * @param[in] s2 sequence s2
 * @param[in] s2Len sequence length of <s2>, strlen(s2)
 * @param[out] result alignment result <score, ndig, alen>
 * @param[in] tbl pre-allocated score table
 * @param[in] del pre-allocated deletion table
 * @param[in] ins pre-allocated insertion table
 */
void affine_gap_align(
        const char *s1, size_t s1Len,
        const char *s2, size_t s2Len,
        cell_t *result, cell_t **tbl, int **del, int **ins);


/**
 * Prints the given row to stdout (for debugging purposes).
 *
 * @param[in] tbl the cell_t table
 * @param[in] i the row to print
 * @param[in] the number of columns in the table
 */
void print_row(cell_t **tbl, int i, int ncol);


/**
 * Asks whether the given cell_t alignment result is an edge, based on the
 * given param parameters.
 *
 * @param[in] result the dynamic programming alignment result
 * @param[in] s1 sequence s1
 * @param[in] s1Len sequence length of <s1>, strlen(s1)
 * @param[in] s2 sequence s2
 * @param[in] s2Len sequence length of <s2>, strlen(s2)
 * @param[in] param the heuristic parameters for edge determination
 * @param[out] sscore self score
 * @param[out] maxLen longer of s1Len and s2Len
 * @return TRUE if this is an edge, FALSE otherwise
 */
int is_edge(const cell_t result,
        const char *s1, size_t s1Len,
        const char *s2, size_t s2Len,
        const param_t param, int *sscore, int *maxLen);

#ifdef __cplusplus
}
#endif 

#endif /* DYNAMIC_H_ */
