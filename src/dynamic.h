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

#include <stddef.h>

#include "param.h"
#include "csequence.h"

/**
 * cell_t for dynamic alignment result
 */
typedef struct {
    int score;  /**< alignment score */
    int ndig;   /**< #(matches) */
    int alen;   /**< alignment length */
} cell_t;


/** callback function type for generic dynamic programming (mis)match score */
typedef int (*match_t)(char one, char two);


/**
 * (Mis)match score based on blosum.
 *
 * Assumes inputs are valid characters from the protein alphabet.
 *
 * @see pg_select_blosum
 *
 * @param[in] one a character from the first sequence
 * @param[in] two a character from the first sequence
 * @return the score based on the blosum substitution
 */
int pg_match_blosum(char one, char two);


/**
 * Calculates the score if the given sequence were aligned with itself.
 *
 * @param[in] s the sequence
 * @param[in] callback the match function callback
 * @return the self score
 */
int pg_self_score(const sequence_t *seq, match_t callback);


/**
 * Calculates the score if the given sequence were aligned with itself.
 *
 * @param[in] s the sequence
 * @param[in] callback the match function callback
 * @return the self score
 */
int pg_self_score_blosum(const sequence_t *seq);


/**
 * Allocates a 2-dimensional cell_t array for sequence alignment.
 *
 * @param[in] nrow number of rows (always 2)
 * @param[in] ncol length of longest sequence to align
 * @return the table
 */
cell_t **pg_alloc_tbl(int nrow, int ncol);


/**
 * Allocates a 2-dimensional int array (for sequence alignment).
 *
 * @param nrow number of rows
 * @param ncol number of columns e.g. length of longest sequence to align
 * @return the table
 */
int **pg_alloc_int(int nrow, int ncol);


/**
 * Frees the 2-dimensional cell_t table. Symmetric to alloc_tbl().
 *
 * @param[in] tbl the cell_t table
 * @param[in] the number of rows (always 2)
 */
void pg_free_tbl(cell_t **tbl, int nrow);


/**
 * Frees the 2-dimensional int table. Symmetric to alloc_int().
 *
 * @param[in] tbl the int table
 * @param[in] the number of rows (always 2)
 */
void pg_free_int(int **tbl, int nrow);


/**
 * Select which BLOSUM matrix to use during pg_affine_gap_align().
 * 
 * @param[in] number the blosum number to use
 */
void pg_select_blosum(int number);


/**
 * Implementation of affine gap pairwise sequence alignment.
 *
 * It is a space efficient version: only two rows are required; also mem for
 * all dynamic tables are allocated ONLY ONCE outside of this function call
 * using alloc_tbl() and alloc_int() and passed as tbl, del, and ins arguments.
 *
 * @param[in] s1 sequence s1
 * @param[in] s2 sequence s2
 * @param[out] result alignment result <score, ndig, alen>
 * @param[in] tbl pre-allocated score table
 * @param[in] del pre-allocated deletion table
 * @param[in] ins pre-allocated insertion table
 * @param[in] open gap penalty
 * @param[in] gap extension penalty
 */
void pg_affine_gap_align(const sequence_t *s1, const sequence_t *s2,
                         cell_t *result, cell_t **tbl, int **del, int **ins,
                         int open, int gap, match_t callback);


/**
 * Implementation of affine gap pairwise protein sequence alignment using
 * blosum.
 *
 * It is a space efficient version: only two rows are required; also mem for
 * all dynamic tables are allocated ONLY ONCE outside of this function call
 * using alloc_tbl() and alloc_int() and passed as tbl, del, and ins arguments.
 *
 * @param[in] s1 sequence s1
 * @param[in] s2 sequence s2
 * @param[out] result alignment result <score, ndig, alen>
 * @param[in] tbl pre-allocated score table
 * @param[in] del pre-allocated deletion table
 * @param[in] ins pre-allocated insertion table
 * @param[in] open gap penalty
 * @param[in] gap extension penalty
 */
void pg_affine_gap_align_blosum(
        const sequence_t *s1, const sequence_t *s2,
        cell_t *result, cell_t **tbl, int **del, int **ins,
        int open, int gap);


/**
 * Implementation of affine gap pairwise protein sequence alignment using
 * blosum.
 *
 * It is a space efficient version: only two rows are required; also mem for
 * all dynamic tables are allocated ONLY ONCE outside of this function call
 * using alloc_tbl() and alloc_int() and passed as tbl, del, and ins arguments.
 *
 * @param[in] s1 sequence s1
 * @param[in] s2 sequence s2
 * @param[out] result alignment result <score, ndig, alen>
 * @param[in] tbl pre-allocated score table
 * @param[in] del pre-allocated deletion table
 * @param[in] ins pre-allocated insertion table
 * @param[in] open gap penalty
 * @param[in] gap extension penalty
 */
void pg_affine_gap_align_blosum2(
        const char *s1, size_t l1,
        const char *s2, size_t l2,
        cell_t *result, cell_t **tbl, int **del, int **ins,
        int open, int gap);


/**
 * Prints the given row to stdout (for debugging purposes).
 *
 * @param[in] tbl the cell_t table
 * @param[in] i the row to print
 * @param[in] the number of columns in the table
 */
void pg_print_row(cell_t **tbl, int i, int ncol);


/**
 * Asks whether the given cell_t alignment result is an edge, based on the
 * given param parameters.
 *
 * @param[in] result the dynamic programming alignment result
 * @param[in] s1 sequence s1
 * @param[in] s2 sequence s2
 * @param[in] param the heuristic parameters for edge determination
 * @param[out] sscore self score
 * @param[out] maxLen longer of s1Len and s2Len
 * @return TRUE if this is an edge, FALSE otherwise
 */
int pg_is_edge(const cell_t result, const sequence_t *s1, const sequence_t *s2,
               const param_t param, int *sscore, size_t *maxLen,
               match_t callback);


/**
 * Asks whether the given cell_t alignment result is an edge, based on the
 * given param parameters.
 *
 * @param[in] result the dynamic programming alignment result
 * @param[in] s1 sequence s1
 * @param[in] s2 sequence s2
 * @param[in] param the heuristic parameters for edge determination
 * @param[out] sscore self score
 * @param[out] maxLen longer of s1Len and s2Len
 * @return TRUE if this is an edge, FALSE otherwise
 */
int pg_is_edge_blosum(const cell_t result,
                      const sequence_t *s1, const sequence_t *s2,
                      const param_t param, int *sscore, size_t *maxLen);

#ifdef __cplusplus
}
#endif

#endif /* DYNAMIC_H_ */
