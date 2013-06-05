/**
 * @file align.h
 *
 * @author andy.cj.wu@gmail.com
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2010 Washington State University. All rights reserved.
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 *
 * Routines for dynamic programming of sequence alignment. These routines
 * strive to be as performance-optimal as possible as well as only using
 * standard C data types rather than higher-level C structs. Higher-level C
 * structs would likely simplify the otherwise long argument lists to many of
 * these routines.
 */
#ifndef _PGRAPH_ALIGN_H_
#define _PGRAPH_ALIGN_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <stddef.h>

/**
 * Dynamic programming sequence alignment table cell and result.
 */
typedef struct {
    int score;      /**< alignment score */
    int matches;    /**< number of matches */
    int length;     /**< alignment length */
} pgraph_cell_t;


/**
 * Allocates a 2-dimensional pgraph_cell_t array for sequence alignment.
 *
 * @see pgraph_free_cell_table().
 *
 * @param[in] nrow number of rows
 * @param[in] ncol length of longest sequence to align
 * @return the table
 */
pgraph_cell_t **pgraph_malloc_cell_table(size_t nrow, size_t ncol);


/**
 * Allocates a 2-dimensional int array for sequence alignment.
 *
 * @see pgraph_free_int_table().
 *
 * @param nrow number of rows
 * @param ncol number of columns e.g. length of longest sequence to align
 * @return the table
 */
int **pgraph_malloc_int_table(size_t nrow, size_t ncol);


/**
 * Frees the 2-dimensional pgraph_cell_t table.
 *
 * @see pgraph_malloc_cell_table().
 *
 * @param[in] table the pgraph_cell_t table
 * @param[in] the number of rows
 */
void pgraph_free_cell_table(pgraph_cell_t **table, size_t nrow);


/**
 * Frees the 2-dimensional int table.
 *
 * @see pgraph_malloc_int_table().
 *
 * @param[in] table the int table
 * @param[in] the number of rows
 */
void pgraph_free_int_table(int **table, size_t nrow);


/**
 * Calculates the score if the given sequence were aligned with itself.
 *
 * Simple scoring where any match results in the same given value.
 *
 * @param[in] seq the sequence
 * @param[in] len the length of the sequence
 * @param[in] match the score of any character match
 * @return the self score
 */
int pgraph_self_score(
        const char * const restrict seq, size_t len,
        const int ** const restrict sub,
        const int * const restrict map, char first);


/**
 * Implementation of affine gap pairwise sequence alignment.
 *
 * It is a space efficient version using only two rows. Also, memory for
 * all dynamic tables are allocated ONLY ONCE outside of this function call
 * using pgraph_malloc_cell_table() and pgraph_malloc_int_table() and passed as
 * tbl, del, and ins arguments.
 *
 * @param[in] s1 character sequence one
 * @param[in] s1_len length of character sequence one
 * @param[in] s2 character sequence two
 * @param[in] s2_len length of character sequence two
 * @param[out] result alignment result
 * @param[in] tbl pre-allocated score table
 * @param[in] del pre-allocated deletion table
 * @param[in] ins pre-allocated insertion table
 * @param[in] open gap penalty
 * @param[in] gap extension penalty
 * @param[in] sub substitution matrix e.g. BLOSUM
 * @param[in] map mapping for substitution matrix
 * @param[in] first (zeroth) character of mapping alphabet e.g. 'A'
 */
void pgraph_affine_gap_align(
        const char * const restrict s1, size_t s1_len,
        const char * const restrict s2, size_t s2_len,
        pgraph_cell_t * const restrict result,
        pgraph_cell_t ** const restrict tbl,
        int ** const restrict del,
        int ** const restrict ins,
        int open, int gap,
        const int ** const restrict sub,
        const int * const restrict map, char first);


/**
 * Asks whether the given pgraph_cell_t alignment result is an edge, based on
 * the given param parameters.
 *
 * @param[in] result the dynamic programming alignment result
 * @param[in] s1 character sequence one
 * @param[in] s1_len length of character sequence one
 * @param[in] s2 character sequence two
 * @param[in] s2_len length of character sequence two
 * @param[in] AOL alignment over longer sequence (heuristic)
 * @param[in] SIM match similarity (heuristic)
 * @param[in] OS optimal score over self score (heuristic)
 * @param[out] self_score self score
 * @param[out] max_len longer of s1_len and s2_len
 * @return 1 if this is an edge, 0 otherwise
 */
int pgraph_is_edge(
        const pgraph_cell_t result,
        const char * const restrict s1, size_t s1_len,
        const char * const restrict s2, size_t s2_len,
        int AOL,
        int SIM,
        int OS,
        int * const restrict self_score,
        size_t * const restrict max_len,
        const int ** const restrict sub,
        const int * const restrict map, char first);

/** alphabet for BLOSUM in the order in which the matrix is defined */
static char PGRAPH_ALPHABET_BLOSUM[] = {
        'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
        'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V',
        'B', 'Z', 'X', 'J', 'O', 'U' };

/** mapping from BLOSUM alphabet to BLOSUM index; use as BLOSUM[map[ch-'A']] */
static int PGRAPH_MAP_BLOSUM[] = {
          0,  20,   4,   3,   6,  13,   7,   8,   9,  23,
         11,  10,  12,   2,  24,  14,   5,   1,  15,  16,
         25,  19,  17,  22,  18,  21 };

static int PGRAPH_ALPHABET_DNA[] = { 'A', 'C', 'G', 'T' };

static int PGRAPH_MAP_DNA[] = {
          0,  -1,   1,  -1,  -1,  -1,   2,  -1,  -1,  -1,
         -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,   3,
         -1,  -1,  -1,  -1,  -1,  -1 };

static int PGRAPH_SUB_DNA[4][4] = {
        /*       A   C   G   T */
        /* A */{ 1,  0,  0,  0},
        /* C */{ 0,  1,  0,  0},
        /* G */{ 0,  0,  1,  0},
        /* T */{ 0,  0,  0,  1}
};

#ifdef __cplusplus
}
#endif

#endif /* _PGRAPH_ALIGN_H_ */
