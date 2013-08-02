/**
 * @file alignment.hpp
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
#ifndef _PGRAPH_ALIGNMENT_H_
#define _PGRAPH_ALIGNMENT_H_

#include <cstddef>

using std::size_t;

namespace pgraph {

/**
 * Dynamic programming sequence alignment table cell and result.
 */
typedef struct {
    int score;      /**< alignment score */
    int matches;    /**< number of matches */
    int length;     /**< alignment length */
} cell_t;

/** callback function type for generic dynamic programming (mis)match score */
typedef int (*match_t)(char one, char two);

/**
 * Allocates a 2-dimensional cell_t array for sequence alignment.
 *
 * @see free_cell_table().
 *
 * @param[in] nrow number of rows
 * @param[in] ncol length of longest sequence to align
 * @return the table
 */
cell_t **allocate_cell_table(size_t nrow, size_t ncol);


/**
 * Allocates a 2-dimensional int array for sequence alignment.
 *
 * @see free_int_table().
 *
 * @param nrow number of rows
 * @param ncol number of columns e.g. length of longest sequence to align
 * @return the table
 */
int **allocate_int_table(size_t nrow, size_t ncol);


/**
 * Frees the 2-dimensional cell_t table.
 *
 * @see allocate_cell_table().
 *
 * @param[in] table the cell_t table
 * @param[in] nrow the number of rows
 */
void free_cell_table(cell_t **table, size_t nrow);


/**
 * Frees the 2-dimensional int table.
 *
 * @see allocate_int_table().
 *
 * @param[in] table the int table
 * @param[in] nrow the number of rows
 */
void free_int_table(int **table, size_t nrow);


/**
 * Calculates the score if the given sequence were aligned with itself.
 *
 * This function takes a substitution matrix 'sub' which is used to provide the
 * scoring criteria. The substitution matrix here is typed as an array of
 * arrays (a ragged array) although the length of each row is uniform. The
 * substitution matrix is indexed using the 'map' and 'first' parameters. The
 * 'first' value is subtracted off of the current character being aligned with
 * itself. For example, index = 'B' - 'A' where 'A' is 'first'. Then, instead
 * of simply sub[first][first], we allow one more level of indirection in case
 * the substitution matrix is not sorted by character. The final result looks
 * like:
 *
 * @code
 * char mychar = 'C';
 * int score = sub[map[mychar - first]][map[mychar - first]];
 * @endcode
 *
 * You would call this funtion like so (using blosum for protein sequences):
 *
 * @code
 * int result = self_score(s1, s1_len, _blosum62, MAP_BLOSUM, 'A');
 * @endcode
 *
 * @param[in] seq the sequence
 * @param[in] len the length of the sequence
 * @param[in] sub substitution matrix
 * @param[in] map index mapping for substitution matrix
 * @param[in] first character offset (see doc above)
 * @return the self score
 */
int self_score(const char * const restrict seq, size_t len,
               const int * const restrict * const restrict sub,
               const int * const restrict map, char first);


/**
 * Calculates the score if the given sequence were aligned with itself.
 *
 * This funtion takes a callback function for calculating the substitution
 * score of a character with itself.
 *
 * @param[in] seq the sequence
 * @param[in] len the length of the sequence
 * @param[in] callback function to calculate any character match/mismatch
 * @return the self score
 */
int self_score(const char * const restrict seq, size_t len, match_t callback);


/**
 * Calculates the score if the given sequence were aligned with itself.
 *
 * This funtion takes a simply match and mismatch score.
 *
 * @param[in] seq the sequence
 * @param[in] len the length of the sequence
 * @param[in] match the static score of a character match
 * @return the self score
 */
int self_score(const char * const restrict seq, size_t len, int match);


/**
 * Calculates the score if the given sequence were aligned with itself.
 *
 * This funtion uses internally one of the blosum substitution matrices set
 * previously by using pgraph::select_blosum(int).
 *
 * @param[in] seq the sequence
 * @param[in] len the length of the sequence
 * @return the self score
 */
int self_score_blosum(const char * const restrict seq, size_t len);


/**
 * Select which BLOSUM matrix to use during affine_gap_align().
 *
 * @param[in] number the blosum number to use
 */
void select_blosum(int number);


/**
 * Implementation of affine gap pairwise sequence alignment.
 *
 * It is a space efficient version using only two rows. Also, memory for
 * all dynamic tables are allocated ONLY ONCE outside of this function call
 * using allocate_cell_table() and allocate_int_table() and passed as
 * tbl, del, and ins arguments.
 *
 * @param[in] s1 character sequence one
 * @param[in] s1_len length of character sequence one
 * @param[in] s2 character sequence two
 * @param[in] s2_len length of character sequence two
 * @param[in] tbl pre-allocated score table
 * @param[in] del pre-allocated deletion table
 * @param[in] ins pre-allocated insertion table
 * @param[in] open gap penalty
 * @param[in] gap extension penalty
 * @param[in] sub substitution matrix e.g. BLOSUM
 * @param[in] map mapping for substitution matrix
 * @param[in] first (zeroth) character of mapping alphabet e.g. 'A'
 * @return alignment result
 */
cell_t affine_gap_align(
        const char * const restrict s1, size_t s1_len,
        const char * const restrict s2, size_t s2_len,
        cell_t ** const restrict tbl,
        int ** const restrict del,
        int ** const restrict ins,
        int open, int gap,
        const int * const restrict * const restrict sub,
        const int * const restrict map, char first);


/**
 * Implementation of affine gap pairwise sequence alignment.
 *
 * It is a space efficient version using only two rows. Also, memory for
 * all dynamic tables are allocated ONLY ONCE outside of this function call
 * using allocate_cell_table() and allocate_int_table() and passed as
 * tbl, del, and ins arguments.
 *
 * @param[in] s1 character sequence one
 * @param[in] s1_len length of character sequence one
 * @param[in] s2 character sequence two
 * @param[in] s2_len length of character sequence two
 * @param[in] tbl pre-allocated score table
 * @param[in] del pre-allocated deletion table
 * @param[in] ins pre-allocated insertion table
 * @param[in] open gap penalty
 * @param[in] gap extension penalty
 * @param[in] match callback function to calculate any character match/mismatch
 * @return alignment result
 */
cell_t affine_gap_align(
        const char * const restrict s1, size_t s1_len,
        const char * const restrict s2, size_t s2_len,
        cell_t ** const restrict tbl,
        int ** const restrict del,
        int ** const restrict ins,
        int open, int gap,
        match_t match);


/**
 * Implementation of affine gap pairwise sequence alignment.
 *
 * It is a space efficient version using only two rows. Also, memory for
 * all dynamic tables are allocated ONLY ONCE outside of this function call
 * using allocate_cell_table() and allocate_int_table() and passed as
 * tbl, del, and ins arguments.
 *
 * @param[in] s1 character sequence one
 * @param[in] s1_len length of character sequence one
 * @param[in] s2 character sequence two
 * @param[in] s2_len length of character sequence two
 * @param[in] tbl pre-allocated score table
 * @param[in] del pre-allocated deletion table
 * @param[in] ins pre-allocated insertion table
 * @param[in] open gap penalty
 * @param[in] gap extension penalty
 * @param[in] match score
 * @param[in] mismatch score
 * @return alignment result
 */
cell_t affine_gap_align(
        const char * const restrict s1, size_t s1_len,
        const char * const restrict s2, size_t s2_len,
        cell_t ** const restrict tbl,
        int ** const restrict del,
        int ** const restrict ins,
        int open, int gap,
        int match, int mismatch);


/**
 * Implementation of affine gap pairwise sequence alignment using blosum.
 *
 * It is a space efficient version using only two rows. Also, memory for
 * all dynamic tables are allocated ONLY ONCE outside of this function call
 * using allocate_cell_table() and allocate_int_table() and passed as
 * tbl, del, and ins arguments.
 *
 * @param[in] s1 character sequence one
 * @param[in] s1_len length of character sequence one
 * @param[in] s2 character sequence two
 * @param[in] s2_len length of character sequence two
 * @param[in] tbl pre-allocated score table
 * @param[in] del pre-allocated deletion table
 * @param[in] ins pre-allocated insertion table
 * @param[in] open gap penalty
 * @param[in] gap extension penalty
 * @return alignment result
 */
cell_t affine_gap_align_blosum(
        const char * const restrict s1, size_t s1_len,
        const char * const restrict s2, size_t s2_len,
        cell_t ** const restrict tbl,
        int ** const restrict del,
        int ** const restrict ins,
        int open, int gap);


/**
 * Asks whether the given cell_t alignment result is an edge, based on
 * the given parameters.
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
 * @param[in] sub substitution matrix
 * @param[in] map index mapping for substitution matrix
 * @param[in] first character offset (see doc above)
 * @return true if this is an edge, false otherwise
 */
bool is_edge(
        const cell_t &result,
        const char * const restrict s1, size_t s1_len,
        const char * const restrict s2, size_t s2_len,
        int AOL,
        int SIM,
        int OS,
        int &self_score,
        size_t &max_len,
        const int * const restrict * const restrict sub,
        const int * const restrict map, char first);


/**
 * Asks whether the given cell_t alignment result is an edge, based on
 * the given parameters.
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
 * @param[in] match callback function to calculate any character match/mismatch
 * @return true if this is an edge, false otherwise
 */
bool is_edge(
        const cell_t &result,
        const char * const restrict s1, size_t s1_len,
        const char * const restrict s2, size_t s2_len,
        int AOL,
        int SIM,
        int OS,
        int &self_score,
        size_t &max_len,
        match_t match);


/**
 * Asks whether the given cell_t alignment result is an edge, based on
 * the given parameters.
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
 * @param[in] match the static score of a character match
 * @return true if this is an edge, false otherwise
 */
bool is_edge(
        const cell_t &result,
        const char * const restrict s1, size_t s1_len,
        const char * const restrict s2, size_t s2_len,
        int AOL,
        int SIM,
        int OS,
        int &self_score,
        size_t &max_len,
        int match);


/**
 * Asks whether the given cell_t alignment result is an edge, based on
 * the given parameters.
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
 * @return true if this is an edge, false otherwise
 */
bool is_edge_blosum(
        const cell_t &result,
        const char * const restrict s1, size_t s1_len,
        const char * const restrict s2, size_t s2_len,
        int AOL,
        int SIM,
        int OS,
        int &self_score,
        size_t &max_len);


/** alphabet for BLOSUM in the order in which the matrix is defined */
static char ALPHABET_BLOSUM[] = {
        'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
        'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V',
        'B', 'Z', 'X', 'J', 'O', 'U' };


/** mapping from BLOSUM alphabet to BLOSUM index; use as BLOSUM[map[ch-'A']] */
static int MAP_BLOSUM[] = {
          0,  20,   4,   3,   6,  13,   7,   8,   9,  23,
         11,  10,  12,   2,  24,  14,   5,   1,  15,  16,
         25,  19,  17,  22,  18,  21 };


static int ALPHABET_DNA[] = { 'A', 'C', 'G', 'T' };


static int MAP_DNA[] = {
          0,  -1,   1,  -1,  -1,  -1,   2,  -1,  -1,  -1,
         -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,   3,
         -1,  -1,  -1,  -1,  -1,  -1 };


static int SUB_DNA[4][4] = {
        /*       A   C   G   T */
        /* A */{ 1,  0,  0,  0},
        /* C */{ 0,  1,  0,  0},
        /* G */{ 0,  0,  1,  0},
        /* T */{ 0,  0,  0,  1}
};

}; /* namespace pgraph */

#endif /* _PGRAPH_ALIGNMENT_H_ */
