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

#define USE_SIMILARITIES 0

namespace pgraph {

/**
 * Dynamic programming sequence alignment table cell and result.
 */
typedef struct {
    int score;      /**< alignment score */
    int matches;    /**< number of matches */
#if USE_SIMILARITIES
    int similarities;/**< number of similarities */
#endif
    int length;     /**< alignment length */
} cell_t;

/**
 * Dynamic programming sequence alignment table cell and result.
 */
typedef struct {
    int score;      /**< alignment score */
    int matches;    /**< number of matches */
#if USE_SIMILARITIES
    int similarities;/**< number of similarities */
#endif
    int length;     /**< alignment length */
    int del;
    int ins;
} tbl_t;

/** callback function type for generic dynamic programming (mis)match score */
typedef int (*match_t)(char one, char two);

/**
 * Allocates a 2-dimensional array for sequence alignment.
 *
 * @see free_table().
 *
 * @param[in] nrow number of rows
 * @param[in] ncol length of longest sequence to align
 * @return the table
 */
template <typename T>
T **allocate_table(size_t nrow, size_t ncol)
{
    size_t i = 0;
    T **table = NULL;

    table = new T*[nrow];
    for (i = 0; i < nrow; i++) {
        /* +1 for dynamic align */
        table[i] = new T[ncol + 1];
    }

    return table;
}


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
 * Allocates a 2-dimensional cell_t array for sequence alignment.
 *
 * @see free_cell_table().
 *
 * @param[in] nrow number of rows
 * @param[in] ncol length of longest sequence to align
 * @return the table
 */
tbl_t **allocate_tbl_table(size_t nrow, size_t ncol);


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
 * Frees the 2-dimensional table.
 *
 * @see allocate_table().
 *
 * @param[in] table the table
 * @param[in] nrow the number of rows
 */
template <typename T>
void free_table(T * const restrict * const restrict table, size_t nrow)
{
    size_t i;
    for (i = 0; i < nrow; i++) {
        delete [] table[i];
    }
    delete [] table;
}


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


/** @name Alignment Functions
 *
 * Implementation of pairwise sequence alignment with affine gap penalty.
 *
 * The align_global_affine() functions perform a global (NW) alignment.
 * The align_semi_affine() functions perform a semi-global alignment.
 * The align_local_affine() functions perform a local (SW) alignment.
 *
 * These functions are all space efficient versions using only two rows of the
 * dynamic programming table.
 *
 * Also, memory for all dynamic tables are may be optionally allocated ONLY
 * ONCE outside of this function call using allocate_cell_table() and
 * allocate_int_table() and passed as tbl, del, and ins arguments. If those
 * arguments are left as NULL, memory will be allocated and subsequently
 * deleted during the execution of this function (not efficient).
 *
 * If you want to use the built-in BLOSUM substitution matrix, which is
 * selected using the select_blosum(int) function, then use the align function
 * which does not explicitly have a paremeter for the matching score.
 *
 * If you want to use a callback function to calcuate the match/mismatch score,
 * use the function with the following parameter.
 * @param[in] match callback function to calculate any character match/mismatch
 *
 * If you want to use a fixed/static/unchanging match and mismatch score, use
 * the funtion with the following two parameters.
 * @param[in] match score
 * @param[in] mismatch score
 *
 * If you want to supply your own substitution matrix (and means of indexing
 * into such a matrix), use the function with the following parameters
 * @param[in] sub substitution matrix e.g. BLOSUM
 * @param[in] map mapping for substitution matrix
 * @param[in] first (zeroth) character of mapping alphabet e.g. 'A'
 *
 * When using your own substitution matrix, the score is calculated like so:
 * @code
 * char c1 = 'A';
 * char c2 = 'B';
 * int score = sub[map[c1 - first]][map[c2 - first]];
 * @endcode
 *
 * Parameters which are common to each align function include the following.
 *
 * @param[in] s1 character sequence one
 * @param[in] s1_len length of character sequence one
 * @param[in] s2 character sequence two
 * @param[in] s2_len length of character sequence two
 * @param[in] open gap penalty
 * @param[in] gap extension penalty
 * @param[in] tbl pre-allocated score table
 * @param[in] del pre-allocated deletion table
 * @param[in] ins pre-allocated insertion table
 * @return alignment result
 */

/** @{ */

cell_t align_global_affine(
        const char * const restrict s1, size_t s1_len,
        const char * const restrict s2, size_t s2_len,
        match_t match,
        int open=-10, int gap=-1,
        cell_t * const restrict * const restrict tbl=NULL,
        int * const restrict * const restrict del=NULL,
        int * const restrict * const restrict ins=NULL);

cell_t align_global_affine(
        const char * const restrict s1, size_t s1_len,
        const char * const restrict s2, size_t s2_len,
        const int * const restrict * const restrict sub,
        const int * const restrict map, char first,
        int open=-10, int gap=-1,
        cell_t * const restrict * const restrict tbl=NULL,
        int * const restrict * const restrict del=NULL,
        int * const restrict * const restrict ins=NULL);

cell_t align_global_affine(
        const char * const restrict s1, size_t s1_len,
        const char * const restrict s2, size_t s2_len,
        int match, int mismatch,
        int open=-10, int gap=-1,
        cell_t * const restrict * const restrict tbl=NULL,
        int * const restrict * const restrict del=NULL,
        int * const restrict * const restrict ins=NULL);

cell_t align_global_affine(
        const char * const restrict s1, size_t s1_len,
        const char * const restrict s2, size_t s2_len,
        int open=-10, int gap=-1,
        cell_t * const restrict * const restrict tbl=NULL,
        int * const restrict * const restrict del=NULL,
        int * const restrict * const restrict ins=NULL);

cell_t align_semi_affine(
        const char * const restrict s1, size_t s1_len,
        const char * const restrict s2, size_t s2_len,
        match_t match,
        int open=-10, int gap=-1,
        cell_t * const restrict * const restrict tbl=NULL,
        int * const restrict * const restrict del=NULL,
        int * const restrict * const restrict ins=NULL);

cell_t align_semi_affine(
        const char * const restrict s1, size_t s1_len,
        const char * const restrict s2, size_t s2_len,
        const int * const restrict * const restrict sub,
        const int * const restrict map, char first,
        int open=-10, int gap=-1,
        cell_t * const restrict * const restrict tbl=NULL,
        int * const restrict * const restrict del=NULL,
        int * const restrict * const restrict ins=NULL);

cell_t align_semi_affine(
        const char * const restrict s1, size_t s1_len,
        const char * const restrict s2, size_t s2_len,
        int match, int mismatch,
        int open=-10, int gap=-1,
        cell_t * const restrict * const restrict tbl=NULL,
        int * const restrict * const restrict del=NULL,
        int * const restrict * const restrict ins=NULL);

cell_t align_semi_affine(
        const char * const restrict s1, size_t s1_len,
        const char * const restrict s2, size_t s2_len,
        int open=-10, int gap=-1,
        cell_t * const restrict * const restrict tbl=NULL,
        int * const restrict * const restrict del=NULL,
        int * const restrict * const restrict ins=NULL);

cell_t align_local_affine(
        const char * const restrict s1, size_t s1_len,
        const char * const restrict s2, size_t s2_len,
        match_t match,
        int open=-10, int gap=-1,
        cell_t * const restrict * const restrict tbl=NULL,
        int * const restrict * const restrict del=NULL,
        int * const restrict * const restrict ins=NULL);

cell_t align_local_affine(
        const char * const restrict s1, size_t s1_len,
        const char * const restrict s2, size_t s2_len,
        const int * const restrict * const restrict sub,
        const int * const restrict map, char first,
        int open=-10, int gap=-1,
        cell_t * const restrict * const restrict tbl=NULL,
        int * const restrict * const restrict del=NULL,
        int * const restrict * const restrict ins=NULL);

cell_t align_local_affine(
        const char * const restrict s1, size_t s1_len,
        const char * const restrict s2, size_t s2_len,
        int match, int mismatch,
        int open=-10, int gap=-1,
        cell_t * const restrict * const restrict tbl=NULL,
        int * const restrict * const restrict del=NULL,
        int * const restrict * const restrict ins=NULL);

cell_t align_local_affine(
        const char * const restrict s1, size_t s1_len,
        const char * const restrict s2, size_t s2_len,
        int open=-10, int gap=-1,
        cell_t * const restrict * const restrict tbl=NULL,
        int * const restrict * const restrict del=NULL,
        int * const restrict * const restrict ins=NULL);

cell_t align_global_affine(
        const char * const restrict s1, size_t s1_len,
        const char * const restrict s2, size_t s2_len,
        match_t match,
        int open=-10, int gap=-1,
        int * const restrict * const restrict scr=NULL,
        int * const restrict * const restrict mat=NULL,
#if USE_SIMILARITIES
        int * const restrict * const restrict sim=NULL,
#endif
        int * const restrict * const restrict len=NULL,
        int * const restrict * const restrict del=NULL,
        int * const restrict * const restrict ins=NULL);

cell_t align_global_affine(
        const char * const restrict s1, size_t s1_len,
        const char * const restrict s2, size_t s2_len,
        const int * const restrict * const restrict sub,
        const int * const restrict map, char first,
        int open=-10, int gap=-1,
        int * const restrict * const restrict scr=NULL,
        int * const restrict * const restrict mat=NULL,
#if USE_SIMILARITIES
        int * const restrict * const restrict sim=NULL,
#endif
        int * const restrict * const restrict len=NULL,
        int * const restrict * const restrict del=NULL,
        int * const restrict * const restrict ins=NULL);

cell_t align_global_affine(
        const char * const restrict s1, size_t s1_len,
        const char * const restrict s2, size_t s2_len,
        int match, int mismatch,
        int open=-10, int gap=-1,
        int * const restrict * const restrict scr=NULL,
        int * const restrict * const restrict mat=NULL,
#if USE_SIMILARITIES
        int * const restrict * const restrict sim=NULL,
#endif
        int * const restrict * const restrict len=NULL,
        int * const restrict * const restrict del=NULL,
        int * const restrict * const restrict ins=NULL);

cell_t align_global_affine(
        const char * const restrict s1, size_t s1_len,
        const char * const restrict s2, size_t s2_len,
        int open=-10, int gap=-1,
        int * const restrict * const restrict scr=NULL,
        int * const restrict * const restrict mat=NULL,
#if USE_SIMILARITIES
        int * const restrict * const restrict sim=NULL,
#endif
        int * const restrict * const restrict len=NULL,
        int * const restrict * const restrict del=NULL,
        int * const restrict * const restrict ins=NULL);

cell_t align_semi_affine(
        const char * const restrict s1, size_t s1_len,
        const char * const restrict s2, size_t s2_len,
        match_t match,
        int open=-10, int gap=-1,
        int * const restrict * const restrict scr=NULL,
        int * const restrict * const restrict mat=NULL,
#if USE_SIMILARITIES
        int * const restrict * const restrict sim=NULL,
#endif
        int * const restrict * const restrict len=NULL,
        int * const restrict * const restrict del=NULL,
        int * const restrict * const restrict ins=NULL);

cell_t align_semi_affine(
        const char * const restrict s1, size_t s1_len,
        const char * const restrict s2, size_t s2_len,
        const int * const restrict * const restrict sub,
        const int * const restrict map, char first,
        int open=-10, int gap=-1,
        int * const restrict * const restrict scr=NULL,
        int * const restrict * const restrict mat=NULL,
#if USE_SIMILARITIES
        int * const restrict * const restrict sim=NULL,
#endif
        int * const restrict * const restrict len=NULL,
        int * const restrict * const restrict del=NULL,
        int * const restrict * const restrict ins=NULL);

cell_t align_semi_affine(
        const char * const restrict s1, size_t s1_len,
        const char * const restrict s2, size_t s2_len,
        int match, int mismatch,
        int open=-10, int gap=-1,
        int * const restrict * const restrict scr=NULL,
        int * const restrict * const restrict mat=NULL,
#if USE_SIMILARITIES
        int * const restrict * const restrict sim=NULL,
#endif
        int * const restrict * const restrict len=NULL,
        int * const restrict * const restrict del=NULL,
        int * const restrict * const restrict ins=NULL);

cell_t align_semi_affine(
        const char * const restrict s1, size_t s1_len,
        const char * const restrict s2, size_t s2_len,
        int open=-10, int gap=-1,
        int * const restrict * const restrict scr=NULL,
        int * const restrict * const restrict mat=NULL,
#if USE_SIMILARITIES
        int * const restrict * const restrict sim=NULL,
#endif
        int * const restrict * const restrict len=NULL,
        int * const restrict * const restrict del=NULL,
        int * const restrict * const restrict ins=NULL);

cell_t align_local_affine(
        const char * const restrict s1, size_t s1_len,
        const char * const restrict s2, size_t s2_len,
        match_t match,
        int open=-10, int gap=-1,
        int * const restrict * const restrict scr=NULL,
        int * const restrict * const restrict mat=NULL,
#if USE_SIMILARITIES
        int * const restrict * const restrict sim=NULL,
#endif
        int * const restrict * const restrict len=NULL,
        int * const restrict * const restrict del=NULL,
        int * const restrict * const restrict ins=NULL);

cell_t align_local_affine(
        const char * const restrict s1, size_t s1_len,
        const char * const restrict s2, size_t s2_len,
        const int * const restrict * const restrict sub,
        const int * const restrict map, char first,
        int open=-10, int gap=-1,
        int * const restrict * const restrict scr=NULL,
        int * const restrict * const restrict mat=NULL,
#if USE_SIMILARITIES
        int * const restrict * const restrict sim=NULL,
#endif
        int * const restrict * const restrict len=NULL,
        int * const restrict * const restrict del=NULL,
        int * const restrict * const restrict ins=NULL);

cell_t align_local_affine(
        const char * const restrict s1, size_t s1_len,
        const char * const restrict s2, size_t s2_len,
        int match, int mismatch,
        int open=-10, int gap=-1,
        int * const restrict * const restrict scr=NULL,
        int * const restrict * const restrict mat=NULL,
#if USE_SIMILARITIES
        int * const restrict * const restrict sim=NULL,
#endif
        int * const restrict * const restrict len=NULL,
        int * const restrict * const restrict del=NULL,
        int * const restrict * const restrict ins=NULL);

cell_t align_local_affine(
        const char * const restrict s1, size_t s1_len,
        const char * const restrict s2, size_t s2_len,
        int open=-10, int gap=-1,
        int * const restrict * const restrict scr=NULL,
        int * const restrict * const restrict mat=NULL,
#if USE_SIMILARITIES
        int * const restrict * const restrict sim=NULL,
#endif
        int * const restrict * const restrict len=NULL,
        int * const restrict * const restrict del=NULL,
        int * const restrict * const restrict ins=NULL);

cell_t align_global_affine(
        const char * const restrict s1, size_t s1_len,
        const char * const restrict s2, size_t s2_len,
        match_t match,
        int open=-10, int gap=-1,
        tbl_t * const restrict * const restrict tbl=NULL);

cell_t align_global_affine(
        const char * const restrict s1, size_t s1_len,
        const char * const restrict s2, size_t s2_len,
        const int * const restrict * const restrict sub,
        const int * const restrict map, char first,
        int open=-10, int gap=-1,
        tbl_t * const restrict * const restrict tbl=NULL);

cell_t align_global_affine(
        const char * const restrict s1, size_t s1_len,
        const char * const restrict s2, size_t s2_len,
        int match, int mismatch,
        int open=-10, int gap=-1,
        tbl_t * const restrict * const restrict tbl=NULL);

cell_t align_global_affine(
        const char * const restrict s1, size_t s1_len,
        const char * const restrict s2, size_t s2_len,
        int open=-10, int gap=-1,
        tbl_t * const restrict * const restrict tbl=NULL);

cell_t align_semi_affine(
        const char * const restrict s1, size_t s1_len,
        const char * const restrict s2, size_t s2_len,
        match_t match,
        int open=-10, int gap=-1,
        tbl_t * const restrict * const restrict tbl=NULL);

cell_t align_semi_affine(
        const char * const restrict s1, size_t s1_len,
        const char * const restrict s2, size_t s2_len,
        const int * const restrict * const restrict sub,
        const int * const restrict map, char first,
        int open=-10, int gap=-1,
        tbl_t * const restrict * const restrict tbl=NULL);

cell_t align_semi_affine(
        const char * const restrict s1, size_t s1_len,
        const char * const restrict s2, size_t s2_len,
        int match, int mismatch,
        int open=-10, int gap=-1,
        tbl_t * const restrict * const restrict tbl=NULL);

cell_t align_semi_affine(
        const char * const restrict s1, size_t s1_len,
        const char * const restrict s2, size_t s2_len,
        int open=-10, int gap=-1,
        tbl_t * const restrict * const restrict tbl=NULL);

cell_t align_local_affine(
        const char * const restrict s1, size_t s1_len,
        const char * const restrict s2, size_t s2_len,
        match_t match,
        int open=-10, int gap=-1,
        tbl_t * const restrict * const restrict tbl=NULL);

cell_t align_local_affine(
        const char * const restrict s1, size_t s1_len,
        const char * const restrict s2, size_t s2_len,
        const int * const restrict * const restrict sub,
        const int * const restrict map, char first,
        int open=-10, int gap=-1,
        tbl_t * const restrict * const restrict tbl=NULL);

cell_t align_local_affine(
        const char * const restrict s1, size_t s1_len,
        const char * const restrict s2, size_t s2_len,
        int match, int mismatch,
        int open=-10, int gap=-1,
        tbl_t * const restrict * const restrict tbl=NULL);

cell_t align_local_affine(
        const char * const restrict s1, size_t s1_len,
        const char * const restrict s2, size_t s2_len,
        int open=-10, int gap=-1,
        tbl_t * const restrict * const restrict tbl=NULL);

/** @} */

/**
 * Implementation of affine gap pairwise sequence alignment using blosum and
 * sse vector instructions.
 *
 * @param[in] s1 character sequence one
 * @param[in] s1_len length of character sequence one
 * @param[in] s2 character sequence two
 * @param[in] s2_len length of character sequence two
 * @param[in] open gap penalty
 * @param[in] gap extension penalty
 * @return alignment result
 */
cell_t align_local_affine_ssw(
        const char * const restrict s1, size_t s1_len,
        const char * const restrict s2, size_t s2_len,
        int open=-10, int gap=-1);


/**
 * Implementation of affine gap pairwise sequence alignment using blosum and
 * sse vector instructions from fasta package.
 *
 * @param[in] s1 character sequence one
 * @param[in] s1_len length of character sequence one
 * @param[in] s2 character sequence two
 * @param[in] s2_len length of character sequence two
 * @param[in] open gap penalty
 * @param[in] gap extension penalty
 * @return alignment result
 */
cell_t align_global_affine_fasta(
        const char * const restrict s1, size_t s1_len,
        const char * const restrict s2, size_t s2_len,
        int open=-10, int gap=-1);


/** @name Alignment Functions
 *
 * Asks whether the given cell_t alignment result is an edge, based on
 * the given parameters.
 *
 * Parameters which are common to each edge function include the following.
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
 *
 * If you want to use the built-in BLOSUM substitution matrix, which is
 * selected using the select_blosum(int) function, then use the align function
 * which does not explicitly have a paremeter for the matching score.
 *
 * If you want to use a callback function to calcuate the match/mismatch score,
 * use the function with the following parameter.
 * @param[in] match callback function to calculate any character match/mismatch
 *
 * If you want to use a fixed/static/unchanging match score, use the funtion
 * with the following parameter.
 * @param[in] match score
 *
 * If you want to supply your own substitution matrix (and means of indexing
 * into such a matrix), use the function with the following parameters
 * @param[in] sub substitution matrix e.g. BLOSUM
 * @param[in] map mapping for substitution matrix
 * @param[in] first (zeroth) character of mapping alphabet e.g. 'A'
 *
 * When using your own substitution matrix, the score is calculated like so:
 * @code
 * char c1 = 'A';
 * char c2 = 'B';
 * int score = sub[map[c1 - first]][map[c2 - first]];
 * @endcode
 *
 * @return true if this is an edge, false otherwise
 */

/** @{ */

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

bool is_edge(
        const cell_t &result,
        const char * const restrict s1, size_t s1_len,
        const char * const restrict s2, size_t s2_len,
        int AOL,
        int SIM,
        int OS,
        int &self_score,
        size_t &max_len);

/** @} */

}; /* namespace pgraph */

#endif /* _PGRAPH_ALIGNMENT_H_ */
