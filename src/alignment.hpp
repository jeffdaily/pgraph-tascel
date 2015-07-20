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

#include "parasail.h"

using std::size_t;

namespace pgraph {

/**
 * Calculates the score if the given sequence were aligned with itself.
 *
 * @param[in] seq the sequence
 * @param[in] len the length of the sequence
 * @param[in] sub substitution matrix
 * @return the self score
 */
int self_score(const char * const restrict seq, size_t len,
               const parasail_matrix_t *matrix);

/** @name Edge Functions
 *
 * Asks whether the given parasail alignment result is an edge, based on
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
        const parasail_result_t *result,
        const char * const restrict s1, size_t s1_len,
        const char * const restrict s2, size_t s2_len,
        int AOL,
        int SIM,
        int OS,
        int &self_score,
        size_t &max_len,
        const parasail_matrix_t *matrix);

/** @} */

}; /* namespace pgraph */

#endif /* _PGRAPH_ALIGNMENT_H_ */
