/**
 * @file alignment.cpp
 *
 * @author andy.cj.wu@gmail.com
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2010 Washington State University. All rights reserved.
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 */
#include "config.h"

#include <cassert>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <iostream>

#include "parasail.h"

#include "alignment.hpp"

using std::cerr;
using std::endl;

namespace pgraph {


int self_score(
        const char * const restrict seq, size_t len,
        const parasail_matrix_t *matrix)
{
    size_t i = 0;
    int score = 0;

    for (i = 0; i < len; i++) {
        int offset = matrix->mapper[seq[i]];
        score += matrix->matrix[offset*matrix->size + offset];
    }

    return score;
}


#define IS_EDGE_ASSERT  \
    assert(s1);         \
    assert(s2);         \
    assert(s1_len);     \
    assert(s2_len);     \

#define IS_EDGE_BODY(SELF_SCORE1, SELF_SCORE2)                              \
    if (s1_len > s2_len) {                                                  \
        max_len = s1_len;                                                   \
        self_score_ = SELF_SCORE1;                                          \
    }                                                                       \
    else {                                                                  \
        max_len = s2_len;                                                   \
        self_score_ = SELF_SCORE2;                                          \
    }                                                                       \
                                                                            \
    if (result->score <= 0) {                                               \
        return false;                                                       \
    }                                                                       \
    else if ((result->length * 100 >= AOL * int(max_len))                   \
             && (result->matches * 100 >= SIM * result->length)             \
             && (result->score * 100 >= OS * self_score_)) {                \
        return true;                                                        \
    }                                                                       \
    else {                                                                  \
        return false;                                                       \
    }


bool is_edge(
        const parasail_result_t *result,
        const char * const restrict s1, size_t s1_len,
        const char * const restrict s2, size_t s2_len,
        int AOL,
        int SIM,
        int OS,
        int &self_score_,
        size_t &max_len,
        const parasail_matrix_t *matrix)
{
    IS_EDGE_ASSERT
    assert(result);
    assert(matrix);

    IS_EDGE_BODY(self_score(s1, s1_len, matrix),
                 self_score(s2, s2_len, matrix))
}

}; /* namespace pgraph */

