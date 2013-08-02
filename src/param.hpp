/**
 * @file param.hpp
 *
 * @author andy.cj.wu@gmail.com
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2010 Washington State University. All rights reserved.
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 */
#ifndef _PGRAPH_PARAM_H_
#define _PGRAPH_PARAM_H_

namespace pgraph {

/**
 * parameters for is_edge test, packed in struct
 */
typedef struct {
    int AOL;            /**< AlignOverLongerSeq */
    int SIM;            /**< MatchSimilarity */
    int OS;             /**< OptimalScoreOverSelfScore */
    int exact_match_len;/**< exact match length cutoff */
    int window_size;    /**< slide window size */
    int open;           /**< open penalty for affine gap alignment */
    int gap;            /**< gap extension penalty for affine gap alignment */
} param_t;


int get_param_int(const char *param_file, const char *key);

void get_params(const char *param_file, param_t *param);

}; /* namespace pgraph */

#endif /* _PGRAPH_PARAM_H_ */
