/**
 * @file param.h
 *
 * @author andy.cj.wu@gmail.com
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2010 Washington State University. All rights reserved.
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 */
#ifndef _PARAM_H_
#define _PARAM_H_

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


int pg_get_param_int(const char *param_file, const char *key);

void pg_get_params(const char *param_file, param_t *param);

#endif /* _PARAM_H_ */
