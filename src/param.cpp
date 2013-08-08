/**
 * @file param.cpp
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2010 Washington State University. All rights reserved.
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 */
#include "config.h"

#include <cassert>

#include "param.hpp"
#include "mpix.hpp"

#define COMMENT '#'


namespace pgraph {

#define CFG_MAX_LINE_LEN 1024
int get_param_int(const char *param_file, const char *key)
{
    FILE *fp = NULL;
    char line[CFG_MAX_LINE_LEN];
    int val;
    char *p = NULL;

    fp = fopen(param_file, "r");
    if (NULL == fopen) {
        perror("get_param_int: fopen");
        exit(EXIT_FAILURE);
    }

    while (fgets(line, CFG_MAX_LINE_LEN, fp)) {
        /* comment line starts with '#' */
        if (strchr(line, COMMENT)) {
            continue;
        }

        /* empty line */
        if (line[0] == '\n') {
            continue;
        }

        /* param line */
        if (strstr(line, key)) {
            p = line;
            while (isspace(*p++));

            /* %*s used to skip the first item */
            sscanf(p, "%*s %d\n", &val);
            fclose(fp);
            return val;
        }
    }

    fclose(fp);
    printf("cannot find param value for key: [%s]\n", key);

    return -1;
}


void get_params(const char *param_file, param_t *param)
{
    assert(param);

    /* read in parameters */
    param->window_size = get_param_int(param_file, "SlideWindowSize");
    param->exact_match_len = get_param_int(param_file, "ExactMatchLen");
    param->AOL = get_param_int(param_file, "AlignOverLongerSeq");
    param->SIM = get_param_int(param_file, "MatchSimilarity");
    param->OS = get_param_int(param_file, "OptimalScoreOverSelfScore");
}

}; /* namespace pgraph */

