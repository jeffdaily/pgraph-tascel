/**
 * @file loadseq.c
 *
 * @author andy.cj.wu@gmail.com
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2010 Washington State University. All rights reserved.
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 */
#include <assert.h>
#include <ctype.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "elib.h"
#include "loadseq.h"
#include "wlib.h"

#define MAX_LINE_LEN 16384

void load_all_sequences(
        const char *fileName,
        size_t sequence_count, sequence_t *sequences,
        size_t *_n_chars, size_t *_max_seq_len)
{
    FILE *fp = NULL;
    char *line = NULL;
    size_t line_max = 0;
    size_t line_len = 0;
    size_t i = 0; 
    size_t max_seq_len = 0;
    size_t n_chars = 0;

    if (LINE_MAX >= MAX_LINE_LEN) {
        /* Use maximum line size of MAX_LINE_LEN.  If LINE_MAX is bigger than
         * our limit, sysconf() can't report a smaller limit. */
        line_max = MAX_LINE_LEN;
    } else {
        long limit = sysconf(_SC_LINE_MAX);
        line_max = (limit < 0 || limit > MAX_LINE_LEN) ?
            MAX_LINE_LEN : (int)limit;
    }
    /* line_max + 1 leaves room for nul byte added by fgets */
    line = emalloc(line_max+1);

    fp = efopen(fileName, "r"); 

    while (fgets(line, line_max+1, fp)) {
        line_len = strlen(line); 
        assert(line_len <= line_max);
        /* add an '$' end to each seq */
        line[line_len-1] = '\0';

        if (line[0] == FASTA_FLAG) {
            sequences[i].gid = estrdup(line);
        }
        else if (isalpha(line[0])) {
            line[line_len-1] = DOLLAR;
            line[line_len] = '\0';
            /* strlen does not include the '\0' */
            sequences[i].strLen = line_len;
            n_chars += line_len;
            if (line_len > max_seq_len) {
                max_seq_len = line_len;
            }
            sequences[i].str = estrdup(line); 
            ++i;
        }
        else {
            warn("empty line in fasta file? will continue!!"); 
        }
    }

    fclose(fp);
    assert(i == sequence_count);

    /* return */
    *_n_chars = n_chars;
    *_max_seq_len = max_seq_len;
}


void print_sequence(sequence_t *sequences, size_t index)
{
    printf("StrLen=%d\n", sequences[index].strLen);
    printf("--------------------\n");
    printf("%s\n%s\n", sequences[index].gid, sequences[index].str);
}


void free_sequences(sequence_t *sequences, size_t size)
{
    size_t i;
    for(i=0; i<size; ++i){
        free(sequences[i].gid);
        free(sequences[i].str);
    } 
}
