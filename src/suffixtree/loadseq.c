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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "elib.h"
#include "loadseq.h"
#include "wlib.h"


void load_all_sequences(
        const char *fileName,
        size_t sequence_count, sequence_t *sequences,
        size_t *_n_chars, size_t *_max_seq_len)
{
    FILE *fp = NULL;
    char line[MAX_LINE_LEN];
    size_t lineLen;
    size_t i = 0; 
    size_t max_seq_len = 0;
    size_t n_chars = 0;

    fp = efopen(fileName, "r"); 

    while (fgets(line, MAX_LINE_LEN, fp)) {
        lineLen = strlen(line); 
        assert(lineLen <= MAX_LINE_LEN);
        /* add an '$' end to each seq */
        line[lineLen-1] = '\0';

        if (line[0] == FASTA_FLAG) {
            sequences[i].gid = estrdup(line);
        }
        else if (isalpha(line[0])) {
            line[lineLen-1] = DOLLAR;
            line[lineLen] = '\0';
            /* strlen does not include the '\0' */
            sequences[i].strLen = lineLen;
            n_chars += lineLen;
            if (lineLen > max_seq_len) {
                max_seq_len = lineLen;
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
