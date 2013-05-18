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
#include "readline.h"
#include "wlib.h"

#define MAX_LINE_LEN 16384

void load_all_sequences(
    const char *fileName,
    size_t sequence_count, sequence_t *sequences,
    size_t *_n_chars, size_t *_max_seq_len)
{
    struct line_reader lr;
    FILE *fp = NULL;
    size_t len = 0;
    char *line = NULL;
    size_t i = 0;
    size_t max_seq_len = 0;
    size_t n_chars = 0;

    fp = efopen(fileName, "r");

    lr_init(&lr, fp);
    while (line = next_line(&lr, &len)) {
        line[len - 1] = '\0'; /* remove newline */
        if (line[0] == FASTA_FLAG) {
            sequences[i].gid = estrdup(line);
        }
        else if (isalpha(line[0])) {
            line[len - 1] = DOLLAR; /* use DOLLAR instead of newline */
            /* strlen does not include the '\0' */
            sequences[i].strLen = len;
            n_chars += len;
            if (len > max_seq_len) {
                max_seq_len = len;
            }
            sequences[i].str = estrdup(line);
            ++i;
        }
        else {
            warn("empty line in fasta file? will continue!!");
        }
    }
    if (!feof(fp)) {
        perror("next_line");
        exit(EXIT_FAILURE);
    }
    lr_free(&lr);

    fclose(fp);
    if (i != sequence_count) {
        printf("i=%zu != sequence_count=%zu\n", i, sequence_count);
    }
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
    for (i = 0; i < size; ++i) {
        free(sequences[i].gid);
        free(sequences[i].str);
    }
}
