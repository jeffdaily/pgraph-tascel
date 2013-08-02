/**
 * @file csequence.c
 *
 * @author andy.cj.wu@gmail.com
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2010 Washington State University. All rights reserved.
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 */
#include <errno.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "csequence.h"


sequences_t* pg_load_fasta(const char *file_name, char delimiter)
{
    char *buffer = NULL;
    FILE *in = NULL;
    long file_size = 0;

    /* preconditions */
    if (NULL == file_name) {
        fprintf(stderr, "pg_load_fasta: NULL == file_name");
        exit(EXIT_FAILURE);
    }

    /* open file */
    in = fopen(file_name, "r");
    if (NULL == in) {
        perror("pg_load_fasta: fopen");
        exit(EXIT_FAILURE);
    }

    /* determine file size */
    if (0 != fseek(in, 0L, SEEK_END)) {
        perror("pg_load_fasta: fseek");
        exit(EXIT_FAILURE);
    }
    file_size = ftell(in);

    /* reset the file */
    errno = 0;
    rewind(in);
    if (0 != errno) {
        perror("pg_load_fasta: rewind");
        exit(EXIT_FAILURE);
    }

    /* allocate file buffer */
    buffer = malloc(file_size+1);
    if (NULL == buffer) {
        perror("pg_load_fasta: malloc file buffer");
        exit(EXIT_FAILURE);
    }

    /* read entire file */
    if (1 != fread(buffer, file_size, 1, in)) {
        fprintf(stderr, "pg_load_fasta: fread fasta file");
        exit(EXIT_FAILURE);
    }
    buffer[file_size] = '\0';

    /* close file */
    if (0 != fclose(in)) {
        perror("pg_load_fasta: fclose");
        exit(EXIT_FAILURE);
    }
    in = NULL;

    return pg_parse_fasta(buffer, file_size, delimiter);
}


sequences_t* pg_parse_fasta(char *buffer, size_t file_size, char delimiter)
{
    sequences_t *retval = NULL;
    sequence_t *sequences = NULL;
    size_t n_sequences = 0;
    size_t sequence_count = 0;
    size_t i = 0;
    size_t r = 0;
    size_t w = 0;
    size_t n_chars = 0;
    size_t max_seq_size = 0;

    /* preconditions */
    if (NULL == buffer) {
        fprintf(stderr, "pg_parse_fasta: NULL == file_name");
        exit(EXIT_FAILURE);
    }

    /* scan the buffer to see how many sequences are contained;
     * a sequence is indicated by a line starting with the '>' character;
     * while we're at it, replace strategic newlines with null and
     * remove/compress the others */
    n_sequences = 0;
    r = 0;
    w = 0;
    while (r<file_size) {
        if (buffer[r] == '>') {
            ++n_sequences;
            r++; /* we don't store the '>' -- we use the space for delim */
            while (buffer[r] != '\n') {
                buffer[w++] = buffer[r++];
            }
            buffer[w++] = '\0';
            r++;
        }
        else {
            while (buffer[r] != '\n') {
                buffer[w++] = buffer[r++];
            }
            r++;
            if (delimiter != '\0') {
                /* We don't add the delimiter if it's null because this throws
                 * off the strlen during later parsing. Users can specify any
                 * other character just fine e.g. '$'. */
                buffer[w++] = delimiter;
            }
            if (buffer[r] == '>') {
                buffer[w++] = '\0';
            }
        }
    }
    while (w<file_size) {
        buffer[w++] = '\0';
    }

    /* allocate storage for the sequences */
    sequences = malloc(n_sequences * sizeof(sequence_t));
    if (NULL == sequences) {
        perror("DP_parse_fasta: calloc");
        exit(EXIT_FAILURE);
    }

    /* now rescan the buffer and index its contents */
    i = 0;
    for (sequence_count=0; sequence_count<n_sequences; ++sequence_count) {
        sequence_t *seq = &sequences[sequence_count];
        seq->gid = &buffer[i];
        i += strlen(&buffer[i])+1;
        seq->str = &buffer[i];
        seq->size = strlen(&buffer[i]);
        i += seq->size+1;
        n_chars += seq->size;
        if (seq->size > max_seq_size) {
            max_seq_size = seq->size;
        }
    }

    /* allocate return structure */
    retval = malloc(sizeof(sequences_t));
    if (NULL == retval) {
        perror("pg_parse_fasta: malloc");
        exit(EXIT_FAILURE);
    }

    /* return values */
    retval->seq = sequences;
    retval->size = sequence_count;
    retval->data = buffer;
    retval->n_chars = n_chars;
    retval->max_seq_size = max_seq_size;

    return retval;
}


void pg_print_sequence(sequence_t *sequences, size_t index)
{
    printf("StrLen=%zu\n", sequences[index].size);
    printf("--------------------\n");
    printf("%s\n%s\n", sequences[index].gid, sequences[index].str);
}


void pg_free_sequences(sequences_t *sequences)
{
    free(sequences->seq);
    free(sequences->data);
    free(sequences);
}
