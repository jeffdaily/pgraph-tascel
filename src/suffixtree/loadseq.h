/**
 * @file loadseq.h
 *
 * @author andy.cj.wu@gmail.com
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2010 Washington State University. All rights reserved.
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 */
#ifndef LOAD_SEQ_H_
#define LOAD_SEQ_H_

#include <stddef.h>

#define MAX_LINE_LEN 100000
#define FASTA_FLAG '>'
#define DOLLAR 'U' /* dollar delimiter */
#define BEGIN 'O'  /* it is oh, not zero */
#define SIGMA 26U  /* size of alphabet, '*'->'J' */


/**
 * Sequence info for fasta file.
 */
typedef struct {
    char *gid;  /**< gid of fasta file */
    char *str;  /**< actual string of fasta file */
    int strLen; /**< string length */
} sequence_t;


/**
 * Opens and parses the given fasta file.
 *
 * Assumes that the given sequence_t array 'sequences' has been preallocated
 * with enough space to hold all of sequences.
 *
 * @param[in] file_name the fasta file
 * @param[in] sequence_count number of sequences in the file; asserted
 * @param[in] sequences preallocated array of sequences to hold parsed sequences
 * @param[out] n_chars number of characters in the fasta file
 * @param[out] max_seq_len longest sequence parsed from fasta file
 */
void load_all_sequences(const char *file_name,
        size_t sequence_count, sequence_t *sequences,
        size_t *n_chars, size_t *max_seq_len);


/**
 * Prints the given sequence at index to stdout.
 *
 * Useful for debugging purposes.
 *
 * @param[in] sequences array of sequences
 * @param[in] index the index within the array to print
 */
void print_sequence(sequence_t *sequences, size_t index);


/**
 * Frees the given array of sequences.
 *
 * The sequence_t has nested allocations which are also free'd.
 *
 * @param[in] sequences array of sequences
 * @param[in] size the number of sequences in the array
 */
void free_sequences(sequence_t *sequences, size_t size);

#endif /* LOAD_SEQ_H_ */
