/**
 * @file dynamic.h
 *
 * @author andy.cj.wu@gmail.com
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2010 Washington State University. All rights reserved.
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 *
 * @todo TODO remove SIGMA and NROW; dynamically assign alphabet
 */
#ifndef SEQUENCE_H_
#define SEQUENCE_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <stddef.h>


/**
 * Sequence info for fasta file.
 */
typedef struct {
    const char *gid; /**< gid of fasta file */
    const char *str; /**< actual string of fasta file */
    size_t size;     /**< string length */
} sequence_t;


/**
 * Array of sequence_t instances and associated metadata.
 */
typedef struct {
    sequence_t *seq;     /**< array of sequence_t instances */
    size_t size;         /**< number of sequence_t instances in this array */
    char *data;          /**< blob of fasta data */
    size_t n_chars;      /** number of characters in the fasta file */
    size_t max_seq_size; /**< longest sequence parsed from fasta file */
} sequences_t;


/**
 * Opens and parses the given fasta file.
 *
 * @param[in] file_name the fasta file
 * @param[in] delimiter character to indicate end of a sequence e.g. '$'
 * @return array of sequences and metadata
 */
sequences_t* pg_load_fasta(const char *file_name, char delimiter);


/**
 * Parses the given fasta buffer.
 *
 * @param[in] file_name the fasta file
 * @param[in] delimiter character to indicate end of a sequence e.g. '$'
 * @return array of sequences and metadata
 */
sequences_t* pg_parse_fasta(char *buffer, size_t size, char delimiter);


/**
 * Free the memory associated with the given sequences_t instance.
 *
 * @param[in] sequences to free
 */
void pg_free_sequences(sequences_t *sequences);


/**
 * Prints the given sequence at index to stdout.
 *
 * Useful for debugging purposes.
 *
 * @param[in] sequences array of sequences
 * @param[in] index the index within the array to print
 */
void pg_print_sequence(sequence_t *sequences, size_t index);

#ifdef __cplusplus
}
#endif

#endif /* SEQUENCE_H_ */
