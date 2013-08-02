/**
 * @file Sequence.hpp
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 */
#ifndef _PGRAPH_SEQUENCE_H_
#define _PGRAPH_SEQUENCE_H_

#include <cassert>
#include <cstddef>
#include <cstring>
#include <iostream>
#include <string>

#include "alignment.hpp"

using std::ostream;
using std::size_t;
using std::string;

namespace pgraph {


/**
 * A biological sequence.
 *
 * Similar to std::string, but immutable and the memory potentially may not be
 * owned by any instance but rather as part of a larger SequenceDatabase
 * instance (or some other data source). If the data is owned, it will be
 * deleted when the Sequence instance is deleted (or goes out of scope).
 *
 * A Sequence is internally represented as a contiguous sequence of characters.
 * Part of those characters may represent a name for the sequence.
 */
class Sequence
{
    public:
        /** Default constructor. */
        Sequence();

        /**
         * Shallow copy constructor.
         *
         * The new instance does not own the data.
         *
         * @param[in] that other Sequence instance
         */
        Sequence(const Sequence &that);

        /**
         * Constructs from character sequence.
         *
         * By default, the new instance does not own the data. The data is
         * assumed to not have a name -- it is simply an unlabelled sequence.
         * It is an error for the data to contain whitespace.
         *
         * The data must be null terminated. 
         *
         * @param[in] data null-terminated character string
         * @param[in] owns whether this instance assumes ownership of data
         */
        Sequence(const char *data, const bool &owns=false);

        /**
         * Constructs from character sequence and length.
         *
         * By default, the new instance does not own the data. The data is
         * assumed to not have a name -- it is simply an unlabelled sequence.
         * It is an error for the data to contain whitespace.
         *
         * The data does not need to be null-terminated.
         *
         * @param[in] data possibly null-terminated character string
         * @param[in] length of data
         * @param[in] owns whether this instance assumes ownership of data
         */
        Sequence(const char *data,
                 const size_t &length,
                 const bool &owns=false);

        /**
         * Constructs from data array and both ID and sequence offsets.
         *
         * By default, the new instance does not own the data.
         * It is an error for the data to contain whitespace.
         * The offset and length for both the ID and the sequence must be
         * specified. The ID offset and length may both be 0 to indicate that
         * no ID is present in the data buffer.
         *
         * The data does not need to be null-terminated.
         *
         * @param[in] data possibly null-terminated character string
         * @param[in] id_offset starting index of ID
         * @param[in] id_length length of ID
         * @param[in] sequence_offset starting index of sequence
         * @param[in] sequence_length length of sequence
         * @param[in] owns whether this instance assumes ownership of data
         */
        Sequence(const char *data,
                 const size_t &id_offset,
                 const size_t &id_length,
                 const size_t &sequence_offset,
                 const size_t &sequence_length,
                 const bool &owns=false);

        /** Destructs instance and deletes data, if owner. */
        virtual ~Sequence();

        /**
         * Retrieves ID and length of ID.
         *
         * @param[out] id the id as a character buffer (not null terminated)
         * @param[out] length of the id character buffer
         */
        void get_id(const char * &id, size_t &length) {
            id = &data[id_offset];
            length = id_length;
        }

        /**
         * Retrieves length of sequence.
         *
         * @return length of the sequence character buffer
         */
        size_t get_sequence_length() {
            return sequence_length;
        }

        /**
         * Retrieves sequence and length of sequence.
         *
         * @param[out] sequence the sequence as a character buffer
         *             (not null * terminated)
         * @param[out] length of the sequence character buffer
         */
        void get_sequence(const char * &sequence, size_t &length) {
            sequence = &data[sequence_offset];
            length = sequence_length;
        }

        /**
         * Retrieves data buffer for both ID and sequence and its size.
         *
         * @param[out] data the buffer backing this Sequence
         * @param[out] length the size of the data buffer
         */
        void get_data(const char * &data, size_t &length) {
            data = this->data;
            if (id_length == 0) {
                /* no ID present in data */
                length = sequence_length;
            }
            else {
                /* ID present in data, probably at the front */
                if (id_offset < sequence_offset) {
                    length = sequence_offset + sequence_length;
                }
                else {
                    assert(id_offset != sequence_offset);
                    length = id_offset + id_length;
                }
            }
        }

        /** affine gap align, returning score, #matches, and align length */
        void align(const Sequence &that, int &score, int &ndig, int &alen);

        /** overload of string cast */
        operator string() const;

        friend ostream &operator << (ostream &os, const Sequence &s);

    private:
        bool is_owner;  /**< whether this instance should delete the data */
        const char *data;   /**< all data, possibly containing ID */
        size_t id_offset;   /**< offset into data for start of ID */
        size_t id_length;   /**< length of ID */
        size_t sequence_offset; /**< offset into data for start of sequence */
        size_t sequence_length; /**< length of sequence */
};


/* for easy use with functions in alignment.hpp */

/**
 * Implementation of affine gap pairwise sequence alignment.
 *
 * It is a space efficient version using only two rows. Also, memory for
 * all dynamic tables are allocated ONLY ONCE outside of this function call
 * using allocate_cell_table() and allocate_int_table() and passed as
 * tbl, del, and ins arguments.
 *
 * @param[in] s1 first Sequence
 * @param[in] s2 second Sequence
 * @param[in] tbl pre-allocated score table
 * @param[in] del pre-allocated deletion table
 * @param[in] ins pre-allocated insertion table
 * @param[in] open gap penalty
 * @param[in] gap extension penalty
 * @param[in] sub substitution matrix e.g. BLOSUM
 * @param[in] map mapping for substitution matrix
 * @param[in] first (zeroth) character of mapping alphabet e.g. 'A'
 * @return alignment result
 */
cell_t affine_gap_align(
        const Sequence &s1,
        const Sequence &s2,
        cell_t ** const restrict tbl,
        int ** const restrict del,
        int ** const restrict ins,
        int open, int gap,
        const int ** const restrict sub,
        const int * const restrict map, char first);


/**
 * Implementation of affine gap pairwise sequence alignment.
 *
 * It is a space efficient version using only two rows. Also, memory for
 * all dynamic tables are allocated ONLY ONCE outside of this function call
 * using allocate_cell_table() and allocate_int_table() and passed as
 * tbl, del, and ins arguments.
 *
 * @param[in] s1 first Sequence
 * @param[in] s2 second Sequence
 * @param[in] tbl pre-allocated score table
 * @param[in] del pre-allocated deletion table
 * @param[in] ins pre-allocated insertion table
 * @param[in] open gap penalty
 * @param[in] gap extension penalty
 * @param[in] match callback function to calculate any character match/mismatch
 * @return alignment result
 */
cell_t affine_gap_align(
        const Sequence &s1,
        const Sequence &s2,
        cell_t ** const restrict tbl,
        int ** const restrict del,
        int ** const restrict ins,
        int open, int gap,
        match_t match);


/**
 * Implementation of affine gap pairwise sequence alignment.
 *
 * It is a space efficient version using only two rows. Also, memory for
 * all dynamic tables are allocated ONLY ONCE outside of this function call
 * using allocate_cell_table() and allocate_int_table() and passed as
 * tbl, del, and ins arguments.
 *
 * @param[in] s1 first Sequence
 * @param[in] s2 second Sequence
 * @param[in] tbl pre-allocated score table
 * @param[in] del pre-allocated deletion table
 * @param[in] ins pre-allocated insertion table
 * @param[in] open gap penalty
 * @param[in] gap extension penalty
 * @param[in] match score
 * @param[in] mismatch score
 * @return alignment result
 */
cell_t affine_gap_align(
        const Sequence &s1,
        const Sequence &s2,
        cell_t ** const restrict tbl,
        int ** const restrict del,
        int ** const restrict ins,
        int open, int gap,
        int match, int mismatch);


/**
 * Implementation of affine gap pairwise sequence alignment using blosum.
 *
 * It is a space efficient version using only two rows. Also, memory for
 * all dynamic tables are allocated ONLY ONCE outside of this function call
 * using allocate_cell_table() and allocate_int_table() and passed as
 * tbl, del, and ins arguments.
 *
 * @param[in] s1 first Sequence
 * @param[in] s2 second Sequence
 * @param[in] tbl pre-allocated score table
 * @param[in] del pre-allocated deletion table
 * @param[in] ins pre-allocated insertion table
 * @param[in] open gap penalty
 * @param[in] gap extension penalty
 * @return alignment result
 */
cell_t affine_gap_align_blosum(
        const Sequence &s1,
        const Sequence &s2,
        cell_t ** const restrict tbl,
        int ** const restrict del,
        int ** const restrict ins,
        int open, int gap);


/**
 * Asks whether the given cell_t alignment result is an edge, based on
 * the given parameters.
 *
 * @param[in] result the dynamic programming alignment result
 * @param[in] s1 first Sequence
 * @param[in] s2 second Sequence
 * @param[in] AOL alignment over longer sequence (heuristic)
 * @param[in] SIM match similarity (heuristic)
 * @param[in] OS optimal score over self score (heuristic)
 * @param[out] self_score self score
 * @param[out] max_len longer of s1_len and s2_len
 * @param[in] sub substitution matrix
 * @param[in] map index mapping for substitution matrix
 * @param[in] first character offset (see doc above)
 * @return true if this is an edge, false otherwise
 */
bool is_edge(
        const cell_t &result,
        const Sequence &s1,
        const Sequence &s2,
        int AOL,
        int SIM,
        int OS,
        int &self_score,
        size_t &max_len,
        const int ** const restrict sub,
        const int * const restrict map, char first);


/**
 * Asks whether the given cell_t alignment result is an edge, based on
 * the given parameters.
 *
 * @param[in] result the dynamic programming alignment result
 * @param[in] s1 first Sequence
 * @param[in] s2 second Sequence
 * @param[in] AOL alignment over longer sequence (heuristic)
 * @param[in] SIM match similarity (heuristic)
 * @param[in] OS optimal score over self score (heuristic)
 * @param[out] self_score self score
 * @param[out] max_len longer of s1_len and s2_len
 * @param[in] match callback function to calculate any character match/mismatch
 * @return true if this is an edge, false otherwise
 */
bool is_edge(
        const cell_t &result,
        const Sequence &s1,
        const Sequence &s2,
        int AOL,
        int SIM,
        int OS,
        int &self_score,
        size_t &max_len,
        match_t match);


/**
 * Asks whether the given cell_t alignment result is an edge, based on
 * the given parameters.
 *
 * @param[in] result the dynamic programming alignment result
 * @param[in] s1 first Sequence
 * @param[in] s2 second Sequence
 * @param[in] AOL alignment over longer sequence (heuristic)
 * @param[in] SIM match similarity (heuristic)
 * @param[in] OS optimal score over self score (heuristic)
 * @param[out] self_score self score
 * @param[out] max_len longer of s1_len and s2_len
 * @param[in] match the static score of a character match
 * @return true if this is an edge, false otherwise
 */
bool is_edge(
        const cell_t &result,
        const Sequence &s1,
        const Sequence &s2,
        int AOL,
        int SIM,
        int OS,
        int &self_score,
        size_t &max_len,
        int match);


/**
 * Asks whether the given cell_t alignment result is an edge, based on
 * the given parameters.
 *
 * @param[in] result the dynamic programming alignment result
 * @param[in] s1 first Sequence
 * @param[in] s2 second Sequence
 * @param[in] AOL alignment over longer sequence (heuristic)
 * @param[in] SIM match similarity (heuristic)
 * @param[in] OS optimal score over self score (heuristic)
 * @param[out] self_score self score
 * @param[out] max_len longer of s1_len and s2_len
 * @return true if this is an edge, false otherwise
 */
bool is_edge_blosum(
        const cell_t &result,
        const Sequence &s1,
        const Sequence &s2,
        int AOL,
        int SIM,
        int OS,
        int &self_score,
        size_t &max_len);


}; /* namespace pgraph */

#endif /* _PGRAPH_SEQUENCE_H_ */
