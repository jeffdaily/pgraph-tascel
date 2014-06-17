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

using ::std::ostream;
using ::std::size_t;
using ::std::string;

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
        virtual ~Sequence() {
            if (is_owner) {
                delete [] buffer;
            }
            buffer = NULL;
        }

        /**
         * Retrieves ID and length of ID -- do not delete.
         *
         * @param[out] id the id as a character buffer (not null terminated)
         * @param[out] length of the id character buffer
         */
        void get_id(const char * &id, size_t &length) const {
            id = &buffer[id_offset];
            length = id_length;
        }

        /**
         * Retrieves copy of ID as a string.
         *
         * @param[out] id the id as a string
         */
        void get_id(string &id) const {
            id.assign(&buffer[id_offset], id_length);
        }

        /**
         * Retrieves length of sequence ID.
         *
         * @return length of the sequence ID character buffer
         */
        size_t get_id_length() const {
            return id_length;
        }

        /**
         * Retrieves length of sequence.
         *
         * @return length of the sequence character buffer
         */
        size_t get_sequence_length() const {
            return sequence_length;
        }

        /**
         * Retrieves length of sequence.
         *
         * @return length of the sequence character buffer
         */
        size_t size() const {
            return sequence_length;
        }

        /**
         * Retrieves pointer to sequence data -- do not delete.
         *
         * As with string::data(), the buffer is not guaranteed to be
         * NULL-terminated.
         *
         * @return length of the sequence character buffer
         */
        const char * data() const {
            return &buffer[sequence_offset];
        }

        /**
         * Retrieves sequence and length of sequence -- do not delete.
         *
         * @param[out] sequence the sequence as a character buffer
         *             (not null * terminated)
         * @param[out] length of the sequence character buffer
         */
        void get_sequence(const char * &sequence, size_t &length) const {
            sequence = &buffer[sequence_offset];
            length = sequence_length;
        }

        /**
         * Retrieves copy of sequence as string.
         *
         * @param[out] sequence the sequence as a string
         */
        void get_sequence(string &sequence) const {
            sequence.assign(&buffer[sequence_offset], sequence_length);
        }

        /**
         * Retrieves data buffer for both ID and sequence and its size.
         *
         * @param[out] data the buffer backing this Sequence
         * @param[out] length the size of the data buffer
         */
        void get_buffer(const char * &data, size_t &length) const {
            data = this->buffer;
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

        /**
         * Retrieves copy of data buffer as a string.
         *
         * @param[out] data copy of the buffer backing this Sequence
         */
        void get_buffer(string &data) const {
            size_t length;
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
            data.assign(this->buffer, length);
        }

        const char & operator[](size_t i) const {
            assert(i < sequence_length);
            return buffer[sequence_offset+i];
        }

        void uses_delimiter(bool value) { has_delimiter = value; }
        bool uses_delimiter() const { return has_delimiter; }

        /** overload of string cast */
        operator string() const {
            assert(NULL != buffer);
            assert(sequence_length > 0);
            return string(&buffer[sequence_offset], sequence_length);
        }

        friend ostream &operator << (ostream &os, const Sequence &s);

    private:
        bool is_owner;  /**< whether this instance should delete the data */
        bool has_delimiter; /**< whether last char is a special delimiter */
        const char *buffer;   /**< all data, possibly containing ID */
        size_t id_offset;   /**< offset into data for start of ID */
        size_t id_length;   /**< length of ID */
        size_t sequence_offset; /**< offset into data for start of sequence */
        size_t sequence_length; /**< length of sequence */
};


/* for easy use with functions in alignment.hpp */

cell_t align_global_affine(
        const Sequence &s1,
        const Sequence &s2,
        match_t match,
        int open=-10, int gap=-1,
        cell_t * const restrict * const restrict tbl=NULL,
        int * const restrict * const restrict del=NULL,
        int * const restrict * const restrict ins=NULL);

cell_t align_global_affine(
        const Sequence &s1,
        const Sequence &s2,
        const int * const restrict * const restrict sub,
        const int * const restrict map, char first,
        int open=-10, int gap=-1,
        cell_t * const restrict * const restrict tbl=NULL,
        int * const restrict * const restrict del=NULL,
        int * const restrict * const restrict ins=NULL);

cell_t align_global_affine(
        const Sequence &s1,
        const Sequence &s2,
        int match, int mismatch,
        int open=-10, int gap=-1,
        cell_t * const restrict * const restrict tbl=NULL,
        int * const restrict * const restrict del=NULL,
        int * const restrict * const restrict ins=NULL);

cell_t align_global_affine(
        const Sequence &s1,
        const Sequence &s2,
        int open=-10, int gap=-1,
        cell_t * const restrict * const restrict tbl=NULL,
        int * const restrict * const restrict del=NULL,
        int * const restrict * const restrict ins=NULL);

cell_t align_global_affine_sse(
        const Sequence &s1,
        const Sequence &s2,
        int open, int gap,
        int * const restrict scr,
        int * const restrict del,
        int * const restrict mch,
        int * const restrict len);

cell_t align_semi_affine(
        const Sequence &s1,
        const Sequence &s2,
        match_t match,
        int open=-10, int gap=-1,
        cell_t * const restrict * const restrict tbl=NULL,
        int * const restrict * const restrict del=NULL,
        int * const restrict * const restrict ins=NULL);

cell_t align_semi_affine(
        const Sequence &s1,
        const Sequence &s2,
        const int * const restrict * const restrict sub,
        const int * const restrict map, char first,
        int open=-10, int gap=-1,
        cell_t * const restrict * const restrict tbl=NULL,
        int * const restrict * const restrict del=NULL,
        int * const restrict * const restrict ins=NULL);

cell_t align_semi_affine(
        const Sequence &s1,
        const Sequence &s2,
        int match, int mismatch,
        int open=-10, int gap=-1,
        cell_t * const restrict * const restrict tbl=NULL,
        int * const restrict * const restrict del=NULL,
        int * const restrict * const restrict ins=NULL);

cell_t align_semi_affine(
        const Sequence &s1,
        const Sequence &s2,
        int open=-10, int gap=-1,
        cell_t * const restrict * const restrict tbl=NULL,
        int * const restrict * const restrict del=NULL,
        int * const restrict * const restrict ins=NULL);

cell_t align_semi_affine_sse(
        const Sequence &s1,
        const Sequence &s2,
        int open, int gap,
        int * const restrict scr,
        int * const restrict del,
        int * const restrict mch,
        int * const restrict len);

cell_t align_local_affine(
        const Sequence &s1,
        const Sequence &s2,
        match_t match,
        int open=-10, int gap=-1,
        cell_t * const restrict * const restrict tbl=NULL,
        int * const restrict * const restrict del=NULL,
        int * const restrict * const restrict ins=NULL);

cell_t align_local_affine(
        const Sequence &s1,
        const Sequence &s2,
        const int * const restrict * const restrict sub,
        const int * const restrict map, char first,
        int open=-10, int gap=-1,
        cell_t * const restrict * const restrict tbl=NULL,
        int * const restrict * const restrict del=NULL,
        int * const restrict * const restrict ins=NULL);

cell_t align_local_affine(
        const Sequence &s1,
        const Sequence &s2,
        int match, int mismatch,
        int open=-10, int gap=-1,
        cell_t * const restrict * const restrict tbl=NULL,
        int * const restrict * const restrict del=NULL,
        int * const restrict * const restrict ins=NULL);

cell_t align_local_affine(
        const Sequence &s1,
        const Sequence &s2,
        int open=-10, int gap=-1,
        cell_t * const restrict * const restrict tbl=NULL,
        int * const restrict * const restrict del=NULL,
        int * const restrict * const restrict ins=NULL);

cell_t align_local_affine_sse(
        const Sequence &s1,
        const Sequence &s2,
        int open, int gap,
        int * const restrict scr,
        int * const restrict del,
        int * const restrict mch,
        int * const restrict len);

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

bool is_edge(
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
