/**
 * @file Sequence.h
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 */
#ifndef SEQUENCE_H_
#define SEQUENCE_H_

#include <iostream>
#include <string>

using std::ostream;
using std::string;


/**
 * A biological sequence.
 *
 * Similar to std::string, but immutable and the memory potentially may not be
 * owned by any instance but rather as part of a larger SequenceDatabase
 * instance (or some other data source).
 */
class Sequence
{
    public:
        /** default */
        Sequence();

        /** shallow copy; new instance does not own the data */
        Sequence(const Sequence &seq);

        /** from character sequence; new instance does not own the data */
        Sequence(const char *seq);

        /** from character sequence; new instance does not own the data */
        Sequence(const char *seq, size_t len);

        /** from character sequence; new instance does not own the data */
        Sequence(const char *id, size_t id_size,
                 const char *data, size_t data_size);

        /** delete data, if owner */
        virtual ~Sequence();

        /** retrieve ID for sequence */
        const char* get_id() const { return id; }

        /** retrieve ID size for sequence */
        size_t get_id_size() const { return id_size; }

        /** retrieve sequence data */
        const char* get_data() const { return data; }

        /** retrieve sequence data size */
        size_t get_data_size() const { return data_size; }

        /** affine gap align, returning score, #matches, and align length */
        void align(const Sequence &that, int &score, int &ndig, int &alen);

        /** overload of string cast */
        operator string() const;

        friend ostream &operator << (ostream &os, const Sequence &s);

    private:
        bool is_owner;      /**< whether this instance should delete the data */
        const char *id;     /**< sequence identifier */
        size_t id_size;     /**< length of sequence identifier */
        const char *data;   /**< character sequence */
        size_t data_size;   /**< length of character sequence */
};


#endif /* SEQUENCE_H_ */
