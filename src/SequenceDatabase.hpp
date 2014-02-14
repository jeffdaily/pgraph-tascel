/**
 * @file SequenceDatabase.hpp
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 */
#ifndef SEQUENCE_DATABASE_H_
#define SEQUENCE_DATABASE_H_

#include <cstddef>

#include "Sequence.hpp"

using std::size_t;

namespace pgraph {


/**
 * Collection of ordered Sequence instances, indexed from 0 to N.
 */
class SequenceDatabase
{
    public:
        /**
         * Creates a sequence database from the given file with the given
         * memory budget.
         */
        SequenceDatabase() {}

        /**
         * Destroys the SequenceDatabase.
         */
        virtual ~SequenceDatabase() {}

        /**
         * Returns how many sequences are stored collectively by every process.
         *
         * @return how many sequences are stored collectively by every process.
         */
        virtual size_t size() const = 0;

        /**
         * Returns how many characters are stored collectively by every process.
         *
         * @return how many characters are stored collectively by every process.
         */
        virtual size_t char_size() const = 0;

        /**
         * Gets the size of the longest sequence in this database.
         *
         * @return the size of the longest sequence in this database.
         */
        virtual size_t longest() const = 0;

        /**
         * Returns a reference to the Sequence based on the global index i.
         *
         * @param[in] i the index based on the global count of sequences
         * @return the Sequence reference (owned by this SequenceDatabase)
         */
        virtual Sequence &get_sequence(size_t i) = 0;

        /**
         * Returns a reference to the Sequence based on the global index i.
         *
         * @param[in] i the index based on the global count of sequences
         * @return the Sequence reference (owned by this SequenceDatabase)
         */
        virtual Sequence &operator[](size_t i) = 0;

};

}; /* namespace pgraph */

#endif /* SEQUENCE_DATABASE_H_ */
