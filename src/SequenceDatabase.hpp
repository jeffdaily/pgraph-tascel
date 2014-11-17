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
#include <map>
#include <set>

#include "Sequence.hpp"

using std::map;
using std::set;
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
        SequenceDatabase(char delimiter) : delimiter(delimiter) {}

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
         * Returns a pointer to the Sequence based on the global index i.
         *
         * @param[in] i the index based on the global count of sequences
         * @return the Sequence (possibly owned by this SequenceDatabase)
         */
        virtual Sequence* get_sequence(size_t i) = 0;

        /**
         * Returns a map of Sequence pointers.
         *
         * @param[in] container of (globally-based) indices of sequences
         * @return the map of ID to Sequence instances
         */
        virtual map<size_t,Sequence*> get_sequences(set<size_t> container) {
            map<size_t,Sequence*> retval;
            for (set<size_t>::const_iterator it=container.begin();
                    it!=container.end(); ++it) {
                retval[*it] = get_sequence(*it);
            }
            return retval;
        }

        /**
         * Returns a vector of Sequence pointers in the given range.
         *
         * @param[in] lo first index of sequence range, inclusive
         * @param[in] hi second index of sequence range, inclusive
         * @return the vector of Sequence instances
         */
        virtual vector<Sequence*> get_sequences(size_t lo, size_t hi) {
            assert(hi>=lo);
            vector<Sequence*> retval(hi-lo+1);
            for (size_t i=lo,r=0; i<=hi; ++i) {
                retval[r++] = get_sequence(i);
            }
            return retval;
        }

        /**
         * Returns length of the given Sequence's data block.
         *
         * @returns length of the given Sequence's data block.
         */
        virtual size_t get_sequence_size(size_t i) = 0;

        /**
         * Returns the delimiter appended to each Sequence, if not the null
         * byte.
         */
        virtual char get_delimiter() const { return delimiter; }

        /**
         * Returns whether the given sequence is local or remote.
         *
         * Used primarly to gather performance statistics.
         */
        virtual bool is_local(size_t i) { return true; }

    protected:
        char delimiter;

};

}; /* namespace pgraph */

#endif /* SEQUENCE_DATABASE_H_ */
