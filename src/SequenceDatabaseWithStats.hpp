/**
 * @file SequenceDatabaseWithStats.hpp
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 */
#ifndef SEQUENCE_DATABASE_WITH_STATS_H_
#define SEQUENCE_DATABASE_WITH_STATS_H_

#include <cstddef>

#include "DbStats.hpp"
#include "SequenceDatabase.hpp"
#include "Sequence.hpp"

using std::size_t;

namespace pgraph {


/**
 * Collection of ordered Sequence instances, indexed from 0 to N.
 */
class SequenceDatabaseWithStats : public SequenceDatabase
{
    public:
        /**
         * Creates a sequence database from the given file with the given
         * memory budget.
         */
        SequenceDatabaseWithStats(SequenceDatabase *db)
            : SequenceDatabase(db->get_delimiter())
            , stats()
            , db(db)
        {}

        /**
         * Destroys the SequenceDatabaseWithStats.
         */
        virtual ~SequenceDatabaseWithStats() {}

        /**
         * Returns how many sequences are stored collectively by every process.
         *
         * @return how many sequences are stored collectively by every process.
         */
        virtual size_t size() const { return db->size(); }

        /**
         * Returns how many characters are stored collectively by every process.
         *
         * @return how many characters are stored collectively by every process.
         */
        virtual size_t char_size() const { return db->char_size(); }

        /**
         * Gets the size of the longest sequence in this database.
         *
         * @return the size of the longest sequence in this database.
         */
        virtual size_t longest() const { return db->longest(); }

        /**
         * Returns a pointer to the Sequence based on the global index i.
         *
         * @param[in] i the index based on the global count of sequences
         * @return the Sequence (possibly owned by this SequenceDatabase)
         */
        virtual Sequence* get_sequence(size_t i) {
            Sequence *sequence = NULL;
            double time;
            if (!db->is_local(i)) {
                time = MPI_Wtime();
            }
            sequence = db->get_sequence(i);
            if (!db->is_local(i)) {
                stats.time.push_back(MPI_Wtime()-time);
                stats.bytes.push_back(sequence->get_id_length()+sequence->size());
            }
            return sequence;
        }

        /**
         * Returns length of the given Sequence's data block.
         *
         * @returns length of the given Sequence's data block.
         */
        virtual size_t get_sequence_size(size_t i) {
            return db->get_sequence_size(i);
        }

        /**
         * Returns the delimiter appended to each Sequence, if not the null
         * byte.
         */
        virtual char get_delimiter() const { return db->get_delimiter(); }

        DbStats stats;

    protected:
        SequenceDatabase *db;

};

}; /* namespace pgraph */

#endif /* SEQUENCE_DATABASE_WITH_STATS_H_ */
