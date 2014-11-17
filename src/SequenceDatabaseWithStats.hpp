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
            bool local = (db->is_local(i));
            if (!local) {
                time = MPI_Wtime();
            }
            sequence = db->get_sequence(i);
            if (!local) {
                stats.time.push_back(MPI_Wtime()-time);
                stats.bytes.push_back(sequence->get_id_length()+sequence->size());
            }
            return sequence;
        }

        /**
         * Returns a map of Sequence pointers.
         *
         * @param[in] container of (globally-based) indices of sequences
         * @return the map of ID to Sequence instances
         */
        virtual map<size_t,Sequence*> get_sequences(set<size_t> container) {
            double time = MPI_Wtime();
            map<size_t,Sequence*> retval = db->get_sequences(container);
            size_t bytes = 0;
            for (set<size_t>::const_iterator it=container.begin();
                    it!=container.end(); ++it) {
                if (!db->is_local(*it)) {
                    Sequence* &sequence = retval[*it];
                    bytes += sequence->get_id_length()+sequence->size();
                }
            }
            stats.time.push_back(MPI_Wtime()-time);
            stats.bytes.push_back(bytes);
            return retval;
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
