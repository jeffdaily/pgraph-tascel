/**
 * @file SequenceDatabaseGArray.hpp
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 */
#ifndef SEQUENCE_DATABASE_GARRAY_H_
#define SEQUENCE_DATABASE_GARRAY_H_

#include <cstddef>
#include <exception>
#include <map>
#include <string>

#include <tascel.h>

#include "Sequence.hpp"
#include "SequenceDatabase.hpp"

using std::exception;
using std::map;
using std::size_t;
using std::string;
using tascel::PthreadMutex;

namespace pgraph {


/**
 * Collection of ordered Sequence instances, indexed from 0 to N.
 */
class SequenceDatabaseGArray : public SequenceDatabase
{
    public:
        /**
         * Creates a sequence database from the given file with the given
         * memory budget.
         */
        SequenceDatabaseGArray(const string &file_name, size_t budget, char delimiter);

        /**
         * Destroys the SequenceDatabaseGArray.
         */
        virtual ~SequenceDatabaseGArray();

        /**
         * Returns how many sequences are stored on this process.
         *
         * @return how many sequences are stored on this process.
         */
        virtual size_t get_local_count() const;

        /**
         * Returns how many characters are stored on this process.
         *
         * @return how many characters are stored on this process.
         */
        virtual size_t get_local_size() const;

        /**
         * Returns how many sequences are stored collectively by every process.
         *
         * @return how many sequences are stored collectively by every process.
         */
        virtual size_t get_global_count() const;

        /**
         * Returns how many complete sets of sequences are stored collectively by every process.
         *
         * @return how many complete sets of sequences are stored collectively by every process.
         */
        virtual size_t get_global_replica_count() const;

        /**
         * Returns the index of this replica, 0-based.
         *
         * @return the index of this replica, 0-based.
         */
        virtual size_t get_global_replica_index() const;

        /**
         * Returns how many characters are stored collectively by every process.
         *
         * @return how many characters are stored collectively by every process.
         */
        virtual size_t get_global_size() const;

        /**
         * Returns a reference to the Sequence based on the global index i.
         *
         * @param[in] i the index based on the global count of sequences
         * @return the Sequence reference (owned by this SequenceDatabaseGArray)
         */
        virtual Sequence &get_sequence(size_t i);

        /**
         * Returns a reference to the Sequence based on the global index i.
         *
         * @param[in] i the index based on the global count of sequences
         * @return the Sequence reference (owned by this SequenceDatabaseGArray)
         */
        virtual Sequence &operator[](size_t i);

        /**
         * Sets the number of threads that might possibly access this db.
         *
         * @param[in] num the number of threads
         */
        virtual void set_num_threads(size_t num);

        /**
         * Gets the size of the largest sequence in this database.
         *
         * @return the size of the largest sequence in this database.
         */
        virtual size_t get_max_length() const;

        /**
         * Computes the alignment between the two sequence IDs.
         *
         * @param[in] i first Sequence index
         * @param[in] j second Sequence index
         * @param[out] score TODO
         * @param[out] ndig TODO
         * @param[out] align TODO
         * @param[in] open penalty
         * @param[in] gap penalty
         * @param[in] tid thread ID, defaulting to 0
         */
        virtual void align(size_t i,
                           size_t j,
                           int &score,
                           int &ndig,
                           int &align,
                           int open=-10,
                           int gap=-1,
                           int tid=0);

        /**
         * Computes the alignment between the two sequence IDs.
         *
         * @param[in] i first Sequence index
         * @param[in] j second Sequence index
         * @param[out] score TODO
         * @param[out] ndig TODO
         * @param[out] align TODO
         * @param[in] open penalty
         * @param[in] gap penalty
         * @param[in] tid thread ID, defaulting to 0
         */
        virtual void align_ssw(size_t i,
                           size_t j,
                           int &score,
                           int &ndig,
                           int &align,
                           int open=-10,
                           int gap=-1,
                           int tid=0);

        /**
         * Computes whether an edge exists between the two sequence IDs.
         *
         * @param[in] i first Sequence index
         * @param[in] j second Sequence index
         * @param[in] score TODO
         * @param[in] ndig TODO
         * @param[in] align TODO
         * @param[in] AOL TODO
         * @param[in] SIM TODO
         * @param[in] OS TODO
         * @param[out] sscore self score
         * @param[out] max_len longest of the two sequences 
         * @return the answer
         */
        virtual bool is_edge(size_t i,
                             size_t j,
                             const int &score,
                             const int &ndig,
                             const int &align,
                             const int &AOL,
                             const int &SIM,
                             const int &OS,
                             int &sscore,
                             size_t &max_len);


    protected:
        void read_and_parse_fasta(const string &file_name);
        void read_and_parse_fasta_himem(const string &file_name);
        void read_and_parse_fasta_lomem(const string &file_name);
        void pack_and_index_fasta(
                char *buffer,
                size_t size,
                char delimiter,
                size_t id,
                size_t &new_size);

    private:
        MPI_Comm comm;      /**< communicator */
        int comm_size;      /**< communicator size */
        int comm_rank;      /**< communicator rank */
        MPI_Offset file_size;/**< fasta file size */
        MPI_Offset read_size;/**< how much this MPI rank reads */
        bool is_replicated; /**< is budget sufficient to hold entire file? */
        size_t budget;      /**< max amount of memory to use */
        char delimiter;     /**< delimiter to separate sequences */
        char *local_data;   /**< memory allocated for local sequences */
        map<size_t, Sequence*> local_cache; /**< TODO */
        tascel::GArray<1> *ga;/**< global array holding fasta file */
        vector<long> counts;/**< number of '>' on each MPI rank */
        long count_total;   /**< total number of sequences */
        vector<long> offsets;/**< index into global array for sequences */
        map<size_t, Sequence*> remote_cache; /**< TODO */
        PthreadMutex mutex; /**< controls access to remote cache */
        long max_seq_size;  /**< biggest sequence */
        cell_t ***tbl;      /**< cell table for alignment */
        int ***del;         /**< del table for alignment */
        int ***ins;         /**< ins table for alignment */
};

}; /* namespace pgraph */

#endif /* SEQUENCE_DATABASE_GARRAY_H_ */
