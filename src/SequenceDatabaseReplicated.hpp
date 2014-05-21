/**
 * @file SequenceDatabaseReplicated.hpp
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 */
#ifndef SEQUENCE_DATABASE_REPLICATED_H_
#define SEQUENCE_DATABASE_REPLICATED_H_

#include <cstddef>
#include <string>
#include <vector>

#include <mpi.h>

#include "Sequence.hpp"
#include "SequenceDatabase.hpp"

using ::std::size_t;
using ::std::string;
using ::std::vector;

namespace pgraph {


/**
 * Collection of ordered Sequence instances, indexed from 0 to N.
 */
class SequenceDatabaseReplicated : public SequenceDatabase
{
    public:

        /**
         * Creates a sequence database from the given file with the given
         * memory budget.
         *
         * @pre !filename.empty()
         * @param[in] filename the file to open (fasta format)
         * @param[in] budget the memory budget, 0 for unlimited
         * @param[in] comm MPI communicator
         * @param[in] num_threads number of threads which migth access this DB
         * @param[in] delimiter used to separate sequences
         */
        SequenceDatabaseReplicated(const string &filename,
                              size_t budget,
                              MPI_Comm comm,
                              char delimiter='\0');

        /**
         * Destroys the SequenceDatabaseReplicated.
         */
        virtual ~SequenceDatabaseReplicated();

        /**
         * Returns how many sequences are stored collectively by every process.
         *
         * @return how many sequences are stored collectively by every process.
         */
        virtual size_t size() const { return local_cache.size(); }

        /**
         * Returns how many characters are stored collectively by every process.
         *
         * @return how many characters are stored collectively by every process.
         */
        virtual size_t char_size() const { return _char_size; }

        /**
         * Gets the size of the largest sequence in this database.
         *
         * @return the size of the largest sequence in this database.
         */
        virtual size_t longest() const { return _longest; }

        /**
         * Returns a reference to the Sequence based on the global index i.
         *
         * @param[in] i the index based on the global count of sequences
         * @return the Sequence reference (owned by this SequenceDatabaseReplicated)
         */
        virtual Sequence &get_sequence(size_t i) {
            if (i > local_cache.size()) {
                std::cout << "get_sequence(" << i << ") failed" << std::endl;
            }
            assert(i <= local_cache.size());
            return *local_cache[i];
        }

        /**
         * Returns a reference to the Sequence based on the global index i.
         *
         * @param[in] i the index based on the global count of sequences
         * @return the Sequence reference (owned by this SequenceDatabaseReplicated)
         */
        virtual Sequence &operator[](size_t i) { return get_sequence(i); }

    private:
        /**
         * Preprocesses a fasta file buffer.
         *
         * We remove all newlines unless it is the newline that terminates the
         * sequence ID. We replace the ID-terminating newline with '#'. We replace
         * the sequence terminating newline with the given delimiter.
         *
         * @param[in] buffer the fasta file buffer
         * @param[in] size of the fasta file buffer
         * @param[in] delimiter to append to each sequence (can be '\0')
         * @param[in] index of first sequence
         * @param[out] new_size of the packed buffer
         */
        void pack_and_index_fasta(char *buffer, size_t size,
                                  size_t index, size_t &new_size);

        MPI_Comm comm;      /**< sub communicator for smaller domain */
        int comm_rank;      /**< communicator rank */
        int comm_size;      /**< communicator size */
        string file_name;   /**< fasta file name */
        char *local_data;   /**< memory allocated for local sequences */
        vector<Sequence*> local_cache; /**< TODO */
        size_t _longest;    /**< longest sequence */
        size_t _char_size;  /**< total character count */
};

}; /* namespace pgraph */

#endif /* SEQUENCE_DATABASE_REPLICATED_H_ */
