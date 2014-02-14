/**
 * @file SequenceDatabaseArmci.hpp
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 */
#ifndef SEQUENCE_DATABASE_ARMCI_H_
#define SEQUENCE_DATABASE_ARMCI_H_

#include <cstddef>
#include <exception>
#include <map>
#include <string>
#include <vector>

/* ANL's armci does not extern "C" inside the header */
extern "C" {
#include <armci.h>
}

#include <tascel.h>

#include "Sequence.hpp"
#include "SequenceDatabase.hpp"

using std::exception;
using std::map;
using std::size_t;
using std::string;
using std::vector;
using tascel::PthreadMutex;

namespace pgraph {


/**
 * Collection of ordered Sequence instances, indexed from 0 to N.
 */
class SequenceDatabaseArmci : public SequenceDatabase
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
        SequenceDatabaseArmci(const string &filename,
                              size_t budget,
                              MPI_Comm comm,
                              char delimiter='$');

        /**
         * Destroys the SequenceDatabaseArmci.
         */
        virtual ~SequenceDatabaseArmci();

        /**
         * Returns how many sequences are stored collectively by every process.
         *
         * @return how many sequences are stored collectively by every process.
         */
        virtual size_t size() const;

        /**
         * Returns how many characters are stored collectively by every process.
         *
         * @return how many characters are stored collectively by every process.
         */
        virtual size_t char_size() const;

        /**
         * Gets the size of the largest sequence in this database.
         *
         * @return the size of the largest sequence in this database.
         */
        virtual size_t longest() const { return max_seq_size; }

        /**
         * Returns a reference to the Sequence based on the global index i.
         *
         * @param[in] i the index based on the global count of sequences
         * @return the Sequence reference (owned by this SequenceDatabaseArmci)
         */
        virtual Sequence &get_sequence(size_t i);

        /**
         * Returns a reference to the Sequence based on the global index i.
         *
         * @param[in] i the index based on the global count of sequences
         * @return the Sequence reference (owned by this SequenceDatabaseArmci)
         */
        virtual Sequence &operator[](size_t i);

    private:
        void read_and_parse_fasta();
        void read_and_parse_fasta_lomem(MPI_File in, MPI_Offset file_size);
        void read_and_parse_fasta_himem(MPI_File in, MPI_Offset file_size);

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
        void pack_and_index_fasta(char *buffer, size_t size, char delimiter,
                                  size_t index, size_t &new_size);

        /**
         * Exchanges metadata associated with each local cache.
         */
        void exchange_local_cache();

        MPI_Comm comm_orig; /**< original communicator */
        int comm_orig_rank; /**< original communicator rank */
        int comm_orig_size; /**< original communicator size */
        MPI_Comm comm;      /**< sub communicator for smaller domain */
        int comm_rank;      /**< communicator rank */
        int comm_size;      /**< communicator size */
        ARMCI_Group armci_group_orig;/**< ARMCI group corresponding to comm_orig */
        ARMCI_Group armci_group;/**< ARMCI group corresponding to comm_orig */
        bool is_replicated; /**< is budget sufficient to hold entire file? */
        size_t budget;      /**< max amount of memory to use */
        char delimiter;     /**< delimiter to separate sequences */
        string file_name;   /**< fasta file name */
        char *local_data;   /**< memory allocated for local sequences */
        map<size_t, Sequence*> local_cache; /**< TODO */
        size_t global_count;/**< total number of sequences */
        size_t global_size;/**< total number of characters */
        size_t local_size;/**< local number of characters */
        vector<int> owners; /**< mapping from seq id to rank owner */
        vector<const void*> addresses; /**< where each sequence data block begins */
        vector<size_t> sizes; /**< length of each sequence data block */
        char **ptr_arr;     /**< TODO */
        map<size_t, Sequence*> remote_cache; /**< TODO */
        PthreadMutex mutex; /**< controls access to remote cache */
        size_t max_seq_size;/**< biggest sequence */
};

}; /* namespace pgraph */

#endif /* SEQUENCE_DATABASE_ARMCI_H_ */
