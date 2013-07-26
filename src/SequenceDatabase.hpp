/**
 * @file SequenceDatabase.h
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 */
#ifndef SEQUENCE_DATABASE_H_
#define SEQUENCE_DATABASE_H_

#include <stdint.h>

#include <exception>
#include <map>
#include <string>
#include <vector>

#include "Sequence.hpp"

using std::exception;
using std::map;
using std::string;
using std::vector;


/**
 * Generic excpetion thrown by the SequenceDatabase.
 */
class SequenceDatabaseException : public exception
{
    public:
        SequenceDatabaseException(const char *file, int line,
                                  const char *function,
                                  const char *message) throw();
        virtual ~SequenceDatabaseException() throw();
        virtual const char *what() const throw();
    private:
        string message;
};

/* for brevity */
#define SEQUENCE_DATABASE_EXCEPTION(MSG) \
    throw SequenceDatabaseException(__FILE__,__LINE__,__func__,(MSG))

/**
 * Collection of ordered Sequence instances, indexed from 0 to N.
 */
class SequenceDatabase
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
         */
        explicit SequenceDatabase(const string &filename, size_t budget,
                                  MPI_Comm comm);

        /**
         * Destroys the SequenceDatabase.
         */
        ~SequenceDatabase();

        /**
         * Returns how many sequences are stored on this process.
         *
         * @return how many sequences are stored on this process.
         */
        size_t get_local_count() const;

        /**
         * Returns how many sequences are stored collectively by every process.
         *
         * @return how many sequences are stored collectively by every process.
         */
        size_t get_global_count() const;

        /**
         * Returns a reference to the Sequence based on the global index i.
         *
         * @param[in] i the index based on the global count of sequences
         * @return the Sequence reference (owned by this SequenceDatabase)
         */
        Sequence &get_sequence(size_t i);

        /**
         * Computes the alignment between the two sequence IDs.
         */
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

        MPI_Comm comm;      /**< communicator */
        int comm_rank;      /**< communicator rank */
        int comm_size;      /**< communicator size */
        bool is_replicated; /**< is budget sufficient to hold entire file? */
        size_t budget;      /**< max amount of memory to use */
        string file_name;   /**< fasta file name */
        char *local_data;   /**< memory allocated for local sequences */
        map<size_t, Sequence*> local_cache; /**< TODO */
        size_t global_count;/**< total number of sequences */
        vector<int> owners; /**< mapping from seq id to rank owner */
        vector<const void*> addresses; /**< where each sequence data block begins */
        vector<size_t> sizes; /**< length of each sequence data block */
        char **ptr_arr;     /**< TODO */
        map<size_t, Sequence*> remote_cache; /**< TODO */
};

#endif /* SEQUENCE_DATABASE_H_ */
