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

#include "Sequence.h"

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
                                  const char *function, const char *message) throw();
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
         */
        explicit SequenceDatabase(const string &filename, size_t budget);

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

    private:
        void read_and_parse_fasta();
        void index_file(vector<MPI_Offset> &index);

        size_t budget;
        string file_name;
        MPI_Offset file_size;
        vector<MPI_Offset> global_file_index;
        vector<int> global_owner;
        vector<char *> global_address;
        char *local_cache;
        MPI_Offset local_cache_size;
        vector<char *> global_cache;
        map<size_t, Sequence> sequences;
        map<size_t, char *> remote_cache;
};

#endif /* SEQUENCE_DATABASE_H_ */
