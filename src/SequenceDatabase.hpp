/**
 * @file SequenceDatabase.h
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 */
#ifndef SEQUENCE_DATABASE_H_
#define SEQUENCE_DATABASE_H_

#include <cstddef>
#include <exception>
#include <string>

#include "Sequence.hpp"

using std::exception;
using std::size_t;
using std::string;


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
        SequenceDatabase();

        /**
         * Destroys the SequenceDatabase.
         */
        virtual ~SequenceDatabase();

        /**
         * Returns how many sequences are stored on this process.
         *
         * @return how many sequences are stored on this process.
         */
        virtual size_t get_local_count() const = 0;

        /**
         * Returns how many sequences are stored collectively by every process.
         *
         * @return how many sequences are stored collectively by every process.
         */
        virtual size_t get_global_count() const = 0;

        /**
         * Returns a reference to the Sequence based on the global index i.
         *
         * @param[in] i the index based on the global count of sequences
         * @return the Sequence reference (owned by this SequenceDatabase)
         */
        virtual Sequence &get_sequence(size_t i) = 0;

        /**
         * Sets the number of threads that might possibly access this db.
         *
         * @param[in] num the number of threads
         */
        virtual void set_num_threads(size_t num) = 0;

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
                           int tid=0) = 0;

};

#endif /* SEQUENCE_DATABASE_H_ */
