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
        virtual const char* what() const throw();
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
         * @param filename the file to open (fasta format)
         * @param budget the memory budget, 0 for unlimited
         */
        explicit SequenceDatabase(const string &filename, size_t budget);

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

    private:
        void read_and_parse_fasta();
        void index_file(vector<MPI_Offset> &index);

        string filename;
        size_t budget;
        vector<char> local_cache;
        map<size_t,Sequence> sequences;
};

#endif /* SEQUENCE_DATABASE_H_ */
