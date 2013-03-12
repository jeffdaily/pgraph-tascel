#include <mpi.h>

#include <string>

#include "SequenceDatabase.h"
#include "mpix.h"

using std::string;


SequenceDatabaseException::SequenceDatabaseException(
        const char *file,
        int line,
        const char *function,
        const char *message) throw()
{

}


SequenceDatabaseException::~SequenceDatabaseException() throw()
{

}


const char* SequenceDatabaseException::what() const throw()
{

}



SequenceDatabase::SequenceDatabase(const string &filename, size_t budget)
    :   filename(filename)
    ,   budget(budget)
{

}


size_t SequenceDatabase::get_local_count() const
{
    /** @todo TODO */
    return 0;
}


size_t SequenceDatabase::get_global_count() const
{
    /** @todo TODO */
    return 0;
}


void SequenceDatabase::read_and_parse_fasta()
{
    MPI_Offset filesize;
    MPI_Offset localsize;
    MPI_Offset start;
    MPI_Offset end;
    MPI_Offset overlap=1024;
    MPI_File in;
    int rank;
    int size;
    char *chunk;
    int ierr;

    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_CHECK(ierr);
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_CHECK(ierr);
    ierr = MPI_File_open(MPI_COMM_WORLD, const_cast<char*>(filename.c_str()),
            MPI_MODE_RDONLY, MPI_INFO_NULL, &in);
    MPI_CHECK(ierr);

    /* figure out who reads what */

    ierr = MPI_File_get_size(in, &filesize);
    MPI_CHECK(ierr);
    if (budget > filesize) {
        start = 0;
        end = filesize;
        localsize = filesize+1;
    }
    else {
        localsize = filesize/size;
        if (localsize > budget) {
            SEQUENCE_DATABASE_EXCEPTION(
                    "sequence memory budget not sufficient");
        }
        start = rank * localsize;
        /* we fudge the margins based on the memory budget specified */
        start = start - ((budget-localsize)/2);
        end   = end   + ((budget-localsize)/2);
        /* except ranks near the end */
        if (end > filesize) end = filesize;
        /* except ranks near the front */
        if (start < 0) start = 0;
        localsize = end - start + 1;
    }


    /* allocate memory */
    chunk = new char[localsize];

    /* everyone reads in their part */
    MPI_File_read_at_all(in, start, chunk, localsize, MPI_CHAR,
            MPI_STATUS_IGNORE);

    /* figure out where first valid sequence begins */
    MPI_Offset first_valid_index = 0;
    MPI_Offset last_valid_index = 0;
    MPI_Offset maybe_last_valid_index = 0;

    for (MPI_Offset i=0; i<localsize; ++i) {
        if (chunk[i] == '>') {
            if (0 == first_valid_index) {
                first_valid_index = i;
            }
            if (i > maybe_last_valid_index) {
                last_valid_index = maybe_last_valid_index;
                maybe_last_valid_index = i;
            }
        }
    }

    mpix_print_sync(MPI_COMM_WORLD, "first_valid_index", first_valid_index);
    mpix_print_sync(MPI_COMM_WORLD, "last_valid_index", last_valid_index);

    return;
}
