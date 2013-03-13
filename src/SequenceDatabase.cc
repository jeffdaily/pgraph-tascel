#include <mpi.h>

#include <cassert>
#include <fstream>
#include <sstream>
#include <string>
#include <utility>

#include "Sequence.h"
#include "SequenceDatabase.h"
#include "mpix.h"

using std::ifstream;
using std::make_pair;
using std::ostringstream;
using std::string;


SequenceDatabaseException::SequenceDatabaseException(
        const char *file,
        int line,
        const char *function,
        const char *message) throw()
{
    ostringstream oss;
    oss << file << ": " << function << ": " << line << ": " << message;
    this->message = oss.str();
}


SequenceDatabaseException::~SequenceDatabaseException() throw()
{

}


const char* SequenceDatabaseException::what() const throw()
{
    return this->message.c_str();
}



SequenceDatabase::SequenceDatabase(const string &filename, size_t budget)
    :   filename(filename)
    ,   budget(budget)
{
    read_and_parse_fasta();
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
    MPI_Offset filesize=0;
    MPI_Offset localsize=0;
    MPI_Offset start=0;
    MPI_Offset end=0;
    MPI_Offset overlap=1024;
    MPI_File in=MPI_FILE_NULL;
    int rank=0;
    int size=0;
    char *chunk=NULL;
    int ierr=0;
    MPI_Datatype offset_type=0;
    vector<MPI_Offset> index;
    unsigned long count;
    size_t sid=0;

    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_CHECK(ierr);
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_CHECK(ierr);

    offset_type = MPI_DATATYPE_NULL;
    if      (sizeof(MPI_Offset) == sizeof(int))  { offset_type = MPI_INT; }
    else if (sizeof(MPI_Offset) == sizeof(long)) { offset_type = MPI_LONG; }
    else if (sizeof(MPI_Offset) == sizeof(long long)) { offset_type =  
        MPI_LONG_LONG; }
    else { MPI_Abort(MPI_COMM_WORLD, 1); }

    if (0 == rank) {
        index_file(index);
        count = index.size();
        ierr = MPI_Bcast(&count, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
        MPI_CHECK(ierr);
        ierr = MPI_Bcast(&index[0], count, offset_type, 0, MPI_COMM_WORLD);
        MPI_CHECK(ierr);
    }
    else {
        ierr = MPI_Bcast(&count, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
        MPI_CHECK(ierr);
        index.resize(count);
        ierr = MPI_Bcast(&index[0], count, offset_type, 0, MPI_COMM_WORLD);
        MPI_CHECK(ierr);
    }

    ierr = MPI_File_open(MPI_COMM_WORLD, const_cast<char*>(filename.c_str()),
            MPI_MODE_RDONLY, MPI_INFO_NULL, &in);
    MPI_CHECK(ierr);

    /* figure out who reads what */

    ierr = MPI_File_get_size(in, &filesize);
    MPI_CHECK(ierr);
    start = 0;
    end = filesize;
    localsize = filesize;
    if (budget <= filesize) {
        size_t i=0;
        size_t j=0;

        localsize = filesize/size;
        if (localsize > budget) {
            SEQUENCE_DATABASE_EXCEPTION(
                    "sequence memory budget not sufficient");
        }
        start = rank * localsize;
        /* find the first offset >= to our start */
        i = 0;
        while (i<count && index[i] < start) {
            ++i;
        }
        /* find the first offset >= to our start+localsize */
        j = i;
        while (j<count && index[j] < start+localsize) {
            ++j;
        }
        /* at this point i and/or j may point past the end of the index */
        if (i >= count) {
            i = count-1; /* reasonable? point to last sequence ID? */
        }
        sid = i;
        start = index[i];
        if (j >= count) {
            end = filesize;
        }
        else {
            end = index[j];
        }
        localsize = end - start;
    }

    /* allocate memory */
    local_cache.resize(localsize);

    /* everyone reads in their part */
    ierr = MPI_File_read_at_all(in, start, &local_cache[0], localsize,
            MPI_CHAR, MPI_STATUS_IGNORE);
    MPI_CHECK(ierr);

    /* everyone can close the file now */
    ierr = MPI_File_close(&in);
    MPI_CHECK(ierr);

    assert(local_cache.front() == '>');
    assert(local_cache.back()  != '>');

    /* parse local cache into Sequence instances */
    {
        size_t label_index=0;
        size_t data_index=0;
        size_t data_length=0;
        for (size_t i=0; i<local_cache.size(); ++i) {
            if (local_cache[i] == '>') {
                label_index = i;
                data_index = i;
                data_length = 0;
            }
            else if (local_cache[i] == '\n') {
                if (label_index == data_index) {
                    data_index = i+1;
                }
                else {
                    data_length = i - data_index;
                    sequences.insert(make_pair(sid++,
                            Sequence(&local_cache[data_index], data_length)));
                }
            }
        }
    }
}


void SequenceDatabase::index_file(vector<MPI_Offset> &index)
{
    string line;
    ifstream fin;

    fin.open(filename.c_str());
    while (fin) {
        long pos = fin.tellg();
        getline(fin, line);
        if (line[0] == '>') {
            index.push_back(pos);
        }
    }

    fin.close();
}

