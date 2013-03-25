#include <mpi.h>

#include <cassert>
#include <fstream>
#include <sstream>
#include <string>
#include <utility>

#include <armci.h>

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



SequenceDatabase::SequenceDatabase(const string &file_name, size_t budget)
    :   budget(budget)
    ,   file_name(file_name)
    ,   file_size(0)
    ,   global_file_index()
    ,   global_owner()
    ,   global_address()
    ,   local_cache(NULL)
    ,   local_cache_size(0)
    ,   global_cache()
    ,   sequences()
    ,   remote_cache()
{
    if (!ARMCI_Initialized()) {
        ARMCI_Init();
    }
    read_and_parse_fasta();
}


SequenceDatabase::~SequenceDatabase()
{
    ARMCI_Free(local_cache);
}


size_t SequenceDatabase::get_local_count() const
{
    return sequences.size();
}


size_t SequenceDatabase::get_global_count() const
{
    return global_file_index.size();
}


void SequenceDatabase::read_and_parse_fasta()
{
    MPI_Offset start=0;
    MPI_Offset end=0;
    MPI_Offset overlap=1024;
    MPI_File in=MPI_FILE_NULL;
    int rank=0;
    int size=0;
    char *chunk=NULL;
    int ierr=0;
    unsigned long count;
    size_t sid=0;

    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_CHECK(ierr);
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_CHECK(ierr);

    if (0 == rank) {
        index_file(global_file_index);
    }
    mpix_bcast(global_file_index);
    count = global_file_index.size();

    ierr = MPI_File_open(MPI_COMM_WORLD, const_cast<char*>(file_name.c_str()),
            MPI_MODE_RDONLY, MPI_INFO_NULL, &in);
    MPI_CHECK(ierr);

    /* figure out who reads what */

    ierr = MPI_File_get_size(in, &file_size);
    MPI_CHECK(ierr);
    start = 0;
    end = file_size;
    local_cache_size = file_size;
    if (budget <= file_size) {
        size_t i=0;
        size_t j=0;

        local_cache_size = file_size/size;
        if (local_cache_size > budget) {
            SEQUENCE_DATABASE_EXCEPTION(
                    "sequence memory budget not sufficient");
        }
        start = rank * local_cache_size;
        /* find the first offset >= to our start */
        i = 0;
        while (i<count && global_file_index[i] < start) {
            ++i;
        }
        /* find the first offset >= to our start+local_cache_size */
        j = i;
        while (j<count && global_file_index[j] < start+local_cache_size) {
            ++j;
        }
        /* at this point i and/or j may point past the end of the index */
        if (i >= count) {
            i = count-1; /* reasonable? point to last sequence ID? */
        }
        sid = i;
        start = global_file_index[i];
        if (j >= count) {
            end = file_size;
        }
        else {
            end = global_file_index[j];
        }
        local_cache_size = end - start;
    }

    /* allocate one-sided memory */
    global_cache.resize(size);
    ARMCI_Malloc((void**)&global_cache[0], local_cache_size);
    local_cache = global_cache[rank];

    /* everyone reads in their part */
    ierr = MPI_File_read_at_all(in, start, local_cache, local_cache_size,
            MPI_CHAR, MPI_STATUS_IGNORE);
    MPI_CHECK(ierr);

    /* everyone can close the file now */
    ierr = MPI_File_close(&in);
    MPI_CHECK(ierr);

    assert(local_cache[0] == '>');
    assert(local_cache[local_cache_size-1] != '>');

    global_owner.assign(global_file_index.size(), 0);
    global_address.assign(global_file_index.size(), 0);
    /* parse local cache into Sequence instances */
    {
        size_t label_index=0;
        size_t data_index=0;
        size_t data_length=0;
        for (size_t i=0; i<local_cache_size; ++i) {
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
                    global_owner[sid] = rank;
                    global_address[sid] = &local_cache[data_index];
                    data_length = i - data_index;
                    sequences.insert(make_pair(sid++,
                            Sequence(&local_cache[data_index], data_length)));
                }
            }
        }
    }
    mpix_print_sync(MPI_COMM_WORLD, "mpix_allreduce int");
    mpix_allreduce(global_owner, MPI_SUM);
    mpix_print_sync(MPI_COMM_WORLD, "mpix_allreduce void*");
    mpix_allreduce(global_address, MPI_SUM);
    mpix_print_sync(MPI_COMM_WORLD, "mpix_allreduce DONE");
}


void SequenceDatabase::index_file(vector<MPI_Offset> &index)
{
    string line;
    ifstream fin;

    fin.open(file_name.c_str());
    while (fin) {
        long pos = fin.tellg();
        getline(fin, line);
        if (line[0] == '>') {
            index.push_back(pos);
        }
    }

    fin.close();
}


Sequence& SequenceDatabase::get_sequence(size_t i)
{
    MPI_Offset size;

    if (!sequences.count(i)) {
        if (i-1 == global_file_index.size()) {
            size = file_size-global_file_index[i];
        }
        else {
            size = global_file_index[i+1]-global_file_index[i];
        }

        cout << "SequenceDatabase::get_sequence("<<i<<")"<<endl;
        cout << "\tsize="<<size<<endl;

        remote_cache[i] = (char*)ARMCI_Malloc_local(size);
        (void)memset(remote_cache[i], '@', size);
        assert(NULL != remote_cache[i]);
        cout << "ARMCI_Get(" << (void*)global_address[i]
            << "," << (void*)remote_cache[i]
            << "," << size
            << "," << global_owner[i]
            << ")" << endl;
        ARMCI_Get(global_address[i], remote_cache[i], size, global_owner[i]);
        sequences.insert(make_pair(i, Sequence(remote_cache[i], size)));
    }

    return sequences[i];
}

