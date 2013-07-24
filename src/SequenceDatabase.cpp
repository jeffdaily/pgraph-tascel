/**
 * @file SequenceDatabase.cc
 *
 * @author jeff.daily@pnnl.gov
 *
 * @todo TODO handle fasta files with multiline sequences
 *
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 */
#include <mpi.h>

#include <cassert>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>
#include <utility>

#include <armci.h>

#include "Sequence.hpp"
#include "SequenceDatabase.hpp"
#include "mpix.hpp"

using std::accumulate;
using std::cerr;
using std::endl;
using std::ifstream;
using std::make_pair;
using std::ostringstream;
using std::string;

#define TAG 2345
#define DELIMITER '$'


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


const char *SequenceDatabaseException::what() const throw()
{
    return this->message.c_str();
}



SequenceDatabase::SequenceDatabase(const string &file_name, size_t budget,
                                   MPI_Comm comm)
    :   comm(comm)
    ,   comm_rank(0)
    ,   comm_size(0)
    ,   is_replicated(false)
    ,   budget(budget)
    ,   file_name(file_name)
    ,   local_data(NULL)
    ,   local_cache()
    ,   owners()
    ,   addresses()
    ,   sizes()
    ,   ptr_arr(NULL)
    ,   remote_cache()
{
    int ierr = 0;

    /* rank and size */
    ierr = MPI_Comm_rank(comm, &comm_rank);
    CHECK_MPI_IERR(ierr, comm_rank, comm);
    ierr = MPI_Comm_size(comm, &comm_size);
    CHECK_MPI_IERR(ierr, comm_rank, comm);

    read_and_parse_fasta();
}


SequenceDatabase::~SequenceDatabase()
{
    if (is_replicated) {
    }
    else {
        ARMCI_Free(local_data);
    }
}


size_t SequenceDatabase::get_local_count() const
{
    return local_cache.size();
}


size_t SequenceDatabase::get_global_count() const
{
    return global_count;
}


void SequenceDatabase::read_and_parse_fasta()
{
    MPI_File in = MPI_FILE_NULL;
    MPI_Offset file_size = 0;
    MPI_Offset budget = MPI_Offset(this->budget);
    int ierr = 0;

    /* open file and get file size */
    ierr = MPI_File_open(comm, const_cast<char *>(file_name.c_str()),
                         MPI_MODE_RDONLY, MPI_INFO_NULL, &in);
    CHECK_MPI_IERR(ierr, comm_rank, comm);
    ierr = MPI_File_get_size(in, &file_size);
    CHECK_MPI_IERR(ierr, comm_rank, comm);

    mpix_print_zero(comm, "memory budget", budget);
    mpix_print_zero(comm, "file size", file_size);

    if (budget >= file_size) {
        is_replicated = true;
        read_and_parse_fasta_himem(in, file_size);
    }
    else {
        is_replicated = false;
        read_and_parse_fasta_lomem(in, file_size);
    }
}


void SequenceDatabase::read_and_parse_fasta_lomem(MPI_File in,
                                                  MPI_Offset file_size)
{
    MPI_Offset start = 0;
    MPI_Offset end = 0;
    MPI_Offset chunk_size = 0;
    MPI_Offset local_size = 0;
    int ierr = 0;
    char *file_buffer = NULL;
    char *extra_buffer = NULL;
    int first_index_mine = -1;
    int first_index_other = 0;
    int local_cache_size = 0;
    size_t new_size = 0;
    vector<int> counts(comm_size, 0);
    int first_id = 0;

    /* initialize ARMCI, if needed */
    if (!ARMCI_Initialized()) {
        ARMCI_Init();
    }

    chunk_size = file_size / MPI_Offset(comm_size);
    if (budget < chunk_size) {
        if (0 == comm_rank) cerr << "insufficient memory budget" << endl;
        MPI_Abort(comm, -1);
    }
    start = comm_rank * chunk_size;
    end = (comm_rank+1) * chunk_size;
    if (comm_size == comm_rank+1) {
        /* last proc gets remainder of file*/
        end = file_size;
    }
    local_size = end-start;

    /* allocate memory to hold file buffer */
    file_buffer = new char[local_size+1];

    /* everyone reads in their part */
    ierr = MPI_File_read_at_all(in, start, file_buffer, local_size,
                                MPI_CHAR, MPI_STATUS_IGNORE);
    CHECK_MPI_IERR(ierr, comm_rank, comm);
    file_buffer[local_size] = '\0';

    /* everyone can close the file now */
    ierr = MPI_File_close(&in);
    CHECK_MPI_IERR(ierr, comm_rank, comm);

    /* locate first '>' indicating start of a sequence ID */
    /* also count how many '>' on each process so we can determine which proc
     * owns which sequences */
    assert(first_index_mine != local_size);
    for (MPI_Offset i=0; i<local_size; ++i) {
        if (file_buffer[i] == '>') {
            ++counts[comm_rank];
            if (first_index_mine < 0) {
                first_index_mine = i;
            }
        }
    }
    mpix_allreduce(counts, MPI_MAX, comm);

    /* send broken chunks left */
    if (comm_size > 0) {
        if (0 == comm_rank) {
            ierr = MPI_Recv(&first_index_other, 1, MPI_INT,
                    1, TAG, comm, MPI_STATUS_IGNORE);
            CHECK_MPI_IERR(ierr, comm_rank, comm);
            extra_buffer = new char[first_index_other+1];
            ierr = MPI_Recv(extra_buffer, first_index_other, MPI_CHAR,
                    1, TAG, comm, MPI_STATUS_IGNORE);
            CHECK_MPI_IERR(ierr, comm_rank, comm);
        }
        else if (comm_size == comm_rank+1) {
            ierr = MPI_Send(&first_index_mine, 1, MPI_INT,
                    comm_rank-1, TAG, comm);
            CHECK_MPI_IERR(ierr, comm_rank, comm);
            extra_buffer = new char[1];
            ierr = MPI_Send(file_buffer, first_index_mine, MPI_CHAR,
                    comm_rank-1, TAG, comm);
            CHECK_MPI_IERR(ierr, comm_rank, comm);
        }
        else {
            ierr = MPI_Sendrecv(&first_index_mine, 1, MPI_INT,
                    comm_rank-1, TAG,
                    &first_index_other, 1, MPI_INT,
                    comm_rank+1, TAG, comm, MPI_STATUS_IGNORE);
            CHECK_MPI_IERR(ierr, comm_rank, comm);
            extra_buffer = new char[first_index_other+1];
            ierr = MPI_Sendrecv(file_buffer, first_index_mine, MPI_CHAR,
                    comm_rank-1, TAG,
                    extra_buffer, first_index_other, MPI_CHAR,
                    comm_rank+1, TAG, comm, MPI_STATUS_IGNORE);
            CHECK_MPI_IERR(ierr, comm_rank, comm);
        }
        extra_buffer[first_index_other] = '\0';
    }
    local_cache_size = local_size - first_index_mine + first_index_other;

    assert(comm == MPI_COMM_WORLD); /* until we handle ARMCI groups */

    /* allocate the ARMCI memory to hold the sequences and copy buffers */
    char **ptr_arr = new char*[comm_size];
    (void)ARMCI_Malloc((void**)ptr_arr, local_cache_size+1);
    local_data = ptr_arr[comm_rank];
    if (comm_size > 0) {
        if (0 == comm_rank) {
            (void)memcpy(ptr_arr[comm_rank],
                         &file_buffer[first_index_mine],
                         local_size-first_index_mine);
            (void)memcpy(&ptr_arr[comm_rank][local_size-first_index_mine],
                         extra_buffer,
                         first_index_other);
            ptr_arr[comm_rank][local_cache_size] = '\0';
        }
        else if (comm_size == comm_rank+1) {
            (void)memcpy(ptr_arr[comm_rank],
                         &file_buffer[first_index_mine],
                         local_size-first_index_mine);
            ptr_arr[comm_rank][local_cache_size] = '\0';
        }
        else {
            (void)memcpy(ptr_arr[comm_rank],
                         &file_buffer[first_index_mine],
                         local_size-first_index_mine);
            (void)memcpy(&ptr_arr[comm_rank][local_size-first_index_mine],
                         extra_buffer,
                         first_index_other);
            ptr_arr[comm_rank][local_cache_size] = '\0';
        }
    }

    /* determine the first sequence index */
    for (int i=0; i<comm_rank; ++i) {
        first_id += counts[i];
    }

    pack_and_index_fasta(ptr_arr[comm_rank], local_cache_size, DELIMITER,
                         first_id, new_size);

    exchange_local_cache();
}


void SequenceDatabase::read_and_parse_fasta_himem(MPI_File in,
                                                  MPI_Offset file_size)
{
    int ierr = 0;
    size_t new_size = 0;

    /* read directly into a local buffer; don't use ARMCI */
    local_data = new char[file_size+2]; /* +2 for last delim and null */
    local_data[file_size] = DELIMITER;
    local_data[file_size+1] = '\0';

    /* everyone reads in their part */
    ierr = MPI_File_read_at_all(in, 0, local_data, file_size,
            MPI_CHAR, MPI_STATUS_IGNORE);
    CHECK_MPI_IERR(ierr, comm_rank, comm);

    /* everyone can close the file now */
    ierr = MPI_File_close(&in);
    CHECK_MPI_IERR(ierr, comm_rank, comm);

    /* pack and index the fasta buffer */
    pack_and_index_fasta(local_data, file_size, DELIMITER, 0, new_size);
    local_data[new_size] = DELIMITER;
    local_data[new_size+1] = '\0';
    assert(!local_cache.empty());
    assert(local_cache.size() > 0);

    global_count = local_cache.size();
}


void SequenceDatabase::pack_and_index_fasta(char *buffer,
                                            size_t size,
                                            char delimiter,
                                            size_t id,
                                            size_t &new_size)
{
    size_t r = 0;           /* read index */
    size_t w = 0;           /* write index */
    size_t last_gt = 0;     /* index of last '>' seen */
    size_t last_hash = 0;   /* index of last '#' inserted */

    new_size = 0;

    /* We skip all newlines unless it is the newline that terminates the
     * sequence ID. We replace the ID-terminating newline with '#'. We replace
     * the sequence terminating newline with the given delimiter. */

    if (buffer[0] != '>') {
        cerr << '[' << comm_rank << "] buffer[0] != '>'" << endl;
        MPI_Abort(comm, -1);
    }

    while (r<size) {
        if (buffer[r] == '>') {
            last_gt = w;
            while (r<size && buffer[r] != '\n') {
                buffer[w++] = buffer[r++];
            }
            last_hash = w;
            buffer[w++] = '#';
            r++;
        }
        else {
            while (r<size && buffer[r] != '\n') {
                buffer[w++] = buffer[r++];
            }
            /* peek at next character, either EOF or '>' */
            if (r<size) {
                if (r+1 >= size || buffer[r+1] == '>') {
                    buffer[w++] = delimiter;
                    local_cache[id++] =
                        new Sequence(&buffer[last_gt],
                                     0,
                                     last_hash - last_gt,
                                     last_hash - last_gt + 1,
                                     w - last_hash - 1);
                    assert(!local_cache.empty());
                    assert(local_cache.size() > 0);
                }
            }
            /* or file wasn't terminated with a newline */
            else {
                /* hopefully space enough to append the delimiter anyway */
                assert(w<r);
                buffer[w++] = delimiter;
                local_cache[id++] =
                    new Sequence(&buffer[last_gt],
                                 0,
                                 last_hash - last_gt,
                                 last_hash - last_gt + 1,
                                 w - last_hash - 1);
                assert(!local_cache.empty());
                assert(local_cache.size() > 0);
                assert(0);
            }
            r++;
        }
    }

    assert(!local_cache.empty());
    assert(local_cache.size() > 0);

    new_size = w;
}



Sequence &SequenceDatabase::get_sequence(size_t i)
{
    if (is_replicated) {
        return *local_cache[i];
    } else {
        if (comm_rank == owners.at(i)) {
            return *local_cache[i];
        }
        else if (remote_cache.count(i)) {
            return *remote_cache[i];
        }
        else {
            Sequence *sequence = NULL;
            char *buffer = NULL;
            char *data = NULL;
            size_t id_len = 0;
            size_t data_len = 0;
            size_t j = 0;

            buffer = (char*)ARMCI_Malloc_local(sizes[i]+2);
            ARMCI_Get((void*)addresses[i], buffer, sizes[i], owners[i]);
            assert('>' == buffer[0]);
            assert('$' == buffer[sizes[i]-1]);
            for (j=0; j<sizes[i]; ++j) {
                if ('#' == buffer[j]) {
                    break;
                }
            }
            assert(j<sizes[i]);
            remote_cache[i] = new Sequence(buffer, 0, j, j+1, sizes[i]-j-1);
            return *remote_cache[i];
        }
    }
}


void SequenceDatabase::exchange_local_cache()
{
    map<size_t, Sequence*>::const_iterator it;

    global_count = local_cache.size();
    mpix_allreduce(global_count, MPI_SUM, comm);

    owners.assign(global_count, 0);
    addresses.assign(global_count, 0);
    sizes.assign(global_count, 0);

    for (it=local_cache.begin(); it!=local_cache.end(); ++it) {
        const char *data = NULL;
        size_t data_size = 0;

        assert(it->first >= 0);
        assert(it->first < global_count);
        it->second->get_data(data, data_size);
        owners[it->first] = comm_rank;
        addresses[it->first] = MPI_Aint(data);
        sizes[it->first] = data_size;
    }

    mpix_allreduce(owners, MPI_MAX, comm);
    mpix_allreduce(addresses, MPI_MAX, comm);
    mpix_allreduce(sizes, MPI_MAX, comm);
}

