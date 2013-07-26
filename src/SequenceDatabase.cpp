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

/* ANL's armci does not extern "C" inside the header */
extern "C" {
#include <armci.h>
}

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

#define DEBUG 0
#define TAG 2345
#define DELIMITER '$'

#if DEBUG || 1
string carr_to_string(
    char* *collection,
    const size_t &size,
    char const *delimiter = ",",
    const string &name = "")
{
    std::ostringstream os;

    if (!name.empty()) {
        os << name << "=";
    }

    if (size == 0) {
        os << "{}";
    }
    else {
        os << "{" << (void*)collection[0];
        for (size_t i=1; i < size; ++i) {
            os << delimiter << (void*)collection[i];
        }
        os << "}";
    }

    return os.str();
}
#endif


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

    ARMCI_Init();

#if defined(ARMCI_TS)
    ARMCI_Group group_world;
    ARMCI_Group_get_world(&group_world);
    ARMCI_Group_rank(&group_world, &comm_rank);
    ARMCI_Group_size(&group_world, &comm_size);
    comm = armci_group_comm(&group_world);
#else
    /* rank and size */
    ierr = MPI_Comm_rank(comm, &comm_rank);
    CHECK_MPI_IERR(ierr, comm_rank, comm);
    ierr = MPI_Comm_size(comm, &comm_size);
    CHECK_MPI_IERR(ierr, comm_rank, comm);
#endif

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

#if 0
    /* initialize ARMCI, if needed */
    if (!ARMCI_Initialized()) {
        ARMCI_Init();
    }
#endif

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
    ptr_arr = new char*[comm_size];
    (void)ARMCI_Malloc((void**)ptr_arr, local_cache_size+1);
    local_data = ptr_arr[comm_rank];
    //mpix_print_sync(comm, "local_data", (void*)local_data);
    //mpix_print_sync(comm, "ptr_arr", carr_to_string(ptr_arr, comm_size));
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

#if DEBUG
    for (size_t p=0; p<comm_size; ++p) {
        if (p == comm_rank) {
            cout << "[" << p << "] ARMCI ptr_arr[me]" << endl;
            cout << ptr_arr[p] << endl;
        }
        ARMCI_Barrier();
    }
#endif

    /* determine the first sequence index */
    for (int i=0; i<comm_rank; ++i) {
        first_id += counts[i];
    }

    pack_and_index_fasta(ptr_arr[comm_rank], local_cache_size,
            DELIMITER, first_id, new_size);

#if DEBUG
    for (size_t p=0; p<comm_size; ++p) {
        if (p == comm_rank) {
            cout << "[" << p << "] ARMCI ptr_arr[me]" << endl;
            cout << (void*)ptr_arr[p] << endl;
            cout << ptr_arr[p] << endl;
        }
        ARMCI_Barrier();
    }
#endif


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
            size_t j = 0;

            buffer = (char*)ARMCI_Malloc_local(sizes[i]+2);
            assert(buffer);
            (void)memset(buffer, 0, sizes[i]+2);
#if 1
            ARMCI_Get((void*)addresses[i], buffer, sizes[i], owners[i]);
#else
                cout << "[" << comm_rank << "] "
                    << addresses[i]
                    << "\t"
                    << addresses[i]
                    << "\t'"
                    << string(addresses[i], sizes[i])
                    << "'" << endl;
#endif
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
    vector<ptrdiff_t> offsets;

    global_count = local_cache.size();
    mpix_allreduce(global_count, MPI_SUM, comm);

    owners.assign(global_count, 0);
    addresses.assign(global_count, 0);
    sizes.assign(global_count, 0);
    offsets.assign(global_count, 0);

    for (it=local_cache.begin(); it!=local_cache.end(); ++it) {
        const char *data = NULL;
        size_t data_size = 0;

        assert(it->first >= 0);
        assert(it->first < global_count);

        it->second->get_data(data, data_size);
        owners[it->first] = comm_rank;
        offsets[it->first] = data - local_data;
        sizes[it->first] = data_size;
        //cout << "[" << comm_rank << "] "
        //    << (void*)data
        //    << "\t"
        //    << MPI_Aint(data)
        //    << "\t'"
        //    << string(data, data_size)
        //    << "'" << endl;
    }

    mpix_allreduce(owners, MPI_SUM, comm);
    mpix_allreduce(offsets, MPI_SUM, comm);
    mpix_allreduce(sizes, MPI_SUM, comm);
    //mpix_print_sync(comm, "owners", vec_to_string(owners));
    //mpix_print_sync(comm, "offsets", vec_to_string(offsets));
    //mpix_print_sync(comm, "sizes", vec_to_string(sizes));

    /* ARMCI is screwey. If the MPI ranks are on the same SMP node, the
     * addresses returned in the ptr_arr are different than the addresses on
     * the procs (due to mmapping the shared memory on each process
     * differently) so we can't use the addresses as we have them. We must make
     * everything relative to each process' own ptr_arr. */
    for (size_t i=0; i<global_count; ++i) {
        addresses[i] = &ptr_arr[owners[i]][offsets[i]];
    }
    //mpix_print_sync(comm, "addresses", vec_to_string(addresses));
}

