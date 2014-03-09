/**
 * @file SequenceDatabaseArmci.cpp
 *
 * @author jeff.daily@pnnl.gov
 *
 * @todo TODO handle fasta files with multiline sequences
 *
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 */
#include "config.h"

#include <mpi.h>

#include <cassert>
#include <cctype>
#include <cstddef>
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

#include <tascel.h>

#include "alignment.hpp"
#include "mpix.hpp"
#include "SequenceDatabaseArmci.hpp"
#include "SequenceDatabase.hpp"
#include "Sequence.hpp"

using std::accumulate;
using std::cerr;
using std::endl;
using std::ifstream;
using std::make_pair;
using std::ostringstream;
using std::size_t;
using std::string;
using tascel::LockGuard;
using tascel::PthreadMutex;

#define TAG 2345

namespace pgraph {


SequenceDatabaseArmci::SequenceDatabaseArmci(
        const string &file_name,
        size_t budget,
        MPI_Comm comm,
        char delimiter)
    :   SequenceDatabase()
    ,   comm_orig(comm)
    ,   comm_orig_rank(0)
    ,   comm_orig_size(0)
    ,   comm(MPI_COMM_NULL)
    ,   comm_rank(0)
    ,   comm_size(0)
    ,   is_replicated(false)
    ,   budget(budget)
    ,   delimiter(delimiter)
    ,   file_name(file_name)
    ,   local_data(NULL)
    ,   local_cache()
    ,   global_count(0)
    ,   global_size(0)
    ,   local_size(0)
    ,   owners()
    ,   addresses()
    ,   sizes()
    ,   ptr_arr(NULL)
    ,   remote_cache()
    ,   mutex()
    ,   max_seq_size(0)
{
    ARMCI_Init();

    /* rank and size */
    comm_orig_rank = mpix::comm_rank(comm_orig);
    comm_orig_size = mpix::comm_size(comm_orig);

    /* create ARMCI parent group -- I hope this works */
    vector<int> pid_list(comm_orig_size, 0);
    pid_list[comm_orig_rank] = comm_orig_rank;
    mpix::allreduce(pid_list, MPI_MAX, comm_orig);
    ARMCI_Group_create(comm_orig_size, &pid_list[0], &armci_group_orig);
    {
        int armci_rank = 0;
        int armci_size = 0;
        ARMCI_Group_rank(&armci_group_orig, &armci_rank);
        ARMCI_Group_size(&armci_group_orig, &armci_size);
        assert(armci_rank == comm_orig_rank);
        assert(armci_size == comm_orig_size);
    }

    read_and_parse_fasta();
}


SequenceDatabaseArmci::~SequenceDatabaseArmci()
{
    map<size_t, Sequence*>::iterator it;

    for (it=local_cache.begin(); it!=local_cache.end(); ++it) {
        delete it->second;
    }

    if (is_replicated) {
        delete [] local_data;
    }
    else {
        ARMCI_Free(local_data);
    }

    ARMCI_Group_free(&armci_group_orig);
}


size_t SequenceDatabaseArmci::size() const
{
    return global_count;
}


size_t SequenceDatabaseArmci::char_size() const
{
    return global_size;
}


void SequenceDatabaseArmci::read_and_parse_fasta()
{
    MPI_File in = MPI_FILE_NULL;
    MPI_Offset file_size = 0;
    MPI_Offset budget = MPI_Offset(this->budget);
    int ierr = 0;

    /* open file and get file size on rank 0*/
    if (0 == comm_orig_rank) {
        ierr = MPI_File_open(MPI_COMM_SELF,
                const_cast<char *>(file_name.c_str()),
                MPI_MODE_RDONLY, MPI_INFO_NULL, &in);
        mpix::check(ierr);
        ierr = MPI_File_get_size(in, &file_size);
        mpix::check(ierr);
        ierr = MPI_File_close(&in);
        mpix::check(ierr);
    }
    mpix::bcast(file_size, 0, comm_orig);

    mpix::print_zero("memory budget", budget, comm_orig);
    mpix::print_zero("file size", file_size, comm_orig);
    
    if (budget >= file_size) {
        is_replicated = true;
        comm = comm_orig;
        comm_rank = comm_orig_rank;
        comm_size = comm_orig_size;
        armci_group = armci_group_orig;
#if 0
        ierr = MPI_File_open(MPI_COMM_SELF,
                const_cast<char *>(file_name.c_str()),
                MPI_MODE_RDONLY, MPI_INFO_NULL, &in);
        mpix::check(ierr);
        read_and_parse_fasta_himem(in, file_size);
#else
        read_and_parse_fasta_himem(MPI_FILE_NULL, file_size);
#endif
    }
    else {
        MPI_Offset ranks_per_domain = 0;
        MPI_Offset domain_count = 0;
        int color = 0;

        is_replicated = false;
        /* split comm_orig based on budget memory domains */
        ranks_per_domain = file_size / budget + (file_size%budget ? 1 : 0);
        mpix::print_zero("ranks_per_domain", ranks_per_domain, comm_orig);
        domain_count = comm_orig_size / ranks_per_domain;
        mpix::print_zero("domain_count", domain_count, comm_orig);
        color = comm_orig_rank / ranks_per_domain;
        assert(color <= domain_count);
        if (color == domain_count) {
            color -= 1;
        }
        //mpix::print_sync("color", color, comm_orig);
        ierr = MPI_Comm_split(comm_orig, color, comm_orig_rank, &comm);
        mpix::check(ierr);
        ierr = MPI_Comm_rank(comm, &comm_rank);
        mpix::check(ierr);
        ierr = MPI_Comm_size(comm, &comm_size);
        mpix::check(ierr);
        vector<int> pid_list(comm_size, 0);
        pid_list[comm_rank] = comm_orig_rank;
        mpix::allreduce(pid_list, MPI_MAX, comm);
        ARMCI_Group_create_child(comm_size, &pid_list[0], &armci_group, &armci_group_orig);
        {
            int armci_rank = 0;
            int armci_size = 0;
            ARMCI_Group_rank(&armci_group, &armci_rank);
            ARMCI_Group_size(&armci_group, &armci_size);
            assert(armci_rank == comm_rank);
            assert(armci_size == comm_size);
        }
        ierr = MPI_File_open(comm,
                const_cast<char *>(file_name.c_str()),
                MPI_MODE_RDONLY, MPI_INFO_NULL, &in);
        mpix::check(ierr);
        read_and_parse_fasta_lomem(in, file_size);
    }
}


void SequenceDatabaseArmci::read_and_parse_fasta_lomem(MPI_File in,
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

    chunk_size = file_size / MPI_Offset(comm_size);
    if (MPI_Offset(budget) < chunk_size) {
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
    mpix::check(ierr);
    file_buffer[local_size] = '\0';

    /* everyone can close the file now */
    ierr = MPI_File_close(&in);
    mpix::check(ierr);

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
    mpix::allreduce(counts, MPI_MAX, comm);

    /* send broken chunks left */
    if (comm_size > 0) {
        if (0 == comm_rank) {
            ierr = MPI_Recv(&first_index_other, 1, MPI_INT,
                    1, TAG, comm, MPI_STATUS_IGNORE);
            mpix::check(ierr);
            extra_buffer = new char[first_index_other+1];
            ierr = MPI_Recv(extra_buffer, first_index_other, MPI_CHAR,
                    1, TAG, comm, MPI_STATUS_IGNORE);
            mpix::check(ierr);
        }
        else if (comm_size == comm_rank+1) {
            ierr = MPI_Send(&first_index_mine, 1, MPI_INT,
                    comm_rank-1, TAG, comm);
            mpix::check(ierr);
            extra_buffer = new char[1];
            ierr = MPI_Send(file_buffer, first_index_mine, MPI_CHAR,
                    comm_rank-1, TAG, comm);
            mpix::check(ierr);
        }
        else {
            ierr = MPI_Sendrecv(&first_index_mine, 1, MPI_INT,
                    comm_rank-1, TAG,
                    &first_index_other, 1, MPI_INT,
                    comm_rank+1, TAG, comm, MPI_STATUS_IGNORE);
            mpix::check(ierr);
            extra_buffer = new char[first_index_other+1];
            ierr = MPI_Sendrecv(file_buffer, first_index_mine, MPI_CHAR,
                    comm_rank-1, TAG,
                    extra_buffer, first_index_other, MPI_CHAR,
                    comm_rank+1, TAG, comm, MPI_STATUS_IGNORE);
            mpix::check(ierr);
        }
        extra_buffer[first_index_other] = '\0';
    }
    local_cache_size = local_size - first_index_mine + first_index_other;

    //assert(comm == MPI_COMM_WORLD); /* until we handle ARMCI groups */

    /* allocate the ARMCI memory to hold the sequences and copy buffers */
    ptr_arr = new char*[comm_size];
    (void)ARMCI_Malloc_group((void**)ptr_arr, local_cache_size+1, &armci_group);
    local_data = ptr_arr[comm_rank];
    //mpix::print_sync("local_data", (void*)local_data, comm);
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
    delete [] file_buffer;

    /* determine the first sequence index */
    for (int i=0; i<comm_rank; ++i) {
        first_id += counts[i];
    }

    pack_and_index_fasta(ptr_arr[comm_rank], local_cache_size,
            delimiter, first_id, new_size);
    mpix::allreduce(global_size, MPI_SUM, comm);
    //mpix::print_sync("global_size", global_size, comm);

    exchange_local_cache();
}


void SequenceDatabaseArmci::read_and_parse_fasta_himem(MPI_File in,
                                                  MPI_Offset file_size)
{
    int ierr = 0;
    size_t new_size = 0;
    MPI_Offset file_size_verify;

    /* read directly into a local buffer; don't use ARMCI */
    local_data = new char[file_size+2]; /* +2 for last delim and null */
    local_data[file_size] = delimiter;
    local_data[file_size+1] = '\0';

    if (in == MPI_FILE_NULL) {
        mpix::read_file_bcast(file_name, local_data, file_size_verify, comm_orig);
        assert(file_size == file_size_verify);
    }
    else {
        /* everyone reads in their part */
        ierr = MPI_File_read_at_all(in, 0, local_data, file_size,
                MPI_CHAR, MPI_STATUS_IGNORE);
        mpix::check(ierr);

        /* everyone can close the file now */
        ierr = MPI_File_close(&in);
        mpix::check(ierr);
    }

    /* pack and index the fasta buffer */
    pack_and_index_fasta(local_data, file_size, delimiter, 0, new_size);
    local_data[new_size] = delimiter;
    local_data[new_size+1] = '\0';
    assert(!local_cache.empty());
    assert(local_cache.size() > 0);

    global_count = local_cache.size();
    //mpix::print_sync("global_size", global_size, comm);
}


void SequenceDatabaseArmci::pack_and_index_fasta(char *buffer,
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
                    size_t id_offset = 0;
                    size_t id_length = last_hash - last_gt;
                    size_t sequence_offset = last_hash - last_gt + 1;
                    size_t sequence_length = w - last_hash - 1;
                    if (delimiter == '\0') {
                        sequence_length -= 1;
                    }
                    local_cache[id++] =
                        new Sequence(&buffer[last_gt],
                                     id_offset,
                                     id_length,
                                     sequence_offset,
                                     sequence_length);
                    assert(!local_cache.empty());
                    assert(local_cache.size() > 0);
                    {
                        size_t l = local_cache[id-1]->get_sequence_length();
                        max_seq_size = l > max_seq_size ? l : max_seq_size;
                        local_size += l;
                    }
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
                {
                    size_t l = local_cache[id-1]->get_sequence_length();
                    max_seq_size = l > max_seq_size ? l : max_seq_size;
                    local_size += l;
                }
                assert(0);
            }
            r++;
        }
    }

    assert(!local_cache.empty());
    assert(local_cache.size() > 0);

    new_size = w;

    mpix::allreduce(max_seq_size, MPI_MAX, comm);
    mpix::print_zero("max_seq_size", max_seq_size, comm);
    //mpix::print_sync("local_size", local_size, comm);
    global_size = local_size;
}



Sequence &SequenceDatabaseArmci::get_sequence(size_t i)
{
    if (is_replicated) {
        return *local_cache[i];
    } else {
        if (comm_rank == owners.at(i)) {
            return *local_cache[i];
        }
        else if (remote_cache.count(i)) {
            LockGuard<PthreadMutex> guard(mutex);
            return *remote_cache[i];
        }
        else {
            Sequence *sequence = NULL;
            char *buffer = NULL;
            size_t j = 0;
            int global_rank = ARMCI_Absolute_id(&armci_group, owners[i]);

#if 1
            {
                LockGuard<PthreadMutex> guard(mutex);
                buffer = (char*)ARMCI_Malloc_local(sizes[i]+2);
            }
#endif
            assert(buffer);
            (void)memset(buffer, 0, sizes[i]+2);
#if 1
            {
                LockGuard<PthreadMutex> guard(mutex);
                ARMCI_Get((void*)addresses[i], buffer, sizes[i], global_rank);
            }
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
            //assert(delimiter == buffer[sizes[i]-1]);
            for (j=0; j<sizes[i]; ++j) {
                if ('#' == buffer[j]) {
                    break;
                }
            }
            assert(j<sizes[i]);
            sequence = new Sequence(buffer, 0, j, j+1, sizes[i]-j-1, true);
            {
                LockGuard<PthreadMutex> guard(mutex);
                remote_cache[i] = sequence;
            }
            return *sequence;
        }
    }
}


Sequence &SequenceDatabaseArmci::operator[](size_t i)
{
    return this->get_sequence(i);
}


void SequenceDatabaseArmci::exchange_local_cache()
{
    map<size_t, Sequence*>::const_iterator it;
    vector<ptrdiff_t> offsets;

    global_count = local_cache.size();
    mpix::allreduce(global_count, MPI_SUM, comm);

    owners.assign(global_count, 0);
    addresses.assign(global_count, 0);
    sizes.assign(global_count, 0);
    offsets.assign(global_count, 0);

    for (it=local_cache.begin(); it!=local_cache.end(); ++it) {
        const char *data = NULL;
        size_t data_size = 0;

        assert(it->first < global_count);

        it->second->get_buffer(data, data_size);
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

    mpix::allreduce(owners, MPI_SUM, comm);
    mpix::allreduce(offsets, MPI_SUM, comm);
    mpix::allreduce(sizes, MPI_SUM, comm);
    //mpix::print_sync("owners", vec_to_string(owners), comm);
    //mpix::print_sync("offsets", vec_to_string(offsets), comm);
    //mpix::print_sync("sizes", vec_to_string(sizes), comm);

    /* ARMCI is screwey. If the MPI ranks are on the same SMP node, the
     * addresses returned in the ptr_arr are different than the addresses on
     * the procs (due to mmapping the shared memory on each process
     * differently) so we can't use the addresses as we have them. We must make
     * everything relative to each process' own ptr_arr. */
    for (size_t i=0; i<global_count; ++i) {
        addresses[i] = &ptr_arr[owners[i]][offsets[i]];
    }
    //mpix::print_sync("addresses", vec_to_string(addresses), comm);
}


}; /* namespace pgraph */

