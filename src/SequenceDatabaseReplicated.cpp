/**
 * @file SequenceDatabaseReplicated.cpp
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

#include <tascel.h>

#include "alignment.hpp"
#include "csequence.h"
#include "Sequence.hpp"
#include "SequenceDatabase.hpp"
#include "SequenceDatabaseReplicated.hpp"
#include "mpix.hpp"

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


SequenceDatabaseReplicated::SequenceDatabaseReplicated(
        const string &file_name,
        size_t budget,
        MPI_Comm comm,
        size_t num_threads,
        char delimiter)
    :   SequenceDatabase()
    ,   comm_orig(comm)
    ,   comm_orig_rank(0)
    ,   comm_orig_size(0)
    ,   comm(MPI_COMM_NULL)
    ,   comm_rank(0)
    ,   comm_size(0)
    ,   budget(budget)
    ,   num_threads(num_threads)
    ,   delimiter(delimiter)
    ,   file_name(file_name)
    ,   local_data(NULL)
    ,   local_cache()
    ,   local_size(0)
    ,   max_seq_size(0)
    ,   tbl(NULL)
    ,   del(NULL)
    ,   ins(NULL)
{
    int ierr = 0;

    /* rank and size */
    ierr = MPI_Comm_rank(comm_orig, &comm_orig_rank);
    MPI_CHECK_IERR(ierr, comm_orig_rank, comm_orig);
    ierr = MPI_Comm_size(comm_orig, &comm_orig_size);
    MPI_CHECK_IERR(ierr, comm_orig_rank, comm_orig);

    read_and_parse_fasta();
}


SequenceDatabaseReplicated::~SequenceDatabaseReplicated()
{
    map<size_t, Sequence*>::iterator it;

    for (it=local_cache.begin(); it!=local_cache.end(); ++it) {
        delete it->second;
    }

    for (size_t worker=0; worker<num_threads; ++worker) {
#if 1
        delete [] tbl[worker][1];
        delete [] tbl[worker][0];
        delete [] tbl[worker];
        delete [] del[worker][1];
        delete [] del[worker][0];
        delete [] del[worker];
        delete [] ins[worker][1];
        delete [] ins[worker][0];
        delete [] ins[worker];
#else
        free_cell_table(2, tbl[worker]);
        free_int_table(2, del[worker]);
        free_int_table(2, ins[worker]);
#endif
    }
    delete [] tbl;
    delete [] del;
    delete [] ins;
    delete [] local_data;
}


size_t SequenceDatabaseReplicated::get_local_count() const
{
    return local_cache.size();
}


size_t SequenceDatabaseReplicated::get_local_size() const
{
    return local_size;
}


size_t SequenceDatabaseReplicated::get_global_count() const
{
    return get_local_count();
}


size_t SequenceDatabaseReplicated::get_global_size() const
{
    return get_local_size();
}


void SequenceDatabaseReplicated::read_and_parse_fasta()
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
        MPI_CHECK_IERR(ierr, comm_orig_rank, comm_orig);
        ierr = MPI_File_get_size(in, &file_size);
        MPI_CHECK_IERR(ierr, comm_orig_rank, comm_orig);
        ierr = MPI_File_close(&in);
        MPI_CHECK_IERR(ierr, comm_orig_rank, comm_orig);
    }
    mpix_bcast(file_size, 0, comm_orig);

    mpix_print_zero("memory budget", budget, comm_orig);
    mpix_print_zero("file size", file_size, comm_orig);
    
    if (budget >= file_size) {
        comm = comm_orig;
        comm_rank = comm_orig_rank;
        comm_size = comm_orig_size;
        ierr = MPI_File_open(MPI_COMM_SELF,
                const_cast<char *>(file_name.c_str()),
                MPI_MODE_RDONLY, MPI_INFO_NULL, &in);
        MPI_CHECK_IERR(ierr, comm_orig_rank, comm_orig);
        read_and_parse_fasta_himem(in, file_size);
    }
    else {
        assert(0);
    }
}


void SequenceDatabaseReplicated::read_and_parse_fasta_himem(MPI_File in,
                                                  MPI_Offset file_size)
{
    int ierr = 0;
    size_t new_size = 0;

    /* read directly into a local buffer */
    local_data = new char[file_size+2]; /* +2 for last delim and null */
    local_data[file_size] = delimiter;
    local_data[file_size+1] = '\0';
    mpix_print_zero("allocated file buffer", comm_orig);

    /* everyone reads in their part */
    ierr = MPI_File_read_at_all(in, 0, local_data, file_size,
            MPI_CHAR, MPI_STATUS_IGNORE);
    MPI_CHECK_IERR(ierr, comm_rank, comm);
    mpix_print_zero("read file", comm_orig);

    /* everyone can close the file now */
    ierr = MPI_File_close(&in);
    MPI_CHECK_IERR(ierr, comm_rank, comm);
    mpix_print_zero("closed file", comm_orig);

    /* pack and index the fasta buffer */
    pack_and_index_fasta(local_data, file_size, delimiter, 0, new_size);
    local_data[new_size] = delimiter;
    local_data[new_size+1] = '\0';
    assert(!local_cache.empty());
    assert(local_cache.size() > 0);
    mpix_print_zero("packed and indexed file", comm_orig);
}


void SequenceDatabaseReplicated::pack_and_index_fasta(char *buffer,
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

    mpix_allreduce(max_seq_size, MPI_MAX, comm);
    mpix_print_zero("max_seq_size", max_seq_size, comm);

    assert(num_threads > 0);
    tbl = new cell_t**[num_threads];
    del = new int**[num_threads];
    ins = new int**[num_threads];
    for (size_t worker=0; worker<num_threads; ++worker) {
#if 1
        size_t max_seq_size_p1 = max_seq_size + 1;
        tbl[worker] = new cell_t*[2];
        tbl[worker][0] = new cell_t[max_seq_size_p1];
        tbl[worker][1] = new cell_t[max_seq_size_p1];
        del[worker] = new int*[2];
        del[worker][0] = new int[max_seq_size_p1];
        del[worker][1] = new int[max_seq_size_p1];
        ins[worker] = new int*[2];
        ins[worker][0] = new int[max_seq_size_p1];
        ins[worker][1] = new int[max_seq_size_p1];
#if USE_MEMSET
        (void)memset(tbl[worker][0], 0, sizeof(cell_t)*max_seq_size_p1);
        (void)memset(tbl[worker][1], 0, sizeof(cell_t)*max_seq_size_p1);
        (void)memset(del[worker][0], 0, sizeof(int)*max_seq_size_p1);
        (void)memset(del[worker][1], 0, sizeof(int)*max_seq_size_p1);
        (void)memset(ins[worker][0], 0, sizeof(int)*max_seq_size_p1);
        (void)memset(ins[worker][1], 0, sizeof(int)*max_seq_size_p1);
#endif
#else
        tbl[worker] = allocate_cell_table(2, max_seq_size);
        del[worker] = allocate_int_table(2, max_seq_size);
        ins[worker] = allocate_int_table(2, max_seq_size);
#endif
    }
}



Sequence &SequenceDatabaseReplicated::get_sequence(size_t i)
{
    return *local_cache[i];
}


Sequence &SequenceDatabaseReplicated::operator[](size_t i)
{
    return this->get_sequence(i);
}


void SequenceDatabaseReplicated::set_num_threads(size_t num)
{
    assert(0); // not used yet
    assert(tbl == NULL);
    assert(del == NULL);
    assert(ins == NULL);

    tbl = new cell_t**[num];
    assert(tbl);
    del = new int**[num];
    assert(del);
    ins = new int**[num];
    assert(ins);

    for (size_t tid=0; tid<num; ++tid) {
        tbl[tid] = NULL;
        del[tid] = NULL;
        ins[tid] = NULL;
    }
}


void SequenceDatabaseReplicated::align(size_t i,
                                  size_t j,
                                  int &score,
                                  int &ndig,
                                  int &alen,
                                  int open,
                                  int gap,
                                  int tid)
{
    Sequence &s1 = get_sequence(i);
    Sequence &s2 = get_sequence(j);
    const char *c1 = NULL;
    const char *c2 = NULL;
    size_t l1 = 0;
    size_t l2 = 0;
    cell_t result;

    assert(max_seq_size > 0);

#if 0
    if (tbl[tid] == NULL) {
        tbl[tid] = new cell_t*[2];
        tbl[tid][0] = new cell_t[max_seq_size];
        tbl[tid][1] = new cell_t[max_seq_size];
        assert(del[tid] == NULL);
        del[tid] = new int*[2];
        del[tid][0] = new int[max_seq_size];
        del[tid][1] = new int[max_seq_size];
        assert(ins[tid] == NULL);
        ins[tid] = new int*[2];
        ins[tid][0] = new int[max_seq_size];
        ins[tid][1] = new int[max_seq_size];
    }
#else
    assert(tbl);
    assert(tbl[tid]);
    assert(tbl[tid][0]);
    assert(tbl[tid][1]);
    assert(del);
    assert(del[tid]);
    assert(del[tid][0]);
    assert(del[tid][1]);
    assert(ins);
    assert(ins[tid]);
    assert(ins[tid][0]);
    assert(ins[tid][1]);
#endif

    s1.get_sequence(c1,l1);
    assert(c1);
    assert(l1);
    s2.get_sequence(c2,l2);
    assert(c2);
    assert(l2);
    result = affine_gap_align_blosum(c1, l1, c2, l2,
            tbl[tid], del[tid], ins[tid], open, gap);

    score = result.score;
    ndig = result.matches;
    alen = result.length;
}


void SequenceDatabaseReplicated::align_ssw(size_t i,
                                      size_t j,
                                      int &score,
                                      int &ndig,
                                      int &alen,
                                      int open,
                                      int gap,
                                      int /*tid*/)
{
    Sequence &s1 = get_sequence(i);
    Sequence &s2 = get_sequence(j);
    const char *c1 = NULL;
    const char *c2 = NULL;
    size_t l1 = 0;
    size_t l2 = 0;
    cell_t result;

    assert(max_seq_size > 0);

    s1.get_sequence(c1,l1);
    assert(c1);
    assert(l1);
    s2.get_sequence(c2,l2);
    assert(c2);
    assert(l2);
    result = affine_gap_align_blosum_ssw(c1, l1, c2, l2, open, gap);

    score = result.score;
    ndig = result.matches;
    alen = result.length;
}


bool SequenceDatabaseReplicated::is_edge(size_t i,
                                    size_t j,
                                    const int &score,
                                    const int &ndig,
                                    const int &alen,
                                    const int &AOL,
                                    const int &SIM,
                                    const int &OS,
                                    int &sscore,
                                    size_t &max_len)
{
    Sequence &s1 = get_sequence(i);
    Sequence &s2 = get_sequence(j);
    const char *c1;
    const char *c2;
    size_t l1;
    size_t l2;
    int nmatch = 0;

    /* DO NOT need to compute the sscore value every time, if the later check
     * failed at the first step, the sscore computation is wasted */
    s1.get_sequence(c1,l1);
    s2.get_sequence(c2,l2);

    assert(l1);
    assert(l2);

    if (l1 > l2) {
        max_len = l1;
        sscore = self_score_blosum(c1,l1);
    }
    else {
        max_len = l2;
        sscore = self_score_blosum(c2,l2);
    }

    nmatch = ndig;

    /* order the condition in strict->loose way, performance perspective
     * comparison using integers, no overflow could happen */
    if (score <= 0) {
        return false;
    }
    else if ((alen * 100 >= AOL * int(max_len))
             && (nmatch * 100 >= SIM * alen)
             && (score * 100 >= OS * sscore)) {
        return true;
    }
    else {
        return false;
    }
}

}; /* namespace pgraph */

