/**
 * @file SequenceDatabaseGArray.cpp
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 */
#include "config.h"

#include <cassert>
#include <cstddef>
#include <exception>
#include <iostream>
#include <numeric>
#include <string>

#include "mpix.hpp"
#include "Sequence.hpp"
#include "SequenceDatabase.hpp"
#include "SequenceDatabaseGArray.hpp"
#include <tascel.h>

using std::cerr;
using std::endl;
using std::exception;
using std::size_t;
using std::string;
using namespace tascel;

namespace pgraph {


SequenceDatabaseGArray::SequenceDatabaseGArray(
        const string &file_name,
        size_t budget,
        char delimiter)
    : SequenceDatabase()
    , comm(MPI_COMM_WORLD)
    , comm_size(0)
    , comm_rank(0)
    , file_size(0)
    , read_size(0)
    , is_replicated(false)
    , budget(budget)
    , delimiter(delimiter)
    , local_data(NULL)
    , local_cache()
    , ga(NULL)
    , counts()
    , count_total(0)
    , offsets()
    , remote_cache()
    , mutex()
    , max_seq_size(0)
    , tbl(NULL)
    , del(NULL)
    , ins(NULL)
{
    int err = MPI_SUCCESS;

    err = MPI_Comm_size(comm, &comm_size);
    assert(MPI_SUCCESS == err);
    err = MPI_Comm_rank(comm, &comm_rank);
    assert(MPI_SUCCESS == err);

    read_and_parse_fasta(file_name);

    assert(NUM_WORKERS > 0);
    tbl = new cell_t**[NUM_WORKERS];
    del = new int**[NUM_WORKERS];
    ins = new int**[NUM_WORKERS];
    for (int worker=0; worker<NUM_WORKERS; ++worker) {
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


void SequenceDatabaseGArray::read_and_parse_fasta(const string &file_name)
{
    MPI_Offset budget = MPI_Offset(this->budget);

    /* determine file size and create GArray */
    file_size = mpix_get_file_size(file_name);

    mpix_print_zero("memory budget", budget, comm);
    mpix_print_zero("file size", file_size, comm);

    if (budget >= file_size) {
        is_replicated = true;
        read_and_parse_fasta_himem(file_name);
    }
    else {
        read_and_parse_fasta_lomem(file_name);
    }
}


void SequenceDatabaseGArray::read_and_parse_fasta_himem(const string &file_name)
{
    MPI_File fh = MPI_FILE_NULL;
    int err = MPI_SUCCESS;
    size_t new_size = 0;

    err = MPI_File_open(MPI_COMM_SELF,
            const_cast<char *>(file_name.c_str()),
            MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    assert(MPI_SUCCESS == err);

    /* read directly into a local buffer */
    local_data = new char[file_size+2]; /* +2 for last delim and null */
    local_data[file_size] = delimiter;
    local_data[file_size+1] = '\0';

    /* everyone reads in their part */
    err = MPI_File_read_at_all(fh, 0, local_data, file_size,
            MPI_CHAR, MPI_STATUS_IGNORE);
    MPI_CHECK_IERR(err, comm_rank, comm);

    /* everyone can close the file now */
    err = MPI_File_close(&fh);
    MPI_CHECK_IERR(err, comm_rank, comm);

    /* pack and index the fasta buffer */
    pack_and_index_fasta(local_data, file_size, delimiter, 0, new_size);
    local_data[new_size] = delimiter;
    local_data[new_size+1] = '\0';
    assert(!local_cache.empty());
    assert(local_cache.size() > 0);

    count_total = local_cache.size();
    mpix_print_sync("global_size", count_total, comm);
}


void SequenceDatabaseGArray::read_and_parse_fasta_lomem(const string &file_name)
{
    int err = MPI_SUCCESS;
    MPI_File fh = MPI_FILE_NULL;
    char *file_buffer = NULL;
    Array<Offset,1> dims;
    Array<Offset,1> blocks;
    Array<int,1> proc_grid;
    ProcDistribution distribution = kBlocked;
    Order order = kRowMajor;
    MPI_Offset lo = 0;
    MPI_Offset hi = 0;
    Box<1> lbox;
    Box<1> gbox;
    GArrayRequest *req1 = GArrayRequest::construct();
    GArrayRequest *req2 = GArrayRequest::construct();
    Dispatcher<NullMutex> dispatcher;

    dims[0] = file_size;
    blocks[0] = (file_size + comm_size - 1) / comm_size;
    proc_grid[0] = comm_size;
    ga = new GArray<1>(dims, blocks, proc_grid, distribution, order, sizeof(char));

    /* determine the [lo,hi) range owned by this process */
    lo = comm_rank * blocks[0];
    hi = lo + blocks[0]; /* exclusive hi, inclusive lo */
    if (comm_rank == comm_size-1) {
        hi = file_size; /* last rank reads to end of file */
    }
    assert (lo >= 0);
    assert (hi <= file_size);
    read_size = hi-lo;

    /* open file, read the piece we own, and close file */
    err = MPI_File_open(comm, const_cast<char *>(file_name.c_str()),
            MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    assert(MPI_SUCCESS == err);
    file_buffer = new char[read_size];
    err = MPI_File_read_at_all(fh, lo, file_buffer, read_size, MPI_CHAR, MPI_STATUS_IGNORE);
    assert(MPI_SUCCESS == err);
    err = MPI_File_close(&fh);
    assert(MPI_SUCCESS == err);

    /* GArray doesn't have ::access(), so we must ::put() data to it. */
    lbox.lo()[0] = 0;
    lbox.hi()[0] = read_size;
    gbox.lo()[0] = lo;
    gbox.hi()[0] = hi;
    ga->put(gbox, file_buffer, lbox, lbox, req1, req2);
    dispatcher.registerCodelet(req1);
    dispatcher.registerCodelet(req2);
    while (!dispatcher.empty()) {
        Codelet* codelet;
        if ((codelet = dispatcher.progress()) != NULL) {
            codelet->execute();
        }
        {
            AmListenObjCodelet<NullMutex>* lcodelet;
            if((lcodelet=theAm().amListeners[0]->progress()) != NULL) {
                lcodelet->execute();
            }
        }
    }
    //ga->dumpPrint<char>();
    delete req1;
    delete req2;

    /* parse file_buffer and count '>' characters */
    counts.resize(comm_size);
    for (int i=0; i<comm_size; ++i) {
        counts[i] = 0;
    }
    for (MPI_Offset i=0; i<read_size; ++i) {
        if (file_buffer[i] == '>') {
            ++counts[comm_rank];
        }
    }
    mpix_allreduce(counts, MPI_SUM);
    count_total = std::accumulate(counts.begin(), counts.end(), 0);

    /* determine the first '>' index */
    long first_sequence_id = 0;
    for (int i=0; i<comm_rank; ++i) {
        first_sequence_id += counts[i];
    }

    /* reparse file_buffer to mark sequence offsets */
    long sequence_id = first_sequence_id;
    offsets.resize(count_total);
    for (int i=0; i<comm_size; ++i) {
        offsets[i] = 0;
    }
    for (MPI_Offset i=0; i<read_size; ++i) {
        if (file_buffer[i] == '>') {
            offsets[sequence_id++] = i+lo;
        }
    }
    mpix_allreduce(offsets, MPI_SUM);

    for (size_t i=1; i<offsets.size(); ++i) {
        long size = offsets[i] - offsets[i-1];
        if (size > max_seq_size) {
            max_seq_size = size;
        }
    }
    if ((file_size - offsets.back()) > max_seq_size) {
        max_seq_size = file_size - offsets.back();
    }

    delete [] file_buffer;
}


void SequenceDatabaseGArray::pack_and_index_fasta(char *buffer,
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
                        long l = local_cache[id-1]->get_sequence_length();
                        max_seq_size = l > max_seq_size ? l : max_seq_size;
                        count_total += l;
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
                    long l = local_cache[id-1]->get_sequence_length();
                    max_seq_size = l > max_seq_size ? l : max_seq_size;
                    count_total += l;
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
}


SequenceDatabaseGArray::~SequenceDatabaseGArray()
{
    map<size_t, Sequence*>::iterator it;

    for (it=remote_cache.begin(); it!=remote_cache.end(); ++it) {
        delete it->second;
    }

    for (int worker=0; worker<NUM_WORKERS; ++worker) {
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

    if (is_replicated) {
        delete [] local_data;
    }
    else {
        delete ga;
    }
}


size_t SequenceDatabaseGArray::get_local_count() const
{
    return counts[comm_rank];
}


size_t SequenceDatabaseGArray::get_local_size() const
{
    return read_size;
}


size_t SequenceDatabaseGArray::get_global_count() const
{
    return count_total;
}


size_t SequenceDatabaseGArray::get_global_replica_count() const
{
    return 1;
}


size_t SequenceDatabaseGArray::get_global_replica_index() const
{
    return 0;
}


size_t SequenceDatabaseGArray::get_global_size() const
{
    return file_size;
}


Sequence &SequenceDatabaseGArray::get_sequence(size_t i)
{
    if (is_replicated) {
        return *local_cache[i];
    }
    else if (remote_cache.count(i)) {
        LockGuard<PthreadMutex> guard(mutex);
        return *remote_cache[i];
    }
    else {
        Sequence *sequence = NULL;
        char *buffer = NULL;
        size_t size = 0;
        Box<1> lbox;
        Box<1> gbox;
        GArrayRequest *req1 = GArrayRequest::construct();
        GArrayRequest *req2 = GArrayRequest::construct();
        Dispatcher<NullMutex> dispatcher;
        bool last_id = (long(i) == count_total-1);

        size = last_id ? file_size - offsets[i] : offsets[i+1] - offsets[i];
        buffer = new char[size+1];
        (void)memset(buffer, 0, size+1);

        /* fetch the sequence data */
        lbox.lo()[0] = 0;
        lbox.hi()[0] = size - 1;
        gbox.lo()[0] = offsets[i];
        gbox.hi()[0] = last_id ? file_size - 1 : offsets[i+1] - 1;
        ga->get(gbox, buffer, lbox, lbox, req1, req2);
        dispatcher.registerCodelet(req1);
        dispatcher.registerCodelet(req2);
        while (!dispatcher.empty()) {
            Codelet* codelet;
            if ((codelet = dispatcher.progress()) != NULL) {
                codelet->execute();
            }
            {
                AmListenObjCodelet<NullMutex>* lcodelet;
                if((lcodelet=theAm().amListeners[0]->progress()) != NULL) {
                    lcodelet->execute();
                }
            }
        }
        delete req1;
        delete req2;

        /* pack fasta */
        {
            size_t r = 0;           /* read index */
            size_t w = 0;           /* write index */
            size_t last_gt = 0;     /* index of last '>' seen */
            size_t last_hash = 0;   /* index of last '#' inserted */

            /* We skip all newlines unless it is the newline that terminates
             * the sequence ID. We replace the ID-terminating newline with '#'.
             * We replace the sequence terminating newline with the given
             * delimiter. */
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
                            size_t sequence_length = w - last_hash - 2;
                            if (delimiter == '\0') {
                                sequence_length -= 1;
                            }
                            sequence =  new Sequence(&buffer[last_gt],
                                        id_offset,
                                        id_length,
                                        sequence_offset,
                                        sequence_length,
                                        true);
                        }
                    }
                    /* or file wasn't terminated with a newline */
                    else {
                        /* hopefully space enough to append the delimiter anyway */
                        buffer[w++] = delimiter;
                        size_t id_offset = 0;
                        size_t id_length = last_hash - last_gt;
                        size_t sequence_offset = last_hash - last_gt + 1;
                        size_t sequence_length = w - last_hash - 2;
                        if (delimiter == '\0') {
                            sequence_length -= 1;
                        }
                        sequence =  new Sequence(&buffer[last_gt],
                                id_offset,
                                id_length,
                                sequence_offset,
                                sequence_length,
                                true);
                    }
                    r++;
                }
            }
        }
        {
            LockGuard<PthreadMutex> guard(mutex);
            if (remote_cache.count(i)) {
                /* another thread already added */
                delete sequence;
                sequence = remote_cache[i];
            }
            else {
                /* we are first to add */
                remote_cache[i] = sequence;
            }
        }
        return *sequence;
    }
}


Sequence &SequenceDatabaseGArray::operator[](size_t i)
{
    return this->get_sequence(i);
}


void SequenceDatabaseGArray::set_num_threads(size_t num)
{
}


size_t SequenceDatabaseGArray::get_max_length() const
{
    return 0;
}


void SequenceDatabaseGArray::align(size_t i, size_t j, int &score, int &ndig, int &alen, int open, int gap, int tid)
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


void SequenceDatabaseGArray::align_ssw(size_t i, size_t j, int &score, int &ndig, int &align, int open, int gap, int tid)
{
}


bool SequenceDatabaseGArray::is_edge(size_t i, size_t j, const int &score, const int &ndig, const int &alen, const int &AOL, const int &SIM, const int &OS, int &sscore, size_t &max_len)
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

};

