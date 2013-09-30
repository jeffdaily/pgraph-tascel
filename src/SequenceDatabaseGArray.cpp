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
#include <string>

#include "mpix.hpp"
#include "Sequence.hpp"
#include "SequenceDatabase.hpp"
#include "SequenceDatabaseGArray.hpp"
#include <tascel.h>

using std::exception;
using std::size_t;
using std::string;
using namespace tascel;

namespace pgraph {


SequenceDatabaseGArray::SequenceDatabaseGArray(const string &filename)
    : SequenceDatabase()
    , file_name(filename)
    , ga(NULL)
{
    Array<Offset,1> dims;
    Array<Offset,1> blocks;
    Array<int,1> proc_grid;
    ProcDistribution distribution = kBlocked;
    Order order = kRowMajor;
    MPI_Offset file_size = mpix_get_file_size(file_name);
    MPI_Offset lo = 0;
    MPI_Offset hi = 0;
    MPI_Offset read_size = 0;
    MPI_Comm comm = MPI_COMM_WORLD;
    int comm_size = 0;
    int comm_rank = 0;
    int err = MPI_SUCCESS;
    MPI_File fh;
    char *file_buffer = NULL;
    Box<1> lbox;
    Box<1> gbox;
    GArrayRequest *req1 = GArrayRequest::construct();
    GArrayRequest *req2 = GArrayRequest::construct();
    Dispatcher<NullMutex> dispatcher;

    err = MPI_Comm_size(comm, &comm_size);
    assert(MPI_SUCCESS == err);
    err = MPI_Comm_rank(comm, &comm_rank);
    assert(MPI_SUCCESS == err);

    dims[0] = file_size;
    blocks[0] = (file_size + comm_size - 1) / comm_size;
    proc_grid[0] = comm_size;
    ga = new GArray<1>(dims, blocks, proc_grid, distribution, order, sizeof(char));

    /* GArray doesn't have ::access(), so we must ::put() data to it. */

    lo = comm_rank * blocks[0];
    hi = lo + blocks[0]; /* exclusive hi, inclusive lo */
    if (comm_rank == comm_size-1) {
        hi = file_size; /* last rank reads to end of file */
    }
    assert (lo >= 0);
    assert (hi <= file_size);
    read_size = hi-lo;

    err = MPI_File_open(comm, const_cast<char *>(file_name.c_str()),
            MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    assert(MPI_SUCCESS == err);

    file_buffer = new char[read_size];
    err = MPI_File_read_at_all(fh, lo, file_buffer, read_size, MPI_CHAR, MPI_STATUS_IGNORE);
    assert(MPI_SUCCESS == err);

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
    ga->dumpPrint<char>();

    delete req1;
    delete req2;
    delete [] file_buffer;
}


SequenceDatabaseGArray::~SequenceDatabaseGArray()
{
    delete ga;
}


size_t SequenceDatabaseGArray::get_local_count() const
{
    return 0;
}


size_t SequenceDatabaseGArray::get_local_size() const
{
    return 0;
}


size_t SequenceDatabaseGArray::get_global_count() const
{
    return 0;
}


size_t SequenceDatabaseGArray::get_global_replica_count() const
{
    return 0;
}


size_t SequenceDatabaseGArray::get_global_replica_index() const
{
    return 0;
}


size_t SequenceDatabaseGArray::get_global_size() const
{
    return 0;
}


Sequence &SequenceDatabaseGArray::get_sequence(size_t i)
{
}


Sequence &SequenceDatabaseGArray::operator[](size_t i)
{
}


void SequenceDatabaseGArray::set_num_threads(size_t num)
{
}


size_t SequenceDatabaseGArray::get_max_length() const
{
    return 0;
}


void SequenceDatabaseGArray::align(size_t i, size_t j, int &score, int &ndig, int &align, int open, int gap, int tid)
{
}


void SequenceDatabaseGArray::align_ssw(size_t i, size_t j, int &score, int &ndig, int &align, int open, int gap, int tid)
{
}


bool SequenceDatabaseGArray::is_edge(size_t i, size_t j, const int &score, const int &ndig, const int &align, const int &AOL, const int &SIM, const int &OS, int &sscore, size_t &max_len)
{
    return false;
}

};

