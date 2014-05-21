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
#include <cstddef>
#include <iostream>
#include <string>

#include "alignment.hpp"
#include "mpix.hpp"
#include "Sequence.hpp"
#include "SequenceDatabase.hpp"
#include "SequenceDatabaseReplicated.hpp"

using ::std::cerr;
using ::std::cout;
using ::std::endl;
using ::std::size_t;
using ::std::string;

namespace pgraph {


SequenceDatabaseReplicated::SequenceDatabaseReplicated(
        const string &file_name,
        size_t budget,
        MPI_Comm comm,
        char delimiter)
    :   SequenceDatabase(delimiter)
    ,   comm(comm)
    ,   comm_rank(0)
    ,   comm_size(0)
    ,   file_name(file_name)
    ,   local_data(NULL)
    ,   local_cache()
    ,   _longest(0)
    ,   _char_size(0)
{
    MPI_Offset file_size = 0;
    MPI_Offset file_size_out = 0;
    size_t new_size = 0;

    /* rank and size */
    comm_rank = mpix::comm_rank(comm);
    comm_size = mpix::comm_size(comm);
    file_size = mpix::get_file_size(file_name, comm);

    if (0 == comm_rank) {
        cout << "sequence file size is " << file_size << endl;
    }
    assert(MPI_Offset(budget) >= file_size);

    /* read directly into a local buffer */
    local_data = new char[file_size+2]; /* +2 for last delim and null */
    local_data[file_size] = delimiter;
    local_data[file_size+1] = '\0';

    mpix::read_file(file_name, local_data, file_size_out, comm);
    assert(file_size == file_size_out);

    /* pack and index the fasta buffer */
    pack_and_index_fasta(local_data, file_size, 0, new_size);
    local_data[new_size] = delimiter;
    local_data[new_size+1] = '\0';
    assert(!local_cache.empty());
    assert(local_cache.size() > 0);

    if (0 == comm_rank) {
        cout << "longest sequence in file is " << _longest << endl;
    }
}


SequenceDatabaseReplicated::~SequenceDatabaseReplicated()
{
    vector<Sequence*>::iterator it;

    for (it=local_cache.begin(); it!=local_cache.end(); ++it) {
        delete *it;
    }

    delete [] local_data;
}


void SequenceDatabaseReplicated::pack_and_index_fasta(char *buffer,
                                            size_t size,
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
            /* either
             * 1) file wasn't terminated with a newline, or
             * 2) newline was the last character, or
             * 3) we have additional characters and we care about '>' */
            if (r == size || r+1 == size || buffer[r+1] == '>') {
                /* in case of (1) or (2), hopefully space for delimiter anyway */
                buffer[w++] = delimiter;
                size_t id_offset = 0;
                size_t id_length = last_hash - last_gt;
                size_t sequence_offset = id_length + 1;
                size_t sequence_length = w - last_hash - 1;
                if (delimiter == '\0') {
                    sequence_length -= 1;
                }
                Sequence *sequence = new Sequence(&buffer[last_gt],
                            id_offset,
                            id_length,
                            sequence_offset,
                            sequence_length);
                sequence->uses_delimiter(delimiter != '\0');
                local_cache.push_back(sequence);
                id++;
                assert(!local_cache.empty());
                size_t l = local_cache[id-1]->get_sequence_length();
                _longest = l > _longest ? l : _longest;
                _char_size += l;
            }
            r++;
        }
    }

    assert(!local_cache.empty());

    new_size = w;
}

}; /* namespace pgraph */

