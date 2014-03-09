/**
 * @file mpix.cpp
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 */
#include "config.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <cassert>
#include <cstddef>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#include <mpi.h>

#include "mpix.hpp"

using std::cout;
using std::endl;
using std::size_t;
using std::string;
using std::vector;


void mpix::print_sync(
        const string &name, const vector<string> &what, MPI_Comm comm)
{
    int rank = mpix::rank(comm);
    int size = mpix::size(comm);

    for (int i = 0; i < size; ++i) {
        if (i == rank) {
            size_t j;
            for (j = 0; j < what.size(); ++j) {
                cout << "[" << rank << "] "
                     << name << "[" << j << "]="
                     << what.at(j) << endl;
            }
        }
        MPI_Barrier(comm);
    }
}


void mpix::print_sync(const string &name, const string &what, MPI_Comm comm)
{
    int rank = mpix::rank(comm);
    int size = mpix::size(comm);

    for (int i = 0; i < size; ++i) {
        if (i == rank) {
            cout << "[" << rank << "] " << name << "=" << what << endl;
        }
        MPI_Barrier(comm);
    }
}


void mpix::print_sync(const string &what, MPI_Comm comm)
{
    int rank = mpix::rank(comm);
    int size = mpix::size(comm);

    for (int i = 0; i < size; ++i) {
        if (i == rank) {
            cout << "[" << rank << "] " << what << endl;
        }
        MPI_Barrier(comm);
    }
}


void mpix::print_zero(
    const string &name, const vector<string> &what, MPI_Comm comm)
{
    int rank = mpix::rank(comm);

    if (0 == rank) {
        size_t j;
        for (j = 0; j < what.size(); ++j) {
            cout << name << "[" << j << "]=" << what.at(j) << endl;
        }
    }
    MPI_Barrier(comm);
}


void mpix::print_zero(const string &name, const string &what, MPI_Comm comm)
{
    int rank = mpix::rank(comm);

    if (0 == rank) {
        cout << "[" << rank << "] " << name << "=" << what << endl;
    }
    MPI_Barrier(comm);
}


void mpix::print_zero(const string &what, MPI_Comm comm)
{
    int rank = mpix::rank(comm);

    if (0 == rank) {
        cout << "[" << rank << "] " << what << endl;
    }
    MPI_Barrier(comm);
}


/**
 * Collectively return size of the given file.
 *
 * All processes within the MPI_Comm instance will receive the file size.
 *
 * @param[in] comm instance
 * @param[in] file_name to open
 *
 * @return the file size
 */
MPI_Offset mpix::get_file_size(const string &file_name, MPI_Comm comm)
{
    MPI_Offset file_size = 0;
    int rank = mpix::rank(comm);

    /* process 0 stats the file to determine its size */
    if (0 == rank) {
        struct stat statbuf;
        int retval;

        retval = stat(file_name.c_str(), &statbuf);
        if (-1 == retval) {
            perror("stat");
            printf("unable to stat file '%s'\n", file_name.c_str());
            MPI_Abort(comm, 1);
        }
        file_size = statbuf.st_size;
    }

    /* the file_size is broadcast to all */
    mpix::bcast(file_size, 0, comm);

    return file_size;
}


/**
 * Collectively read a file.
 *
 * All processes within the MPI_Comm instance will receive the entire contents
 * of the file.
 *
 * @param[in] file_name to open
 * @param[out] file_buffer to store file contents
 * @param[out] file_size of the given file
 * @param[out] chunk_size max size to read at one time
 * @param[in] comm instance
 */
void mpix::read_file(
    const string &file_name, char *&file_buffer,
    MPI_Offset &file_size, long chunk_size, MPI_Comm comm)
{
    mpix::read_file_mpiio(file_name, file_buffer, file_size, chunk_size, comm);
    //mpix::read_file_bcast(file_name, file_buffer, file_size, chunk_size, comm);
}


/**
 * Collectively read an entire file using MPI_File_read_all().
 *
 * All processes within the MPI_Comm instance will receive the entire contents
 * of the file.
 *
 * @param[in] file_name to open
 * @param[out] file_buffer to store file contents
 * @param[out] file_size of the given file
 * @param[out] chunk_size max size to read at one time
 * @param[in] comm instance
 */
void mpix::read_file_mpiio(
    const string &file_name, char *&file_buffer,
    MPI_Offset &file_size, long chunk_size, MPI_Comm comm)
{
    int rank = mpix::rank(comm);
    MPI_Status status;
    MPI_File fh;

    /* allocate a buffer for the file, of the entire size */
    file_size = mpix::get_file_size(file_name, comm);
    if (NULL == file_buffer) {
        file_buffer = new char[file_size];
    }

    if (file_size > chunk_size) {
        long offset = 0;
        int count = 0;

        /* read file contents in chunks */
        MPI_CHECK(MPI_File_open(comm, const_cast<char *>(file_name.c_str()),
                                MPI_MODE_RDONLY | MPI_MODE_UNIQUE_OPEN,
                                MPI_INFO_NULL, &fh));
        while (offset < file_size) {
            long message_size = chunk_size;
            if (offset + chunk_size > file_size) {
                message_size = file_size % chunk_size;
            }
            if (0 == rank) {
                printf("reading chunk %ld->%ld message_size=%ld\n",
                       offset, offset + message_size, message_size);
            }
            MPI_CHECK(MPI_File_read_at_all(fh, offset, file_buffer + offset,
                                           message_size, MPI_CHAR, &status));
            MPI_CHECK(MPI_Get_count(&status, MPI_CHAR, &count));
            assert(count == message_size);
            offset += count;
        }
    }
    else {
        int count = 0;

        /* all procs read the entire file */
        MPI_CHECK(MPI_File_open(comm, const_cast<char *>(file_name.c_str()),
                                /*MPI_MODE_RDONLY|MPI_MODE_UNIQUE_OPEN,*/
                                MPI_MODE_RDONLY,
                                MPI_INFO_NULL, &fh));
        MPI_CHECK(MPI_File_read_at_all(fh, 0, file_buffer,
                                       file_size, MPI_CHAR, &status));
        MPI_CHECK(MPI_Get_count(&status, MPI_CHAR, &count));
        assert(count == file_size);
        MPI_CHECK(MPI_File_close(&fh));
    }
}


/**
 * Collectively read a file using process 0 and bcast.
 *
 * All processes within the MPI_Comm instance will receive the entire contents
 * of the file.
 *
 * @param[in] file_name to open
 * @param[out] file_buffer to store file contents
 * @param[out] file_size of the given file
 * @param[out] chunk_size max size to read at one time
 * @param[in] comm instance
 */
void mpix::read_file_bcast(
    const string &file_name, char *&file_buffer,
    MPI_Offset &file_size, long chunk_size, MPI_Comm comm)
{
    int rank = mpix::rank(comm);

    /* allocate a buffer for the file, of the entire size */
    file_size = mpix::get_file_size(file_name, comm);
    if (NULL == file_buffer) {
        file_buffer = new char[file_size];
    }

    if (0 == rank) {
        FILE *file = NULL;
        size_t read_count = 0;

        file = fopen(file_name.c_str(), "r");
        if (NULL == file) {
            perror("fopen");
            printf("unable to open file on process 0\n");
            MPI_Abort(comm, 1);
        }

        read_count = fread(file_buffer, file_size, 1, file);
        if (0 == read_count) {
            printf("unable to read file on process 0\n");
            MPI_Abort(comm, 1);
        }

    }

    if (file_size > chunk_size) {
        /* bcast file contents in chunks */
        long offset = 0;
        while (offset < file_size) {
            long message_size = chunk_size;
            if (offset + chunk_size > file_size) {
                message_size = file_size % chunk_size;
            }
            if (0 == rank) {
                printf("broadcasting chunk %ld->%ld message_size=%ld\n",
                       offset, offset + message_size, message_size);
            }
            MPI_CHECK(MPI_Bcast(&file_buffer[offset],
                                message_size, MPI_CHAR, 0, comm));
            offset += chunk_size;
        }
    }
    else {
        /* bcast entire file contents */
        MPI_CHECK(MPI_Bcast(file_buffer, file_size, MPI_CHAR, 0, comm));
    }
}
