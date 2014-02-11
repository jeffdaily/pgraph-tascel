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


void mpix_bcast(string &object, int root, MPI_Comm comm)
{
    int rank = mpix_rank(comm);
    int size = 0;

    if (rank == root) {
        size = int(object.size());
        mpix_bcast(size, root, comm);
        MPI_CHECK_C(MPI_Bcast(const_cast<char*>(object.data()),
                    size, MPI_CHAR, root, comm));
    }
    else {
        char *data = NULL;

        mpix_bcast(size, root, comm);
        data = new char[size];
        MPI_CHECK_C(MPI_Bcast(data, size, MPI_CHAR, root, comm));
        object.assign(data, size);
        delete [] data;
    }
}


void mpix_bcast(vector<string> &object, int root, MPI_Comm comm)
{
    int rank = mpix_rank(comm);
    int size = 0;

    size = object.size();
    mpix_bcast(size, root, comm);

    if (rank != root) {
        object.resize(size);
    }

    for (vector<string>::iterator it=object.begin(); it!=object.end(); ++it) {
        mpix_bcast(*it, root, comm);
    }
}


/* MPI standard does not guarantee all procs receive argc and argv */
void mpix_bcast_argv(
        int argc, char **argv, vector<string> &all_argv, MPI_Comm comm)
{
    int rank = mpix_rank(comm);;

    if (0 == rank) {
        MPI_CHECK(MPI_Bcast(&argc, 1, MPI_INT, 0, comm));
        for (int i = 0; i < argc; ++i) {
            int length = strlen(argv[i]) + 1;
            MPI_CHECK(MPI_Bcast(&length, 1, MPI_INT, 0, comm));
            MPI_CHECK(MPI_Bcast(argv[i], length, MPI_CHAR, 0, comm));
            all_argv.push_back(argv[i]);
        }
    }
    else {
        int all_argc;
        MPI_CHECK(MPI_Bcast(&all_argc, 1, MPI_INT, 0, comm));
        for (int i = 0; i < all_argc; ++i) {
            int length;
            char buffer[ARG_LEN_MAX];
            MPI_CHECK(MPI_Bcast(&length, 1, MPI_INT, 0, comm));
            MPI_CHECK(MPI_Bcast(buffer, length, MPI_CHAR, 0, comm));
            all_argv.push_back(buffer);
        }
    }
}


void mpix_print_sync(
        const string &name, const vector<string> &what, MPI_Comm comm)
{
    int rank = mpix_rank(comm);
    int size = mpix_size(comm);

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


void mpix_print_sync(const string &name, const string &what, MPI_Comm comm)
{
    int rank = mpix_rank(comm);
    int size = mpix_size(comm);

    for (int i = 0; i < size; ++i) {
        if (i == rank) {
            cout << "[" << rank << "] " << name << "=" << what << endl;
        }
        MPI_Barrier(comm);
    }
}


void mpix_print_sync(const string &what, MPI_Comm comm)
{
    int rank = mpix_rank(comm);
    int size = mpix_size(comm);

    for (int i = 0; i < size; ++i) {
        if (i == rank) {
            cout << "[" << rank << "] " << what << endl;
        }
        MPI_Barrier(comm);
    }
}


void mpix_print_zero(
    const string &name, const vector<string> &what, MPI_Comm comm)
{
    int rank = mpix_rank(comm);

    if (0 == rank) {
        size_t j;
        for (j = 0; j < what.size(); ++j) {
            cout << name << "[" << j << "]=" << what.at(j) << endl;
        }
    }
    MPI_Barrier(comm);
}


void mpix_print_zero(const string &name, const string &what, MPI_Comm comm)
{
    int rank = mpix_rank(comm);

    if (0 == rank) {
        cout << "[" << rank << "] " << name << "=" << what << endl;
    }
    MPI_Barrier(comm);
}


void mpix_print_zero(const string &what, MPI_Comm comm)
{
    int rank = mpix_rank(comm);

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
MPI_Offset mpix_get_file_size(const string &file_name, MPI_Comm comm)
{
    MPI_Offset file_size = 0;
    int rank = mpix_rank(comm);

    /* process 0 stats the file to determine its size */
    if (0 == rank) {
        struct stat statbuf;
        int retval;

        retval = stat(file_name.c_str(), &statbuf);
        if (-1 == retval) {
            perror("stat");
            printf("unable to stat sequence file\n");
            MPI_Abort(comm, 1);
        }
        file_size = statbuf.st_size;
    }

    /* the file_size is broadcast to all */
    mpix_bcast(file_size, 0, comm);
#if 0
    if (0 == rank) {
        printf("file_size=%ld\n", file_size);
    }
#endif
    mpix_print_sync("file_size", file_size, comm);

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
void mpix_read_file(
    const string &file_name, char *&file_buffer,
    MPI_Offset &file_size, long chunk_size, MPI_Comm comm)
{
    mpix_read_file_mpiio(file_name, file_buffer, file_size, chunk_size, comm);
    //mpix_read_file_bcast(file_name, file_buffer, file_size, chunk_size, comm);
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
void mpix_read_file_mpiio(
    const string &file_name, char *&file_buffer,
    MPI_Offset &file_size, long chunk_size, MPI_Comm comm)
{
    int rank = mpix_rank(comm);
    MPI_Status status;
    MPI_File fh;

    /* allocate a buffer for the file, of the entire size */
    file_size = mpix_get_file_size(file_name, comm);
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
void mpix_read_file_bcast(
    const string &file_name, char *&file_buffer,
    MPI_Offset &file_size, long chunk_size, MPI_Comm comm)
{
    int rank = mpix_rank(comm);

    /* allocate a buffer for the file, of the entire size */
    file_size = mpix_get_file_size(file_name, comm);
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
