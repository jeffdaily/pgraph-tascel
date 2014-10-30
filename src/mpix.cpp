/**
 * @file mpix-inl.hpp
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 *
 * Inline implementations of templated and/or overloaded MPI functions.
 */

#include <mpi.h>

#include <sys/stat.h>

#include <cassert>
#include <cstddef>
#include <iostream>
#include <map>
#include <numeric>
#include <sstream>
#include <string>
#include <typeinfo>
#include <utility>
#include <vector>

#include "mpix.hpp"
#include "mpix_helper.hpp"

using ::std::accumulate;
using ::std::cerr;
using ::std::cout;
using ::std::endl;
using ::std::map;
using ::std::ostringstream;
using ::std::pair;
using ::std::size_t;
using ::std::string;
using ::std::vector;

namespace mpix {

void init(int &argc, char **&argv)
{
    check(MPI_Init(&argc,&argv));
}


void init_thread(int &argc, char **&argv, int requested)
{
    int provided;
    check(MPI_Init_thread(&argc, &argv, requested, &provided));
    if (provided < requested) {
        cerr << "[" << comm_rank(MPI_COMM_WORLD) << "] MPI ERROR"
            << ": " << "MPI_Init_thread provided < requested"
            << ": " << "provided=" << provided
            << ": " << "requested=" << requested
            << endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
}


void finalize()
{
    check(MPI_Finalize());
}


void check(int errorcode)
{
    if (MPI_SUCCESS != errorcode) {
        cerr << "[" << comm_rank(MPI_COMM_WORLD) << "] MPI ERROR"
            << ": " << errorcode
            << ": " << error_class(errorcode)
            << ": " << error_string(errorcode)
            << endl;
        MPI_Abort(MPI_COMM_WORLD, errorcode);
    }
}


int error_class(int errorcode)
{
    int eclass;
    MPI_Error_class(errorcode, &eclass);
    return eclass;
}


string error_string(int errorcode)
{
    string result;
    char estring[MPI_MAX_ERROR_STRING];
    int estring_length;

    MPI_Error_string(errorcode, estring, &estring_length);
    result.assign(estring, estring_length);

    return result;
}


MPI_Comm comm_dup(MPI_Comm orig)
{
    MPI_Comm dup;
    check(MPI_Comm_dup(orig, &dup));
    return dup;
}


void comm_free(MPI_Comm &comm)
{
    check(MPI_Comm_free(&comm));
}


int comm_rank(MPI_Comm comm)
{
    int result = 0;
    check(MPI_Comm_rank(comm, &result));
    return result;
}


int comm_size(MPI_Comm comm)
{
    int result = 0;
    check(MPI_Comm_size(comm, &result));
    return result;
}


void barrier(MPI_Comm comm)
{
    check(MPI_Barrier(comm));
}


void type_contiguous(int count, MPI_Datatype oldtype, MPI_Datatype &newtype)
{
    check(MPI_Type_contiguous(count, oldtype, &newtype));
}


void type_create_struct(int count, int blocklengths[], MPI_Aint displacements[], MPI_Datatype types[], MPI_Datatype &newtype)
{
    check(MPI_Type_create_struct(count, blocklengths, displacements, types, &newtype));
}


void type_commit(MPI_Datatype &type)
{
    check(MPI_Type_commit(&type));
}


void type_free(MPI_Datatype &type)
{
    check(MPI_Type_free(&type));
}


#define PAIR(A,B) pair<A,B>

#define MPIX_GET_MPI_DATATYPE_IMPL(T,M) \
MPI_Datatype get_mpi_datatype(T object) { return M; }
MPIX_GET_MPI_DATATYPE_IMPL(char,                   MPI_CHAR)
MPIX_GET_MPI_DATATYPE_IMPL(short,                  MPI_SHORT)
MPIX_GET_MPI_DATATYPE_IMPL(int,                    MPI_INT)
MPIX_GET_MPI_DATATYPE_IMPL(long,                   MPI_LONG)
MPIX_GET_MPI_DATATYPE_IMPL(long long,              MPI_LONG_LONG)
MPIX_GET_MPI_DATATYPE_IMPL(signed char,            MPI_BYTE)
MPIX_GET_MPI_DATATYPE_IMPL(unsigned char,          MPI_UNSIGNED_CHAR)
MPIX_GET_MPI_DATATYPE_IMPL(unsigned short,         MPI_UNSIGNED_SHORT)
MPIX_GET_MPI_DATATYPE_IMPL(unsigned int,           MPI_UNSIGNED)
MPIX_GET_MPI_DATATYPE_IMPL(unsigned long,          MPI_UNSIGNED_LONG)
MPIX_GET_MPI_DATATYPE_IMPL(unsigned long long,     MPI_UNSIGNED_LONG_LONG)
MPIX_GET_MPI_DATATYPE_IMPL(float,                  MPI_FLOAT)
MPIX_GET_MPI_DATATYPE_IMPL(double,                 MPI_DOUBLE)
MPIX_GET_MPI_DATATYPE_IMPL(long double,            MPI_LONG_DOUBLE)
MPIX_GET_MPI_DATATYPE_IMPL(PAIR(short, int),       MPI_SHORT_INT)
MPIX_GET_MPI_DATATYPE_IMPL(PAIR(int, int),         MPI_2INT)
MPIX_GET_MPI_DATATYPE_IMPL(PAIR(long, int),        MPI_LONG_INT)
MPIX_GET_MPI_DATATYPE_IMPL(PAIR(float, int),       MPI_FLOAT_INT)
MPIX_GET_MPI_DATATYPE_IMPL(PAIR(double, int),      MPI_DOUBLE_INT)
MPIX_GET_MPI_DATATYPE_IMPL(PAIR(long double, int), MPI_LONG_DOUBLE_INT)

#define MPIX_BCAST_IMPL_ALLP(A,B) \
MPIX_BCAST_IMPL_ONE(PAIR(A,B))    \
MPIX_BCAST_IMPL_ARR(PAIR(A,B))    \
MPIX_BCAST_IMPL_VEC(PAIR(A,B))

MPIX_BCAST_IMPL_ALL(char)
MPIX_BCAST_IMPL_ALL(short)
MPIX_BCAST_IMPL_ALL(int)
MPIX_BCAST_IMPL_ALL(long)
MPIX_BCAST_IMPL_ALL(long long)
MPIX_BCAST_IMPL_ALL(signed char)
MPIX_BCAST_IMPL_ALL(unsigned char)
MPIX_BCAST_IMPL_ALL(unsigned short)
MPIX_BCAST_IMPL_ALL(unsigned int)
MPIX_BCAST_IMPL_ALL(unsigned long)
MPIX_BCAST_IMPL_ALL(unsigned long long)
MPIX_BCAST_IMPL_ALL(float)
MPIX_BCAST_IMPL_ALL(double)
MPIX_BCAST_IMPL_ALL(long double)
MPIX_BCAST_IMPL_ALLP(short, int)
MPIX_BCAST_IMPL_ALLP(int, int)
MPIX_BCAST_IMPL_ALLP(long, int)
MPIX_BCAST_IMPL_ALLP(float, int)
MPIX_BCAST_IMPL_ALLP(double, int)
MPIX_BCAST_IMPL_ALLP(long double, int)

void bcast(string &object, int root, MPI_Comm comm)
{
    int size = int(object.size());

    bcast(size, root, comm);
    if (comm_rank(comm) == root) {
        bcast(const_cast<char*>(object.data()), size, root, comm);
    }
    else {
        char *data = new char[size];
        bcast(data, size, root, comm);
        object.assign(data, size);
        delete [] data;
    }
}


void bcast(vector<string> &object, int root, MPI_Comm comm)
{
    int size = int(object.size());

    bcast(size, root, comm);
    if (comm_rank(comm) != root) {
        object.resize(size);
    }
    for (vector<string>::iterator it=object.begin(); it!=object.end(); ++it) {
        bcast(*it, root, comm);
    }
}


vector<string> bcast(int argc, char **argv, MPI_Comm comm)
{
    vector<string> result(argc);

    if (0 == comm_rank(comm)) {
        for (int i = 0; i < argc; ++i) {
            result[i] = argv[i];
        }
    }
    bcast(result, 0, comm);

    return result;
}


#define MPIX_REDUCE_IMPL_ALLP(A,B) \
MPIX_REDUCE_IMPL_ONE(PAIR(A,B))    \
MPIX_REDUCE_IMPL_ARR(PAIR(A,B))    \
MPIX_REDUCE_IMPL_VEC(PAIR(A,B))

MPIX_REDUCE_IMPL_ALL(char)
MPIX_REDUCE_IMPL_ALL(short)
MPIX_REDUCE_IMPL_ALL(int)
MPIX_REDUCE_IMPL_ALL(long)
MPIX_REDUCE_IMPL_ALL(long long)
MPIX_REDUCE_IMPL_ALL(signed char)
MPIX_REDUCE_IMPL_ALL(unsigned char)
MPIX_REDUCE_IMPL_ALL(unsigned short)
MPIX_REDUCE_IMPL_ALL(unsigned int)
MPIX_REDUCE_IMPL_ALL(unsigned long)
MPIX_REDUCE_IMPL_ALL(unsigned long long)
MPIX_REDUCE_IMPL_ALL(float)
MPIX_REDUCE_IMPL_ALL(double)
MPIX_REDUCE_IMPL_ALL(long double)
MPIX_REDUCE_IMPL_ALLP(short, int)
MPIX_REDUCE_IMPL_ALLP(int, int)
MPIX_REDUCE_IMPL_ALLP(long, int)
MPIX_REDUCE_IMPL_ALLP(float, int)
MPIX_REDUCE_IMPL_ALLP(double, int)
MPIX_REDUCE_IMPL_ALLP(long double, int)


#define MPIX_ALLREDUCE_IMPL_ALLP(A,B) \
MPIX_ALLREDUCE_IMPL_ONE(PAIR(A,B))    \
MPIX_ALLREDUCE_IMPL_ARR(PAIR(A,B))    \
MPIX_ALLREDUCE_IMPL_VEC(PAIR(A,B))

MPIX_ALLREDUCE_IMPL_ALL(char)
MPIX_ALLREDUCE_IMPL_ALL(short)
MPIX_ALLREDUCE_IMPL_ALL(int)
MPIX_ALLREDUCE_IMPL_ALL(long)
MPIX_ALLREDUCE_IMPL_ALL(long long)
MPIX_ALLREDUCE_IMPL_ALL(signed char)
MPIX_ALLREDUCE_IMPL_ALL(unsigned char)
MPIX_ALLREDUCE_IMPL_ALL(unsigned short)
MPIX_ALLREDUCE_IMPL_ALL(unsigned int)
MPIX_ALLREDUCE_IMPL_ALL(unsigned long)
MPIX_ALLREDUCE_IMPL_ALL(unsigned long long)
MPIX_ALLREDUCE_IMPL_ALL(float)
MPIX_ALLREDUCE_IMPL_ALL(double)
MPIX_ALLREDUCE_IMPL_ALL(long double)
MPIX_ALLREDUCE_IMPL_ALLP(short, int)
MPIX_ALLREDUCE_IMPL_ALLP(int, int)
MPIX_ALLREDUCE_IMPL_ALLP(long, int)
MPIX_ALLREDUCE_IMPL_ALLP(float, int)
MPIX_ALLREDUCE_IMPL_ALLP(double, int)
MPIX_ALLREDUCE_IMPL_ALLP(long double, int)


#define MPIX_ALLTOALL_IMPL_ALLP(A,B)    \
MPIX_ALLTOALL_IMPL_SENDRECV(PAIR(A,B))  \
MPIX_ALLTOALL_IMPL_SEND_RECV(PAIR(A,B))

MPIX_ALLTOALL_IMPL_ALL(char)
MPIX_ALLTOALL_IMPL_ALL(short)
MPIX_ALLTOALL_IMPL_ALL(int)
MPIX_ALLTOALL_IMPL_ALL(long)
MPIX_ALLTOALL_IMPL_ALL(long long)
MPIX_ALLTOALL_IMPL_ALL(signed char)
MPIX_ALLTOALL_IMPL_ALL(unsigned char)
MPIX_ALLTOALL_IMPL_ALL(unsigned short)
MPIX_ALLTOALL_IMPL_ALL(unsigned int)
MPIX_ALLTOALL_IMPL_ALL(unsigned long)
MPIX_ALLTOALL_IMPL_ALL(unsigned long long)
MPIX_ALLTOALL_IMPL_ALL(float)
MPIX_ALLTOALL_IMPL_ALL(double)
MPIX_ALLTOALL_IMPL_ALL(long double)
MPIX_ALLTOALL_IMPL_ALLP(short, int)
MPIX_ALLTOALL_IMPL_ALLP(int, int)
MPIX_ALLTOALL_IMPL_ALLP(long, int)
MPIX_ALLTOALL_IMPL_ALLP(float, int)
MPIX_ALLTOALL_IMPL_ALLP(double, int)
MPIX_ALLTOALL_IMPL_ALLP(long double, int)


#define MPIX_GATHER_IMPL_ALLP(A,B) \
MPIX_GATHER_IMPL_ONE(PAIR(A,B))    \
MPIX_GATHER_IMPL_ARR(PAIR(A,B))    \
MPIX_GATHER_IMPL_VEC(PAIR(A,B))

MPIX_GATHER_IMPL_ALL(char)
MPIX_GATHER_IMPL_ALL(short)
MPIX_GATHER_IMPL_ALL(int)
MPIX_GATHER_IMPL_ALL(long)
MPIX_GATHER_IMPL_ALL(long long)
MPIX_GATHER_IMPL_ALL(signed char)
MPIX_GATHER_IMPL_ALL(unsigned char)
MPIX_GATHER_IMPL_ALL(unsigned short)
MPIX_GATHER_IMPL_ALL(unsigned int)
MPIX_GATHER_IMPL_ALL(unsigned long)
MPIX_GATHER_IMPL_ALL(unsigned long long)
MPIX_GATHER_IMPL_ALL(float)
MPIX_GATHER_IMPL_ALL(double)
MPIX_GATHER_IMPL_ALL(long double)
MPIX_GATHER_IMPL_ALLP(short, int)
MPIX_GATHER_IMPL_ALLP(int, int)
MPIX_GATHER_IMPL_ALLP(long, int)
MPIX_GATHER_IMPL_ALLP(float, int)
MPIX_GATHER_IMPL_ALLP(double, int)
MPIX_GATHER_IMPL_ALLP(long double, int)


vector<string> gather(string &object, int root, MPI_Comm comm)
{
    vector<string> result;
    int sendcount = int(object.size());
    vector<int> recvcounts = gather(sendcount, root, comm);

    if (comm_rank(comm) == root) {
        vector<char> recvbuf(accumulate(recvcounts.begin(), recvcounts.end(), 0));
        vector<int> displs(recvcounts.size());
        displs[0] = 0;
        for (size_t i=1; i<displs.size(); ++i) {
            displs[i] = displs[i-1] + recvcounts[i-1];
        }
        check(MPI_Gatherv(const_cast<char*>(object.data()), sendcount, MPI_CHAR,
                    &recvbuf[0], &recvcounts[0], &displs[0],
                    MPI_CHAR, root, comm));
        result.resize(comm_size(comm));
        for (size_t i=0; i<displs.size(); ++i) {
            result[i].assign(&recvbuf[displs[0]], recvcounts[0]);
        }
    }
    else {
        check(MPI_Gatherv(const_cast<char*>(object.data()), sendcount, MPI_CHAR,
                    NULL, NULL, NULL,
                    MPI_CHAR, root, comm));
    }

    return result;
}


MPIX_PRINT_SYNC_IMPL_ALL(char)
MPIX_PRINT_SYNC_IMPL_ALL(short)
MPIX_PRINT_SYNC_IMPL_ALL(int)
MPIX_PRINT_SYNC_IMPL_ALL(long)
MPIX_PRINT_SYNC_IMPL_ALL(long long)
MPIX_PRINT_SYNC_IMPL_ALL(signed char)
MPIX_PRINT_SYNC_IMPL_ALL(unsigned char)
MPIX_PRINT_SYNC_IMPL_ALL(unsigned short)
MPIX_PRINT_SYNC_IMPL_ALL(unsigned int)
MPIX_PRINT_SYNC_IMPL_ALL(unsigned long)
MPIX_PRINT_SYNC_IMPL_ALL(unsigned long long)
MPIX_PRINT_SYNC_IMPL_ALL(float)
MPIX_PRINT_SYNC_IMPL_ALL(double)
MPIX_PRINT_SYNC_IMPL_ALL(long double)
//MPIX_PRINT_SYNC_IMPL_ALLP(short, int)
//MPIX_PRINT_SYNC_IMPL_ALLP(int, int)
//MPIX_PRINT_SYNC_IMPL_ALLP(long, int)
//MPIX_PRINT_SYNC_IMPL_ALLP(float, int)
//MPIX_PRINT_SYNC_IMPL_ALLP(double, int)
//MPIX_PRINT_SYNC_IMPL_ALLP(long double, int)


void print_zero(const string &name, MPI_Comm comm)
{
    barrier(comm);
    if (0 == comm_rank(comm)) {
        cout << name << endl;
    }
    barrier(comm);
}


MPIX_PRINT_ZERO_IMPL_ALL(char)
MPIX_PRINT_ZERO_IMPL_ALL(short)
MPIX_PRINT_ZERO_IMPL_ALL(int)
MPIX_PRINT_ZERO_IMPL_ALL(long)
MPIX_PRINT_ZERO_IMPL_ALL(long long)
MPIX_PRINT_ZERO_IMPL_ALL(signed char)
MPIX_PRINT_ZERO_IMPL_ALL(unsigned char)
MPIX_PRINT_ZERO_IMPL_ALL(unsigned short)
MPIX_PRINT_ZERO_IMPL_ALL(unsigned int)
MPIX_PRINT_ZERO_IMPL_ALL(unsigned long)
MPIX_PRINT_ZERO_IMPL_ALL(unsigned long long)
MPIX_PRINT_ZERO_IMPL_ALL(float)
MPIX_PRINT_ZERO_IMPL_ALL(double)
MPIX_PRINT_ZERO_IMPL_ALL(long double)
//MPIX_PRINT_ZERO_IMPL_ALLP(short, int)
//MPIX_PRINT_ZERO_IMPL_ALLP(int, int)
//MPIX_PRINT_ZERO_IMPL_ALLP(long, int)
//MPIX_PRINT_ZERO_IMPL_ALLP(float, int)
//MPIX_PRINT_ZERO_IMPL_ALLP(double, int)
//MPIX_PRINT_ZERO_IMPL_ALLP(long double, int)


/* file reading */
MPI_Offset get_file_size(const string &file_name, MPI_Comm comm)
{
    MPI_Offset file_size = 0;

    /* process 0 stats the file to determine its size */
    if (0 == comm_rank(comm)) {
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
    bcast(file_size, 0, comm);

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
 * @param[in] comm instance
 */
void read_file(
        const string &file_name,
        char *&file_buffer,
        MPI_Offset &file_size,
        MPI_Comm comm)
{
    //read_file_mpiio(file_name, file_buffer, file_size, comm);
    read_file_bcast(file_name, file_buffer, file_size, comm);
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
 * @param[in] comm instance
 */
void read_file_bcast(
    const string &file_name,
    char *&file_buffer,
    MPI_Offset &file_size,
    MPI_Comm comm)
{
    long chunk_size = 1073741824;
    int rank = comm_rank(comm);

    /* allocate a buffer for the file, of the entire size */
    file_size = get_file_size(file_name, comm);
    if (NULL == file_buffer) {
        file_buffer = new char[file_size];
    }

    if (0 == rank) {
        int retval = 0;
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

        retval = fclose(file);
        if (0 != retval) {
            perror("fclose");
            printf("unable to close file on process 0\n");
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
            check(MPI_Bcast(&file_buffer[offset],
                                message_size, MPI_CHAR, 0, comm));
            offset += chunk_size;
        }
    }
    else {
        /* bcast entire file contents */
        check(MPI_Bcast(file_buffer, file_size, MPI_CHAR, 0, comm));
    }
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
 * @param[in] comm instance
 */
void read_file_mpiio(
        const string &file_name,
        char *&file_buffer,
        MPI_Offset &file_size,
        MPI_Comm comm)
{
    long chunk_size = 1073741824;
    int rank = comm_rank(comm);
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
        check(MPI_File_open(comm, const_cast<char *>(file_name.c_str()),
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
            check(MPI_File_read_at_all(fh, offset, file_buffer + offset,
                                           message_size, MPI_CHAR, &status));
            check(MPI_Get_count(&status, MPI_CHAR, &count));
            assert(count == message_size);
            offset += count;
        }
    }
    else {
        int count = 0;

        /* all procs read the entire file */
        check(MPI_File_open(comm, const_cast<char *>(file_name.c_str()),
                                /*MPI_MODE_RDONLY|MPI_MODE_UNIQUE_OPEN,*/
                                MPI_MODE_RDONLY,
                                MPI_INFO_NULL, &fh));
        check(MPI_File_read_at_all(fh, 0, file_buffer,
                                       file_size, MPI_CHAR, &status));
        check(MPI_Get_count(&status, MPI_CHAR, &count));
        assert(count == file_size);
        check(MPI_File_close(&fh));
    }
}

} /* namespace mpix */

