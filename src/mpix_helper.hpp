/**
 * @file mpix_helper.hpp
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 *
 * Inline implementations of templated and/or overloaded MPI functions.
 */

#include <mpi.h>

#include <sys/stat.h>

#include <iostream>
#include <string>
#include <vector>

#include "mpix.hpp"

using ::std::cerr;
using ::std::cout;
using ::std::endl;
using ::std::string;
using ::std::vector;

namespace mpix {

#define MPIX_BCAST_IMPL_ONE(T)                          \
void bcast(T &object, int root, MPI_Comm comm)          \
{                                                       \
    MPI_Datatype datatype = get_mpi_datatype(object);   \
    check(MPI_Bcast(&object, 1, datatype, root, comm)); \
}

#define MPIX_BCAST_IMPL_ARR(T)                            \
void bcast(T *object, int size, int root, MPI_Comm comm)  \
{                                                         \
    MPI_Datatype datatype = get_mpi_datatype(object[0]);  \
    check(MPI_Bcast(object, size, datatype, root, comm)); \
}

#define MPIX_BCAST_IMPL_VEC(T)                                         \
void bcast(vector<T> &object, int root, MPI_Comm comm)                 \
{                                                                      \
    MPI_Datatype datatype = get_mpi_datatype(object[0]);               \
    check(MPI_Bcast(&object[0], object.size(), datatype, root, comm)); \
}

#define MPIX_BCAST_IMPL_ALL(T) \
MPIX_BCAST_IMPL_ONE(T)         \
MPIX_BCAST_IMPL_ARR(T)         \
MPIX_BCAST_IMPL_VEC(T)


#define MPIX_REDUCE_IMPL_ONE(T)                                                \
void reduce(T &object, MPI_Op op, int root, MPI_Comm comm)                     \
{                                                                              \
    MPI_Datatype datatype = get_mpi_datatype(object);                          \
                                                                               \
    if (root == comm_rank(comm)) {                                             \
        check(MPI_Reduce(MPI_IN_PLACE, &object, 1, datatype, op, root, comm)); \
    }                                                                          \
    else {                                                                     \
        check(MPI_Reduce(&object, NULL, 1, datatype, op, root, comm));         \
    }                                                                          \
}


#define MPIX_REDUCE_IMPL_ARR(T)                                                  \
void reduce(T *object, int size, MPI_Op op, int root, MPI_Comm comm)             \
{                                                                                \
    MPI_Datatype datatype = get_mpi_datatype(object[0]);                         \
                                                                                 \
    if (root == comm_rank(comm)) {                                               \
        check(MPI_Reduce(MPI_IN_PLACE, object, size, datatype, op, root, comm)); \
    }                                                                            \
    else {                                                                       \
        check(MPI_Reduce(object, NULL, size, datatype, op, root, comm));         \
    }                                                                            \
}


#define MPIX_REDUCE_IMPL_VEC(T)                                     \
void reduce(vector<T> &object, MPI_Op op, int root, MPI_Comm comm)  \
{                                                                   \
    MPI_Datatype datatype = get_mpi_datatype(object[0]);            \
                                                                    \
    if (root == comm_rank(comm)) {                                  \
        check(MPI_Reduce(MPI_IN_PLACE, &object[0], object.size(),   \
                    datatype, op, root, comm));                     \
    }                                                               \
    else {                                                          \
        check(MPI_Reduce(&object[0], NULL, object.size(),           \
                    datatype, op, root, comm));                     \
    }                                                               \
}


#define MPIX_REDUCE_IMPL_ALL(T) \
MPIX_REDUCE_IMPL_ONE(T)         \
MPIX_REDUCE_IMPL_ARR(T)         \
MPIX_REDUCE_IMPL_VEC(T)


#define MPIX_ALLREDUCE_IMPL_ONE(T)                                      \
void allreduce(T &object, MPI_Op op, MPI_Comm comm)                     \
{                                                                       \
    MPI_Datatype datatype = get_mpi_datatype(object);                   \
    check(MPI_Allreduce(MPI_IN_PLACE, &object, 1, datatype, op, comm)); \
}


#define MPIX_ALLREDUCE_IMPL_ARR(T)                                        \
void allreduce(T *object, int size, MPI_Op op, MPI_Comm comm)             \
{                                                                         \
    MPI_Datatype datatype = get_mpi_datatype(object[0]);                  \
    check(MPI_Allreduce(MPI_IN_PLACE, object, size, datatype, op, comm)); \
}


#define MPIX_ALLREDUCE_IMPL_VEC(T)                                                     \
void allreduce(vector<T> &object, MPI_Op op, MPI_Comm comm)                            \
{                                                                                      \
    MPI_Datatype datatype = get_mpi_datatype(object[0]);                               \
    check(MPI_Allreduce(MPI_IN_PLACE, &object[0], object.size(), datatype, op, comm)); \
}


#define MPIX_ALLREDUCE_IMPL_ALL(T) \
MPIX_ALLREDUCE_IMPL_ONE(T)         \
MPIX_ALLREDUCE_IMPL_ARR(T)         \
MPIX_ALLREDUCE_IMPL_VEC(T)


/* unfortunately, support of MPI_IN_PLACE for MPI_Alltoall is
 * extremely lacking, e.g., OpenMPI, Cray. */
#define MPIX_ALLTOALL_IMPL_SENDRECV(T)                              \
void alltoall(vector<T> &object, MPI_Comm comm)                     \
{                                                                   \
    int size_comm = comm_size(comm);                                \
    int count = 0;                                                  \
    MPI_Datatype datatype = get_mpi_datatype(object[0]);            \
    vector<T> object_copy(object);                                  \
                                                                    \
    if (object.size() % size_comm != 0) {                           \
        cerr << "[" << comm_rank(MPI_COMM_WORLD) << "] MPI ERROR"   \
            << ": " << "mpix::alltoall has incorrect buffer size"   \
            << endl;                                                \
        MPI_Abort(MPI_COMM_WORLD, -1);                              \
    }                                                               \
    count = object.size() / size_comm;                              \
    check(MPI_Alltoall(&object_copy[0], count, datatype,            \
                &object[0], count, datatype, comm));                \
}


/* unfortunately, support of MPI_IN_PLACE for MPI_Alltoall is
 * extremely lacking, e.g., OpenMPI, Cray. */
#define MPIX_ALLTOALL_IMPL_SEND_RECV(T)                                 \
void alltoall(vector<T> &sendbuf, vector<T> &recvbuf, MPI_Comm comm)    \
{                                                                       \
    int size_comm = comm_size(comm);                                    \
    int count = 0;                                                      \
    MPI_Datatype datatype = get_mpi_datatype(sendbuf[0]);               \
                                                                        \
    if (sendbuf.size() % size_comm != 0) {                              \
        cerr << "[" << comm_rank(MPI_COMM_WORLD) << "] MPI ERROR"       \
            << ": " << "mpix::alltoall has incorrect buffer size"       \
            << endl;                                                    \
        MPI_Abort(MPI_COMM_WORLD, -1);                                  \
    }                                                                   \
    if (sendbuf.size() != recvbuf.size()) {                             \
        cerr << "[" << comm_rank(MPI_COMM_WORLD) << "] MPI ERROR"       \
            << ": " << "mpix::alltoall mismatched buffer size"          \
            << endl;                                                    \
        MPI_Abort(MPI_COMM_WORLD, -1);                                  \
    }                                                                   \
    count = sendbuf.size() / size_comm;                                 \
    check(MPI_Alltoall(&sendbuf[0], count, datatype,                    \
                &recvbuf[0], count, datatype, comm));                   \
}


#define MPIX_ALLTOALL_IMPL_ALL(T)   \
MPIX_ALLTOALL_IMPL_SENDRECV(T)      \
MPIX_ALLTOALL_IMPL_SEND_RECV(T)


#define MPIX_GATHER_IMPL_ONE(T)                             \
vector<T> gather(T object, int root, MPI_Comm comm)         \
{                                                           \
    vector<T> result;                                       \
    MPI_Datatype datatype = get_mpi_datatype(object);       \
                                                            \
    if (comm_rank(comm) == root) {                          \
        result.resize(comm_size(comm));                     \
        check(MPI_Gather(&object, 1, datatype,              \
                    &result[0], 1, datatype, root, comm));  \
    }                                                       \
    else {                                                  \
        check(MPI_Gather(&object, 1, datatype,              \
                    NULL, 0, datatype, root, comm));        \
    }                                                       \
                                                            \
    return result;                                          \
}


#define MPIX_GATHER_IMPL_ARR(T)                                 \
vector<T> gather(T *object, int size, int root, MPI_Comm comm)  \
{                                                               \
    vector<T> result;                                           \
    MPI_Datatype datatype = get_mpi_datatype(*object);          \
                                                                \
    if (comm_rank(comm) == root) {                              \
        result.resize(comm_size(comm) * size);                  \
        check(MPI_Gather(object, size, datatype,                \
                    &result[0], size, datatype, root, comm));   \
    }                                                           \
    else {                                                      \
        check(MPI_Gather(object, size, datatype,                \
                    NULL, size, datatype, root, comm));         \
    }                                                           \
                                                                \
    return result;                                              \
}


#define MPIX_GATHER_IMPL_NTV(T)                                                                 \
void gather(T *sendbuf, int sendcount, T *recvbuf, int recvcount, int root, MPI_Comm comm)      \
{                                                                                               \
    MPI_Datatype datatype = get_mpi_datatype(*sendbuf);                                         \
    check(MPI_Gather(sendbuf, sendcount, datatype, recvbuf, recvcount, datatype, root, comm));  \
}


#define MPIX_GATHER_IMPL_VEC(T)                                 \
vector<T> gather(vector<T> &object, int root, MPI_Comm comm)    \
{                                                               \
    vector<T> result;                                           \
    MPI_Datatype datatype = get_mpi_datatype(object[0]);        \
    int size = int(object.size());                              \
                                                                \
    if (comm_rank(comm) == root) {                              \
        result.resize(comm_size(comm) * size);                  \
        check(MPI_Gather(&object[0], size, datatype,            \
                    &result[0], size, datatype, root, comm));   \
    }                                                           \
    else {                                                      \
        check(MPI_Gather(&object[0], size, datatype,            \
                    NULL, 0, datatype, root, comm));            \
    }                                                           \
                                                                \
    return result;                                              \
}

#define MPIX_GATHER_IMPL_ALL(T) \
MPIX_GATHER_IMPL_ONE(T)         \
MPIX_GATHER_IMPL_ARR(T)         \
MPIX_GATHER_IMPL_VEC(T)         \
MPIX_GATHER_IMPL_NTV(T)


#define MPIX_PRINT_SYNC_IMPL_ONE(T)                                         \
void print_sync(const string &name, T what, MPI_Comm comm)                  \
{                                                                           \
    vector<T> all_what = gather(what, 0, comm);                             \
                                                                            \
    if (0 == comm_rank(comm)) {                                             \
        for (int i = 0, size = comm_size(comm); i < size; ++i) {            \
            cout << "[" << i << "] " << name << "=" << all_what[i] << endl; \
        }                                                                   \
    }                                                                       \
                                                                            \
    barrier(comm);                                                          \
}


#define MPIX_PRINT_SYNC_IMPL_ARR(T)                                       \
void print_sync(const string &name, T *what, int size_, MPI_Comm comm)    \
{                                                                         \
    vector<T> all_what = gather(what, size_, 0, comm);                    \
                                                                          \
    if (0 == comm_rank(comm)) {                                           \
        for (int i = 0, size = comm_size(comm); i < size; ++i) {          \
            cout << "[" << i << "] " << name << "={";                     \
            cout << all_what[i*size_];                                    \
            for (int j = 1; j < size_; ++j) {                             \
                cout << "," << all_what[i*size_ + j];                     \
            }                                                             \
            cout << "}" << endl;                                          \
        }                                                                 \
    }                                                                     \
                                                                          \
    barrier(comm);                                                        \
}


#define MPIX_PRINT_SYNC_IMPL_VEC(T)                                       \
void print_sync(const string &name, vector<T> &what, MPI_Comm comm)       \
{                                                                         \
    vector<T> all_what = gather(what, 0, comm);                           \
                                                                          \
    if (0 == comm_rank(comm)) {                                           \
        for (int i = 0, size = comm_size(comm); i < size; ++i) {          \
            cout << "[" << i << "] " << name << "={";                     \
            cout << all_what[i*what.size()];                              \
            for (int j = 1; j < what.size(); ++j) {                       \
                cout << "," << all_what[i*what.size() + j];               \
            }                                                             \
            cout << "}" << endl;                                          \
        }                                                                 \
    }                                                                     \
                                                                          \
    barrier(comm);                                                        \
}


#define MPIX_PRINT_SYNC_IMPL_ALL(T) \
MPIX_PRINT_SYNC_IMPL_ONE(T)         \
MPIX_PRINT_SYNC_IMPL_ARR(T)         \
MPIX_PRINT_SYNC_IMPL_VEC(T)


#define MPIX_PRINT_ZERO_IMPL_ONE(T)                           \
void print_zero(const string &name, T what, MPI_Comm comm)    \
{                                                             \
    barrier(comm);                                            \
    if (0 == comm_rank(comm)) {                               \
        cout << name << "=" << what << endl;                  \
    }                                                         \
    barrier(comm);                                            \
}


#define MPIX_PRINT_ZERO_IMPL_ARR(T)                                   \
void print_zero(const string &name, T *what, int size, MPI_Comm comm) \
{                                                                     \
    barrier(comm);                                                    \
    if (0 == comm_rank(comm)) {                                       \
        cout << name << "={";                                         \
        cout << what[0];                                              \
        for (int j = 1; j < size; ++j) {                              \
            cout << "," << what[j];                                   \
        }                                                             \
        cout << "}" << endl;                                          \
    }                                                                 \
                                                                      \
    barrier(comm);                                                    \
}


#define MPIX_PRINT_ZERO_IMPL_VEC(T)                                   \
void print_zero(const string &name, vector<T> &what, MPI_Comm comm)   \
{                                                                     \
    barrier(comm);                                                    \
    if (0 == comm_rank(comm)) {                                       \
        cout << name << "={";                                         \
        cout << what[0];                                              \
        for (int j = 1; j < what.size(); ++j) {                       \
            cout << "," << what[j];                                   \
        }                                                             \
        cout << "}" << endl;                                          \
    }                                                                 \
                                                                      \
    barrier(comm);                                                    \
}


#define MPIX_PRINT_ZERO_IMPL_ALL(T) \
MPIX_PRINT_ZERO_IMPL_ONE(T)         \
MPIX_PRINT_ZERO_IMPL_ARR(T)         \
MPIX_PRINT_ZERO_IMPL_VEC(T)


#define MPIX_DECL_ALL(T)                                                \
void bcast(T &object, int root, MPI_Comm comm);                         \
void bcast(T *object, int size, int root, MPI_Comm comm);               \
void bcast(vector<T> &object, int root, MPI_Comm comm);                 \
void reduce(T &object, MPI_Op op, int root, MPI_Comm comm);             \
void reduce(T *object, int size, MPI_Op op, int root, MPI_Comm comm);   \
void reduce(vector<T> &object, MPI_Op op, int root, MPI_Comm comm);     \
void allreduce(T &object, MPI_Op op, MPI_Comm comm);                    \
void allreduce(T *object, int size, MPI_Op op, MPI_Comm comm);          \
void allreduce(vector<T> &object, MPI_Op op, MPI_Comm comm);            \
void alltoall(vector<T> &object, MPI_Comm comm);                        \
void alltoall(vector<T> &sendbuf, vector<T> &recvbuf, MPI_Comm comm);   \
vector<T> gather(T object, int root, MPI_Comm comm);                    \
vector<T> gather(T *object, int size, int root, MPI_Comm comm);         \
void      gather(T *sendbuf, int sendcount, T *recvbuf, int recvcount, int root, MPI_Comm comm); \
vector<T> gather(vector<T> &object, int root, MPI_Comm comm);           \
void print_sync(const string &name, T what, MPI_Comm comm);             \
void print_sync(const string &name, T *what, int size_, MPI_Comm comm); \
void print_sync(const string &name, vector<T> &what, MPI_Comm comm);    \
void print_zero(const string &name, T what, MPI_Comm comm);             \
void print_zero(const string &name, T *what, int size, MPI_Comm comm);  \
void print_zero(const string &name, vector<T> &what, MPI_Comm comm);


#define MPIX_IMPL_ALL(T)    \
MPIX_BCAST_IMPL_ALL(T)      \
MPIX_REDUCE_IMPL_ALL(T)     \
MPIX_ALLREDUCE_IMPL_ALL(T)  \
MPIX_ALLTOALL_IMPL_ALL(T)   \
MPIX_GATHER_IMPL_ALL(T)     \
MPIX_PRINT_SYNC_IMPL_ALL(T) \
MPIX_PRINT_ZERO_IMPL_ALL(T)

} /* namespace mpix */

