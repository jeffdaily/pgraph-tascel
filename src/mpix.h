/**
 * @file mpix.h
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 *
 * Templated and/or overloaded MPI functions with intelligent default
 * parameters.
 */
#ifndef _MPIX_H_
#define _MPIX_H_

#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <mpi.h>

using std::cout;
using std::endl;
using std::ostringstream;
using std::pair;
using std::string;
using std::vector;


#define ARG_LEN_MAX 1024

#define MPI_CHECK(what) do {                              \
    int __err;                                            \
    __err = what;                                         \
    if (MPI_SUCCESS != __err) {                           \
        printf("[%d] FAILED FILE=%s LINE=%d:" #what "\n", \
                rank, __FILE__, __LINE__);                \
        MPI_Abort(MPI_COMM_WORLD, -1);                    \
    }                                                     \
} while (0)

/* for a collective function and we don't have a 'rank' */
#define MPI_CHECK_C(what) do {                            \
    int __err;                                            \
    __err = what;                                         \
    if (MPI_SUCCESS != __err) {                           \
        printf("[C] FAILED FILE=%s LINE=%d:" #what "\n",  \
                __FILE__, __LINE__);                      \
        MPI_Abort(MPI_COMM_WORLD, -1);                    \
    }                                                     \
} while (0)

#if 0
/**
 * Return the MPI_Datatype associated with the given parameter.
 */
template <typename T>
inline MPI_Datatype
mpix_get_mpi_datatype(const T &x)
{
    printf("FAILED FILE=%s LINE=%d: unknown datatype\n", __FILE__, __LINE__);
    MPI_Abort(MPI_COMM_WORLD, -1);
}
#endif
/* helper macro and macro invocations to implement the known types */
#define MPIX_GET_MPI_DATATYPE_IMPL(CTYPE,MTYPE)         \
inline MPI_Datatype                                     \
mpix_get_mpi_datatype(CTYPE&) { return MTYPE; }
#define MPIX_GET_MPI_DATATYPE_IMPL_CONST(CTYPE,MTYPE)   \
inline MPI_Datatype                                     \
mpix_get_mpi_datatype(const CTYPE&) { return MTYPE; }
MPIX_GET_MPI_DATATYPE_IMPL(char, MPI_CHAR);
MPIX_GET_MPI_DATATYPE_IMPL(short, MPI_SHORT);
MPIX_GET_MPI_DATATYPE_IMPL(int, MPI_INT);
MPIX_GET_MPI_DATATYPE_IMPL(long, MPI_LONG);
#if defined(MPI_LONG_LONG_INT) || (defined(MPI_VERSION) && MPI_VERSION >= 2)
MPIX_GET_MPI_DATATYPE_IMPL(long long, MPI_LONG_LONG_INT);
#endif
MPIX_GET_MPI_DATATYPE_IMPL(float, MPI_FLOAT);
MPIX_GET_MPI_DATATYPE_IMPL(double, MPI_DOUBLE);
MPIX_GET_MPI_DATATYPE_IMPL(long double, MPI_LONG_DOUBLE);
MPIX_GET_MPI_DATATYPE_IMPL(unsigned char, MPI_UNSIGNED_CHAR);
MPIX_GET_MPI_DATATYPE_IMPL(unsigned short, MPI_UNSIGNED_SHORT);
MPIX_GET_MPI_DATATYPE_IMPL(unsigned, MPI_UNSIGNED);
MPIX_GET_MPI_DATATYPE_IMPL(unsigned long, MPI_UNSIGNED_LONG);
#if defined(MPI_UNSIGNED_LONG_LONG) || (defined(MPI_VERSION) && MPI_VERSION >= 2)
MPIX_GET_MPI_DATATYPE_IMPL(unsigned long long, MPI_UNSIGNED_LONG_LONG);
#endif
#define MPIX_LIST2(A,B) A,B
MPIX_GET_MPI_DATATYPE_IMPL(pair<MPIX_LIST2(float, int)>, MPI_FLOAT_INT);
MPIX_GET_MPI_DATATYPE_IMPL(pair<MPIX_LIST2(double, int)>, MPI_DOUBLE_INT);
MPIX_GET_MPI_DATATYPE_IMPL(pair<MPIX_LIST2(long double, int)>, MPI_LONG_DOUBLE_INT);
MPIX_GET_MPI_DATATYPE_IMPL(pair<MPIX_LIST2(long, int)>, MPI_LONG_INT);
MPIX_GET_MPI_DATATYPE_IMPL(pair<MPIX_LIST2(short, int)>, MPI_SHORT_INT);
MPIX_GET_MPI_DATATYPE_IMPL(pair<MPIX_LIST2(int, int)>, MPI_2INT);
#undef MPIX_LIST2
MPIX_GET_MPI_DATATYPE_IMPL(void *, MPI_UNSIGNED_LONG);
MPIX_GET_MPI_DATATYPE_IMPL(char *, MPI_UNSIGNED_LONG);

MPIX_GET_MPI_DATATYPE_IMPL_CONST(char, MPI_CHAR);
MPIX_GET_MPI_DATATYPE_IMPL_CONST(short, MPI_SHORT);
MPIX_GET_MPI_DATATYPE_IMPL_CONST(int, MPI_INT);
MPIX_GET_MPI_DATATYPE_IMPL_CONST(long, MPI_LONG);
#if defined(MPI_LONG_LONG_INT) || (defined(MPI_VERSION) && MPI_VERSION >= 2)
MPIX_GET_MPI_DATATYPE_IMPL_CONST(long long, MPI_LONG_LONG_INT);
#endif
MPIX_GET_MPI_DATATYPE_IMPL_CONST(float, MPI_FLOAT);
MPIX_GET_MPI_DATATYPE_IMPL_CONST(double, MPI_DOUBLE);
MPIX_GET_MPI_DATATYPE_IMPL_CONST(long double, MPI_LONG_DOUBLE);
MPIX_GET_MPI_DATATYPE_IMPL_CONST(unsigned char, MPI_UNSIGNED_CHAR);
MPIX_GET_MPI_DATATYPE_IMPL_CONST(unsigned short, MPI_UNSIGNED_SHORT);
MPIX_GET_MPI_DATATYPE_IMPL_CONST(unsigned, MPI_UNSIGNED);
MPIX_GET_MPI_DATATYPE_IMPL_CONST(unsigned long, MPI_UNSIGNED_LONG);
#if defined(MPI_UNSIGNED_LONG_LONG) || (defined(MPI_VERSION) && MPI_VERSION >= 2)
MPIX_GET_MPI_DATATYPE_IMPL_CONST(unsigned long long, MPI_UNSIGNED_LONG_LONG);
#endif
#define MPIX_LIST2(A,B) A,B
MPIX_GET_MPI_DATATYPE_IMPL_CONST(pair<MPIX_LIST2(float, int)>, MPI_FLOAT_INT);
MPIX_GET_MPI_DATATYPE_IMPL_CONST(pair<MPIX_LIST2(double, int)>, MPI_DOUBLE_INT);
MPIX_GET_MPI_DATATYPE_IMPL_CONST(pair<MPIX_LIST2(long double, int)>, MPI_LONG_DOUBLE_INT);
MPIX_GET_MPI_DATATYPE_IMPL_CONST(pair<MPIX_LIST2(long, int)>, MPI_LONG_INT);
MPIX_GET_MPI_DATATYPE_IMPL_CONST(pair<MPIX_LIST2(short, int)>, MPI_SHORT_INT);
MPIX_GET_MPI_DATATYPE_IMPL_CONST(pair<MPIX_LIST2(int, int)>, MPI_2INT);
#undef MPIX_LIST2
MPIX_GET_MPI_DATATYPE_IMPL_CONST(void *, MPI_UNSIGNED_LONG);
MPIX_GET_MPI_DATATYPE_IMPL_CONST(char *, MPI_UNSIGNED_LONG);

#if 0
inline MPI_Datatype
mpix_get_mpi_datatype(const void *&x)
{
    if (sizeof(void *) == sizeof(int)) {
        return MPI_INT;
    }
    else if (sizeof(void *) == sizeof(long)) {
        return MPI_LONG;
    }
#if defined(MPI_LONG_LONG_INT) || (defined(MPI_VERSION) && MPI_VERSION >= 2)
    else if (sizeof(void *) == sizeof(long long)) {
        return MPI_LONG_LONG_INT;
    }
#endif
    else {
        printf("FAILED FILE=%s LINE=%d: unknown pointer datatype\n",
               __FILE__, __LINE__);
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
}
#endif

template <class T>
void mpix_bcast(T &object, int root = 0, MPI_Comm comm = MPI_COMM_WORLD)
{
    MPI_Datatype datatype = mpix_get_mpi_datatype(object);
    MPI_CHECK_C(MPI_Bcast(&object, 1, datatype, root, comm));
}

template <class T>
void mpix_bcast(vector<T> &object, int root = 0, MPI_Comm comm = MPI_COMM_WORLD)
{
    typedef typename vector<T>::size_type size_type;

    int rank = 0;
    size_type size = 0;
    MPI_Datatype datatype_size = MPI_DATATYPE_NULL;
    MPI_Datatype datatype_obj  = MPI_DATATYPE_NULL;

    MPI_CHECK(MPI_Comm_rank(comm, &rank));

    size = object.size();
    datatype_size = mpix_get_mpi_datatype(size);
    datatype_obj  = mpix_get_mpi_datatype(T());

    MPI_CHECK(MPI_Bcast(&size, 1, datatype_size, root, comm));

    if (rank != root) {
        object.resize(size);
    }

    MPI_CHECK(MPI_Bcast(&object[0], size, datatype_obj, root, comm));
}

template <class T>
void mpix_allreduce(T &object, MPI_Op op, MPI_Comm comm = MPI_COMM_WORLD)
{
    MPI_Datatype datatype = mpix_get_mpi_datatype(object);
    MPI_CHECK_C(MPI_Allreduce(MPI_IN_PLACE, &object, 1, datatype, op, comm));
}

template <class T>
void mpix_allreduce(vector<T> &object, MPI_Op op, MPI_Comm comm = MPI_COMM_WORLD)
{
    MPI_Datatype datatype = mpix_get_mpi_datatype(object[0]);
    MPI_CHECK_C(MPI_Allreduce(MPI_IN_PLACE, &object[0], object.size(),
                              datatype, op, comm));
}

/* overloads of mpix_print_sync */
void mpix_print_sync(MPI_Comm comm, const string &name, const vector<string> &what);
void mpix_print_sync(MPI_Comm comm, const string &name, const string &what);
void mpix_print_sync(MPI_Comm comm, const string &what);
/* template to capture remaining generic cases */
template <class T>
void mpix_print_sync(MPI_Comm comm, const string &name, const T &what)
{
    int rank;
    int size;
    T what_copy = what;
    T *all_what = NULL;

    MPI_CHECK(MPI_Comm_rank(comm, &rank));
    MPI_CHECK(MPI_Comm_size(comm, &size));

    all_what = new T[size];
    MPI_CHECK(MPI_Gather(&what_copy, sizeof(T), MPI_CHAR, all_what, sizeof(T),
                         MPI_CHAR, 0, comm));
    if (0 == rank) {
        for (int i = 0; i < size; ++i) {
            cout << "[" << i << "] " << name << "=" << all_what[i] << endl;
        }
    }
    delete [] all_what;
    MPI_Barrier(comm);
}


void mpix_bcast_argv(int argc, char **argv,
                     vector<string> &all_argv, MPI_Comm = MPI_COMM_WORLD);

MPI_Offset mpix_get_file_size(
    const string &file_name, MPI_Comm comm = MPI_COMM_WORLD);

void mpix_read_file(MPI_Comm comm, const string &file_name, char *&file_buffer, MPI_Offset &file_size, long chunk_size = 1073741824);
void mpix_read_file_bcast(MPI_Comm comm, const string &file_name, char *&file_buffer, MPI_Offset &file_size, long chunk_size = 1073741824);
void mpix_read_file_mpiio(MPI_Comm comm, const string &file_name, char *&file_buffer, MPI_Offset &file_size, long chunk_size = 1073741824);

template <class Collection>
string vec_to_string(
    const Collection &collection,
    char const *delimiter = ",",
    const string &name = "")
{
    typedef typename Collection::const_iterator iter;
    std::ostringstream os;
    iter beg = collection.begin();
    iter end = collection.end();

    if (!name.empty()) {
        os << name << "=";
    }

    if (beg == end) {
        os << "{}";
    }
    else {
        os << "{" << *(beg++);
        for (; beg != end; ++beg) {
            os << delimiter << *beg;
        }
        os << "}";
    }

    return os.str();
}

#endif /* _MPIX_H_ */
