/**
 * @file mpix-inl.hpp
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 *
 * Inline implementations of templated and/or overloaded MPI functions.
 */
#ifndef _MPIX_INL_H_
#define _MPIX_INL_H_

#include <sys/stat.h>

#include <cstddef>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <mpi.h>

using ::std::accumulate;
using ::std::cerr;
using ::std::cout;
using ::std::endl;
using ::std::ostringstream;
using ::std::pair;
using ::std::size_t;
using ::std::string;
using ::std::vector;

namespace mpix {

inline MPI_Datatype type_contiguous(int count, MPI_Datatype oldtype)
{
    MPI_Datatype newtype;
    check(MPI_Type_contiguous(count, oldtype, &newtype));
    return newtype;
}

inline MPI_Datatype type_create_struct(int count, const int blocklengths[], const MPI_Aint displacements[], const MPI_Datatype types[])
{
    MPI_Datatype newtype;
    check(MPI_Type_create_struct(count, blocklengths, displacements, types, &newtype));
    return newtype;
}

inline void type_commit(MPI_Datatype &type)
{
    check(MPI_Type_commit(&type));
}

template <typename T>
inline MPI_Datatype build_mpi_datatype(const T&)
{
    return MPI_DATATYPE_NULL;
}

template <typename T>
inline MPI_Datatype get_mpi_datatype(const T& object)
{
    static MPI_Datatype type = build_mpi_datatype(object);
    assert(type != MPI_DATATYPE_NULL);
    return type;
}

/* helper macro and macro invocations to implement the known types */
#define MPIX_GET_MPI_DATATYPE_IMPL(CTYPE,MTYPE)         \
template <>                                             \
inline MPI_Datatype                                     \
get_mpi_datatype<CTYPE>(const CTYPE&) {                 \
    return MTYPE;                                       \
}

MPIX_GET_MPI_DATATYPE_IMPL(char, MPI_CHAR);
MPIX_GET_MPI_DATATYPE_IMPL(short, MPI_SHORT);
MPIX_GET_MPI_DATATYPE_IMPL(int, MPI_INT);
MPIX_GET_MPI_DATATYPE_IMPL(long, MPI_LONG);
MPIX_GET_MPI_DATATYPE_IMPL(float, MPI_FLOAT);
MPIX_GET_MPI_DATATYPE_IMPL(double, MPI_DOUBLE);
MPIX_GET_MPI_DATATYPE_IMPL(long double, MPI_LONG_DOUBLE);
MPIX_GET_MPI_DATATYPE_IMPL(unsigned char, MPI_UNSIGNED_CHAR);
MPIX_GET_MPI_DATATYPE_IMPL(unsigned short, MPI_UNSIGNED_SHORT);
MPIX_GET_MPI_DATATYPE_IMPL(unsigned, MPI_UNSIGNED);
MPIX_GET_MPI_DATATYPE_IMPL(unsigned long, MPI_UNSIGNED_LONG);
#define MPIX_PAIR(A,B) A,B
MPIX_GET_MPI_DATATYPE_IMPL(pair<MPIX_PAIR(float, int)>, MPI_FLOAT_INT);
MPIX_GET_MPI_DATATYPE_IMPL(pair<MPIX_PAIR(double, int)>, MPI_DOUBLE_INT);
MPIX_GET_MPI_DATATYPE_IMPL(pair<MPIX_PAIR(long double, int)>, MPI_LONG_DOUBLE_INT);
MPIX_GET_MPI_DATATYPE_IMPL(pair<MPIX_PAIR(long, int)>, MPI_LONG_INT);
MPIX_GET_MPI_DATATYPE_IMPL(pair<MPIX_PAIR(short, int)>, MPI_SHORT_INT);
MPIX_GET_MPI_DATATYPE_IMPL(pair<MPIX_PAIR(int, int)>, MPI_2INT);
#undef MPIX_PAIR

#if defined(MPI_LONG_LONG_INT) || (defined(MPI_VERSION) && MPI_VERSION >= 2)
MPIX_GET_MPI_DATATYPE_IMPL(long long, MPI_LONG_LONG_INT);
#endif

#if defined(MPI_UNSIGNED_LONG_LONG) || (defined(MPI_VERSION) && MPI_VERSION >= 2)
MPIX_GET_MPI_DATATYPE_IMPL(unsigned long long, MPI_UNSIGNED_LONG_LONG);
#endif

#if SIZEOF_BOOL == SIZEOF_CHAR
MPIX_GET_MPI_DATATYPE_IMPL(bool, MPI_CHAR);
#elif SIZEOF_BOOL == SIZEOF_SHORT
MPIX_GET_MPI_DATATYPE_IMPL(bool, MPI_SHORT);
#elif SIZEOF_BOOL == SIZEOF_INT
MPIX_GET_MPI_DATATYPE_IMPL(bool, MPI_INT);
#else
#error Cannot find MPI datatype for boolean
#endif


inline void init(int &argc, char **&argv)
{
    check(MPI_Init(&argc,&argv));
}


inline void init_thread(int &argc, char **&argv, int requested)
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


inline void finalize()
{
    check(MPI_Finalize());
}


inline MPI_Comm comm_dup(MPI_Comm orig)
{
    MPI_Comm dup;
    check(MPI_Comm_dup(orig, &dup));
    return dup;
}


inline void comm_free(MPI_Comm &comm)
{
    check(MPI_Comm_free(&comm));
}


inline void check(int errorcode)
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


inline int error_class(int errorcode)
{
    int eclass;
    MPI_Error_class(errorcode, &eclass);
    return eclass;
}


inline string error_string(int errorcode)
{
    string result;
    char estring[MPI_MAX_ERROR_STRING];
    int estring_length;

    MPI_Error_string(errorcode, estring, &estring_length);
    result.assign(estring, estring_length);

    return result;
}


inline int comm_rank(MPI_Comm comm)
{
    int result = 0;
    check(MPI_Comm_rank(comm, &result));
    return result;
}


inline int comm_size(MPI_Comm comm)
{
    int result = 0;
    check(MPI_Comm_size(comm, &result));
    return result;
}


inline void barrier(MPI_Comm comm)
{
    check(MPI_Barrier(comm));
}


template <typename T>
inline void bcast(T &object, int root, MPI_Comm comm)
{
    MPI_Datatype datatype = get_mpi_datatype(object);
    check(MPI_Bcast(&object, 1, datatype, root, comm));
}


template <typename T>
inline void bcast(vector<T> &object, int root, MPI_Comm comm)
{
    typedef typename vector<T>::size_type size_type;
    size_type size = object.size();
    MPI_Datatype datatype = get_mpi_datatype(object[0]);

    bcast(size, root, comm);
    if (comm_rank(comm) != root) {
        object.resize(size);
    }
    check(MPI_Bcast(&object[0], size, datatype, root, comm));
}


template <typename T>
inline void bcast(T *object, int size, int root, MPI_Comm comm)
{
    MPI_Datatype datatype = get_mpi_datatype(*object);
    check(MPI_Bcast(object, size, datatype, root, comm));
}


template <>
inline void bcast<string>(string &object, int root, MPI_Comm comm)
{
    typedef string::size_type size_type;
    size_type size = object.size();

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


template <>
inline void bcast<string>(vector<string> &object, int root, MPI_Comm comm)
{
    typedef vector<string>::size_type size_type;
    size_type size = object.size();

    bcast(size, root, comm);
    if (comm_rank(comm) != root) {
        object.resize(size);
    }
    for (vector<string>::iterator it=object.begin(); it!=object.end(); ++it) {
        bcast(*it, root, comm);
    }
}


inline vector<string> bcast(int argc, char **argv, MPI_Comm comm)
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


/* reduce */
template <class T>
inline void reduce(T &object, MPI_Op op, int root, MPI_Comm comm)
{
    MPI_Datatype datatype = get_mpi_datatype(object);

    if (root == comm_rank(comm)) {
        check(MPI_Reduce(MPI_IN_PLACE, &object, 1, datatype, op, root, comm));
    }
    else {
        check(MPI_Reduce(&object, NULL, 1, datatype, op, root, comm));
    }
}


template <typename T>
inline void reduce(vector<T> &object, MPI_Op op, int root, MPI_Comm comm)
{
    MPI_Datatype datatype = get_mpi_datatype(object[0]);

    if (root == comm_rank(comm)) {
        check(MPI_Reduce(MPI_IN_PLACE, &object[0], object.size(),
                    datatype, op, root, comm));
    }
    else {
        check(MPI_Reduce(&object[0], NULL, object.size(),
                    datatype, op, root, comm));
    }
}


template <typename T>
inline void reduce(T *object, int size, MPI_Op op, int root, MPI_Comm comm)
{
    MPI_Datatype datatype = get_mpi_datatype(*object);

    if (root == comm_rank(comm)) {
        check(MPI_Reduce(MPI_IN_PLACE, object, size, datatype, op, root, comm));
    }
    else {
        check(MPI_Reduce(object, NULL, size, datatype, op, root, comm));
    }
}



/* all reduce */
template <typename T>
inline void allreduce(T &object, MPI_Op op, MPI_Comm comm)
{
    MPI_Datatype datatype = get_mpi_datatype(object);
    check(MPI_Allreduce(MPI_IN_PLACE, &object, 1, datatype, op, comm));
}


template <typename T>
inline void allreduce(vector<T> &object, MPI_Op op, MPI_Comm comm)
{
    MPI_Datatype datatype = get_mpi_datatype(object[0]);
    check(MPI_Allreduce(MPI_IN_PLACE, &object[0], object.size(),
                datatype, op, comm));
}


template <typename T>
inline void allreduce(T *object, int size, MPI_Op op, MPI_Comm comm)
{
    MPI_Datatype datatype = get_mpi_datatype(*object);
    check(MPI_Allreduce(MPI_IN_PLACE, object, size, datatype, op, comm));
}


/* all to all */
template <typename T>
inline void alltoall(vector<T> &object, MPI_Comm comm)
{
    /* unfortunately, support of MPI_IN_PLACE for MPI_Alltoall is
     * extremely lacking, e.g., OpenMPI, Cray. */

    int size = comm_size(comm);
    int count = 0;
    MPI_Datatype datatype = get_mpi_datatype(object[0]);
    vector<T> object_copy(object);

    if (object.size() % size != 0) {
        cerr << "[" << comm_rank(MPI_COMM_WORLD) << "] MPI ERROR"
            << ": " << "mpix::alltoall has incorrect buffer size"
            << endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    count = object.size() / size;
    check(MPI_Alltoall(&object_copy[0], count, datatype,
                &object[0], count, datatype, comm));
}


/* all to all */
template <typename T>
inline void alltoall(vector<T> &sendbuf, vector<T> &recvbuf, MPI_Comm comm)
{
    /* unfortunately, support of MPI_IN_PLACE for MPI_Alltoall is
     * extremely lacking, e.g., OpenMPI, Cray. */

    int size = comm_size(comm);
    int count = 0;
    MPI_Datatype datatype = get_mpi_datatype(sendbuf[0]);

    if (sendbuf.size() % size != 0) {
        cerr << "[" << comm_rank(MPI_COMM_WORLD) << "] MPI ERROR"
            << ": " << "mpix::alltoall has incorrect buffer size"
            << endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    if (sendbuf.size() != recvbuf.size()) {
        cerr << "[" << comm_rank(MPI_COMM_WORLD) << "] MPI ERROR"
            << ": " << "mpix::alltoall mismatched buffer size"
            << endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    count = sendbuf.size() / size;
    check(MPI_Alltoall(&sendbuf[0], count, datatype,
                &recvbuf[0], count, datatype, comm));
}


/* gather */
template <typename T>
inline vector<T> gather(const T &object, int root, MPI_Comm comm)
{
    vector<T> result;
    MPI_Datatype datatype = get_mpi_datatype(object);

    if (comm_rank(comm) == root) {
        result.resize(comm_size(comm));
        check(MPI_Gather(&object, 1, datatype,
                    &result[0], 1, datatype, root, comm));
    }
    else {
        check(MPI_Gather(&object, 1, datatype,
                    NULL, 0, datatype, root, comm));
    }

    return result;
}


template <typename T>
inline vector<T> gather(const vector<T> &object, int root, MPI_Comm comm)
{
    vector<T> result;
    MPI_Datatype datatype = get_mpi_datatype(object[0]);
    int size = int(object.size());

    if (comm_rank(comm) == root) {
        result.resize(comm_size(comm) * size);
        check(MPI_Gather(&object[0], size, datatype,
                    &result[0], size, datatype, root, comm));
    }
    else {
        check(MPI_Gather(&object[0], size, datatype,
                    NULL, 0, datatype, root, comm));
    }

    return result;
}


template <typename T>
inline vector<T> gather(const T *object, int size, int root, MPI_Comm comm)
{
    vector<T> result;
    MPI_Datatype datatype = get_mpi_datatype(*object);

    if (comm_rank(comm) == root) {
        result.resize(comm_size(comm) * size);
        check(MPI_Gather(object, size, datatype,
                    &result[0], size, datatype, root, comm));
    }
    else {
        check(MPI_Gather(object, size, datatype,
                    NULL, size, datatype, root, comm));
    }

    return result;
}


template <>
inline vector<string> gather<string>(const string &object, int root, MPI_Comm comm)
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


template <typename T>
inline void gather(const T *sendbuf, int size, T *recvbuf, int root, MPI_Comm comm)
{
    MPI_Datatype datatype = get_mpi_datatype(sendbuf[0]);

    check(MPI_Gather(sendbuf, size, datatype,
                recvbuf, size, datatype, root, comm));
}


template <typename T>
inline void gather(const vector<T> &sendbuf, vector<T> &recvbuf, int root, MPI_Comm comm)
{
    MPI_Datatype datatype = get_mpi_datatype(sendbuf[0]);
    int size = int(sendbuf.size());

    if (comm_rank(comm) == root) {
        recvbuf.resize(comm_size(comm) * size);
    }
    else {
        recvbuf.clear();
    }

    check(MPI_Gather(&sendbuf[0], size, datatype,
                &recvbuf[0], size, datatype, root, comm));
}


/* synchronous printing */
template <class T>
inline void print_sync(const string &name, const T &what, MPI_Comm comm)
{
    vector<T> all_what = gather(what, 0, comm);

    if (0 == comm_rank(comm)) {
        for (int i = 0, size = comm_size(comm); i < size; ++i) {
            cout << "[" << i << "] " << name << "=" << all_what[i] << endl;
        }
    }

    barrier(comm);
}


template <typename T>
inline void print_sync(const string &name, const vector<T> &what, MPI_Comm comm)
{
    vector<T> all_what = gather(what, 0, comm);

    if (0 == comm_rank(comm)) {
        for (int i = 0, size = comm_size(comm); i < size; ++i) {
            cout << "[" << i << "] " << name << "={";
            cout << all_what[i*what.size()];
            for (int j = 1; j < what.size(); ++j) {
                cout << "," << all_what[i*what.size() + j];
            }
            cout << "}" << endl;
        }
    }

    barrier(comm);
}


template <typename T>
inline void print_sync(const string &name, const T *what, int size_, MPI_Comm comm)
{
    vector<T> all_what = gather(what, size_, 0, comm);

    if (0 == comm_rank(comm)) {
        for (int i = 0, size = comm_size(comm); i < size; ++i) {
            cout << "[" << i << "] " << name << "={";
            cout << all_what[i*size_];
            for (int j = 1; j < size_; ++j) {
                cout << "," << all_what[i*size_ + j];
            }
            cout << "}" << endl;
        }
    }

    barrier(comm);
}


template <class T>
inline void print_zero(const string &name, const T &what, MPI_Comm comm)
{
    if (0 == comm_rank(comm)) {
        cout << name << "=" << what << endl;
    }
    barrier(comm);
}




template <typename T>
inline void print_zero(const string &name, const vector<T> &what, MPI_Comm comm)
{
    if (0 == comm_rank(comm)) {
        cout << name << "={";
        cout << what[0];
        for (int j = 1; j < what.size(); ++j) {
            cout << "," << what[j];
        }
        cout << "}" << endl;
    }

    barrier(comm);
}


template <typename T>
inline void print_zero(const string &name, const T *what, int size, MPI_Comm comm)
{
    if (0 == comm_rank(comm)) {
        cout << name << "={";
        cout << what[0];
        for (int j = 1; j < size; ++j) {
            cout << "," << what[j];
        }
        cout << "}" << endl;
    }

    barrier(comm);
}


/* file reading */
inline MPI_Offset get_file_size(const string &file_name, MPI_Comm comm)
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
inline void read_file(
        const string &file_name,
        char *&file_buffer,
        MPI_Offset &file_size,
        MPI_Comm comm)
{
    read_file_mpiio(file_name, file_buffer, file_size, comm);
    //read_file_bcast(file_name, file_buffer, file_size, comm);
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
inline void read_file_bcast(
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
inline void read_file_mpiio(
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

#endif /* _MPIX_INL_H_ */
