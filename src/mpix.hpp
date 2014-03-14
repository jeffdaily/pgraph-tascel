/**
 * @file mpix.hpp
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 *
 * Templated and/or overloaded MPI functions.
 */
#ifndef _MPIX_H_
#define _MPIX_H_

#include <cstddef>
#include <map>
#include <string>
#include <vector>

#include <mpi.h>

using ::std::map;
using ::std::size_t;
using ::std::string;
using ::std::vector;

namespace mpix {

inline void init(int &argc, char **&argv);
inline void init_thread(int &argc, char **&argv, int requested);
inline void finalize();

inline void check(int errorcode);
inline int error_class(int errorcode);
inline string error_string(int errorcode);

inline MPI_Comm comm_dup(MPI_Comm comm);
inline void comm_free(MPI_Comm &comm);
inline int comm_rank(MPI_Comm comm);
inline int comm_size(MPI_Comm comm);
inline void barrier(MPI_Comm comm);

/* broadcast */
template <typename T> inline void           bcast(T &object, int root, MPI_Comm comm);
template <typename T> inline void           bcast(vector<T> &object, int root, MPI_Comm comm);
template <typename T> inline void           bcast(T *object, int size, int root, MPI_Comm comm);
template <>           inline void           bcast<string>(string &object, int root, MPI_Comm comm);
template <>           inline void           bcast<string>(vector<string> &object, int root, MPI_Comm comm);
                      inline vector<string> bcast(int argc, char **argv, MPI_Comm);

/* reduce */
template <typename T> inline void reduce(T &object, MPI_Op op, int root, MPI_Comm comm);
template <typename T> inline void reduce(vector<T> &object, MPI_Op op, int root, MPI_Comm comm);
template <typename T> inline void reduce(T *object, int size, MPI_Op op, int root, MPI_Comm comm);

/* all reduce */
template <typename T> inline void allreduce(T &object, MPI_Op op, MPI_Comm comm);
template <typename T> inline void allreduce(vector<T> &object, MPI_Op op, MPI_Comm comm);
template <typename T> inline void allreduce(T *object, int size, MPI_Op op, MPI_Comm comm);

/* all to all */
template <typename T> inline void alltoall(vector<T> &sendrecv, MPI_Comm comm);
template <typename T> inline void alltoall(vector<T> &sendbuf, vector<T> &recvbuf, MPI_Comm comm);

/* gather */
template <typename T> inline vector<T>      gather(const T &object, int root, MPI_Comm comm);
template <typename T> inline vector<T>      gather(const vector<T> &object, int root, MPI_Comm comm);
template <typename T> inline vector<T>      gather(const T *object, int size, int root, MPI_Comm comm);
template <>           inline vector<string> gather<string>(const string &object, int root, MPI_Comm comm);

template <typename T> inline void           gather(const T *sendbuf, int size, T *recvbuf, int root, MPI_Comm comm);
template <typename T> inline void           gather(const vector<T> &sendbuf, vector<T> &recvbuf, int root, MPI_Comm comm);

/* synchronous printing */
template <typename T> inline void print_sync(const string &label, const T &object, MPI_Comm comm);
template <typename T> inline void print_sync(const string &label, const vector<T> &object, MPI_Comm comm);
template <typename T> inline void print_sync(const string &label, const T *object, int size, MPI_Comm comm);

template <typename T> inline void print_zero(const string &label, const T &object, MPI_Comm comm);
template <typename T> inline void print_zero(const string &label, const vector<T> &object, MPI_Comm comm);
template <typename T> inline void print_zero(const string &label, const T *object, int size, MPI_Comm comm);

/* file reading */
inline MPI_Offset get_file_size(const string &file_name, MPI_Comm comm);
inline void read_file(const string &file_name, char *&file_buffer, MPI_Offset &file_size, MPI_Comm comm);
inline void read_file_bcast(const string &file_name, char *&file_buffer, MPI_Offset &file_size, MPI_Comm comm);
inline void read_file_mpiio(const string &file_name, char *&file_buffer, MPI_Offset &file_size, MPI_Comm comm);

/* data types */
inline MPI_Datatype type_contiguous(int count, MPI_Datatype oldtype);
inline MPI_Datatype type_create_struct(int count, const int blocklengths[], const MPI_Aint displacements[], const MPI_Datatype types[]);
inline void type_commit(MPI_Datatype &type);
inline void type_free(MPI_Datatype &type);
template <typename T> inline MPI_Datatype build_mpi_datatype(const T& object);
template <typename T> inline MPI_Datatype get_mpi_datatype(const T& object);
inline map<string,MPI_Datatype>& get_custom_mpi_datatypes();
inline void add_custom_mpi_datatype(const string &name, MPI_Datatype type);


} /* namespace mpix */

#include "mpix-inl.hpp"

#endif /* _MPIX_H_ */
