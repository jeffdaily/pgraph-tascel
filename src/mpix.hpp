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

#include <mpi.h>

#include <cstddef>
#include <map>
#include <string>
#include <utility>
#include <vector>

using ::std::map;
using ::std::pair;
using ::std::size_t;
using ::std::string;
using ::std::vector;

namespace mpix {

void init(int &argc, char **&argv);
void init_thread(int &argc, char **&argv, int requested);
void finalize();

void check(int errorcode);
int error_class(int errorcode);
string error_string(int errorcode);

MPI_Comm comm_dup(MPI_Comm comm);
void comm_free(MPI_Comm &comm);
int comm_rank(MPI_Comm comm);
int comm_size(MPI_Comm comm);
void barrier(MPI_Comm comm);

/* data types */
void type_contiguous(int count, MPI_Datatype oldtype, MPI_Datatype &newtype);
void type_create_struct(int count, int blocklengths[], MPI_Aint displacements[], MPI_Datatype types[], MPI_Datatype &newtype);
void type_commit(MPI_Datatype &type);
void type_free(MPI_Datatype &type);

MPI_Datatype get_mpi_datatype(char               object);
MPI_Datatype get_mpi_datatype(short              object);
MPI_Datatype get_mpi_datatype(int                object);
MPI_Datatype get_mpi_datatype(long               object);
MPI_Datatype get_mpi_datatype(long long          object);
MPI_Datatype get_mpi_datatype(signed char        object);
MPI_Datatype get_mpi_datatype(unsigned char      object);
MPI_Datatype get_mpi_datatype(unsigned short     object);
MPI_Datatype get_mpi_datatype(unsigned int       object);
MPI_Datatype get_mpi_datatype(unsigned long      object);
MPI_Datatype get_mpi_datatype(unsigned long long object);
MPI_Datatype get_mpi_datatype(float              object);
MPI_Datatype get_mpi_datatype(double             object);
MPI_Datatype get_mpi_datatype(long double        object);
MPI_Datatype get_mpi_datatype(pair<float,       int> object);
MPI_Datatype get_mpi_datatype(pair<double,      int> object);
MPI_Datatype get_mpi_datatype(pair<long,        int> object);
MPI_Datatype get_mpi_datatype(pair<short,       int> object);
MPI_Datatype get_mpi_datatype(pair<int,         int> object);
MPI_Datatype get_mpi_datatype(pair<long double, int> object);

/* broadcast */
void bcast(char               & object, int root, MPI_Comm comm);
void bcast(short              & object, int root, MPI_Comm comm);
void bcast(int                & object, int root, MPI_Comm comm);
void bcast(long               & object, int root, MPI_Comm comm);
void bcast(long long          & object, int root, MPI_Comm comm);
void bcast(signed char        & object, int root, MPI_Comm comm);
void bcast(unsigned char      & object, int root, MPI_Comm comm);
void bcast(unsigned short     & object, int root, MPI_Comm comm);
void bcast(unsigned int       & object, int root, MPI_Comm comm);
void bcast(unsigned long      & object, int root, MPI_Comm comm);
void bcast(unsigned long long & object, int root, MPI_Comm comm);
void bcast(float              & object, int root, MPI_Comm comm);
void bcast(double             & object, int root, MPI_Comm comm);
void bcast(long double        & object, int root, MPI_Comm comm);
void bcast(string             & object, int root, MPI_Comm comm);

void bcast(char               * object, int count, int root, MPI_Comm comm);
void bcast(short              * object, int count, int root, MPI_Comm comm);
void bcast(int                * object, int count, int root, MPI_Comm comm);
void bcast(long               * object, int count, int root, MPI_Comm comm);
void bcast(long long          * object, int count, int root, MPI_Comm comm);
void bcast(signed char        * object, int count, int root, MPI_Comm comm);
void bcast(unsigned char      * object, int count, int root, MPI_Comm comm);
void bcast(unsigned short     * object, int count, int root, MPI_Comm comm);
void bcast(unsigned int       * object, int count, int root, MPI_Comm comm);
void bcast(unsigned long      * object, int count, int root, MPI_Comm comm);
void bcast(unsigned long long * object, int count, int root, MPI_Comm comm);
void bcast(float              * object, int count, int root, MPI_Comm comm);
void bcast(double             * object, int count, int root, MPI_Comm comm);
void bcast(long double        * object, int count, int root, MPI_Comm comm);

void bcast(vector<char>               & object, int root, MPI_Comm comm);
void bcast(vector<short>              & object, int root, MPI_Comm comm);
void bcast(vector<int>                & object, int root, MPI_Comm comm);
void bcast(vector<long>               & object, int root, MPI_Comm comm);
void bcast(vector<long long>          & object, int root, MPI_Comm comm);
void bcast(vector<signed char>        & object, int root, MPI_Comm comm);
void bcast(vector<unsigned char>      & object, int root, MPI_Comm comm);
void bcast(vector<unsigned short>     & object, int root, MPI_Comm comm);
void bcast(vector<unsigned int>       & object, int root, MPI_Comm comm);
void bcast(vector<unsigned long>      & object, int root, MPI_Comm comm);
void bcast(vector<unsigned long long> & object, int root, MPI_Comm comm);
void bcast(vector<float>              & object, int root, MPI_Comm comm);
void bcast(vector<double>             & object, int root, MPI_Comm comm);
void bcast(vector<long double>        & object, int root, MPI_Comm comm);
void bcast(vector<string>             & object, int root, MPI_Comm comm);

vector<string> bcast(int argc, char **argv, MPI_Comm comm);

/* reduce */
void reduce(char               & object, MPI_Op op, int root, MPI_Comm comm);
void reduce(short              & object, MPI_Op op, int root, MPI_Comm comm);
void reduce(int                & object, MPI_Op op, int root, MPI_Comm comm);
void reduce(long               & object, MPI_Op op, int root, MPI_Comm comm);
void reduce(long long          & object, MPI_Op op, int root, MPI_Comm comm);
void reduce(signed char        & object, MPI_Op op, int root, MPI_Comm comm);
void reduce(unsigned char      & object, MPI_Op op, int root, MPI_Comm comm);
void reduce(unsigned short     & object, MPI_Op op, int root, MPI_Comm comm);
void reduce(unsigned int       & object, MPI_Op op, int root, MPI_Comm comm);
void reduce(unsigned long      & object, MPI_Op op, int root, MPI_Comm comm);
void reduce(unsigned long long & object, MPI_Op op, int root, MPI_Comm comm);
void reduce(float              & object, MPI_Op op, int root, MPI_Comm comm);
void reduce(double             & object, MPI_Op op, int root, MPI_Comm comm);
void reduce(long double        & object, MPI_Op op, int root, MPI_Comm comm);

void reduce(char               * object, int count, MPI_Op op, int root, MPI_Comm comm);
void reduce(short              * object, int count, MPI_Op op, int root, MPI_Comm comm);
void reduce(int                * object, int count, MPI_Op op, int root, MPI_Comm comm);
void reduce(long               * object, int count, MPI_Op op, int root, MPI_Comm comm);
void reduce(long long          * object, int count, MPI_Op op, int root, MPI_Comm comm);
void reduce(signed char        * object, int count, MPI_Op op, int root, MPI_Comm comm);
void reduce(unsigned char      * object, int count, MPI_Op op, int root, MPI_Comm comm);
void reduce(unsigned short     * object, int count, MPI_Op op, int root, MPI_Comm comm);
void reduce(unsigned int       * object, int count, MPI_Op op, int root, MPI_Comm comm);
void reduce(unsigned long      * object, int count, MPI_Op op, int root, MPI_Comm comm);
void reduce(unsigned long long * object, int count, MPI_Op op, int root, MPI_Comm comm);
void reduce(float              * object, int count, MPI_Op op, int root, MPI_Comm comm);
void reduce(double             * object, int count, MPI_Op op, int root, MPI_Comm comm);
void reduce(long double        * object, int count, MPI_Op op, int root, MPI_Comm comm);

void reduce(vector<char>               & object, MPI_Op op, int root, MPI_Comm comm);
void reduce(vector<short>              & object, MPI_Op op, int root, MPI_Comm comm);
void reduce(vector<int>                & object, MPI_Op op, int root, MPI_Comm comm);
void reduce(vector<long>               & object, MPI_Op op, int root, MPI_Comm comm);
void reduce(vector<long long>          & object, MPI_Op op, int root, MPI_Comm comm);
void reduce(vector<signed char>        & object, MPI_Op op, int root, MPI_Comm comm);
void reduce(vector<unsigned char>      & object, MPI_Op op, int root, MPI_Comm comm);
void reduce(vector<unsigned short>     & object, MPI_Op op, int root, MPI_Comm comm);
void reduce(vector<unsigned int>       & object, MPI_Op op, int root, MPI_Comm comm);
void reduce(vector<unsigned long>      & object, MPI_Op op, int root, MPI_Comm comm);
void reduce(vector<unsigned long long> & object, MPI_Op op, int root, MPI_Comm comm);
void reduce(vector<float>              & object, MPI_Op op, int root, MPI_Comm comm);
void reduce(vector<double>             & object, MPI_Op op, int root, MPI_Comm comm);
void reduce(vector<long double>        & object, MPI_Op op, int root, MPI_Comm comm);

/* all reduce */
void allreduce(char               & object, MPI_Op op, MPI_Comm comm);
void allreduce(short              & object, MPI_Op op, MPI_Comm comm);
void allreduce(int                & object, MPI_Op op, MPI_Comm comm);
void allreduce(long               & object, MPI_Op op, MPI_Comm comm);
void allreduce(long long          & object, MPI_Op op, MPI_Comm comm);
void allreduce(signed char        & object, MPI_Op op, MPI_Comm comm);
void allreduce(unsigned char      & object, MPI_Op op, MPI_Comm comm);
void allreduce(unsigned short     & object, MPI_Op op, MPI_Comm comm);
void allreduce(unsigned int       & object, MPI_Op op, MPI_Comm comm);
void allreduce(unsigned long      & object, MPI_Op op, MPI_Comm comm);
void allreduce(unsigned long long & object, MPI_Op op, MPI_Comm comm);
void allreduce(float              & object, MPI_Op op, MPI_Comm comm);
void allreduce(double             & object, MPI_Op op, MPI_Comm comm);
void allreduce(long double        & object, MPI_Op op, MPI_Comm comm);

void allreduce(char               * object, int count, MPI_Op op, MPI_Comm comm);
void allreduce(short              * object, int count, MPI_Op op, MPI_Comm comm);
void allreduce(int                * object, int count, MPI_Op op, MPI_Comm comm);
void allreduce(long               * object, int count, MPI_Op op, MPI_Comm comm);
void allreduce(long long          * object, int count, MPI_Op op, MPI_Comm comm);
void allreduce(signed char        * object, int count, MPI_Op op, MPI_Comm comm);
void allreduce(unsigned char      * object, int count, MPI_Op op, MPI_Comm comm);
void allreduce(unsigned short     * object, int count, MPI_Op op, MPI_Comm comm);
void allreduce(unsigned int       * object, int count, MPI_Op op, MPI_Comm comm);
void allreduce(unsigned long      * object, int count, MPI_Op op, MPI_Comm comm);
void allreduce(unsigned long long * object, int count, MPI_Op op, MPI_Comm comm);
void allreduce(float              * object, int count, MPI_Op op, MPI_Comm comm);
void allreduce(double             * object, int count, MPI_Op op, MPI_Comm comm);
void allreduce(long double        * object, int count, MPI_Op op, MPI_Comm comm);

void allreduce(vector<char>               & object, MPI_Op op, MPI_Comm comm);
void allreduce(vector<short>              & object, MPI_Op op, MPI_Comm comm);
void allreduce(vector<int>                & object, MPI_Op op, MPI_Comm comm);
void allreduce(vector<long>               & object, MPI_Op op, MPI_Comm comm);
void allreduce(vector<long long>          & object, MPI_Op op, MPI_Comm comm);
void allreduce(vector<signed char>        & object, MPI_Op op, MPI_Comm comm);
void allreduce(vector<unsigned char>      & object, MPI_Op op, MPI_Comm comm);
void allreduce(vector<unsigned short>     & object, MPI_Op op, MPI_Comm comm);
void allreduce(vector<unsigned int>       & object, MPI_Op op, MPI_Comm comm);
void allreduce(vector<unsigned long>      & object, MPI_Op op, MPI_Comm comm);
void allreduce(vector<unsigned long long> & object, MPI_Op op, MPI_Comm comm);
void allreduce(vector<float>              & object, MPI_Op op, MPI_Comm comm);
void allreduce(vector<double>             & object, MPI_Op op, MPI_Comm comm);
void allreduce(vector<long double>        & object, MPI_Op op, MPI_Comm comm);

/* all to all */
void alltoall(vector<char>               & sendrecv, MPI_Comm comm);
void alltoall(vector<short>              & sendrecv, MPI_Comm comm);
void alltoall(vector<int>                & sendrecv, MPI_Comm comm);
void alltoall(vector<long>               & sendrecv, MPI_Comm comm);
void alltoall(vector<long long>          & sendrecv, MPI_Comm comm);
void alltoall(vector<signed char>        & sendrecv, MPI_Comm comm);
void alltoall(vector<unsigned char>      & sendrecv, MPI_Comm comm);
void alltoall(vector<unsigned short>     & sendrecv, MPI_Comm comm);
void alltoall(vector<unsigned int>       & sendrecv, MPI_Comm comm);
void alltoall(vector<unsigned long>      & sendrecv, MPI_Comm comm);
void alltoall(vector<unsigned long long> & sendrecv, MPI_Comm comm);
void alltoall(vector<float>              & sendrecv, MPI_Comm comm);
void alltoall(vector<double>             & sendrecv, MPI_Comm comm);
void alltoall(vector<long double>        & sendrecv, MPI_Comm comm);

void alltoall(vector<char>               & send, vector<char>               & recv, MPI_Comm comm);
void alltoall(vector<short>              & send, vector<short>              & recv, MPI_Comm comm);
void alltoall(vector<int>                & send, vector<int>                & recv, MPI_Comm comm);
void alltoall(vector<long>               & send, vector<long>               & recv, MPI_Comm comm);
void alltoall(vector<long long>          & send, vector<long long>          & recv, MPI_Comm comm);
void alltoall(vector<signed char>        & send, vector<signed char>        & recv, MPI_Comm comm);
void alltoall(vector<unsigned char>      & send, vector<unsigned char>      & recv, MPI_Comm comm);
void alltoall(vector<unsigned short>     & send, vector<unsigned short>     & recv, MPI_Comm comm);
void alltoall(vector<unsigned int>       & send, vector<unsigned int>       & recv, MPI_Comm comm);
void alltoall(vector<unsigned long>      & send, vector<unsigned long>      & recv, MPI_Comm comm);
void alltoall(vector<unsigned long long> & send, vector<unsigned long long> & recv, MPI_Comm comm);
void alltoall(vector<float>              & send, vector<float>              & recv, MPI_Comm comm);
void alltoall(vector<double>             & send, vector<double>             & recv, MPI_Comm comm);
void alltoall(vector<long double>        & send, vector<long double>        & recv, MPI_Comm comm);

/* gather */
vector<char>               gather(char               object, int root, MPI_Comm comm);
vector<short>              gather(short              object, int root, MPI_Comm comm);
vector<int>                gather(int                object, int root, MPI_Comm comm);
vector<long>               gather(long               object, int root, MPI_Comm comm);
vector<long long>          gather(long long          object, int root, MPI_Comm comm);
vector<signed char>        gather(signed char        object, int root, MPI_Comm comm);
vector<unsigned char>      gather(unsigned char      object, int root, MPI_Comm comm);
vector<unsigned short>     gather(unsigned short     object, int root, MPI_Comm comm);
vector<unsigned int>       gather(unsigned int       object, int root, MPI_Comm comm);
vector<unsigned long>      gather(unsigned long      object, int root, MPI_Comm comm);
vector<unsigned long long> gather(unsigned long long object, int root, MPI_Comm comm);
vector<float>              gather(float              object, int root, MPI_Comm comm);
vector<double>             gather(double             object, int root, MPI_Comm comm);
vector<long double>        gather(long double        object, int root, MPI_Comm comm);

void                       gather(char               * sendbuf, int sendcount, char               * recvbuf, int recvcount, int root, MPI_Comm comm);
void                       gather(short              * sendbuf, int sendcount, short              * recvbuf, int recvcount, int root, MPI_Comm comm);
void                       gather(int                * sendbuf, int sendcount, int                * recvbuf, int recvcount, int root, MPI_Comm comm);
void                       gather(long               * sendbuf, int sendcount, long               * recvbuf, int recvcount, int root, MPI_Comm comm);
void                       gather(long long          * sendbuf, int sendcount, long long          * recvbuf, int recvcount, int root, MPI_Comm comm);
void                       gather(signed char        * sendbuf, int sendcount, signed char        * recvbuf, int recvcount, int root, MPI_Comm comm);
void                       gather(unsigned char      * sendbuf, int sendcount, unsigned char      * recvbuf, int recvcount, int root, MPI_Comm comm);
void                       gather(unsigned short     * sendbuf, int sendcount, unsigned short     * recvbuf, int recvcount, int root, MPI_Comm comm);
void                       gather(unsigned int       * sendbuf, int sendcount, unsigned int       * recvbuf, int recvcount, int root, MPI_Comm comm);
void                       gather(unsigned long      * sendbuf, int sendcount, unsigned long      * recvbuf, int recvcount, int root, MPI_Comm comm);
void                       gather(unsigned long long * sendbuf, int sendcount, unsigned long long * recvbuf, int recvcount, int root, MPI_Comm comm);
void                       gather(float              * sendbuf, int sendcount, float              * recvbuf, int recvcount, int root, MPI_Comm comm);
void                       gather(double             * sendbuf, int sendcount, double             * recvbuf, int recvcount, int root, MPI_Comm comm);
void                       gather(long double        * sendbuf, int sendcount, long double        * recvbuf, int recvcount, int root, MPI_Comm comm);

vector<char>               gather(char               * object, int count, int root, MPI_Comm comm);
vector<short>              gather(short              * object, int count, int root, MPI_Comm comm);
vector<int>                gather(int                * object, int count, int root, MPI_Comm comm);
vector<long>               gather(long               * object, int count, int root, MPI_Comm comm);
vector<long long>          gather(long long          * object, int count, int root, MPI_Comm comm);
vector<signed char>        gather(signed char        * object, int count, int root, MPI_Comm comm);
vector<unsigned char>      gather(unsigned char      * object, int count, int root, MPI_Comm comm);
vector<unsigned short>     gather(unsigned short     * object, int count, int root, MPI_Comm comm);
vector<unsigned int>       gather(unsigned int       * object, int count, int root, MPI_Comm comm);
vector<unsigned long>      gather(unsigned long      * object, int count, int root, MPI_Comm comm);
vector<unsigned long long> gather(unsigned long long * object, int count, int root, MPI_Comm comm);
vector<float>              gather(float              * object, int count, int root, MPI_Comm comm);
vector<double>             gather(double             * object, int count, int root, MPI_Comm comm);
vector<long double>        gather(long double        * object, int count, int root, MPI_Comm comm);

vector<char>               gather(vector<char>               & object, int root, MPI_Comm comm);
vector<short>              gather(vector<short>              & object, int root, MPI_Comm comm);
vector<int>                gather(vector<int>                & object, int root, MPI_Comm comm);
vector<long>               gather(vector<long>               & object, int root, MPI_Comm comm);
vector<long long>          gather(vector<long long>          & object, int root, MPI_Comm comm);
vector<signed char>        gather(vector<signed char>        & object, int root, MPI_Comm comm);
vector<unsigned char>      gather(vector<unsigned char>      & object, int root, MPI_Comm comm);
vector<unsigned short>     gather(vector<unsigned short>     & object, int root, MPI_Comm comm);
vector<unsigned int>       gather(vector<unsigned int>       & object, int root, MPI_Comm comm);
vector<unsigned long>      gather(vector<unsigned long>      & object, int root, MPI_Comm comm);
vector<unsigned long long> gather(vector<unsigned long long> & object, int root, MPI_Comm comm);
vector<float>              gather(vector<float>              & object, int root, MPI_Comm comm);
vector<double>             gather(vector<double>             & object, int root, MPI_Comm comm);
vector<long double>        gather(vector<long double>        & object, int root, MPI_Comm comm);

vector<string> gather(string &object, int root, MPI_Comm comm);

/* synchronous printing */
void print_sync(const string &label, char               object, MPI_Comm comm);
void print_sync(const string &label, short              object, MPI_Comm comm);
void print_sync(const string &label, int                object, MPI_Comm comm);
void print_sync(const string &label, long               object, MPI_Comm comm);
void print_sync(const string &label, long long          object, MPI_Comm comm);
void print_sync(const string &label, signed char        object, MPI_Comm comm);
void print_sync(const string &label, unsigned char      object, MPI_Comm comm);
void print_sync(const string &label, unsigned short     object, MPI_Comm comm);
void print_sync(const string &label, unsigned int       object, MPI_Comm comm);
void print_sync(const string &label, unsigned long      object, MPI_Comm comm);
void print_sync(const string &label, unsigned long long object, MPI_Comm comm);
void print_sync(const string &label, float              object, MPI_Comm comm);
void print_sync(const string &label, double             object, MPI_Comm comm);
void print_sync(const string &label, long double        object, MPI_Comm comm);

void print_sync(const string &label, vector<char>               & object, MPI_Comm comm);
void print_sync(const string &label, vector<short>              & object, MPI_Comm comm);
void print_sync(const string &label, vector<int>                & object, MPI_Comm comm);
void print_sync(const string &label, vector<long>               & object, MPI_Comm comm);
void print_sync(const string &label, vector<long long>          & object, MPI_Comm comm);
void print_sync(const string &label, vector<signed char>        & object, MPI_Comm comm);
void print_sync(const string &label, vector<unsigned char>      & object, MPI_Comm comm);
void print_sync(const string &label, vector<unsigned short>     & object, MPI_Comm comm);
void print_sync(const string &label, vector<unsigned int>       & object, MPI_Comm comm);
void print_sync(const string &label, vector<unsigned long>      & object, MPI_Comm comm);
void print_sync(const string &label, vector<unsigned long long> & object, MPI_Comm comm);
void print_sync(const string &label, vector<float>              & object, MPI_Comm comm);
void print_sync(const string &label, vector<double>             & object, MPI_Comm comm);
void print_sync(const string &label, vector<long double>        & object, MPI_Comm comm);

void print_sync(const string &label, char               * object, int size, MPI_Comm comm);
void print_sync(const string &label, short              * object, int size, MPI_Comm comm);
void print_sync(const string &label, int                * object, int size, MPI_Comm comm);
void print_sync(const string &label, long               * object, int size, MPI_Comm comm);
void print_sync(const string &label, long long          * object, int size, MPI_Comm comm);
void print_sync(const string &label, signed char        * object, int size, MPI_Comm comm);
void print_sync(const string &label, unsigned char      * object, int size, MPI_Comm comm);
void print_sync(const string &label, unsigned short     * object, int size, MPI_Comm comm);
void print_sync(const string &label, unsigned int       * object, int size, MPI_Comm comm);
void print_sync(const string &label, unsigned long      * object, int size, MPI_Comm comm);
void print_sync(const string &label, unsigned long long * object, int size, MPI_Comm comm);
void print_sync(const string &label, float              * object, int size, MPI_Comm comm);
void print_sync(const string &label, double             * object, int size, MPI_Comm comm);
void print_sync(const string &label, long double        * object, int size, MPI_Comm comm);

void print_zero(const string &label, MPI_Comm comm);

void print_zero(const string &label, char               object, MPI_Comm comm);
void print_zero(const string &label, short              object, MPI_Comm comm);
void print_zero(const string &label, int                object, MPI_Comm comm);
void print_zero(const string &label, long               object, MPI_Comm comm);
void print_zero(const string &label, long long          object, MPI_Comm comm);
void print_zero(const string &label, signed char        object, MPI_Comm comm);
void print_zero(const string &label, unsigned char      object, MPI_Comm comm);
void print_zero(const string &label, unsigned short     object, MPI_Comm comm);
void print_zero(const string &label, unsigned int       object, MPI_Comm comm);
void print_zero(const string &label, unsigned long      object, MPI_Comm comm);
void print_zero(const string &label, unsigned long long object, MPI_Comm comm);
void print_zero(const string &label, float              object, MPI_Comm comm);
void print_zero(const string &label, double             object, MPI_Comm comm);
void print_zero(const string &label, long double        object, MPI_Comm comm);

void print_zero(const string &label, vector<char>              & object, MPI_Comm comm);
void print_zero(const string &label, vector<short>             & object, MPI_Comm comm);
void print_zero(const string &label, vector<int>               & object, MPI_Comm comm);
void print_zero(const string &label, vector<long>              & object, MPI_Comm comm);
void print_zero(const string &label, vector<long long>         & object, MPI_Comm comm);
void print_zero(const string &label, vector<signed char>       & object, MPI_Comm comm);
void print_zero(const string &label, vector<unsigned char>     & object, MPI_Comm comm);
void print_zero(const string &label, vector<unsigned short>    & object, MPI_Comm comm);
void print_zero(const string &label, vector<unsigned int>      & object, MPI_Comm comm);
void print_zero(const string &label, vector<unsigned long>     & object, MPI_Comm comm);
void print_zero(const string &label, vector<unsigned long long>& object, MPI_Comm comm);
void print_zero(const string &label, vector<float>             & object, MPI_Comm comm);
void print_zero(const string &label, vector<double>            & object, MPI_Comm comm);
void print_zero(const string &label, vector<long double>       & object, MPI_Comm comm);

void print_zero(const string &label, char               * object, int size, MPI_Comm comm);
void print_zero(const string &label, short              * object, int size, MPI_Comm comm);
void print_zero(const string &label, int                * object, int size, MPI_Comm comm);
void print_zero(const string &label, long               * object, int size, MPI_Comm comm);
void print_zero(const string &label, long long          * object, int size, MPI_Comm comm);
void print_zero(const string &label, signed char        * object, int size, MPI_Comm comm);
void print_zero(const string &label, unsigned char      * object, int size, MPI_Comm comm);
void print_zero(const string &label, unsigned short     * object, int size, MPI_Comm comm);
void print_zero(const string &label, unsigned int       * object, int size, MPI_Comm comm);
void print_zero(const string &label, unsigned long      * object, int size, MPI_Comm comm);
void print_zero(const string &label, unsigned long long * object, int size, MPI_Comm comm);
void print_zero(const string &label, float              * object, int size, MPI_Comm comm);
void print_zero(const string &label, double             * object, int size, MPI_Comm comm);
void print_zero(const string &label, long double        * object, int size, MPI_Comm comm);

/* file reading */
MPI_Offset get_file_size(const string &file_name, MPI_Comm comm);
void read_file(const string &file_name, char *&file_buffer, MPI_Offset &file_size, MPI_Comm comm);
void read_file_bcast(const string &file_name, char *&file_buffer, MPI_Offset &file_size, MPI_Comm comm);
void read_file_mpiio(const string &file_name, char *&file_buffer, MPI_Offset &file_size, MPI_Comm comm);

} /* namespace mpix */

#endif /* _MPIX_H_ */
