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
MPI_Datatype type_contiguous(int count, MPI_Datatype oldtype);
MPI_Datatype type_create_struct(int count, int blocklengths[], MPI_Aint displacements[], MPI_Datatype types[]);
void type_commit(MPI_Datatype &type);
void type_free(MPI_Datatype &type);
MPI_Datatype get_mpi_datatype(          bool object);
MPI_Datatype get_mpi_datatype(          char object);
MPI_Datatype get_mpi_datatype(   signed char object);
MPI_Datatype get_mpi_datatype( unsigned char object);
MPI_Datatype get_mpi_datatype(         short object);
MPI_Datatype get_mpi_datatype(           int object);
MPI_Datatype get_mpi_datatype(          long object);
MPI_Datatype get_mpi_datatype(     long long object);
MPI_Datatype get_mpi_datatype(unsigned short object);
MPI_Datatype get_mpi_datatype(  unsigned int object);
MPI_Datatype get_mpi_datatype( unsigned long object);
MPI_Datatype get_mpi_datatype(         float object);
MPI_Datatype get_mpi_datatype(        double object);
MPI_Datatype get_mpi_datatype(   long double object);
MPI_Datatype get_mpi_datatype(pair<float,       int> object);
MPI_Datatype get_mpi_datatype(pair<double,      int> object);
MPI_Datatype get_mpi_datatype(pair<long,        int> object);
MPI_Datatype get_mpi_datatype(pair<short,       int> object);
MPI_Datatype get_mpi_datatype(pair<int,         int> object);
MPI_Datatype get_mpi_datatype(pair<long double, int> object);

/* broadcast */
void bcast(              char & object, int root, MPI_Comm comm);
void bcast(       signed char & object, int root, MPI_Comm comm);
void bcast(     unsigned char & object, int root, MPI_Comm comm);
void bcast(             short & object, int root, MPI_Comm comm);
void bcast(               int & object, int root, MPI_Comm comm);
void bcast(              long & object, int root, MPI_Comm comm);
void bcast(         long long & object, int root, MPI_Comm comm);
void bcast(    unsigned short & object, int root, MPI_Comm comm);
void bcast(      unsigned int & object, int root, MPI_Comm comm);
void bcast(     unsigned long & object, int root, MPI_Comm comm);
void bcast(unsigned long long & object, int root, MPI_Comm comm);
void bcast(            string & object, int root, MPI_Comm comm);

void bcast(              char * object, int count, int root, MPI_Comm comm);
void bcast(       signed char * object, int count, int root, MPI_Comm comm);
void bcast(     unsigned char * object, int count, int root, MPI_Comm comm);
void bcast(             short * object, int count, int root, MPI_Comm comm);
void bcast(               int * object, int count, int root, MPI_Comm comm);
void bcast(              long * object, int count, int root, MPI_Comm comm);
void bcast(         long long * object, int count, int root, MPI_Comm comm);
void bcast(    unsigned short * object, int count, int root, MPI_Comm comm);
void bcast(      unsigned int * object, int count, int root, MPI_Comm comm);
void bcast(     unsigned long * object, int count, int root, MPI_Comm comm);
void bcast(unsigned long long * object, int count, int root, MPI_Comm comm);
void bcast(            string * object, int count, int root, MPI_Comm comm);

void bcast(vector<              char> & object, int root, MPI_Comm comm);
void bcast(vector<       signed char> & object, int root, MPI_Comm comm);
void bcast(vector<     unsigned char> & object, int root, MPI_Comm comm);
void bcast(vector<             short> & object, int root, MPI_Comm comm);
void bcast(vector<               int> & object, int root, MPI_Comm comm);
void bcast(vector<              long> & object, int root, MPI_Comm comm);
void bcast(vector<         long long> & object, int root, MPI_Comm comm);
void bcast(vector<    unsigned short> & object, int root, MPI_Comm comm);
void bcast(vector<      unsigned int> & object, int root, MPI_Comm comm);
void bcast(vector<     unsigned long> & object, int root, MPI_Comm comm);
void bcast(vector<unsigned long long> & object, int root, MPI_Comm comm);
void bcast(vector<            string> & object, int root, MPI_Comm comm);

vector<string> bcast(int argc, char **argv, MPI_Comm comm);

/* reduce */
void reduce(              char & object, MPI_Op op, int root, MPI_Comm comm);
void reduce(       signed char & object, MPI_Op op, int root, MPI_Comm comm);
void reduce(     unsigned char & object, MPI_Op op, int root, MPI_Comm comm);
void reduce(             short & object, MPI_Op op, int root, MPI_Comm comm);
void reduce(               int & object, MPI_Op op, int root, MPI_Comm comm);
void reduce(              long & object, MPI_Op op, int root, MPI_Comm comm);
void reduce(         long long & object, MPI_Op op, int root, MPI_Comm comm);
void reduce(    unsigned short & object, MPI_Op op, int root, MPI_Comm comm);
void reduce(      unsigned int & object, MPI_Op op, int root, MPI_Comm comm);
void reduce(     unsigned long & object, MPI_Op op, int root, MPI_Comm comm);
void reduce(unsigned long long & object, MPI_Op op, int root, MPI_Comm comm);

void reduce(              char * object, int count, MPI_Op op, int root, MPI_Comm comm);
void reduce(       signed char * object, int count, MPI_Op op, int root, MPI_Comm comm);
void reduce(     unsigned char * object, int count, MPI_Op op, int root, MPI_Comm comm);
void reduce(             short * object, int count, MPI_Op op, int root, MPI_Comm comm);
void reduce(               int * object, int count, MPI_Op op, int root, MPI_Comm comm);
void reduce(              long * object, int count, MPI_Op op, int root, MPI_Comm comm);
void reduce(         long long * object, int count, MPI_Op op, int root, MPI_Comm comm);
void reduce(    unsigned short * object, int count, MPI_Op op, int root, MPI_Comm comm);
void reduce(      unsigned int * object, int count, MPI_Op op, int root, MPI_Comm comm);
void reduce(     unsigned long * object, int count, MPI_Op op, int root, MPI_Comm comm);
void reduce(unsigned long long * object, int count, MPI_Op op, int root, MPI_Comm comm);

void reduce(vector<              char> & object, MPI_Op op, int root, MPI_Comm comm);
void reduce(vector<       signed char> & object, MPI_Op op, int root, MPI_Comm comm);
void reduce(vector<     unsigned char> & object, MPI_Op op, int root, MPI_Comm comm);
void reduce(vector<             short> & object, MPI_Op op, int root, MPI_Comm comm);
void reduce(vector<               int> & object, MPI_Op op, int root, MPI_Comm comm);
void reduce(vector<              long> & object, MPI_Op op, int root, MPI_Comm comm);
void reduce(vector<         long long> & object, MPI_Op op, int root, MPI_Comm comm);
void reduce(vector<    unsigned short> & object, MPI_Op op, int root, MPI_Comm comm);
void reduce(vector<      unsigned int> & object, MPI_Op op, int root, MPI_Comm comm);
void reduce(vector<     unsigned long> & object, MPI_Op op, int root, MPI_Comm comm);
void reduce(vector<unsigned long long> & object, MPI_Op op, int root, MPI_Comm comm);

/* all reduce */
void allreduce(              char & object, MPI_Op op, MPI_Comm comm);
void allreduce(       signed char & object, MPI_Op op, MPI_Comm comm);
void allreduce(     unsigned char & object, MPI_Op op, MPI_Comm comm);
void allreduce(             short & object, MPI_Op op, MPI_Comm comm);
void allreduce(               int & object, MPI_Op op, MPI_Comm comm);
void allreduce(              long & object, MPI_Op op, MPI_Comm comm);
void allreduce(         long long & object, MPI_Op op, MPI_Comm comm);
void allreduce(    unsigned short & object, MPI_Op op, MPI_Comm comm);
void allreduce(      unsigned int & object, MPI_Op op, MPI_Comm comm);
void allreduce(     unsigned long & object, MPI_Op op, MPI_Comm comm);
void allreduce(unsigned long long & object, MPI_Op op, MPI_Comm comm);

void allreduce(              char * object, int count, MPI_Op op, MPI_Comm comm);
void allreduce(       signed char * object, int count, MPI_Op op, MPI_Comm comm);
void allreduce(     unsigned char * object, int count, MPI_Op op, MPI_Comm comm);
void allreduce(             short * object, int count, MPI_Op op, MPI_Comm comm);
void allreduce(               int * object, int count, MPI_Op op, MPI_Comm comm);
void allreduce(              long * object, int count, MPI_Op op, MPI_Comm comm);
void allreduce(         long long * object, int count, MPI_Op op, MPI_Comm comm);
void allreduce(    unsigned short * object, int count, MPI_Op op, MPI_Comm comm);
void allreduce(      unsigned int * object, int count, MPI_Op op, MPI_Comm comm);
void allreduce(     unsigned long * object, int count, MPI_Op op, MPI_Comm comm);
void allreduce(unsigned long long * object, int count, MPI_Op op, MPI_Comm comm);

void allreduce(vector<              char> & object, MPI_Op op, MPI_Comm comm);
void allreduce(vector<       signed char> & object, MPI_Op op, MPI_Comm comm);
void allreduce(vector<     unsigned char> & object, MPI_Op op, MPI_Comm comm);
void allreduce(vector<             short> & object, MPI_Op op, MPI_Comm comm);
void allreduce(vector<               int> & object, MPI_Op op, MPI_Comm comm);
void allreduce(vector<              long> & object, MPI_Op op, MPI_Comm comm);
void allreduce(vector<         long long> & object, MPI_Op op, MPI_Comm comm);
void allreduce(vector<    unsigned short> & object, MPI_Op op, MPI_Comm comm);
void allreduce(vector<      unsigned int> & object, MPI_Op op, MPI_Comm comm);
void allreduce(vector<     unsigned long> & object, MPI_Op op, MPI_Comm comm);
void allreduce(vector<unsigned long long> & object, MPI_Op op, MPI_Comm comm);

/* all to all */
void alltoall(vector<              char> & sendrecv, MPI_Comm comm);
void alltoall(vector<       signed char> & sendrecv, MPI_Comm comm);
void alltoall(vector<     unsigned char> & sendrecv, MPI_Comm comm);
void alltoall(vector<             short> & sendrecv, MPI_Comm comm);
void alltoall(vector<               int> & sendrecv, MPI_Comm comm);
void alltoall(vector<              long> & sendrecv, MPI_Comm comm);
void alltoall(vector<         long long> & sendrecv, MPI_Comm comm);
void alltoall(vector<    unsigned short> & sendrecv, MPI_Comm comm);
void alltoall(vector<      unsigned int> & sendrecv, MPI_Comm comm);
void alltoall(vector<     unsigned long> & sendrecv, MPI_Comm comm);
void alltoall(vector<unsigned long long> & sendrecv, MPI_Comm comm);

void alltoall(vector<              char> & send, vector<              char> & recv, MPI_Comm comm);
void alltoall(vector<       signed char> & send, vector<       signed char> & recv, MPI_Comm comm);
void alltoall(vector<     unsigned char> & send, vector<     unsigned char> & recv, MPI_Comm comm);
void alltoall(vector<             short> & send, vector<             short> & recv, MPI_Comm comm);
void alltoall(vector<               int> & send, vector<               int> & recv, MPI_Comm comm);
void alltoall(vector<              long> & send, vector<              long> & recv, MPI_Comm comm);
void alltoall(vector<         long long> & send, vector<         long long> & recv, MPI_Comm comm);
void alltoall(vector<    unsigned short> & send, vector<    unsigned short> & recv, MPI_Comm comm);
void alltoall(vector<      unsigned int> & send, vector<      unsigned int> & recv, MPI_Comm comm);
void alltoall(vector<     unsigned long> & send, vector<     unsigned long> & recv, MPI_Comm comm);
void alltoall(vector<unsigned long long> & send, vector<unsigned long long> & recv, MPI_Comm comm);

/* gather */
vector<              char> gather(              char & object, int root, MPI_Comm comm);
vector<       signed char> gather(       signed char & object, int root, MPI_Comm comm);
vector<     unsigned char> gather(     unsigned char & object, int root, MPI_Comm comm);
vector<             short> gather(             short & object, int root, MPI_Comm comm);
vector<               int> gather(               int & object, int root, MPI_Comm comm);
vector<              long> gather(              long & object, int root, MPI_Comm comm);
vector<         long long> gather(         long long & object, int root, MPI_Comm comm);
vector<    unsigned short> gather(    unsigned short & object, int root, MPI_Comm comm);
vector<      unsigned int> gather(      unsigned int & object, int root, MPI_Comm comm);
vector<     unsigned long> gather(     unsigned long & object, int root, MPI_Comm comm);
vector<unsigned long long> gather(unsigned long long & object, int root, MPI_Comm comm);

vector<              char> gather(              char * object, int count, int root, MPI_Comm comm);
vector<       signed char> gather(       signed char * object, int count, int root, MPI_Comm comm);
vector<     unsigned char> gather(     unsigned char * object, int count, int root, MPI_Comm comm);
vector<             short> gather(             short * object, int count, int root, MPI_Comm comm);
vector<               int> gather(               int * object, int count, int root, MPI_Comm comm);
vector<              long> gather(              long * object, int count, int root, MPI_Comm comm);
vector<         long long> gather(         long long * object, int count, int root, MPI_Comm comm);
vector<    unsigned short> gather(    unsigned short * object, int count, int root, MPI_Comm comm);
vector<      unsigned int> gather(      unsigned int * object, int count, int root, MPI_Comm comm);
vector<     unsigned long> gather(     unsigned long * object, int count, int root, MPI_Comm comm);
vector<unsigned long long> gather(unsigned long long * object, int count, int root, MPI_Comm comm);

vector<              char> gather(vector<              char> & object, int root, MPI_Comm comm);
vector<       signed char> gather(vector<       signed char> & object, int root, MPI_Comm comm);
vector<     unsigned char> gather(vector<     unsigned char> & object, int root, MPI_Comm comm);
vector<             short> gather(vector<             short> & object, int root, MPI_Comm comm);
vector<               int> gather(vector<               int> & object, int root, MPI_Comm comm);
vector<              long> gather(vector<              long> & object, int root, MPI_Comm comm);
vector<         long long> gather(vector<         long long> & object, int root, MPI_Comm comm);
vector<    unsigned short> gather(vector<    unsigned short> & object, int root, MPI_Comm comm);
vector<      unsigned int> gather(vector<      unsigned int> & object, int root, MPI_Comm comm);
vector<     unsigned long> gather(vector<     unsigned long> & object, int root, MPI_Comm comm);
vector<unsigned long long> gather(vector<unsigned long long> & object, int root, MPI_Comm comm);

/* synchronous printing */
void print_sync(const string &label, const               char & object, MPI_Comm comm);
void print_sync(const string &label, const        signed char & object, MPI_Comm comm);
void print_sync(const string &label, const      unsigned char & object, MPI_Comm comm);
void print_sync(const string &label, const              short & object, MPI_Comm comm);
void print_sync(const string &label, const                int & object, MPI_Comm comm);
void print_sync(const string &label, const               long & object, MPI_Comm comm);
void print_sync(const string &label, const          long long & object, MPI_Comm comm);
void print_sync(const string &label, const     unsigned short & object, MPI_Comm comm);
void print_sync(const string &label, const       unsigned int & object, MPI_Comm comm);
void print_sync(const string &label, const      unsigned long & object, MPI_Comm comm);
void print_sync(const string &label, const unsigned long long & object, MPI_Comm comm);

void print_sync(const string &label, const vector<              char> & object, MPI_Comm comm);
void print_sync(const string &label, const vector<       signed char> & object, MPI_Comm comm);
void print_sync(const string &label, const vector<     unsigned char> & object, MPI_Comm comm);
void print_sync(const string &label, const vector<             short> & object, MPI_Comm comm);
void print_sync(const string &label, const vector<               int> & object, MPI_Comm comm);
void print_sync(const string &label, const vector<              long> & object, MPI_Comm comm);
void print_sync(const string &label, const vector<         long long> & object, MPI_Comm comm);
void print_sync(const string &label, const vector<    unsigned short> & object, MPI_Comm comm);
void print_sync(const string &label, const vector<      unsigned int> & object, MPI_Comm comm);
void print_sync(const string &label, const vector<     unsigned long> & object, MPI_Comm comm);
void print_sync(const string &label, const vector<unsigned long long> & object, MPI_Comm comm);

void print_sync(const string &label, const               char * object, int size, MPI_Comm comm);
void print_sync(const string &label, const        signed char * object, int size, MPI_Comm comm);
void print_sync(const string &label, const      unsigned char * object, int size, MPI_Comm comm);
void print_sync(const string &label, const              short * object, int size, MPI_Comm comm);
void print_sync(const string &label, const                int * object, int size, MPI_Comm comm);
void print_sync(const string &label, const               long * object, int size, MPI_Comm comm);
void print_sync(const string &label, const          long long * object, int size, MPI_Comm comm);
void print_sync(const string &label, const     unsigned short * object, int size, MPI_Comm comm);
void print_sync(const string &label, const       unsigned int * object, int size, MPI_Comm comm);
void print_sync(const string &label, const      unsigned long * object, int size, MPI_Comm comm);
void print_sync(const string &label, const unsigned long long * object, int size, MPI_Comm comm);

void print_zero(const string &label, MPI_Comm comm);

void print_zero(const string &label, const               char & object, MPI_Comm comm);
void print_zero(const string &label, const        signed char & object, MPI_Comm comm);
void print_zero(const string &label, const      unsigned char & object, MPI_Comm comm);
void print_zero(const string &label, const              short & object, MPI_Comm comm);
void print_zero(const string &label, const                int & object, MPI_Comm comm);
void print_zero(const string &label, const               long & object, MPI_Comm comm);
void print_zero(const string &label, const          long long & object, MPI_Comm comm);
void print_zero(const string &label, const     unsigned short & object, MPI_Comm comm);
void print_zero(const string &label, const       unsigned int & object, MPI_Comm comm);
void print_zero(const string &label, const      unsigned long & object, MPI_Comm comm);
void print_zero(const string &label, const unsigned long long & object, MPI_Comm comm);

void print_zero(const string &label, const vector<              char> & object, MPI_Comm comm);
void print_zero(const string &label, const vector<       signed char> & object, MPI_Comm comm);
void print_zero(const string &label, const vector<     unsigned char> & object, MPI_Comm comm);
void print_zero(const string &label, const vector<             short> & object, MPI_Comm comm);
void print_zero(const string &label, const vector<               int> & object, MPI_Comm comm);
void print_zero(const string &label, const vector<              long> & object, MPI_Comm comm);
void print_zero(const string &label, const vector<         long long> & object, MPI_Comm comm);
void print_zero(const string &label, const vector<    unsigned short> & object, MPI_Comm comm);
void print_zero(const string &label, const vector<      unsigned int> & object, MPI_Comm comm);
void print_zero(const string &label, const vector<     unsigned long> & object, MPI_Comm comm);
void print_zero(const string &label, const vector<unsigned long long> & object, MPI_Comm comm);

void print_zero(const string &label, const               char * object, int size, MPI_Comm comm);
void print_zero(const string &label, const        signed char * object, int size, MPI_Comm comm);
void print_zero(const string &label, const      unsigned char * object, int size, MPI_Comm comm);
void print_zero(const string &label, const              short * object, int size, MPI_Comm comm);
void print_zero(const string &label, const                int * object, int size, MPI_Comm comm);
void print_zero(const string &label, const               long * object, int size, MPI_Comm comm);
void print_zero(const string &label, const          long long * object, int size, MPI_Comm comm);
void print_zero(const string &label, const     unsigned short * object, int size, MPI_Comm comm);
void print_zero(const string &label, const       unsigned int * object, int size, MPI_Comm comm);
void print_zero(const string &label, const      unsigned long * object, int size, MPI_Comm comm);
void print_zero(const string &label, const unsigned long long * object, int size, MPI_Comm comm);

/* file reading */
MPI_Offset get_file_size(const string &file_name, MPI_Comm comm);
void read_file(const string &file_name, char *&file_buffer, MPI_Offset &file_size, MPI_Comm comm);
void read_file_bcast(const string &file_name, char *&file_buffer, MPI_Offset &file_size, MPI_Comm comm);
void read_file_mpiio(const string &file_name, char *&file_buffer, MPI_Offset &file_size, MPI_Comm comm);

} /* namespace mpix */

#endif /* _MPIX_H_ */
