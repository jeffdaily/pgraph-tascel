/**
 * @author jeff.daily@pnnl.gov
 *
 * A first attempt at creating a forest of tries for a given protein dataset.
 *
 * We generate the tree index as a forest of emerging at a specified
 * depth <= W, so that the individual subtrees can be independently traversed
 * in parallel to generate promising pairs.
 *
 * This involves binning the suffix prefixes for each input sequence.
 */
#include "config.h"

#include <mpi.h>

#include <stdint.h>

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <strstream>
#include <utility>
#include <vector>

#include "mpix.h"
#include "timer.h"

#define ARG_LEN_MAX 1024
#define W 4

using namespace std;
typedef uint64_t idx_t;
#define MAX(A,B) (A) > (B) ? (A) : (B)

int rank;
int nprocs;
vector<string> sequences;
#define VEC_BINS 0
#if VEC_BINS
vector<vector<pair<idx_t,idx_t> > > bins(456976); // 26**4
#else
typedef struct {
    int bin_id;
    int seq_id;
    int offset;
} test_t;
bool cmp_test_t(const test_t &i, const test_t &j) {
    return i.bin_id < j.bin_id;
}
vector<test_t> bins; // 26**4
#endif
static int pow_26[] = { 1, 26, 676, 17576 };

static inline int prefix_hash(const string &prefix)
{
    int hash = 0;
    int ordinal=0;
    string::const_reverse_iterator rit;
    for (rit=prefix.rbegin(); rit<prefix.rend(); ++rit) {
        hash += pow_26[ordinal] * ((*rit)-'A');
        ++ordinal;
    }
    return hash;
}

static inline int prefix_hash(const string &seq, const idx_t &j)
{
    int hash = 0;
    int ordinal=W-1;
    for (idx_t i=j; i<j+W; ++i) {
        hash += pow_26[ordinal] * (seq[i]-'A');
        --ordinal;
    }
    return hash;
}


int main(int argc, char **argv)
{
    MPI_Comm comm = MPI_COMM_NULL;
    vector<string> all_argv;
    long file_size = -1;
    char *file_buffer = NULL;
    MPI_File fh;
    MPI_Status status;
    long seg_count = 0;
    idx_t max_seq_len = 0;
    idx_t tot_seq_len = 0;

    timer_init();

    /* initialize MPI */
    MPI_CHECK(MPI_Init(&argc, &argv));
    MPI_CHECK(MPI_Comm_rank(MPI_COMM_WORLD, &rank));
    MPI_CHECK(MPI_Comm_size(MPI_COMM_WORLD, &nprocs));
    MPI_CHECK(MPI_Comm_dup(MPI_COMM_WORLD, &comm));

    /* MPI standard does not guarantee all procs receive argc and argv */
    if (0 == rank) {
        MPI_CHECK(MPI_Bcast(&argc, 1, MPI_INT, 0, comm));
        for (int i=0; i<argc; ++i) {
            int length = strlen(argv[i])+1;
            MPI_CHECK(MPI_Bcast(&length, 1, MPI_INT, 0, comm));
            MPI_CHECK(MPI_Bcast(argv[i], length, MPI_CHAR, 0, comm));
            all_argv.push_back(argv[i]);
        }
    } else {
        int all_argc;
        MPI_CHECK(MPI_Bcast(&all_argc, 1, MPI_INT, 0, comm));
        for (int i=0; i<all_argc; ++i) {
            int length;
            char buffer[ARG_LEN_MAX];
            MPI_CHECK(MPI_Bcast(&length, 1, MPI_INT, 0, comm));
            MPI_CHECK(MPI_Bcast(buffer, length, MPI_CHAR, 0, comm)); all_argv.push_back(buffer);
        }
    }

#if DEBUG
    /* print the command line arguments */
    for (int i=0; i<nprocs; ++i) {
        if (i == rank) {
            int j;
            for (j=0; j<all_argv.size(); ++j) {
                printf("[%d] argv[%d]=%s\n", rank, j, all_argv[j].c_str());
            }
        }
        MPI_Barrier(comm);
    }
#endif

    /* sanity check that we got the correct number of arguments */
    if (all_argv.size() <= 1 || all_argv.size() >= 3) {
        if (0 == rank) {
            if (all_argv.size() <= 1) {
                printf("missing input file\n");
            }
            else if (all_argv.size() >= 3) {
                printf("too many arguments\n");
            }
        }
        MPI_Comm_free(&comm);
        MPI_Finalize();
        return 1;
    }

    /* process 0 open the file locally to determine its size */
    if (0 == rank) {
        FILE *file = fopen(all_argv[1].c_str(), "r");
        if (NULL == file) {
            printf("unable to open file on process 0\n");
            MPI_Abort(comm, 1);
        }
        if (0 != fseek(file, 0, SEEK_END)) {
            printf("unable to seek to end of file on process 0\n");
            MPI_Abort(comm, 1);
        }
        file_size = ftell(file);
        if (-1 == file_size) {
            printf("unable to get size of file on process 0\n");
            MPI_Abort(comm, 1);
        }
        if (0 != fclose(file)) {
            printf("unable to get close file on process 0\n");
            MPI_Abort(comm, 1);
        }
    }

    /* the file_size is broadcast to all */
    MPI_CHECK(MPI_Bcast(&file_size, 1, MPI_LONG, 0, comm));
    printf("[%d] file_size=%ld\n", rank, file_size);
    /* allocate a buffer for the file, of the entire size */
    /* TODO: this is not memory efficient since we allocate a buffer to read
     * the entire input file and then parse the buffer into a vector of
     * strings, essentially doubling the memory requirement */
    file_buffer = new char[file_size];

    /* all procs read the entire file */
    MPI_CHECK(MPI_File_open(comm, const_cast<char*>(all_argv[1].c_str()),
                MPI_MODE_RDONLY|MPI_MODE_UNIQUE_OPEN,
                MPI_INFO_NULL, &fh));
    MPI_CHECK(MPI_File_read_all(fh, file_buffer, file_size, MPI_CHAR, &status));
    MPI_CHECK(MPI_File_close(&fh));

    /* each process counts how many '>' characters are in the file_buffer */
    for (int i=0; i<file_size; ++i) {
        if (file_buffer[i] == '>') {
            ++seg_count;
        }
    }
#if 1
    /* print the seg_count on each process */
    for (int i=0; i<nprocs; ++i) {
        if (i == rank) {
            printf("[%d] seg_count=%ld\n", rank, seg_count);
        }
        MPI_Barrier(comm);
    }
#endif

    /* TODO declare these at the top */
    /* each process indexes the file_buffer */
    istrstream input_stream(const_cast<const char*>(file_buffer), file_size);
    string line;
    string sequence;
    sequences.reserve(seg_count);
    while (getline(input_stream, line)) {
        if (line[0] == '>') {
            if (!sequence.empty()) {
                sequences.push_back(sequence);
                max_seq_len = MAX(max_seq_len, sequence.size());
                tot_seq_len += sequence.size();
            }
            sequence.clear();
            continue;
        }
        sequence += line;
    }
    /* add the last sequence in the file since we wouldn't encounter another
     * '>' character but rather an EOF */
    if (!sequence.empty()) {
        sequences.push_back(sequence);
        max_seq_len = MAX(max_seq_len, sequence.size());
        tot_seq_len += sequence.size();
    }
    sequence.clear();
#if 1
    /* print the seg_count on each process */
    for (int i=0; i<nprocs; ++i) {
        if (i == rank) {
            printf("[%d] sequences.size()=%ld sequences.capacity()=%ld\n",
                    rank, long(sequences.size()), long(sequences.capacity()));
        }
        MPI_Barrier(comm);
    }
    /* print the max_seq_len on each process */
    for (int i=0; i<nprocs; ++i) {
        if (i == rank) {
            printf("[%d] max_seq_len=%lu\n",
                    rank, (unsigned long)(max_seq_len));
            printf("[%d] tot_seq_len=%lu\n",
                    rank, (unsigned long)(tot_seq_len));
        }
        MPI_Barrier(comm);
    }
#endif
#if DEBUG
    /* print the first sequence on each process */
    for (int i=0; i<nprocs; ++i) {
        if (i == rank) {
            printf("[%d] sequences[0]=%s\n", rank, sequences[0].c_str());
        }
        MPI_Barrier(comm);
    }
    /* print the last sequence on each process */
    for (int i=0; i<nprocs; ++i) {
        if (i == rank) {
            printf("[%d] sequences.back()=%s\n",
                    rank, sequences.back().c_str());
        }
        MPI_Barrier(comm);
    }
#endif
#if !VEC_BINS
    bins.reserve(tot_seq_len);
#endif

    /* iterate over the sequences, binning their suffixes */
    long long t = timer_start();
    long long ti = timer_start();
    idx_t ilimit=sequences.size();
    for (idx_t i=0; i<ilimit; ++i) {
        if (i % 10000 == 0) {
            long long elapsed = timer_end(ti);
            ti = timer_start();
            cout << i << "/" << ilimit << "\t" << elapsed << endl;
        }
        string& seq = sequences[i];
        if (seq.size() < W) {
            continue;
        }
        for (idx_t j=0,jlimit=seq.size()-W; j<jlimit; ++j) {
#if VEC_BINS
            int hash = prefix_hash(seq,j);
            vector<pair<idx_t,idx_t> > &bin = bins[hash];
            bin.push_back(make_pair(i,j));
#else
            test_t t = {prefix_hash(seq,j), i, j};
            bins.push_back(t);
#endif
        }
    }
#if !VEC_BINS
    //sort(bins.begin(), bins.end(), cmp_test_t);
#endif
    t = timer_end(t);
    cout << "timer " << t << endl;
    long long capacity=0;
    long long size=0;
    long long zeros=0;
#if VEC_BINS
    for (idx_t i=0; i<bins.size(); ++i) {
        capacity+=bins[i].capacity();
        size+=bins[i].size();
        if (bins[i].size() == 0) {
            ++zeros;
        }
    }
#else
    capacity=bins.capacity();
    size=bins.size();
#endif
    cout << "capacity=" << capacity << endl;
    cout << "size=" << size << endl;
    cout << "zeros=" << zeros << endl;

    /* clean up */
    delete [] file_buffer;

    MPI_Comm_free(&comm);
    MPI_Finalize();

    return 0;
}
