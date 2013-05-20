/**
 * @author jeff.daily@pnnl.gov
 *
 * Serially perform sequence alignment. No work stealing. No threading.
 */
#include "config.h"

#include <sys/time.h>

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <strstream>
#include <vector>

#include "AlignStats.hpp"
#include "combinations.h"
#include "dynamic.h"

using namespace std;

#if USE_CLOCK_GETTTIME
double timespec_diff(double time_past=0.0)
{
    struct timespec current;
    clock_gettime(CLOCK_MONOTONIC, &current);
    return (current.tv_sec * 1000000000.0) + current.tv_nsec - time_past;
}
#else
double timespec_diff(double time_past=0.0)
{
    struct timeval current;
    gettimeofday(&current, NULL);
    //return (current.tv_sec * 1000000.0) + current.tv_usec - time_past;
    return current.tv_sec + (current.tv_usec/1000000.0) - time_past;
}
#endif

#define DEBUG 0

#ifndef OUTPUT_EDGES
#define OUTPUT_EDGES 1
#endif
//#define SEP "\t"
#define SEP ","
#define CACHE_RESULTS 1
#define ALL_RESULTS 1

#if OUTPUT_EDGES
class EdgeResult {
    public:
        unsigned long id1;
        unsigned long id2;
        double a;
        double b;
        double c;
#if ALL_RESULTS
        bool is_edge;
#endif

        EdgeResult(
                unsigned long id1, unsigned long id2,
                double a, double b, double c
#if ALL_RESULTS
                ,bool is_edge
#endif
                )
            : id1(id1)
            , id2(id2)
            , a(a)
            , b(b)
            , c(c)
#if ALL_RESULTS
            , is_edge(is_edge)
#endif
        {}

        friend ostream& operator << (ostream &os, const EdgeResult &edge) {
            os << edge.id1
                << SEP << edge.id2
                << SEP << edge.a
                << SEP << edge.b
                << SEP << edge.c;
#if ALL_RESULTS
            if (edge.is_edge) {
                os << SEP << "Edge";
            }
            else {
                os << SEP << "NotAnEdge";
            }
#endif
            return os;
        }
};
#endif


cell_t **tbl;
int **del;
int **ins;
vector<string> sequences;
#if OUTPUT_EDGES
#if CACHE_RESULTS
vector<EdgeResult> edge_results;
#else
ofstream edges;
#endif
#endif
AlignStats stats;


static void alignment_task(unsigned long seq_id[2])
{
    cell_t result;
    is_edge_param_t param;
    int is_edge_answer = 0;
    double t = 0;
    unsigned long i = 0;
    int sscore = 0;
    int maxLen = 0;

    result.score = 0;
    result.ndig = 0;
    result.alen = 0;
    param.AOL = 8;
    param.SIM = 4;
    param.OS = 3;

    t = timespec_diff();
    affine_gap_align(
            sequences[seq_id[0]].c_str(), sequences[seq_id[0]].size(),
            sequences[seq_id[1]].c_str(), sequences[seq_id[1]].size(),
            &result, tbl, del, ins);
    is_edge_answer = is_edge(result,
            sequences[seq_id[0]].c_str(), sequences[seq_id[0]].size(),
            sequences[seq_id[1]].c_str(), sequences[seq_id[1]].size(),
            param, &sscore, &maxLen);
    ++stats.align_counts;

    if (is_edge_answer || ALL_RESULTS)
    {
#if DEBUG
        cout << ": aligned " << seq_id[0] << " " << seq_id[1]
             << ": (score,ndig,alen)=("
             << result.score << ","
             << result.ndig << ","
             << result.alen << ")"
             << ": edge? " << is_edge_answer << endl;
#endif
#if OUTPUT_EDGES
#if CACHE_RESULTS
        assert(maxLen);
        assert(result.alen);
        edge_results.push_back(EdgeResult(
                    seq_id[0], seq_id[1],
#if 1
                    1.0*result.alen/maxLen,
                    1.0*result.ndig/result.alen,
                    1.0*result.score/sscore
#else
                    result.alen,
                    result.ndig,
                    result.score
#endif
#if ALL_RESULTS
                    ,is_edge_answer
#endif
                    ));
#else
        edges << seq_id[0] << SEP << seq_id[1] << SEP
            << 1.0*result.alen/maxLen << SEP
            << 1.0*result.ndig/result.alen << SEP
            << 1.0*result.score/sscore
#if ALL_RESULTS
            << SEP << is_edge_answer
#endif
            << endl;
#endif
#endif
        if (is_edge_answer) {
            ++stats.edge_counts;
        }
    }
    t = timespec_diff(t);
    stats.align_times_tot += t;
    stats.calc_min(t);
    stats.calc_max(t);
}


int main(int argc, char **argv)
{
    long file_size = -1;
    char *file_buffer = NULL;
    long seg_count = 0;
    size_t max_seq_len = 0;
    unsigned long nCk;
    unsigned long nalignments;
    bool use_combinations = true;
    int fixed_size = 0;

    /* initialize dynamic code */
    init_map(SIGMA);

#if DEBUG
    /* print the command line arguments */
    for (int j=0; j<argc; ++j) {
        printf("argv[%d]=%s\n", j, argv[j]);
    }
#endif

    /* sanity check that we got the correct number of arguments */
    if (argc <= 1 || argc >= 4) {
        if (argc <= 1) {
            printf("missing input file\n");
        }
        else if (argc >= 4) {
            printf("too many arguments\n");
        }
        printf("usage: brute_force_serial sequence_file fixed\n");
        return 1;
    }

    if (argc >= 3) {
        fixed_size = atoi(argv[2]);
        use_combinations = false;
    }

    /* process 0 open the file locally to determine its size */
    FILE *file = fopen(argv[1], "r");
    if (NULL == file) {
        printf("unable to open file\n");
        return 1;
    }
    if (0 != fseek(file, 0, SEEK_END)) {
        printf("unable to seek to end of file\n");
        return 1;
    }
    file_size = ftell(file);
    if (-1 == file_size) {
        printf("unable to get size of file\n");
        return 1;
    }

    /* allocate a buffer for the file, of the entire size */
    /* TODO: this is not memory efficient since we allocate a buffer to read
     * the entire input file and then parse the buffer into a vector of
     * strings, essentially doubling the memory requirement */
    file_buffer = new char[file_size];

    rewind(file);

    if (0 == fread(file_buffer, file_size, 1, file)) {
        printf("unable to read entire file\n");
        return 1;
    }

    if (0 != fclose(file)) {
        printf("unable to get close file\n");
        return 1;
    }

    /* count how many '>' characters are in the file_buffer */
    for (int i=0; i<file_size; ++i) {
        if (file_buffer[i] == '>') {
            ++seg_count;
        }
    }
#if DEBUG
    /* print the seg_count */
    printf("seg_count=%ld\n", seg_count);
#endif

    /* TODO declare these at the top */
    /* index the file_buffer */
    istrstream input_stream(const_cast<const char*>(file_buffer), file_size);
    string line;
    string sequence;
    sequences.reserve(seg_count);
    while (getline(input_stream, line)) {
        if (line[0] == '>') {
            if (!sequence.empty()) {
                sequences.push_back(sequence);
                max_seq_len = max(max_seq_len, sequence.size());
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
        max_seq_len = max(max_seq_len, sequence.size());
    }
    sequence.clear();
    delete [] file_buffer;
#if DEBUG
    /* print the seg_count */
    printf("sequences.size()=%ld sequences.capacity()=%ld\n",
            long(sequences.size()), long(sequences.capacity()));
    /* print the max_seq_len */
    printf("max_seq_len=%ld\n", long(max_seq_len));
#endif
#if DEBUG
    /* print the first sequence */
    printf("sequences[0]=%s\n", sequences[0].c_str());
    /* print the last sequence */
    printf("sequences.back()=%s\n", sequences.back().c_str());
#endif
    /* how many combinations of sequences are there? */
    nCk = binomial_coefficient(sequences.size(), 2);

    if (use_combinations) {
        nalignments = nCk;
        printf("brute force %lu C 2 has %lu combinations\n",
                sequences.size(), nCk);
    }
    else {
        printf("aligning the first %lu against all %lu for %lu alignments\n",
                fixed_size, sequences.size(), sequences.size()*fixed_size);
        nalignments = sequences.size() * fixed_size;
    }

    /* some more dynamic initialization */
    assert(NROW == 2);
    tbl = alloc_tbl(NROW, max_seq_len);
    del = alloc_int(NROW, max_seq_len);
    ins = alloc_int(NROW, max_seq_len);

#if OUTPUT_EDGES
#if CACHE_RESULTS
    edge_results.reserve(nalignments);
#else
    edges.open("edges.txt");
#endif
#endif
    /* perform the alignments */
    double mytimer = timespec_diff();
    unsigned long seq_id[2];
    if (use_combinations) {
        k_combination(0, 2, seq_id);
        for (unsigned long i=0; i<nalignments; ++i) {
            alignment_task(seq_id);
            next_combination(2, seq_id);
        }
    }
    else {
        printf("fixed_size=%d\n",fixed_size);
        for (unsigned long i=0; i<sequences.size(); ++i) {
            for (unsigned long j=0; j<fixed_size; ++j) {
                seq_id[0] = i;
                seq_id[1] = j;
                alignment_task(seq_id);
            }
        }
    }
    mytimer = timespec_diff(mytimer);
    cout << "mytimer=" << mytimer << endl;
    cout << stats.getHeader() << endl;      
    cout << stats << endl;
    if (use_combinations) {
#if OUTPUT_EDGES
#if CACHE_RESULTS
        ofstream edge_out("edges.txt");
        for (size_t i=0,limit=edge_results.size(); i<limit; ++i) {
            edge_out << edge_results[i] << endl;
        }
        edge_out.close();
#else
        edges.close();
#endif
#endif
    }
    else {
        unsigned long e=0;
        ofstream out("edges_fixed.csv");
        out << "id";
        for (unsigned long i=0; i<fixed_size; i++) {
            out << "," << i << "_alen/maxlen"
                << "," << i << "_ndig/alen"
                << "," << i << "_score/sscore";
        }
        out << ",class" << endl;
        for (unsigned long i=0; i<sequences.size(); ++i) {
            out << i;
            bool is_homologous = false;
            for (unsigned long j=0; j<fixed_size; ++j) {
                EdgeResult r = edge_results[e++];
                out << "," << r.a
                    << "," << r.b
                    << "," << r.c;
                is_homologous = is_homologous || r.is_edge;
            }
            if (is_homologous) {
                out << "," << "Homologous" << endl;
            }
            else {
                out << "," << "Not" << endl;
            }
        }
        out.close();
    }
    /* clean up */
    free_tbl(tbl, NROW);
    free_int(del, NROW);
    free_int(ins, NROW);

    return 0;
}
