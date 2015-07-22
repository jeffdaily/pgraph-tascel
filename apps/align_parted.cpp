/**
 * @author jeff.daily@pnnl.gov
 *
 * Alignment of an input dataset using work stealing and combinations of
 * input partitions.
 */
#include "config.h"

/* 3rd party headers */
#include <mpi.h>

#include <tascel.h>
#include <tascel/UniformTaskCollection.h>
#include <tascel/UniformTaskCollectionSplit.h>
#if THREADED
#include "tascelx.hpp"
#endif

#include <parasail.h>
#include <parasail/io.h>

/* C/C++ headers */
#include <algorithm>
#include <cctype>
#include <cfloat>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <set>
#include <stack>
#include <utility>
#include <vector>

/* pgraph contrib headers */
#include "sais.h"

/* pgraph headers */
#include "AlignStats.hpp"
#include "alignment.hpp"
#include "combinations.h"
#include "EdgeResult.hpp"
#include "Bootstrap.hpp"
#include "mpix.hpp"
#include "mpix_types.hpp"
#include "Parameters.hpp"
#include "SuffixArrayStats.hpp"

using namespace ::std;
using namespace ::tascel;
using namespace ::pgraph;

typedef pair<int,int> Pair;

typedef set<Pair> PairSet;

typedef struct {
    size_t id1;
    size_t id2;
    char type;
} task_description_t;

typedef struct {
    int rank;
    int nprocs;
    UniformTaskCollection **utcs;
    AlignStats *stats_align;
    SuffixArrayStats *stats_sa;
    const char *sequences;
    long n_sequences;
    unsigned long *SID;
    vector<unsigned long> *BEG;
    vector<unsigned long> *END;
    char sentinal;
    vector<EdgeResult> *edge_results;
    Parameters *parameters;
    parasail_function_t *aligner;
    const parasail_matrix_t *matrix;
} local_data_t;

struct quad {
    int lcp;
    int lb;
    int rb;
    vector<quad> children;

    quad()
        : lcp(0), lb(0), rb(INT_MAX), children() {}
    quad(int lcp, int lb, int rb)
        : lcp(lcp), lb(lb), rb(rb), children() {}
    quad(int lcp, int lb, int rb, vector<quad> children)
        : lcp(lcp), lb(lb), rb(rb), children(children) {}

    bool empty() { return rb == INT_MAX; }
};

static int inner_main(int argc, char **argv);

static void pair_check(
        UniformTaskCollection *utc,
        unsigned long &count_generated,
        PairSet &pairs,
        const int &i,
        const int &j,
        const int * const restrict SA,
        const unsigned char * const restrict BWT,
        const unsigned long * const restrict SID,
        unsigned long sid_crossover,
        const char &sentinal);

static void process(
        UniformTaskCollection *utc,
        unsigned long &count_generated,
        PairSet &pairs,
        const quad &q,
        const int * const restrict SA,
        const unsigned char * const restrict BWT,
        const unsigned long * const restrict SID,
        unsigned long sid_crossover,
        const char &sentinal,
        const int &cutoff);

static void SA_filter(
        UniformTaskCollection *utc,
        unsigned long *SID,
        const char *T,
        long n,
        char sentinal,
        unsigned long sid_crossover,
        int cutoff,
        SuffixArrayStats &stats_sa);

static int trank(int thd);

static string get_edges_filename(int rank);

static bool length_filter(size_t s1Len, size_t s2Len, size_t cutOff);

static void dual_task(
        UniformTaskCollection * utc,
        void *bigd, int bigd_len,
        void *pldata, int pldata_len,
        vector<void *> data_bufs,
        int thd);

static void alignment_task(
        UniformTaskCollection * /*utc*/,
        void *bigd, int /*bigd_len*/,
        void *pldata, int /*pldata_len*/,
        vector<void *> /*data_bufs*/,
        int thd);

static void sa_task(
        UniformTaskCollection * /*utc*/,
        void *bigd, int /*bigd_len*/,
        void *pldata, int /*pldata_len*/,
        vector<void *> /*data_bufs*/,
        int thd);


int main(int argc, char **argv)
{
    int retval;

    try {
        retval = inner_main(argc, argv);
    }
    catch (const std::bad_alloc &ba) {
        MPI_Abort(MPI_COMM_WORLD, -1);
        ::std::cerr << "bad_alloc: " << ba.what() << ::std::endl;
        retval = -1;
    }

    return retval;
}


static int inner_main(int argc, char **argv)
{
    double time_main = 0.0;
    double time = 0.0;
    int rank = 0;
    int nprocs = 0;
    vector<string> all_argv;
    UniformTaskCollection **utcs = NULL;
    AlignStats *stats_align = NULL;
    SuffixArrayStats *stats_sa = NULL;
    vector<EdgeResult> *edge_results = NULL;
    Parameters *parameters = NULL;
    local_data_t *local_data = NULL;
    char *file_buffer = NULL;
    char *packed_buffer = NULL;
    long packed_size = 0;
    MPI_Offset file_size = 0;
    unsigned long sid = 0;
    unsigned long *SID = NULL;
    vector<unsigned long> BEG;
    vector<unsigned long> END;
    char sentinal = 0;
    int cutoff = 7;

    /* init pgraph, which inits MPI line */
    pgraph::initialize(argc, argv);
    rank = mpix::comm_rank(pgraph::comm);
    nprocs = mpix::comm_size(pgraph::comm);
    time_main = MPI_Wtime();

    /* initialize tascel */
    TascelConfig::initialize(NUM_WORKERS_DEFAULT, pgraph::comm);
#if THREADED
    allocate_threads();
    initialize_server_thread();
#endif

    /* initialize global data */
    utcs = new UniformTaskCollection*[NUM_WORKERS];
    stats_align = new AlignStats[NUM_WORKERS];
    stats_sa = new SuffixArrayStats[NUM_WORKERS];
    edge_results = new vector<EdgeResult>[NUM_WORKERS];
    parameters = new Parameters;
    local_data = new local_data_t;
    local_data->rank = rank;
    local_data->nprocs = nprocs;
    local_data->utcs = utcs;
    local_data->stats_align = stats_align;
    local_data->stats_sa = stats_sa;
    local_data->edge_results = edge_results;
    local_data->parameters = parameters;

    /* MPI standard does not guarantee all procs receive argc and arg */
    all_argv = mpix::bcast(argc, argv, pgraph::comm);

    /* print the command line arguments */
    if (0 == rank) {
        ostringstream header;
        header.fill('-');
        header << left << setw(79) << "--- Command Line ";
        cout << header.str() << endl;
        for (size_t j=0; j<all_argv.size(); ++j) {
            cout << "argv[" << j << "]='" << all_argv[j] << "'" << endl;
        }
    }

    /* sanity check that we got the correct number of arguments */
    if (all_argv.size() <= 1 || all_argv.size() >= 4) {
        if (0 == rank) {
            if (all_argv.size() <= 1) {
                cout << "missing input file" << endl;
            }
            else if (all_argv.size() >= 4) {
                cout << "too many arguments" << endl;
            }
            cout << "usage: align sequence_file <config_file>" << endl;
        }
        TascelConfig::finalize();
        pgraph::finalize();
        return 1;
    }
    else if (all_argv.size() >= 3) {
        parameters->parse(all_argv[2].c_str(), pgraph::comm);
    }
    else if (all_argv.size() >= 2) {
        /* do nothing */
    }

    /* print parameters */
    if (0 == rank) {
        ostringstream header;
        header.fill('-');
        header << left << setw(79) << "-- Parameters ";
        cout << header.str() << endl;
        cout << *parameters << endl;
        cout << endl;
    }
    cutoff = parameters->exact_match_length;

    /* lookup function and matrix */
    local_data->aligner = parasail_lookup_function(parameters->function.c_str());
    local_data->matrix = parasail_matrix_lookup(parameters->matrix.c_str());
    if (NULL == local_data->aligner) {
        cout << "specified function not found" << endl;
        TascelConfig::finalize();
        pgraph::finalize();
        return 1;
    }

    time = MPI_Wtime();
    mpix::read_file(all_argv[1], file_buffer, file_size, pgraph::comm);
    time = MPI_Wtime() - time;
    if (0 == rank) {
        cout << "time sequence db open " << time << endl;
    }

    /* pack, then scan packed buffer to build various indexes */
    time = MPI_Wtime();
    packed_buffer = parasail_pack_buffer(file_buffer, file_size, &packed_size);
    /* done with original file buffer */
    delete [] file_buffer;
    /* allocate sequence ID lookup */
    SID = new unsigned long[packed_size];
    /* determine sentinal and actual end of file (last char) */
    {
        int off = 0;
        while (!isgraph(packed_buffer[packed_size-off])) {
            ++off;
        }
        sentinal = packed_buffer[packed_size-off];
        packed_size = packed_size - off + 1;
    }
    /* scan packed_buffer from left to build sequence ID and end index */
    sid = 0;
    BEG.push_back(0);
    for (unsigned long i=0; i<packed_size; ++i) {
        SID[i] = sid;
        if (packed_buffer[i] == sentinal) {
            END.push_back(i);
            BEG.push_back(i+1);
            ++sid;
        }
    }
    assert(0 != sid);
    BEG.pop_back();
    assert(BEG.size() == END.size());
    time = MPI_Wtime() - time;
    if (0 == rank) {
        cout << "number of sequences: " << sid << endl;;
        cout << "time pack and index db " << time << endl;
    }
    local_data->sequences = packed_buffer;
    local_data->n_sequences = sid;
    local_data->SID = SID;
    local_data->BEG = &BEG;
    local_data->END = &END;
    local_data->sentinal = sentinal;

    /* how many combinations of sequences are there? */
    unsigned long ntasks = binomial_coefficient(sid, 2);
    if (0 == rank) {
        cout << "brute force "
            << sid
            << " choose 2 has "
            << ntasks
            << " combinations"
            << endl;
    }

    while (sid % parameters->sa_block_size == 1) {
        if (0 == rank) {
            cout << "sa_block_size parameter left a remainder of 1; increasing by 1" << endl;
            parameters->sa_block_size += 1;
        }
    }

    unsigned long global_num_workers = nprocs*NUM_WORKERS;
    unsigned long tasks_per_worker = ntasks / global_num_workers;
    unsigned long max_tasks_per_worker = ntasks / global_num_workers;
    max_tasks_per_worker += ntasks % global_num_workers;
    max_tasks_per_worker *= 100;
    if (0 == rank) {
        printf("global_num_workers=%lu\n", global_num_workers);
        printf("tasks_per_worker=%lu\n", tasks_per_worker);
        printf("max_tasks_per_worker=%lu\n", max_tasks_per_worker);
    }
    MPI_Barrier(pgraph::comm);

    /* populate tasks */
    vector<double> populate_time(NUM_WORKERS, 0.0);
    vector<double> populate_count(NUM_WORKERS, 0.0);
    {
        TslFuncRegTbl *frt = new TslFuncRegTbl;
        TslFunc tf = frt->add(dual_task);
        TaskCollProps props;
        props.functions(tf, frt)
            .taskSize(sizeof(task_description_t))
            .maxTasks(parameters->memory_worker / sizeof(task_description_t))
            .localData(local_data, sizeof(void*));
        int parts = (sid + parameters->sa_block_size - 1) / parameters->sa_block_size;
        int tiles = parts*(parts-1)/2 + parts;
        if (0 == rank) {
            printf("sequences split into %d parts for a total of %d tiles\n",
                    parts, tiles);
        }
        int crank = 0;
        int wsize = nprocs*NUM_WORKERS;
        task_description_t desc;
        map<int,int> wrank_to_worker;
        for (int worker=0; worker<NUM_WORKERS; ++worker) {
            populate_time[worker] = MPI_Wtime();
            utcs[worker] = new UniformTaskCollectionSplit(props, worker);
            wrank_to_worker[trank(worker)] = worker;
        }
        for (int i=0; i<parts; ++i) {
            for (int j=i; j<parts; ++j) {
                desc.id1 = i;
                desc.id2 = j;
                desc.type = 0;
                if (wrank_to_worker.count(crank)) {
                    int worker = wrank_to_worker[crank];
                    utcs[worker]->addTask(&desc, sizeof(desc));
                    populate_count[worker] += 1;
                }
                ++crank;
                crank = crank % wsize;
            }
        }
        for (int worker=0; worker<NUM_WORKERS; ++worker) {
            populate_time[worker] = MPI_Wtime() - populate_time[worker];
        }
    }

    if (parameters->print_stats) {
        if (0 == rank) {
            cout << "utcs capacity = " << utcs[0]->capacity() << endl;
        }
    }

    if (parameters->print_stats) {
        populate_time = mpix::gather(populate_time, 0, pgraph::comm);
        populate_count = mpix::gather(populate_count, 0, pgraph::comm);
        if (0 == rank) {
            Stats times;
            Stats counts;
            int p = cout.precision();
            char f = cout.fill();
            cout << setfill('-');
            cout << left;
            cout << setw(79) << "--- Task Population Statistics" << endl;
            cout << setfill(f);
            cout << setprecision(2);
            cout << right;
            cout << setw(5) << "pid";
            cout << setw(15) << "populate_time";
            cout << setw(15) << "populate_count";
            cout << endl;
            for (size_t i=0; i<populate_time.size(); ++i) {
                times.push_back(populate_time[i]);
                counts.push_back(populate_count[i]);
                cout << setw(5) << std::right << i
                    << setw(15) << fixed << populate_time[i]
                    << setw(15) << fixed << populate_count[i]
                    << endl;
            }
            cout << string(79, '=') << endl;
            cout << "       " << Stats::header() << endl;
            cout << "Times  " << times << endl;
            cout << "Counts " << counts << endl;
            cout.precision(p);
        }
    }
    populate_time.clear();
    populate_count.clear();

    amBarrier();

    {
#if THREADED
        initialize_worker_threads(utcs);
#endif
        double mytimer = MPI_Wtime();
        utcs[0]->process();
        mytimer = MPI_Wtime() - mytimer;
        if (0 == trank(0)) {
            cout << "mytimer=" << mytimer << endl;
        }
#if THREADED
        finalize_worker_threads();
        finalize_server_thread();
        deallocate_threads();
#endif
    }

    amBarrier();
    MPI_Barrier(pgraph::comm);

    if (parameters->print_stats) {
        vector<AlignStats> rstats = mpix::gather(stats_align, NUM_WORKERS, 0, pgraph::comm);
        /* synchronously print alignment stats all from process 0 */
        if (0 == rank) {
            Stats edge_counts;
            Stats align_counts;
            Stats align_skipped;
            Stats time_align;
            Stats time_kcomb;
            Stats time_total;
            Stats work;
            Stats work_skipped;
            ostringstream header;
            int p = cout.precision();

            header.fill('-');
            header << left << setw(79) << "--- Align Stats ";
            Stats::width(11);
            cout << header.str() << endl;
            cout << setprecision(2);
            cout << right << setw(5) << "pid" << AlignStats::header() << endl;
            for(int i=0; i<nprocs*NUM_WORKERS; i++) {
                edge_counts.push_back(rstats[i].edge_counts);
                align_counts.push_back(rstats[i].align_counts);
                align_skipped.push_back(rstats[i].align_skipped);
                time_align.push_back(rstats[i].time_align.sum());
                time_kcomb.push_back(rstats[i].time_kcomb);
                time_total.push_back(rstats[i].time_total);
                work.push_back(rstats[i].work);
                work_skipped.push_back(rstats[i].work_skipped);
                cout << right << setw(5) << i << rstats[i] << endl;
            }
            Stats::width(21);
            cout << setprecision(1);
            cout << string(79, '=') << endl;
            cout << "           " << Stats::header() << endl;
            cout << "      Edges" << edge_counts << endl;
            cout << " Alignments" << align_counts << endl;
            cout << "  AlignSkip" << align_skipped << endl;
            cout << "     TAlign" << time_align << endl;
            cout << "     TTotal" << time_total << endl;
            cout << "       Work" << work << endl;
            cout << "WorkSkipped" << work_skipped << endl;
            cout << string(79, '-') << endl;
            cout.precision(p);
        }
    }

    if (parameters->print_stats) {
        vector<SuffixArrayStats> rstats = mpix::gather(stats_sa, NUM_WORKERS, 0, pgraph::comm);
        /* synchronously print tree stats all from process 0 */
        if (0 == rank) {
            SuffixArrayStats cumulative;
            Stats arrays_per_worker;
            Stats times_per_worker;
            ostringstream header;
            int p = cout.precision();

            header.fill('-');
            header << left << setw(79) << "--- Suffix Array Stats ";
            cout << header.str() << endl;
            cout << setprecision(2);
            Stats::width(13);
            for(int i=0; i<nprocs*NUM_WORKERS; i++) {
                cout << right << setw(5) << i;
                cout << right << setw(14) << "name";
                cout << Stats::header() << endl;
                cumulative += rstats[i];
                arrays_per_worker.push_back(rstats[i].arrays);
                times_per_worker.push_back(
                        rstats[i].time_build.sum()+
                        rstats[i].time_process.sum());
                cout << rstats[i] << endl;
            }
            cout << string(79, '=') << endl;
            cout << right << setw(5) << "TOTAL";
            cout << right << setw(14) << "name";
            cout << Stats::header() << endl;
            cout << cumulative;
            cout << right << setw(19) << "ArraysPerWorker" << arrays_per_worker << endl;
            cout << right << setw(19) << "TimesPerWorker" << times_per_worker << endl;
            cout << "first array" << setw(25) << cumulative.time_first << endl;
            cout << " last array" << setw(25) << cumulative.time_last << endl;
            cout << "       diff" << setw(25) << cumulative.time_last - cumulative.time_first << endl;
            cout.precision(p);
            cout << string(79, '-') << endl;
        }
    }

    if (parameters->print_stats) {
        vector<StealingStats> stt(NUM_WORKERS);
        vector<StealingStats> rstt(NUM_WORKERS*nprocs);

        for(int i=0; i<NUM_WORKERS; i++) {
            stt[i] = utcs[i]->getStats();
        }
        rstt = mpix::gather(stt, 0, pgraph::comm);
        /* synchronously print stealing stats all from process 0 */
        if (0 == rank) {
            StealingStats tstt;
            cout<<" pid"<<rstt[0].formatString()<<endl;
            for(int i=0; i<nprocs*NUM_WORKERS; i++) {
                tstt += rstt[i];
                cout<<std::setw(4)<<std::right<<i<<rstt[i]<<endl;
            }
            cout << string(79, '=') << endl;
            cout<<"TOT "<<tstt<<endl;
        }
    }

    amBarrier();
    MPI_Barrier(pgraph::comm);

    if (parameters->output_to_disk) {
        double time = MPI_Wtime();
        ofstream edge_out(get_edges_filename(rank).c_str());
        for (int worker=0; worker<NUM_WORKERS; ++worker) {
            for (size_t i=0,limit=edge_results[worker].size(); i<limit; ++i) {
                edge_out << edge_results[worker][i] << endl;
            }
        }
        edge_out.close();
        time = MPI_Wtime() - time;
        if (0 == rank) {
            ostringstream header;
            header.fill('-');
            header << left << setw(79) << "--- Disk Output Stats ";
            Stats::width(11);
            cout << header.str() << endl;
            cout << Stats::header() << endl;
            vector<double> times = mpix::gather(time, 0, pgraph::comm);
            Stats ts;
            for (int i=0; i<nprocs; ++i ) {
                ts.push_back(times[i]);
            }
            cout << ts << endl;
        }
        else {
            (void)mpix::gather(time, 0, pgraph::comm);
        }
    }

    /* clean up */
    for (int worker=0; worker<NUM_WORKERS; ++worker) {
        delete utcs[worker];
    }
    delete [] utcs;
    delete [] stats_align;
    delete [] edge_results;
    delete parameters;
    free(packed_buffer);
    delete [] SID;
    delete local_data;

    time_main = MPI_Wtime() - time_main;
    if (0 == rank) {
        cout << "time_main " << time_main << " seconds" << endl;
    }
    TascelConfig::finalize();
    pgraph::finalize();

    return 0;
}

static void pair_check(
        UniformTaskCollection *utc,
        unsigned long &count_generated,
        PairSet &pairs,
        const int &i,
        const int &j,
        const int * const restrict SA,
        const unsigned char * const restrict BWT,
        const unsigned long * const restrict SID,
        unsigned long sid_crossover,
        const char &sentinal)
{
    const int &sidi = SID[SA[i]];
    const int &sidj = SID[SA[j]];
    pair<PairSet::iterator, bool> result;
    result.second = false;
    if (BWT[i] != BWT[j] || BWT[i] == sentinal) {
        if (0 == sid_crossover) {
            if (sidi != sidj) {
                ++count_generated;
                if (sidi < sidj) {
                    result = pairs.insert(make_pair(sidi,sidj));
                }
                else {
                    result = pairs.insert(make_pair(sidj,sidi));
                }
            }
        }
        else {
            if ((sidi < sid_crossover && sidj >= sid_crossover)
                    || (sidj < sid_crossover && sidi >= sid_crossover)) {
                ++count_generated;
                if (sidi < sidj) {
                    result = pairs.insert(make_pair(sidi,sidj));
                }
                else {
                    result = pairs.insert(make_pair(sidj,sidi));
                }
            }
        }
    }
    if (result.second) {
        if (result.first->first >= result.first->second) {
            printf("uh oh %d %d\n", result.first->first, result.first->second);
        }
        assert(result.first->first < result.first->second);
        task_description_t desc;
        desc.id1 = result.first->first;
        desc.id2 = result.first->second;
        desc.type = 1;
        utc->addTask(&desc, sizeof(desc));
    }
}

/* try to reduce number of duplicate pairs generated */
/* we observe that l-intervals (i.e. internal nodes) always have at
 * least two children, but these children could be singleton
 * l-intervals, e.g., [i..j]=[1..1], in addition to l-intervals with
 * non-singleton ranges/quads. For each l-interval, we take the cross
 * product of its child l-intervals. Naively, we could take the cross
 * product of the entire lb/rb range of the l-interval, but this
 * generates too many duplicate pairs. Instead, the complexity should be
 * bounded by the number of exact matches...
 */
static void process(
        UniformTaskCollection *utc,
        unsigned long &count_generated,
        PairSet &pairs,
        const quad &q,
        const int * const restrict SA,
        const unsigned char * const restrict BWT,
        const unsigned long * const restrict SID,
        unsigned long sid_crossover,
        const char &sentinal,
        const int &cutoff)
{
    const int n_children = q.children.size();
    int child_index = 0;

    if (q.lcp < cutoff) return;

    if (n_children) {
        for (int i=q.lb; i<=q.rb; ++i) {
            int j = i+1;
            if (child_index < n_children) {
                if (i >= q.children[child_index].lb) {
                    j = q.children[child_index].rb+1;
                    if (i >= q.children[child_index].rb) {
                        ++child_index;
                    }
                }
            }
            for (/*nope*/; j<=q.rb; ++j) {
                pair_check(utc, count_generated, pairs, i, j, SA, BWT, SID, sid_crossover, sentinal);
            }
        }
    }
    else {
        for (int i=q.lb; i<=q.rb; ++i) {
            for (int j=i+1; j<=q.rb; ++j) {
                pair_check(utc, count_generated, pairs, i, j, SA, BWT, SID, sid_crossover, sentinal);
            }
        }
    }
}

/* SID and T arrays are of size n. */
static void SA_filter(
        UniformTaskCollection *utc,
        unsigned long *SID,
        const char *T,
        long n,
        char sentinal,
        unsigned long sid_crossover,
        int cutoff,
        SuffixArrayStats &stats_sa)
{
    int rank = mpix::comm_rank(pgraph::comm);
    int nprocs = mpix::comm_size(pgraph::comm);
    unsigned long *SID_local = NULL;
    int *SA = NULL;
    int *LCP = NULL;
    unsigned char *BWT = NULL;
    double start = 0;
    double finish = 0;
    int i = 0;
    unsigned long sid = 0;
    vector<unsigned long> BEG;
    vector<unsigned long> END;
    double time = 0.0;
    unsigned long count_possible = 0;
    unsigned long count_generated = 0;
    unsigned long sid_crossover_local = 0;
    PairSet pairs;
    double time_build = 0.0;
    double time_process = 0.0;

    if (stats_sa.time_first == 0.0) {
        stats_sa.time_first = MPI_Wtime();
    }

    time_build = MPI_Wtime();
    time = MPI_Wtime();
    SID_local = new unsigned long[n+1];
    /* scan T from left to build sequence ID and end index */
    sid = 0;
    BEG.push_back(0);
    for (unsigned long i=0; i<n; ++i) {
        SID_local[i] = sid;
        if (T[i] == sentinal) {
            END.push_back(i);
            BEG.push_back(i+1);
            if (0 == sid_crossover_local && SID[i] == sid_crossover) {
                sid_crossover_local = sid;
            }
            ++sid;
        }
    }
    assert(0 != sid);
    BEG.pop_back();
    assert(BEG.size() == END.size());
    time = MPI_Wtime() - time;
#if 0
    if (0 == rank) {
        cout << "number of local sequences: " << sid << endl;;
        cout << "time index local db " << time << endl;
    }
#endif

    /* Allocate memory for enhanced SA. */
    SA = new int[n+1]; /* +1 for LCP */
    LCP = new int[n+1]; /* +1 for lcp tree */
    BWT = new unsigned char[n+1];
    if((SA == NULL) || (LCP == NULL) || (BWT == NULL))
    {
        cerr << "Cannot allocate ESA memory." << endl;
        exit(EXIT_FAILURE);
    }

    /* Construct the suffix and LCP arrays.
     * The following sais routine is from Fischer, with bugs fixed. */
    start = parasail_time();
    if(sais((const unsigned char *)T, SA, LCP, (int)n) != 0) {
        cerr << "Cannot allocate memory." << endl;
        exit(EXIT_FAILURE);
    }
    finish = parasail_time();
    //cout << "induced SA time: " << finish-start << " seconds" << endl;

    /* construct naive BWT: */
    start = parasail_time();
    for (i = 0; i < n; ++i) {
        BWT[i] = (SA[i] > 0) ? T[SA[i]-1] : sentinal;
    }
    finish = parasail_time();
    //cout << "naive BWT time: " << finish-start << " seconds" << endl;

    /* "fix" the LCP array to clamp LCP's that are too long */
    start = parasail_time();
    for (i = 0; i < n; ++i) {
        int len = END[SID_local[SA[i]]] - SA[i]; /* don't include sentinal */
        if (LCP[i] > len) LCP[i] = len;
    }
    finish = parasail_time();
    //cout << "clamp LCP time: " << finish-start << " seconds" << endl;

    stats_sa.time_build.push_back(MPI_Wtime() - time_build);
    time_process = MPI_Wtime();

    /* The GSA we create will put all sentinals either at the beginning
     * or end of the SA. We don't want to count all of the terminals,
     * nor do we want to process them in our bottom-up traversal. */
    /* do the sentinals appear at the beginning or end of SA? */
    int bup_start = 1;
    int bup_stop = n;
    if (T[SA[0]] == sentinal) {
        /* sentinals at beginning */
        bup_start = sid+1;
        bup_stop = n;
    }
    else if (T[SA[n-1]] == sentinal) {
        /* sentinals at end */
        bup_start = 1;
        bup_stop = n-sid;
    }
    else {
        cerr << "sentinals not found at beginning or end of SA" << endl;
        exit(EXIT_FAILURE);
    }

    /* DFS of enhanced SA, from Abouelhoda et al */
    start = parasail_time();
    count_generated = 0;
    LCP[n] = 0; /* doesn't really exist, but for the root */
    {
        stack<quad> the_stack;
        quad last_interval;
        the_stack.push(quad());
        for (i = bup_start; i <= bup_stop; ++i) {
            int lb = i - 1;
            while (LCP[i] < the_stack.top().lcp) {
                the_stack.top().rb = i - 1;
                last_interval = the_stack.top();
                the_stack.pop();
                process(utc, count_generated, pairs, last_interval, SA, BWT, SID, sid_crossover, sentinal, cutoff);
                lb = last_interval.lb;
                if (LCP[i] <= the_stack.top().lcp) {
                    last_interval.children.clear();
                    the_stack.top().children.push_back(last_interval);
                    last_interval = quad();
                }
            }
            if (LCP[i] > the_stack.top().lcp) {
                if (!last_interval.empty()) {
                    last_interval.children.clear();
                    the_stack.push(quad(LCP[i],lb,INT_MAX,vector<quad>(1, last_interval)));
                    last_interval = quad();
                }
                else {
                    the_stack.push(quad(LCP[i],lb,INT_MAX));
                }
            }
        }
        the_stack.top().rb = bup_stop - 1;
        process(utc, count_generated, pairs, the_stack.top(), SA, BWT, SID, sid_crossover, sentinal, cutoff);
    }
    finish = parasail_time();
    stats_sa.time_process.push_back(MPI_Wtime() - time_process);
    if (0 == sid_crossover) {
        count_possible = ((unsigned long)sid)*((unsigned long)sid-1)/2;
    }
    else {
        count_possible = (sid-sid_crossover_local)*sid_crossover_local;
    }
#if 0
    cout << "ESA time: " << finish-start << " seconds" << endl;
    cout << "possible pairs: " << count_possible << endl;
    cout << "generated pairs: " << count_generated << endl;
    cout << "unique pairs: " << pairs.size() << endl;
#endif

    stats_sa.arrays++;
    stats_sa.suffixes.push_back(n);
    stats_sa.pairs.push_back(count_generated);
    stats_sa.time_last = MPI_Wtime();

    /* Deallocate memory. */
    delete [] SID_local;
    delete [] SA;
    delete [] LCP;
    delete [] BWT;
}

static int trank(int thd)
{
    return (theTwoSided().getProcRank().toInt() * NUM_WORKERS) + thd;
}

static string get_edges_filename(int rank)
{
    ostringstream str;
    str << "edges." << rank << ".txt";
    return str.str();
}

static void dual_task(
        UniformTaskCollection * utc,
        void *bigd, int bigd_len,
        void *pldata, int pldata_len,
        vector<void *> data_bufs,
        int thd)
{
    task_description_t *desc = (task_description_t*)bigd; 
    if (desc->type == 0) {
        sa_task(utc, bigd, bigd_len, pldata, pldata_len, data_bufs, thd);
    }
    else if (desc->type == 1) {
        alignment_task(utc, bigd, bigd_len, pldata, pldata_len, data_bufs, thd);
    }
    else {
        cerr << "bad task type (" << desc->type << ")" << endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
}

static void alignment_task(
        UniformTaskCollection * /*utc*/,
        void *bigd, int /*bigd_len*/,
        void *pldata, int /*pldata_len*/,
        vector<void *> /*data_bufs*/,
        int thd)
{
    task_description_t *desc = (task_description_t*)bigd; 
    local_data_t *local_data = (local_data_t*)pldata; 
    bool is_edge_answer = false;
    double t = 0;
    double tt = 0;
    int sscore;
    size_t max_len;

    AlignStats *stats = local_data->stats_align;
    const char *sequences = local_data->sequences;
    vector<unsigned long> &BEG = *(local_data->BEG);
    vector<unsigned long> &END = *(local_data->END);
    vector<EdgeResult> *edge_results = local_data->edge_results;
    parasail_function_t *aligner = local_data->aligner;
    const parasail_matrix_t *matrix = local_data->matrix;
    Parameters *parameters = local_data->parameters;

    int open = parameters->open;
    int gap = parameters->gap;
    int AOL = parameters->AOL;
    int SIM = parameters->SIM;
    int OS = parameters->OS;
    bool do_alignment = parameters->perform_alignments;

    tt = MPI_Wtime();

    size_t i = desc->id1;
    size_t j = desc->id2;
    int i_beg = BEG[i];
    int i_end = END[i];
    int s1Len = i_end-i_beg;
    int j_beg = BEG[j];
    int j_end = END[j];
    int s2Len = j_end-j_beg;

    const char * c1 = &sequences[i_beg];
    const char * c2 = &sequences[j_beg];

    if (parameters->use_length_filter) {
        do_alignment &= length_filter(s1Len, s2Len, AOL*SIM/100);
    }

    if (do_alignment)
    {
        stats[thd].work += s1Len * s2Len;
        ++stats[thd].align_counts;
        t = MPI_Wtime();
        parasail_result_t *result;
        result = aligner(
                c1, s1Len, c2, s2Len, -open, -gap, matrix);
        is_edge_answer = is_edge(
                result, c1, s1Len, c2, s2Len, AOL, SIM, OS, sscore, max_len, matrix);

        if (parameters->output_to_disk
                && (is_edge_answer || parameters->output_all))
        {
            edge_results[thd].push_back(
                    EdgeResult(
                        i, j,
                        1.0*result->length/max_len,
                        1.0*result->matches/result->length,
                        1.0*result->score/sscore,
                        is_edge_answer)
                    );
        }
        parasail_result_free(result);
        if (is_edge_answer) {
            ++stats[thd].edge_counts;
        }
        t = MPI_Wtime() - t;
        stats[thd].time_align.push_back(t);
    }
    else {
        stats[thd].work_skipped += s1Len * s2Len;
        stats[thd].align_skipped += 1;
    }

    tt = MPI_Wtime() - tt;
    stats[thd].time_total += tt;
}

static void sa_task(
        UniformTaskCollection * utc,
        void *bigd, int /*bigd_len*/,
        void *pldata, int /*pldata_len*/,
        vector<void *> /*data_bufs*/,
        int thd)
{
    task_description_t *desc = (task_description_t*)bigd; 
    local_data_t *local_data = (local_data_t*)pldata; 
    int block_size = local_data->parameters->sa_block_size;
    int id1 = desc->id1;
    int id2 = desc->id2;
    int id1_beg = id1 * block_size;
    int id2_beg = id2 * block_size;
    int id1_end = id1_beg + block_size - 1;
    int id2_end = id2_beg + block_size - 1;
    int beg1 = (*local_data->BEG)[id1_beg];
    int beg2 = (*local_data->BEG)[id2_beg];
    int end1 = 0;
    int end2 = 0;
    assert(id1 <= id2);
    if (id1_end > local_data->n_sequences) {
        id1_end = local_data->n_sequences - 1;
    }
    end1 = (*local_data->END)[id1_end];
    if (id2_end > local_data->n_sequences) {
        id2_end = local_data->n_sequences - 1;
    }
    end2 = (*local_data->END)[id2_end];
    int len1 = end1 - beg1 + 1;
    int len2 = end2 - beg2 + 1;
    char *sequences = NULL;
    unsigned long *SID = NULL;
    int cutoff = local_data->parameters->exact_match_length;
    int sid_crossover = 0;
    SuffixArrayStats *stats_sa = local_data->stats_sa;

    if (id1 == id2) {
        sequences = new char[len1+1];
        SID = new unsigned long[len1];
        copy(&local_data->sequences[beg1],
             &local_data->sequences[beg1+len1],
             &sequences[0]);
        copy(&local_data->SID[beg1],
             &local_data->SID[beg1+len1],
             &SID[0]);
        sequences[len1] = '\0';
        SA_filter(utc, SID, sequences, len1, local_data->sentinal, sid_crossover, cutoff, stats_sa[thd]);
    }
    else {
        sequences = new char[len1+len2+1];
        SID = new unsigned long[len1+len2];
        copy(&local_data->sequences[beg1],
             &local_data->sequences[beg1+len1],
             &sequences[0]);
        copy(&local_data->SID[beg1],
             &local_data->SID[beg1+len1],
             &SID[0]);
        copy(&local_data->sequences[beg2],
             &local_data->sequences[beg2+len2],
             &sequences[len1]);
        copy(&local_data->SID[beg2],
             &local_data->SID[beg2+len2],
             &SID[len1]);
        sid_crossover = id2_beg;
        sequences[len1+len2] = '\0';
        SA_filter(utc, SID, sequences, len1+len2, local_data->sentinal, sid_crossover, cutoff, stats_sa[thd]);
    }

    delete [] sequences;
    delete [] SID;
}

static bool length_filter(size_t s1Len, size_t s2Len, size_t cutOff)
{
    bool result = true;

    if (s1Len <= s2Len) {
        if (100 * s1Len < cutOff * s2Len) {
            result = false;
        }
    }
    else {
        if (100 * s2Len < cutOff * s1Len) {
            result = false;
        }
    }

    return result;
}

