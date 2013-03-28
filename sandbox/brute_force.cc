/**
 * @author jeff.daily@pnnl.gov
 *
 * A first attempt at a brute force alignment of an input dataset using work
 * stealing. Each MPI task reads the input file.
 */
#include "config.h"

#include <mpi.h>
#include <tascel.h>
#include <tascel/UniformTaskCollection.h>
#include <tascel/UniformTaskCollSplitHybrid.h>

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <strstream>
#include <vector>

#include "AlignStats.h"
#include "combinations.h"
#include "dynamic.h"
#include "mpix.h"

using namespace std;
using namespace tascel;

#ifndef MULTIPLE_PAIRS_PER_TASK
#define MULTIPLE_PAIRS_PER_TASK 0
#endif

/* MULTIPLE_PAIRS_PER_TASK with ROUND_ROBIN_SEQUENCE_PAIRS not supported */
#if MULTIPLE_PAIRS_PER_TASK
#undef ROUND_ROBIN_SEQUENCE_PAIRS
#define ROUND_ROBIN_SEQUENCE_PAIRS 0
#endif

#ifndef ROUND_ROBIN_SEQUENCE_PAIRS
#define ROUND_ROBIN_SEQUENCE_PAIRS 0
#endif

#ifndef SORT_TASKS_LOCALLY
#define SORT_TASKS_LOCALLY 0
#endif

#ifndef SORT_SEQUENCES
#define SORT_SEQUENCES 0
#endif

#ifndef DUMP_COSTS
#define DUMP_COSTS 0
#endif

#ifndef OUTPUT_EDGES
#define OUTPUT_EDGES 0
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
                << SEP << edge.c
#if ALL_RESULTS
                << SEP << edge.is_edge
#endif
                ;
            return os;
        }
};
#endif

int rank = 0;
int nprocs = 0;
cell_t ***tbl = 0;
int ***del = 0;
int ***ins = 0;
UniformTaskCollSplitHybrid** utcs = 0;
AlignStats *stats = 0;
static pthread_t *threadHandles = 0;
static unsigned *threadRanks = 0;
vector<string> sequences;
#if DUMP_COSTS
ofstream *out = 0;
#endif
#if OUTPUT_EDGES
#if CACHE_RESULTS
vector<EdgeResult> *edge_results = 0;
#else
ofstream *edges = 0;
#endif
#endif
// Synchronization for worker threads
pthread_barrier_t workersStart, workersEnd;
// Synchronization for server thread
pthread_barrier_t serverStart, serverEnd;
volatile bool serverEnabled = true;
#if MULTIPLE_PAIRS_PER_TASK
unsigned long combinations_per_task=0;
#endif


bool longest_string_first(const string &i, const string &j) {
    return i.size() > j.size();
}


void *serverThd(void *args)
{
#if defined(SET_AFFINITY)
    cpu_set_t cpuset;
    pthread_t thread;
    CPU_ZERO(&cpuset);
    thread = pthread_self();

    int rank = theTwoSided().getProcRank().toInt();

    CPU_SET(NUM_WORKERS, &cpuset);

    //int ret = pthread_setaffinity_np(thread, sizeof(cpu_set_t), &cpuset);
    int ret = sched_setaffinity(0, sizeof(cpu_set_t), &cpuset);
#endif

    // When enabled execute any active messages that arrive
    while (1) {
        pthread_barrier_wait(&serverStart);
        while (serverEnabled) {
            AmListenObjCodelet<NullMutex>* lcodelet;
            if((lcodelet=theAm().amListeners[0]->progress()) != NULL) {
                lcodelet->execute();
            }
            Codelet* codelet;
            if ((codelet = serverDispatcher.progress()) != NULL) {
                codelet->execute();
                // Assume that codelet was an AmRequest and needs to be freed
                delete reinterpret_cast<AmRequest*>(codelet);
            }
        }
        pthread_barrier_wait(&serverEnd);
    }
    return NULL;
}

void *workerThd(void *args)
{
    const unsigned threadRank = *(unsigned*)args;

#if defined(SET_AFFINITY)
    cpu_set_t cpuset;
    pthread_t thread;
    CPU_ZERO(&cpuset);
    thread = pthread_self();

    int rank = theTwoSided().getProcRank().toInt();

    CPU_SET(threadRank, &cpuset);

    int ret = pthread_setaffinity_np(thread, sizeof(cpu_set_t), &cpuset);
#endif

    while (1) {
        pthread_barrier_wait(&workersStart);
        utcs[threadRank]->process(threadRank);
        pthread_barrier_wait(&workersEnd);
    }

    return NULL;
}

void amBarrier()
{
    int epoch = pgrp->signalBarrier();
    while(!pgrp->testBarrier(epoch)) {
        AmListenObjCodelet<NullMutex>* codelet;
        if((codelet=theAm().amListeners[0]->progress()) != NULL) {
            codelet->execute();
        }
    }
}

void amBarrierThd()
{
    int epoch = pgrp->signalBarrier();
    while(!pgrp->testBarrier(epoch)) { }
}

static int trank(int thd)
{
    return (theTwoSided().getProcRank().toInt() * NUM_WORKERS) + thd;
}

#if DUMP_COSTS
static string get_filename(int thd)
{
    ostringstream str;
    str << "cost." << trank(thd) << ".txt";
    return str.str();
}
#endif


#if OUTPUT_EDGES
#if CACHE_RESULTS
static string get_edges_filename(int rank)
{
    ostringstream str;
    str << "edges." << rank << ".txt";
    return str.str();
}
#else
static string get_edges_filename(int thd)
{
    ostringstream str;
    str << "edges." << trank(thd) << ".txt";
    return str.str();
}
#endif
#endif


#if MULTIPLE_PAIRS_PER_TASK
typedef struct {
    unsigned long id;
} task_description;
#else
typedef struct {
    unsigned long id1;
    unsigned long id2;
} task_description;

#if SORT_TASKS_LOCALLY
class SortedTask {
    public:
        unsigned long _id1;
        unsigned long _id2;
        unsigned long cost;
        SortedTask()
            : _id1(0)
            , _id2(0)
            , cost(0)
        { }
        SortedTask(unsigned long id1, unsigned long id2)
            : _id1(id1)
            , _id2(id2)
            , cost(sequences[id1].size() * sequences[id2].size())
        { }
        void assign(unsigned long id1, unsigned long id2) {
            _id1 = id1;
            _id2 = id2;
            cost = sequences[id1].size() * sequences[id2].size();
        }
        bool operator < (const SortedTask &other) const {
            return cost > other.cost;
        }
};
#endif

#endif

static void alignment_task(
        UniformTaskCollection *utc,
        void *_bigd, int bigd_len,
        void *pldata, int pldata_len,
        vector<void *> data_bufs, int thd) {
    task_description *desc = (task_description*)_bigd;
#if MULTIPLE_PAIRS_PER_TASK
    unsigned long task_id = desc->id;
    unsigned long combo_start_index = task_id * combinations_per_task;
#endif
    unsigned long seq_id[2];
    cell_t result;
    is_edge_param_t param;
    int is_edge_answer = 0;
    double t = 0;
    unsigned long i;
    int sscore;
    int maxLen;

#if MULTIPLE_PAIRS_PER_TASK
    k_combination(combo_start_index, 2, seq_id);
    for (i=0; i<combinations_per_task; ++i)
#else
    seq_id[0] = desc->id1;
    seq_id[1] = desc->id2;
#endif
    {
#if MULTIPLE_PAIRS_PER_TASK
        /* the last task may run out of sequences */
        if (seq_id[0] >= sequences.size() || seq_id[1] >= sequences.size()) {
            break;
        }
#endif
        t = MPI_Wtime();
        affine_gap_align(
                sequences[seq_id[0]].c_str(), sequences[seq_id[0]].size(),
                sequences[seq_id[1]].c_str(), sequences[seq_id[1]].size(),
                &result, tbl[thd], del[thd], ins[thd]);
        param.AOL = 8;
        param.SIM = 4;
        param.OS = 3;
        is_edge_answer = is_edge(result,
                sequences[seq_id[0]].c_str(), sequences[seq_id[0]].size(),
                sequences[seq_id[1]].c_str(), sequences[seq_id[1]].size(),
                param, &sscore, &maxLen);
        ++stats[thd].align_counts;

        if (is_edge_answer || ALL_RESULTS)
        {
#if DEBUG
            cout << trank(thd)
                << ": aligned " << seq_id[0] << " " << seq_id[1]
                << ": (score,ndig,alen)=("
                << result.score << ","
                << result.ndig << ","
                << result.alen << ")"
                << ": edge? " << is_edge_answer << endl;
#endif
#if OUTPUT_EDGES
#if CACHE_RESULTS
            edge_results[thd].push_back(EdgeResult(
                        seq_id[0], seq_id[1],
#if 0
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
            edges[thd] << seq_id[0] << SEP << seq_id[1] << SEP
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
                ++stats[thd].edge_counts;
            }
        }
        t = MPI_Wtime() - t;
#if DUMP_COSTS
        out[thd]
            << seq_id[0] << SEP
            << seq_id[1] << SEP
            << sequences[seq_id[0]].size() << SEP
            << sequences[seq_id[1]].size() << SEP
            << sequences[seq_id[0]].size()*sequences[seq_id[1]].size() << SEP
            << t << endl;
#endif
        stats[thd].align_times_tot += t;
        stats[thd].calc_min(t);
        stats[thd].calc_max(t);
#if MULTIPLE_PAIRS_PER_TASK
        next_combination(2, seq_id);
#endif
    }
}

#if ROUND_ROBIN_SEQUENCE_PAIRS
unsigned long populate_tasks_rr(
        unsigned long ntasks, unsigned long tasks_per_worker, int worker)
{
    task_description desc;
    int nworkers = nprocs * NUM_WORKERS;
    int wrank = trank(worker);
    unsigned long count = 0;
    unsigned long seq_id[2];
    unsigned long nseq = sequences.size();
    unsigned long remainder = ntasks % (nprocs*NUM_WORKERS);

    /* distribute remainder tasks to first 'remainder' workers */
    if (wrank < long(remainder)) {
        ++tasks_per_worker;
    }

    init_combination(2, seq_id); /* seq_id = {0,1} */
    inc_combination2(wrank, seq_id); /* start seq_id at this worker's rank */
#if SORT_TASKS_LOCALLY
    vector<SortedTask> tasks(tasks_per_worker);
    while (count < tasks_per_worker) {
        //assert(seq_id[1] < nseq);
        tasks[count++].assign(seq_id[0], seq_id[1]);
        inc_combination2(nworkers, seq_id);
    }
    sort(tasks.begin(), tasks.end());
    for (unsigned long i=0; i<count; ++i) {
        desc.id1 = tasks[i]._id1;
        desc.id2 = tasks[i]._id2;
#if DEBUG
        cout << wrank << " added " << desc.id1 << "," << desc.id2 << " cost=" << tasks[i].cost << endl;
#endif
        utcs[worker]->addTask(&desc, sizeof(desc));
    }
#else
    while (count < tasks_per_worker) {
        //assert(seq_id[1] < nseq);
        desc.id1 = seq_id[0];
        desc.id2 = seq_id[1];
#if DEBUG
        cout << wrank << " added " << seq_id[0] << "," << seq_id[1] << endl;
#endif
        utcs[worker]->addTask(&desc, sizeof(desc));
        count++;
        inc_combination2(nworkers, seq_id);
    }
#endif

    return count;
}

#else

unsigned long populate_tasks(
        unsigned long ntasks, unsigned long tasks_per_worker, int worker)
{
    task_description desc;
    int wrank = trank(worker);
    unsigned long count = 0;
    unsigned long i;
    unsigned long lower_limit = wrank*tasks_per_worker;
    unsigned long upper_limit = lower_limit + tasks_per_worker;
    unsigned long remainder = ntasks % (nprocs*NUM_WORKERS);

#if MULTIPLE_PAIRS_PER_TASK
#else
    unsigned long seq_id[2];
    k_combination(lower_limit, 2, seq_id);
#endif
    for (i=lower_limit; i<upper_limit; ++i) {
#if MULTIPLE_PAIRS_PER_TASK
        desc.id = i;
#else
        desc.id1 = seq_id[0];
        desc.id2 = seq_id[1];
#if DEBUG
        cout << wrank << " added " << seq_id[0] << "," << seq_id[1] << endl;
#endif
        next_combination(2, seq_id);
#endif
        utcs[worker]->addTask(&desc, sizeof(desc));
        count++;
    }
    /* if I'm the last worker, add the remainder of the tasks */
    if (wrank == nprocs*NUM_WORKERS-1) {
        for (/*ignore*/; i<upper_limit+remainder; ++i) {
            count++;
#if MULTIPLE_PAIRS_PER_TASK
            desc.id = i;
#else
            desc.id1 = seq_id[0];
            desc.id2 = seq_id[1];
            next_combination(2, seq_id);
#endif
            utcs[worker]->addTask(&desc, sizeof(desc));
        }
    }

    return count;
}

#endif

int main(int argc, char **argv)
{
    MPI_Comm comm = MPI_COMM_NULL;
    int provided;
    vector<string> all_argv;
    long file_size = -1;
    char *file_buffer = NULL;
    MPI_File fh;
    MPI_Status status;
    long seg_count = 0;
    size_t max_seq_len = 0;
    unsigned long nCk;

    /* initialize MPI */
#if defined(THREADED)
    MPI_CHECK(MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided));
    assert(provided == MPI_THREAD_MULTIPLE);
#else
    MPI_CHECK(MPI_Init(&argc, &argv));
#endif
    MPI_CHECK(MPI_Comm_dup(MPI_COMM_WORLD, &comm));
    MPI_CHECK(MPI_Comm_rank(comm, &rank));
    MPI_CHECK(MPI_Comm_size(comm, &nprocs));

    /* initialize tascel */
    TascelConfig::initialize(NUM_WORKERS_DEFAULT, comm);
    tbl = new cell_t**[NUM_WORKERS];
    del = new int**[NUM_WORKERS];
    ins = new int **[NUM_WORKERS];
    utcs = new UniformTaskCollSplitHybrid*[NUM_WORKERS];
    stats = new AlignStats[NUM_WORKERS];
    threadHandles = new pthread_t[NUM_WORKERS + NUM_SERVERS];
    threadRanks = new unsigned[NUM_WORKERS + NUM_SERVERS];
#if DUMP_COSTS
    out = new ofstream[NUM_WORKERS];
#endif
#if OUTPUT_EDGES
#if CACHE_RESULTS
    edge_results = new vector<EdgeResult>[NUM_WORKERS];
#else
    edges = new ofstream[NUM_WORKERS];
#endif
#endif
    for (int worker=0; worker<NUM_WORKERS; ++worker) {
        threadRanks[worker] = worker;
    }

    /* initialize dynamic code */
    init_map(SIGMA);

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
            MPI_CHECK(MPI_Bcast(buffer, length, MPI_CHAR, 0, comm));
            all_argv.push_back(buffer);
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
    if (all_argv.size() <= 1 || all_argv.size() >= 4) {
        if (0 == rank) {
            if (all_argv.size() <= 1) {
                printf("missing input file\n");
            }
            else if (all_argv.size() >= 4) {
                printf("too many arguments\n");
            }
            printf("usage: brute_force sequence_file [combinations_per_task]\n");
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
    if (0 == trank(0)) {
        printf("file_size=%ld\n", file_size);
    }
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
#if DEBUG
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
#if SORT_SEQUENCES
    double *sort_times = new double[nprocs];
    double  sort_time = MPI_Wtime();
    /* sort the sequences, largest first */
    sort(sequences.begin(), sequences.end(), longest_string_first);
    sort_time = MPI_Wtime() - sort_time;
    MPI_CHECK(MPI_Gather(&sort_time, 1, MPI_DOUBLE, sort_times, 1, MPI_DOUBLE, 0, comm));
    if (0 == rank) {
        double tally = 0;
        cout << "rank sort_time" << endl;
        for(int i=0; i<nprocs; i++) {
            tally += sort_times[i];
            cout << std::setw(4) << std::right << i
                << setw(14) << fixed << sort_times[i] << endl;
        }
        cout << "==============================================" << endl;
        cout << setw(4) << right << "T" << setw(14) << fixed << tally << endl;
    }
    delete [] sort_times;
#endif
#if DEBUG
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
            printf("[%d] max_seq_len=%ld\n", rank, long(max_seq_len));
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
    /* how many combinations of sequences are there? */
    nCk = binomial_coefficient(sequences.size(), 2);
    if (0 == trank(0)) {
        printf("brute force %lu C 2 has %lu combinations\n",
                sequences.size(), nCk);
    }

#if MULTIPLE_PAIRS_PER_TASK
    /* set the combinations_per_task variable */
    if (all_argv.size() >= 3) {
        combinations_per_task = strtoll(all_argv[2].c_str(), NULL, 10);
    }
    else {
        combinations_per_task = sequences.size();
    }
    double selectivity = 1.0;
    if (all_argv.size() == 4) {
      selectivity  = fabs(min(1.0,atof(all_argv[3].c_str())));
    }
    unsigned long nalignments = (long)(0.5+selectivity*nCk);
    unsigned long ntasks = nalignments / combinations_per_task;
    if (nalignments % combinations_per_task != 0) {
        ntasks += 1;
    }
    unsigned long global_num_workers = nprocs*NUM_WORKERS;
    unsigned long tasks_per_worker = ntasks / global_num_workers;
    unsigned long max_tasks_per_worker = ntasks / global_num_workers;
    max_tasks_per_worker += ntasks % global_num_workers;
#else
    double selectivity = 1.0;
    if (all_argv.size() == 3) {
      selectivity  = fabs(min(1.0,atof(all_argv[2].c_str())));
    }
    unsigned long nalignments = (long)(0.5+selectivity*nCk);
    unsigned long ntasks = nalignments;
    unsigned long global_num_workers = nprocs*NUM_WORKERS;
    unsigned long tasks_per_worker = ntasks / global_num_workers;
    unsigned long max_tasks_per_worker = ntasks / global_num_workers;
    max_tasks_per_worker += ntasks % global_num_workers;
#endif
    max_tasks_per_worker *= 10;
    if (0 == trank(0)) {
#if MULTIPLE_PAIRS_PER_TASK
        printf("combinations_per_task=%lu\n", combinations_per_task);
#endif
        printf("selectivity=%lf\n", selectivity);
        printf("nalignments=%lu\n", nalignments);
        printf("ntasks=%lu\n", ntasks);
        printf("global_num_workers=%lu\n", global_num_workers);
        printf("tasks_per_worker=%lu\n", tasks_per_worker);
        printf("max_tasks_per_worker=%lu\n", max_tasks_per_worker);
    }
    MPI_Barrier(comm);

    /* some more dynamic initialization */
    assert(NROW == 2);
    for (int worker=0; worker<NUM_WORKERS; ++worker) {
        tbl[worker] = alloc_tbl(NROW, max_seq_len);
        del[worker] = alloc_int(NROW, max_seq_len);
        ins[worker] = alloc_int(NROW, max_seq_len);
    }

    /* the tascel part */
    double populate_times[NUM_WORKERS];
    unsigned long count[NUM_WORKERS];
    unsigned long global_count = 0;
    for (int worker=0; worker<NUM_WORKERS; ++worker)
    {
#if DUMP_COSTS
        out[worker].open(get_filename(worker).c_str());
#endif
#if OUTPUT_EDGES
#if CACHE_RESULTS
        edge_results[worker].reserve(tasks_per_worker);
#else
        edges[worker].open(get_edges_filename(worker).c_str());
#endif
#endif
        UniformTaskCollSplitHybrid*& utc = utcs[worker];
        TslFuncRegTbl *frt = new TslFuncRegTbl();
        TslFunc tf = frt->add(alignment_task);
        TaskCollProps props;
        props.functions(tf, frt)
            .taskSize(sizeof(task_description))
            .maxTasks(max_tasks_per_worker);
        utc = new UniformTaskCollSplitHybrid(props, worker);

        /* add some tasks */
        populate_times[worker] = MPI_Wtime();
#if ROUND_ROBIN_SEQUENCE_PAIRS
        count[worker] = populate_tasks_rr(ntasks, tasks_per_worker, worker);
#else
        count[worker] = populate_tasks(ntasks, tasks_per_worker, worker);
#endif
        populate_times[worker] = MPI_Wtime() - populate_times[worker];
    }
#if DEBUG
    double *g_populate_times = new double[nprocs*NUM_WORKERS];
    MPI_CHECK(MPI_Gather(populate_times, NUM_WORKERS, MPI_DOUBLE,
                g_populate_times, NUM_WORKERS, MPI_DOUBLE, 0, comm));
    if (0 == rank) {
        double tally = 0;
        cout << " pid populate_time" << endl;
        for(unsigned i=0; i<nprocs*NUM_WORKERS; i++) {
            tally += g_populate_times[i];
            cout << std::setw(4) << std::right << i
                << setw(14) << fixed << g_populate_times[i] << endl;
        }
        cout << "==============================================" << endl;
        cout << setw(4) << right << "T" << setw(14) << fixed << tally << endl;
    }
    delete [] g_populate_times;
#endif
#if DEBUG
    unsigned long *task_counts = new unsigned long[nprocs*NUM_WORKERS];
    MPI_CHECK(MPI_Gather(count, NUM_WORKERS, MPI_UNSIGNED_LONG,
                task_counts, NUM_WORKERS, MPI_UNSIGNED_LONG, 0, comm));
    if (0 == rank) {
        unsigned long tally = 0;
        cout << " pid        ntasks" << endl;
        for(unsigned i=0; i<nprocs*NUM_WORKERS; i++) {
            tally += task_counts[i];
            cout << std::setw(4) << std::right << i
                << setw(14) << task_counts[i] << endl;
        }
        cout << "==============================================" << endl;
        cout << setw(4) << right << "T" << setw(14) << tally << endl;
        if (tally != ntasks) {
            cout << "tally != ntasks\t" << tally << " != " << ntasks << endl;
        }
    }
    delete [] task_counts;
#endif

    amBarrier();

#if defined(THREADED)
#if defined(SET_AFFINITY)
    cpu_set_t cpuset;
    pthread_t thread;
    CPU_ZERO(&cpuset);
    thread = pthread_self();
    CPU_SET(0, &cpuset);
    pthread_setaffinity_np(thread, sizeof(cpu_set_t), &cpuset);
#endif

    pthread_barrier_init(&workersStart, 0, NUM_WORKERS);
    pthread_barrier_init(&workersEnd, 0, NUM_WORKERS);
    pthread_barrier_init(&serverStart, 0, NUM_SERVERS + 1);
    pthread_barrier_init(&serverEnd, 0, NUM_SERVERS + 1);
    asm("mfence");
    for (unsigned i = 1; i < NUM_WORKERS; ++i) {
      pthread_create(&threadHandles[i], NULL, workerThd, &threadRanks[i]);
    }
    pthread_create(&threadHandles[NUM_WORKERS], NULL,
           serverThd, &threadRanks[NUM_WORKERS]);
#endif

#if defined(THREADED)
    serverEnabled = true;
    pthread_barrier_wait(&serverStart);
    MPI_Barrier(comm);
    pthread_barrier_wait(&workersStart);
#endif

    double mytimer = MPI_Wtime();
    utcs[0]->process(0);
    mytimer = MPI_Wtime() - mytimer;
    if (0 == trank(0)) {
        cout << "mytimer=" << mytimer << endl;
    }

#if defined(THREADED)
    pthread_barrier_wait(&workersEnd);
    amBarrierThd();

    serverEnabled = false;
    asm("mfence");
    pthread_barrier_wait(&serverEnd);

    amBarrier();
#endif

#if 0
    for (int i = 0; i < NUM_WORKERS; i++) {
        utcs[i]->restore();
    }
#endif

    amBarrier();
    MPI_Barrier(comm);

    AlignStats * rstats = new AlignStats[NUM_WORKERS*nprocs];
    MPI_Gather(stats, sizeof(AlignStats)*NUM_WORKERS, MPI_CHAR, 
	       rstats, sizeof(AlignStats)*NUM_WORKERS, MPI_CHAR, 
	       0, comm);

    /* synchronously print alignment stats all from process 0 */
    if (0 == rank) {
        AlignStats totals;
        cout << " pid" << rstats[0].getHeader() << endl;      
        for(int i=0; i<nprocs*NUM_WORKERS; i++) {
            totals += rstats[i];
            cout << std::setw(4) << std::right << i << rstats[i] << endl;
        }
        cout << "==============================================" << endl;
        cout << setw(4) << right << "TOT" << totals << endl;
    }
    delete [] rstats;
    rstats=NULL;

    StealingStats stt[NUM_WORKERS];
    for(int i=0; i<NUM_WORKERS; i++) {
      stt[i] = utcs[i]->getStats();
    }

    StealingStats * rstt = new StealingStats[NUM_WORKERS*nprocs];
    MPI_Gather(stt, sizeof(StealingStats)*NUM_WORKERS, MPI_CHAR, 
	       rstt, sizeof(StealingStats)*NUM_WORKERS, MPI_CHAR, 
	       0, comm);

    /* synchronously print stealing stats all from process 0 */
    if (0 == rank) {
      StealingStats tstt;
      cout<<" pid "<<rstt[0].formatString()<<endl;      
      for(int i=0; i<nprocs*NUM_WORKERS; i++) {
	tstt += rstt[i];
	cout<<std::setw(4)<<std::right<<i<<rstt[i]<<endl;
      }
      cout<<"=============================================="<<endl;
      cout<<"TOT "<<tstt<<endl;
    }
    delete [] rstt;
    rstt=NULL;

    amBarrier();
    MPI_Barrier(comm);

#if OUTPUT_EDGES
#if CACHE_RESULTS
    ofstream edge_out(get_edges_filename(rank).c_str());
#endif
#endif
    /* clean up */
    for (int worker=0; worker<NUM_WORKERS; ++worker) {
        free_tbl(tbl[worker], NROW);
        free_int(del[worker], NROW);
        free_int(ins[worker], NROW);
        delete utcs[worker];
#if DUMP_COSTS
        out[worker].close();
#endif
#if OUTPUT_EDGES
#if CACHE_RESULTS
        for (size_t i=0,limit=edge_results[worker].size(); i<limit; ++i) {
            edge_out << edge_results[worker][i] << endl;
        }
#else
        edges[worker].close();
#endif
#endif
    }
#if OUTPUT_EDGES
#if CACHE_RESULTS
    edge_out.close();
#endif
#endif

    TascelConfig::finalize();
    MPI_Comm_free(&comm);
    MPI_Finalize();

    return 0;
}
