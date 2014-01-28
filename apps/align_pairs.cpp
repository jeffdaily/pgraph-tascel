/**
 * @author jeff.daily@pnnl.gov
 *
 * Alignment of an input dataset using work stealing. Each MPI task reads the
 * input file.
 */
#include "config.h"

#include <mpi.h>
#include <tascel.h>
#include <tascel/UniformTaskCollection.h>
#include <tascel/UniformTaskCollSplitHybrid.h>

#include <algorithm>
#include <cassert>
#include <cctype>
#include <cfloat>
#include <cmath>
#include <csignal>
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

#include "AlignStats.hpp"
#include "constants.h"
#include "combinations.h"
#include "alignment.hpp"
#include "mpix.hpp"
#include "Parameters.hpp"
#include "tascelx.hpp"
#include "SequenceDatabase.hpp"
#ifdef USE_GARRAY
#include "SequenceDatabaseGArray.hpp"
#elif HAVE_ARMCI
#include "SequenceDatabaseArmci.hpp"
#else
#include "SequenceDatabaseReplicated.hpp"
#endif

using namespace std;
using namespace tascel;
using namespace pgraph;

#define SEP ","
#define ALL_RESULTS 0
#define OUTPUT_EDGES 1


//#define PAUSE_ON_ERROR
#ifdef PAUSE_ON_ERROR
static void (*SigSegvOrig)(int) = NULL;

static void SigSegvHandler(int sig)
{
    fprintf(stderr,"(%d): Segmentation Violation ... pausing\n",
            getpid() );
    pause();
    if (SigSegvOrig == SIG_DFL) {
        signal(sig, SIG_DFL);
    }
    else if (SigSegvOrig == SIG_IGN) {
    }
    else {
        SigSegvOrig(sig);
    }
}

static void TrapSigSegv()
{
    if ((SigSegvOrig=signal(SIGSEGV, SigSegvHandler)) == SIG_ERR) {
        fprintf(stderr, "TrapSigSegv: error from signal setting SIGSEGV");
        exit(EXIT_FAILURE);
    }
}
#endif


class EdgeResult {
    public:
        unsigned long id1;
        unsigned long id2;
        double a;
        double b;
        double c;
        bool is_edge;

        EdgeResult(
                unsigned long id1, unsigned long id2,
                double a, double b, double c, bool is_edge=true)
            : id1(id1)
            , id2(id2)
            , a(a)
            , b(b)
            , c(c)
            , is_edge(is_edge)
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

int rank = 0;
int nprocs = 0;
UniformTaskCollSplitHybrid** utcs = 0;
AlignStats *stats = 0;
SequenceDatabase *sequences = 0;
vector<EdgeResult> *edge_results = 0;
Parameters parameters;

#if defined(THREADED)
static pthread_t *threadHandles = 0;
static unsigned *threadRanks = 0;
// Synchronization for worker threads
pthread_barrier_t workersStart, workersEnd;
// Synchronization for server thread
pthread_barrier_t serverStart, serverEnd;
volatile bool serverEnabled = true;
#endif


static int trank(int thd)
{
    return (theTwoSided().getProcRank().toInt() * NUM_WORKERS) + thd;
}

static string get_pairs_filename(int rank, const string &dir)
{
    ostringstream str;
    str << dir << "/pairs." << rank << ".txt";
    return str.str();
}


static string get_edges_filename(int rank)
{
    ostringstream str;
    str << "edges." << rank << ".txt";
    return str.str();
}


typedef struct {
    unsigned long id1;
    unsigned long id2;
} task_description;


static size_t parse_memory_budget(const string& value)
{
    long budget = 0;
    char budget_multiplier = 0;
    istringstream iss(value);

    if (isdigit(*value.rbegin())) {
        iss >> budget;
    }
    else {
        iss >> budget >> budget_multiplier;
    }

    if (budget <= 0) {
        cerr << "memory budget must be positive real number" << endl;
        assert(budget > 0);
    }

    if (budget_multiplier == 'b' || budget_multiplier == 'B') {
        budget *= 1; /* byte */
    }
    else if (budget_multiplier == 'k' || budget_multiplier == 'K') {
        budget *= 1024; /* kilobyte */
    }
    else if (budget_multiplier == 'm' || budget_multiplier == 'M') {
        budget *= 1048576; /* megabyte */
    }
    else if (budget_multiplier == 'g' || budget_multiplier == 'G') {
        budget *= 1073741824; /* gigabyte */
    }
    else if (budget_multiplier != 0) {
        cerr << "unrecognized size multiplier" << endl;
        assert(0);
    }

    assert(budget > 0);
    return size_t(budget);
}


static void alignment_task(
        UniformTaskCollection *utc,
        void *_bigd, int bigd_len,
        void *pldata, int pldata_len,
        vector<void *> data_bufs, int thd) {
    task_description *desc = (task_description*)_bigd;
    unsigned long seq_id[2];
    cell_t result;
    bool is_edge_answer = false;
    double t = 0;
    double tt = 0;
    int sscore;
    size_t max_len;

    int open = -10;
    int gap = -1;
    int AOL = 8;
    int SIM = 4;
    int OS = 3;

    tt = MPI_Wtime();

    seq_id[0] = desc->id1;
    seq_id[1] = desc->id2;
    unsigned long s1Len = (*sequences)[seq_id[0]].get_sequence_length();
    unsigned long s2Len = (*sequences)[seq_id[1]].get_sequence_length();
#ifdef LENGTH_FILTER
    int cutOff = AOL * SIM;
    if ((s1Len <= s2Len && (100 * s1Len < cutOff * s2Len))
            || (s2Len < s1Len && (100 * s2Len < cutOff * s1Len))) {
        stats[thd].work_skipped += s1Len * s2Len;
        ++stats[thd].align_skipped;
    }
    else
#endif
    {
        stats[thd].work += s1Len * s2Len;
        ++stats[thd].align_counts;
        t = MPI_Wtime();
#if defined(NOALIGN)
#else
#if USE_SSW
        sequences->align_ssw(seq_id[0], seq_id[1], result.score, result.matches, result.length, open, gap, thd);
#else
        sequences->align(seq_id[0], seq_id[1], result.score, result.matches, result.length, open, gap, thd);
#endif
        is_edge_answer = sequences->is_edge(
                seq_id[0],
                seq_id[1],
                result.score, result.matches, result.length,
                AOL, SIM, OS,
                sscore, max_len);

        if (is_edge_answer || ALL_RESULTS)
        {
#if DEBUG
            cout << trank(thd)
                << ": aligned " << seq_id[0] << " " << seq_id[1]
                << ": (score,ndig,alen)=("
                << result.score << ","
                << result.matches << ","
                << result.length << ")"
                << ": edge? " << is_edge_answer << endl;
#endif
            edge_results[thd].push_back(EdgeResult(
                        seq_id[0], seq_id[1],
#if 1
                        1.0*result.length/max_len,
                        1.0*result.matches/result.length,
                        1.0*result.score/sscore
#else
                        result.length,
                        result.matches,
                        result.score
#endif
#if ALL_RESULTS
                        ,is_edge_answer
#endif
                        ));
            if (is_edge_answer) {
                ++stats[thd].edge_counts;
            }
        }
#endif
        t = MPI_Wtime() - t;
        stats[thd].align_times_tot += t;
        stats[thd].calc_min(t);
        stats[thd].calc_max(t);
    }

    tt = MPI_Wtime() - tt;
    stats[thd].total_times_tot += tt;
}


unsigned long populate_tasks(MPI_Comm comm, string filename)
{
    unsigned long count = 0;
    task_description desc;
#if 1
    ifstream ifs(filename.c_str());
    string line;

    int worker = 0;
    while (getline(ifs, line)) {
        istringstream iss(line);
        size_t id1;
        size_t id2;
        char ignore;
        iss >> id1 >> ignore >> id2;
        assert(ignore == ',');
        desc.id1 = id1;
        desc.id2 = id2;
        utcs[worker]->addTask(&desc, sizeof(desc));
        //cout << trank(worker) << " added " << id1 << "," << id2 << endl;
        worker = (worker + 1) % NUM_WORKERS;
        count++;
    }
#else
    FILE *file = fopen(filename.c_str(), "r");
    if (NULL == file) {
        perror("fopen");
        MPI_Abort(comm, 1);
    }
    fseek(file, 0L, SEEK_END);
    size_t file_size = ftell(file);
    fseek(file, 0L, SEEK_SET);
    char *file_buffer = new char[file_size+1];
    size_t read_count = fread(file_buffer, file_size, 1, file);
    if (0 == read_count) {
        printf("unable to read pairs file\n");
        MPI_Abort(comm, 1);
    }

    mpix_print_sync("read pairs file, now to parse...", comm);

    istringstream iss(file_buffer);
    size_t id1;
    size_t id2;
    char ignore1;
    //char ignore2;
    int worker = 0;
    while (iss >> id1 >> ignore1 >> id2)
        assert(ignore1 == ',');
        //assert(ignore2 == '\n');
        desc.id1 = id1;
        desc.id2 = id2;
        utcs[worker]->addTask(&desc, sizeof(desc));
        cout << trank(worker) << " added " << id1 << "," << id2 << endl;
        worker = (worker + 1) % NUM_WORKERS;
        count++;
#endif

    return count;
}


int main(int argc, char **argv)
{
    MPI_Comm comm = MPI_COMM_NULL;
    vector<string> all_argv;
    unsigned long nCk;

#ifdef PAUSE_ON_ERROR
    TrapSigSegv();
#endif

    /* initialize MPI */
#if defined(THREADED)
    {
        int provided;
        MPI_CHECK(MPI_Init_thread(&argc, &argv,
                    MPI_THREAD_MULTIPLE, &provided));
        assert(provided == MPI_THREAD_MULTIPLE);
    }
#else
    MPI_CHECK(MPI_Init(&argc, &argv));
#endif
    MPI_CHECK(MPI_Comm_dup(MPI_COMM_WORLD, &comm));
    MPI_CHECK(MPI_Comm_rank(comm, &rank));
    MPI_CHECK(MPI_Comm_size(comm, &nprocs));

    /* initialize tascel */
    TascelConfig::initialize(NUM_WORKERS_DEFAULT, comm);
    utcs = new UniformTaskCollSplitHybrid*[NUM_WORKERS];
    stats = new AlignStats[NUM_WORKERS];
    edge_results = new vector<EdgeResult>[NUM_WORKERS];
#if defined(THREADED)
    threadHandles = new pthread_t[NUM_WORKERS + NUM_SERVERS];
    threadRanks = new unsigned[NUM_WORKERS + NUM_SERVERS];
    for (int worker=0; worker<NUM_WORKERS; ++worker) {
        threadRanks[worker] = worker;
    }
#endif

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
            for (size_t j=0; j<all_argv.size(); ++j) {
                printf("[%d] argv[%zd]=%s\n", rank, j, all_argv[j].c_str());
            }
        }
        MPI_Barrier(comm);
    }
#endif

    /* sanity check that we got the correct number of arguments */
    if (all_argv.size() <= 2 || all_argv.size() >= 6) {
        if (0 == rank) {
            if (all_argv.size() <= 1) {
                printf("missing input file\n");
            }
            else if (all_argv.size() <= 2) {
                printf("missing memory budget\n");
            }
            else if (all_argv.size() <= 3) {
                printf("missing parameterst\n");
            }
            else if (all_argv.size() <= 4) {
                printf("missing pairs directory\n");
            }
            else if (all_argv.size() >= 6) {
                printf("too many arguments\n");
            }
            printf("usage: align sequence_file memory_budget\n");
        }
        MPI_Comm_free(&comm);
        MPI_Finalize();
        return 1;
    }
    else if (all_argv.size() >= 4) {
        parameters.parse(all_argv[3].c_str(), comm);
    }
    if (0 == trank(0)) {
        printf("----------------------------------------------\n");
        printf("%-20s: %d\n", "slide size", parameters.window_size);
        printf("%-20s: %d\n", "exactMatch len", parameters.exact_match_len);
        printf("%-20s: %d\n", "AlignOverLongerSeq", parameters.AOL);
        printf("%-20s: %d\n", "MatchSimilarity", parameters.SIM);
        printf("%-20s: %d\n", "OptimalScoreOverSelfScore", parameters.OS);
        printf("----------------------------------------------\n");
    }

#ifdef USE_GARRAY
    sequences = new SequenceDatabaseGArray(all_argv[1],
            parse_memory_budget(all_argv[2].c_str()), '\0');
#elif HAVE_ARMCI
    sequences = new SequenceDatabaseArmci(all_argv[1],
            parse_memory_budget(all_argv[2].c_str()), comm, NUM_WORKERS, '\0');
#else
    sequences = new SequenceDatabaseReplicated(all_argv[1],
            parse_memory_budget(all_argv[2].c_str()), comm, NUM_WORKERS, '\0');
#endif

    /* how many combinations of sequences are there? */
    nCk = binomial_coefficient(sequences->get_global_count(), 2);
    if (0 == trank(0)) {
        printf("brute force %lu C 2 has %lu combinations\n",
                sequences->get_global_count(), nCk);
    }

    unsigned long ntasks = nCk;
    unsigned long global_num_workers = nprocs*NUM_WORKERS;
    unsigned long tasks_per_worker = ntasks / global_num_workers;
    unsigned long max_tasks_per_worker = ntasks / global_num_workers;
    max_tasks_per_worker += ntasks % global_num_workers;
    max_tasks_per_worker *= 10;
#if 1
    unsigned long GB_2 = 536870912U;
    unsigned long GB_4 = 268435456U;
    max_tasks_per_worker = std::min(max_tasks_per_worker, GB_2/sizeof(task_description));
#else
    max_tasks_per_worker = max_tasks_per_worker * 0.001; /* approx. selectivity */
#endif
    if (0 == trank(0)) {
        printf("ntasks=%lu\n", ntasks);
        printf("global_num_workers=%lu\n", global_num_workers);
        printf("tasks_per_worker=%lu\n", tasks_per_worker);
        printf("max_tasks_per_worker=%lu\n", max_tasks_per_worker);
    }
    MPI_Barrier(comm);

    /* the tascel part */
    for (int worker=0; worker<NUM_WORKERS; ++worker)
    {
        //edge_results[worker].reserve(tasks_per_worker);
        UniformTaskCollSplitHybrid*& utc = utcs[worker];
        TslFuncRegTbl *frt = new TslFuncRegTbl();
        TslFunc tf = frt->add(alignment_task);
        TaskCollProps props;
        props.functions(tf, frt)
            .taskSize(sizeof(task_description))
            .maxTasks(max_tasks_per_worker);
        utc = new UniformTaskCollSplitHybrid(props, worker);
    }
    /* add some tasks */
    double populate_time = MPI_Wtime();
    unsigned long count = populate_tasks(comm, get_pairs_filename(rank, all_argv[4]));
    populate_time = MPI_Wtime() - populate_time;
#if DEBUG || 1
    double *g_populate_times = new double[nprocs];
    MPI_CHECK(MPI_Gather(&populate_time, 1, MPI_DOUBLE,
                g_populate_times, 1, MPI_DOUBLE, 0, comm));
    if (0 == rank) {
        double tally = 0;
        cout << " pid populate_time" << endl;
        for(int i=0; i<nprocs; i++) {
            tally += g_populate_times[i];
            cout << std::setw(4) << std::right << i
                << setw(14) << fixed << g_populate_times[i] << endl;
        }
        cout << "==============================================" << endl;
        cout << setw(4) << right << "T" << setw(14) << fixed << tally << endl;
    }
    delete [] g_populate_times;
#endif
#if DEBUG || 1
    unsigned long *task_counts = new unsigned long[nprocs];
    MPI_CHECK(MPI_Gather(&count, 1, MPI_UNSIGNED_LONG,
                task_counts, 1, MPI_UNSIGNED_LONG, 0, comm));
    if (0 == rank) {
        unsigned long tally = 0;
        cout << " pid        ntasks" << endl;
        for(int i=0; i<nprocs; i++) {
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
    set_affinity(0);
    pthread_barrier_init(&workersStart, 0, NUM_WORKERS);
    pthread_barrier_init(&workersEnd, 0, NUM_WORKERS);
    pthread_barrier_init(&serverStart, 0, NUM_SERVERS + 1);
    pthread_barrier_init(&serverEnd, 0, NUM_SERVERS + 1);
    MFENCE
    for (int i = 1; i < NUM_WORKERS; ++i) {
        worker_thread_args *args = new worker_thread_args(
                threadRanks[i], utcs[i], &workersStart, &workersEnd);;
        pthread_create(&threadHandles[i], NULL, worker_thread, args);
    }
    {
        server_thread_args *args = new server_thread_args(
                &serverEnabled, &serverStart, &serverEnd);
        pthread_create(&threadHandles[NUM_WORKERS], NULL,
           server_thread, args);
    }
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
    MFENCE
    pthread_barrier_wait(&serverEnd);

    amBarrier();
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
    delete [] stats;

    StealingStats *stt = new StealingStats[NUM_WORKERS];
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
    delete [] stt;

    amBarrier();
    MPI_Barrier(comm);

#if OUTPUT_EDGES
    ofstream edge_out(get_edges_filename(rank).c_str());
#endif
    /* clean up */
    for (int worker=0; worker<NUM_WORKERS; ++worker) {
        delete utcs[worker];
#if OUTPUT_EDGES
        for (size_t i=0,limit=edge_results[worker].size(); i<limit; ++i) {
            edge_out << edge_results[worker][i] << endl;
        }
#endif
    }
#if OUTPUT_EDGES
    edge_out.close();
#endif

    delete [] utcs;
    delete [] edge_results;
    delete sequences;

    TascelConfig::finalize();
    MPI_Comm_free(&comm);
    MPI_Finalize();

    return 0;
}
