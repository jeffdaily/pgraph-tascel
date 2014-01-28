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
#include <tascel/tmpi.h>

#include <algorithm>
#include <cassert>
#include <cctype>
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

#include "AlignStats.hpp"
#include "TreeStats.hpp"
#include "constants.h"
#include "combinations.h"
#include "alignment.hpp"
#include "mpix.hpp"
#include "tascelx.hpp"
#include "SequenceDatabase.hpp"
#include "SequenceDatabaseArmci.hpp"
#include "SequenceDatabaseReplicated.hpp"
#include "SuffixBuckets.hpp"
#include "SuffixTree.hpp"

#define DEBUG 1
#define USE_CALLBACK 1

using namespace std;
using namespace tascel;
using namespace pgraph;


int rank = 0;
int nprocs = 0;
UniformTaskCollSplitHybrid** utcs = 0;
TreeStats *treestats = 0;
SequenceDatabase *sequences = 0;
vector<pair<size_t,size_t> > *pair_results = 0;
double *timeout = 0;
Parameters parameters;
SuffixBuckets2 *suffix_buckets = NULL;
#if defined(THREADED)
static pthread_t *threadHandles = 0;
static unsigned *threadRanks = 0;
// Synchronization for worker threads
pthread_barrier_t workersStart, workersEnd;
// Synchronization for server thread
pthread_barrier_t serverStart, serverEnd;
volatile bool serverEnabled = true;
#endif
long ntasks = 0;
long **counter = NULL;
PthreadMutex mutex;


static int trank(int thd)
{
    return (theTwoSided().getProcRank().toInt() * NUM_WORKERS) + thd;
}


static string get_pairs_filename(int rank)
{
    ostringstream str;
    str << "pairs." << rank << ".txt";
    return str.str();
}


typedef struct {
    unsigned long id;
} task_desc_tree;


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


#if 0
static size_t powz(size_t base, size_t n)
{
    size_t p = 1;

    for(/*empty*/; n > 0; --n) {
        assert(p < SIZE_MAX/base);
        p *= base;
    }

    return p;
}


static string entry_string(size_t index, size_t k)
{
    string retval = "";
    for (unsigned long i = k-1; i >= 0; --i) {
        size_t tmp = powz(SIGMA,i);
        size_t quo = index / tmp;
        size_t rem = index % tmp;
        assert(quo < SIGMA);
        retval += char('A'+char(quo));
    }
    return retval;
}
#endif


#define CALLBACK_PAIRS 0
struct Callback {
    int worker;
    vector<pair<size_t,size_t> > pairs;

    Callback(int worker)
        : worker(worker)
        , pairs()
    {}

    bool operator()(pair<size_t,size_t> the_pair) {
        bool retval = false;
        double time = MPI_Wtime();
#if CALLBACK_PAIRS
        pairs.push_back(the_pair);
#else
        pair_results[worker].push_back(the_pair);
#endif
        if (time - timeout[worker] > 30.0) { /* seconds */
            retval = true;
        }
        return retval;
    }
};



static void tree_task(
        UniformTaskCollection *utc,
        void *_bigd, int bigd_len,
        void *pldata, int pldata_len,
        vector<void *> data_bufs, int worker) {
    task_desc_tree *tree_desc = (task_desc_tree*)_bigd;
    unsigned long i = tree_desc->id;
#if DEBUG && 0
    cout << trank(worker) << " tree_task id=" << i << endl;
#endif
#if USE_SET
    set<pair<size_t,size_t> > local_pairs;
#endif
    /* suffix tree construction & processing */
    Bucket *bucket = suffix_buckets->get(i);
    if (NULL != bucket->suffixes) {
        double atimer = MPI_Wtime();
        SuffixTree *tree = new SuffixTree(
                sequences, bucket, parameters);
        treestats[worker].size_internal += tree->get_size_internal();
        treestats[worker].fanout += tree->get_avg_fanout();
        treestats[worker].avgdepth += tree->get_avg_depth();
        treestats[worker].deepest += tree->get_deepest();
        treestats[worker].suffix_avg_length += tree->get_suffix_avg_length();
        treestats[worker].suffix_max_length += tree->get_suffix_max_length();
        treestats[worker].time_build += MPI_Wtime() - atimer;
        treestats[worker].times += MPI_Wtime() - atimer;
        treestats[worker].trees += 1;
        treestats[worker].suffixes += bucket->size;
        atimer = MPI_Wtime();
#if USE_SET
        tree->generate_pairs(local_pairs);
#elif USE_CALLBACK
        timeout[worker] = MPI_Wtime();
        Callback callback(worker);
        bool ret = tree->generate_pairs_cb(callback);
        if (ret) {
            cout << trank(worker)
                << " tree_task id=" << i
                //<< " st=" << entry_string(i,parameters.window_size)
                << " aborted" << endl;
        }
#if CALLBACK_PAIRS
        else {
            pair_results[worker].insert(pair_results[worker].end(),
                    callback.pairs.begin(), callback.pairs.end());
        }
#endif
#else
        tree->generate_pairs(pair_results[worker]);
#endif
        treestats[worker].time_process += MPI_Wtime() - atimer;
        treestats[worker].times += MPI_Wtime() - atimer;
        delete tree;
        suffix_buckets->rem(i, bucket);
    }
    
#if 0
    pair_results[worker].insert(pair_results[worker].end(),
            local_pairs.begin(), local_pairs.end());
#endif

#if ARMCI_COUNTER
    long count;
    {
        LockGuard<PthreadMutex> guard(mutex);
        (void)ARMCI_Rmw(ARMCI_FETCH_AND_ADD_LONG, &count, counter[0], 1, 0);
    }

    if (0 == trank(worker)) {
        cout << count << "/" << ntasks << endl;
    }
#endif
}


ostream& operator << (ostream &os, const pair<size_t,size_t> &p)
{
    os << p.first << "," << p.second;
    return os;
}


int main(int argc, char **argv)
{
    MPI_Comm comm = MPI_COMM_NULL;
    vector<string> all_argv;
    unsigned long nCk;
    time_t t1 = 0;                  /* start timer */
    time_t t2 = 0;                  /* stop timer */

    try {

    /* initialize MPI */
#if defined(THREADED)
    {
        int provided;
        MPI_CHECK(TMPI_Init_thread(&argc, &argv,
                    MPI_THREAD_MULTIPLE, &provided));
        assert(provided == MPI_THREAD_MULTIPLE);
    }
#else
    MPI_CHECK(TMPI_Init(&argc, &argv));
#endif
    MPI_CHECK(MPI_Comm_dup(MPI_COMM_WORLD, &comm));
    MPI_CHECK(MPI_Comm_rank(comm, &rank));
    MPI_CHECK(MPI_Comm_size(comm, &nprocs));

    /* initialize tascel */
    double totaltime = MPI_Wtime();
    TascelConfig::initialize(NUM_WORKERS_DEFAULT, comm);
    utcs = new UniformTaskCollSplitHybrid*[NUM_WORKERS];
    treestats = new TreeStats[NUM_WORKERS];
    pair_results = new vector<pair<size_t,size_t> >[NUM_WORKERS];
    timeout = new double[NUM_WORKERS];
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

#if DEBUG && 0
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
    if (all_argv.size() <= 2 || all_argv.size() >= 5) {
        if (0 == rank) {
            if (all_argv.size() <= 1) {
                printf("missing input file\n");
            }
            else if (all_argv.size() <= 2) {
                printf("missing memory budget\n");
            }
            else if (all_argv.size() >= 5) {
                printf("too many arguments\n");
            }
            printf("usage: align sequence_file memory_budget param\n");
        }
        MPI_Comm_free(&comm);
        MPI_Finalize();
        return 1;
    }
    else if (all_argv.size() >= 4) {
        parameters.parse(all_argv[3].c_str(), comm);
    }

    char hostname[256];
    int hostname_retval = gethostname(hostname, 256);
    assert(0 == hostname_retval);
    mpix_print_sync("hostname", string(hostname), comm);

#define GB (1073741824U)
#if 1
    sequences = new SequenceDatabaseArmci(all_argv[1],
            parse_memory_budget(all_argv[2].c_str()), comm, NUM_WORKERS, DOLLAR);
#else
    sequences = new SequenceDatabaseReplicated(all_argv[1],
            parse_memory_budget(all_argv[2].c_str()), comm, NUM_WORKERS, DOLLAR);
#endif

    (void) time(&t1);
    suffix_buckets = new SuffixBuckets2(sequences, parameters, comm);
    mpix_print_sync("after SuffixBuckets2 ctor", comm);
    (void) time(&t2);
    if (0 == trank(0)) {
        printf("Bucketing finished in <%lld> secs\n", (long long)(t2-t1));
        fflush(stdout);
    }

    unsigned long global_num_workers = nprocs*NUM_WORKERS;
    unsigned long max_tasks_per_worker = suffix_buckets->n_buckets/global_num_workers + 1;
    unsigned long max_tasks_per_worker_orig = max_tasks_per_worker;
    unsigned long GB_4 = 268435456U;
    unsigned long GB_2 = 536870912U;
    //max_tasks_per_worker = std::min(max_tasks_per_worker, GB_4/sizeof(task_desc_tree));
    max_tasks_per_worker = GB/sizeof(task_desc_tree);
    if (0 == trank(0)) {
        printf("ntasks=%lu\n", suffix_buckets->n_buckets);
        printf("global_num_workers=%lu\n", global_num_workers);
        printf("     max_tasks_per_worker=%lu\n", max_tasks_per_worker);
        printf("max_tasks_per_worker_orig=%lu\n", max_tasks_per_worker_orig);
        printf("                     GB_2=%lu\n", GB_2);
        fflush(stdout);
    }
    if (0 == trank(0)) {
        printf("----------------------------------------------\n");
        printf("%-20s: %d\n", "slide size", parameters.window_size);
        printf("%-20s: %d\n", "exactMatch len", parameters.exact_match_len);
        printf("%-20s: %d\n", "AlignOverLongerSeq", parameters.AOL);
        printf("%-20s: %d\n", "MatchSimilarity", parameters.SIM);
        printf("%-20s: %d\n", "OptimalScoreOverSelfScore", parameters.OS);
        printf("----------------------------------------------\n");
        fflush(stdout);
    }

    MPI_Barrier(comm);

    /* the tascel part */
    for (int worker=0; worker<NUM_WORKERS; ++worker)
    {
        UniformTaskCollSplitHybrid*& utc = utcs[worker];
        TslFuncRegTbl *frt = new TslFuncRegTbl();
        TslFunc tf = frt->add(tree_task);
        TaskCollProps props;
        props.functions(tf, frt)
            .taskSize(sizeof(task_desc_tree))
            .maxTasks(max_tasks_per_worker);
        utc = new UniformTaskCollSplitHybrid(props, worker);
    }

    MPI_Barrier(comm);

    if (0 == trank(0)) {
        printf("created UTCs\n");
        fflush(stdout);
    }

    int worker = 0;
    int task_count = 0;
    for (size_t j=rank; j<suffix_buckets->n_buckets; j+=nprocs) {
        int bucket_index = j / nprocs;
        if (suffix_buckets->buckets[bucket_index].size > 0) {
            task_desc_tree desc;
            desc.id = j;
#if DEBUG && 0
            cout << trank(worker) << " added " << desc.id << endl;
#endif
            utcs[worker]->addTask(&desc, sizeof(task_desc_tree));
            worker = (worker + 1) % NUM_WORKERS;
            task_count += 1;
        }
    }

    MPI_Barrier(comm);

    mpix_print_sync("done adding tasks", rank, comm);
    int *gtask_counts = new int[nprocs];
    MPI_Gather(&task_count, 1, MPI_INT, gtask_counts, 1, MPI_INT, 0, comm);
    if (0 == rank) {
        cout << setw(4) << right << "pid " << "task-counts" << endl;
        for(int i=0; i<nprocs; i++) {
            ntasks += gtask_counts[i];
            cout << std::setw(4) << right << i << " " << gtask_counts[i] << endl;
        }
    }
    delete [] gtask_counts;

    /* this could slow things down but we need some way of debugging
     * global progress ... */
#if ARMCI_COUNTER
    long bytes = 0;
    int retval = 0;
    counter = new long*[nprocs];
    bytes = (0 == rank) ? sizeof(long) : 0;
    retval = ARMCI_Malloc((void**)counter, bytes);
    assert(0 == retval);
    ARMCI_Barrier();
    if (0 == rank) {
        counter[0][0] = 0; /* init counter */
    }
    ARMCI_Barrier();
#endif

    MPI_Barrier(comm);

    amBarrier();

    MPI_Barrier(comm);

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

#if 1
    size_t *gremote_buckets = new size_t[nprocs];
    MPI_Gather(&suffix_buckets->count_remote_buckets, sizeof(size_t), MPI_CHAR,
            gremote_buckets, sizeof(size_t), MPI_CHAR, 0, comm);
    if (0 == rank) {
        cout << setw(4) << right << "pid " << "remote-buckets" << endl;
        for(int i=0; i<nprocs; i++) {
            cout << std::setw(4) << right << i << " " << gremote_buckets[i] << endl;
        }
    }
    delete [] gremote_buckets;
    size_t *gremote_suffixes = new size_t[nprocs];
    MPI_Gather(&suffix_buckets->count_remote_suffixes, sizeof(size_t), MPI_CHAR,
            gremote_suffixes, sizeof(size_t), MPI_CHAR, 0, comm);
    if (0 == rank) {
        cout << setw(4) << right << "pid " << "remote-suffixes" << endl;
        for(int i=0; i<nprocs; i++) {
            cout << std::setw(4) << right << i << " " << gremote_suffixes[i] << endl;
        }
    }
    delete [] gremote_suffixes;
#endif

#if 1
    TreeStats * gtreestats = new TreeStats[NUM_WORKERS*nprocs];
    MPI_Gather(treestats, sizeof(TreeStats)*NUM_WORKERS, MPI_CHAR, 
	       gtreestats, sizeof(TreeStats)*NUM_WORKERS, MPI_CHAR, 
	       0, comm);

    /* synchronously print tree stats all from process 0 */
    if (0 == rank) {
        TreeStats totals;
        TreeStats mins = gtreestats[0];
        TreeStats maxs = gtreestats[0];
        cout << setw(4) << right << "pid " << gtreestats[0].getHeader() << endl;
        for(int i=0; i<nprocs*NUM_WORKERS; i++) {
            totals += gtreestats[i];
            mins < gtreestats[i];
            maxs > gtreestats[i];
            cout << std::setw(4) << right << i << " " << gtreestats[i] << endl;
        }
        cout << "==============================================" << endl;
        cout << setw(4) << right << "TOT " << totals << endl;
        cout << setw(4) << right << "MIN " << mins << endl;
        cout << setw(4) << right << "MAX " << maxs << endl;
        cout << setw(4) << right << "AVG " << (totals/(nprocs*NUM_WORKERS)) << endl;
        cout << setw(4) << right << "STD " << (totals.stddev(nprocs*NUM_WORKERS, gtreestats)) << endl;
        cout << "==============================================" << endl;
    }
    delete [] gtreestats;
    gtreestats=NULL;
#endif

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
            cout<<std::setw(4)<<std::right<<i<<" "<<rstt[i]<<endl;
        }
        cout<<"=============================================="<<endl;
        cout<<"TOT "<<tstt<<endl;
    }
    delete [] rstt;
    rstt=NULL;
    delete [] stt;

    amBarrier();
    MPI_Barrier(comm);

    double write_time = MPI_Wtime();
    ofstream pair_out(get_pairs_filename(rank).c_str());
    /* clean up */
    for (int worker=0; worker<NUM_WORKERS; ++worker) {
        delete utcs[worker];
#if SORT_THEM || 1
        sort( pair_results[worker].begin(), pair_results[worker].end() );
        pair_results[worker].erase( unique( pair_results[worker].begin(), pair_results[worker].end() ), pair_results[worker].end() );
#endif
        for (size_t i=0,limit=pair_results[worker].size(); i<limit; ++i) {
            pair_out << pair_results[worker][i] << endl;
        }
    }
    pair_out.close();
    amBarrier();
    MPI_Barrier(comm);
    write_time = MPI_Wtime() - write_time;
    if (0 == trank(0)) {
        cout << "write time = " << write_time << endl;
    }

    delete suffix_buckets;
    delete [] utcs;
    delete [] pair_results;
    delete sequences;
#if ARMCI_COUNTER
    retval = ARMCI_Free(counter[nprocs]);
    delete [] counter;
#endif

    totaltime = MPI_Wtime() - totaltime;
    if (0 == trank(0)) {
        cout << "totaltime = " << totaltime << endl;
    }
    TascelConfig::finalize();
    MPI_Comm_free(&comm);
    MPI_Finalize();

    } catch (const std::exception& ex) {
        cout << rank << "\tstd::exception\t" << ex.what() << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    } catch (const std::string& ex) {
        cout << rank << "\tstring exception\t" << ex << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    } catch (...) {
        cout << rank << "\tunknown exception\t" << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    return 0;
}
