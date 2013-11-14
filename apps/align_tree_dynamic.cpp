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
#include "SuffixBuckets.hpp"
#include "SuffixTree.hpp"

using namespace std;
using namespace tascel;
using namespace pgraph;

#define SEP ","
#define ALL_RESULTS 1


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
                << SEP << edge.is_edge
                ;
            return os;
        }
};

int rank = 0;
int nprocs = 0;
UniformTaskCollSplitHybrid** utcs = 0;
AlignStats *stats = 0;
TreeStats *treestats = 0;
SequenceDatabase *sequences = 0;
#if OUTPUT_EDGES
vector<EdgeResult> *edge_results = 0;
#endif
Parameters parameters;
SuffixBuckets *suffix_buckets = NULL;

#if defined(GLOBAL_DUPLICATES)
set<pair<size_t,size_t> > *pairs = NULL;
#endif

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

static string get_edges_filename(int rank)
{
    ostringstream str;
    str << "edges." << rank << ".txt";
    return str.str();
}


typedef struct {
    unsigned long id1;
    unsigned long id2;
} task_desc_align;


typedef struct {
    unsigned long id1;
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


static void alignment_task(
        UniformTaskCollection *utc,
        void *_bigd, int bigd_len,
        void *pldata, int pldata_len,
        vector<void *> data_bufs, int thd) {
    task_desc_align *desc = (task_desc_align*)_bigd;
    unsigned long seq_id[2];
    cell_t result;
    bool is_edge_answer = false;
    double t = 0;
    int sscore;
    size_t max_len;

    int open = -10;
    int gap = -1;
    int AOL = 8;
    int SIM = 4;
    int OS = 3;

    seq_id[0] = desc->id1;
    seq_id[1] = desc->id2;
    {
        t = MPI_Wtime();
        unsigned long s1Len = (*sequences)[seq_id[0]].get_sequence_length();
        unsigned long s2Len = (*sequences)[seq_id[1]].get_sequence_length();
        stats[thd].work += s1Len * s2Len;
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
        ++stats[thd].align_counts;

        if (is_edge_answer || ALL_RESULTS)
        {
#if DEBUG
            cout << trank(thd)
                << ": aligned " << seq_id[0] << " " << seq_id[1]
                << ": (score,matches,alen)=("
                << result.score << ","
                << result.matches << ","
                << result.length << ")"
                << ": edge? " << is_edge_answer << endl;
#endif
#if OUTPUT_EDGES
            edge_results[thd].push_back(EdgeResult(
                        seq_id[0], seq_id[1],
#if 0
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
#endif
            if (is_edge_answer) {
                ++stats[thd].edge_counts;
            }
        }
        t = MPI_Wtime() - t;
        stats[thd].align_times_tot += t;
        stats[thd].calc_min(t);
        stats[thd].calc_max(t);
    }
}


struct Callback {
    int worker;
    set<pair<size_t,size_t> > local_pairs;

    Callback(int worker)
        : worker(worker)
        , local_pairs()
    {}

    void operator()(pair<size_t,size_t> the_pair) {
        bool proceed = false;
#if defined(GLOBAL_DUPLICATES)
        if (pairs[worker].count(the_pair) == 0) {
            pairs[worker].insert(the_pair);
            proceed = true;
        }
#else
        if (local_pairs.count(the_pair) == 0) {
            local_pairs.insert(the_pair);
            proceed = true;
        }
#endif
        if (proceed) {
            task_desc_align desc;
            desc.id1 = the_pair.first;
            desc.id2 = the_pair.second;
#if DEBUG
            cout << trank(worker) << " added " << desc.id1 << "," << desc.id2 << endl;
#endif
            utcs[worker]->addTask(&desc, sizeof(desc));
        }
    }
};


static void tree_task(
        UniformTaskCollection *utc,
        void *_bigd, int bigd_len,
        void *pldata, int pldata_len,
        vector<void *> data_bufs, int worker) {
    task_desc_tree *tree_desc = (task_desc_tree*)_bigd;
    unsigned long i = tree_desc->id1;
    set<pair<size_t,size_t> > local_pairs;
    assert(i < suffix_buckets->buckets_size);
    /* suffix tree construction & processing */
    if (NULL != suffix_buckets->buckets[i].suffixes) {
        double atimer = MPI_Wtime();
        SuffixTree *tree = new SuffixTree(
                sequences, &(suffix_buckets->buckets[i]), parameters);
        treestats[worker].size_internal += tree->get_size_internal();
        treestats[worker].fanout += tree->get_avg_fanout();
        treestats[worker].avgdepth += tree->get_avg_depth();
        treestats[worker].deepest += tree->get_deepest();
        treestats[worker].suffix_avg_length += tree->get_suffix_avg_length();
        treestats[worker].suffix_max_length += tree->get_suffix_max_length();
        treestats[worker].time_build += MPI_Wtime() - atimer;
        atimer = MPI_Wtime();
#if defined(CALLBACK)
        tree->generate_pairs_cb(Callback(worker));
#else
        tree->generate_pairs(local_pairs);
#endif
        treestats[worker].time_process += MPI_Wtime() - atimer;
        delete tree;
    }
    
    task_desc_align desc;
    for (set<pair<size_t,size_t> >::iterator it=local_pairs.begin();
            it!=local_pairs.end(); ++it) {
#if defined(GLOBAL_DUPLICATES)
        if (pairs[worker].count(*it) == 0) {
            pairs[worker].insert(*it);
        }
        else {
            continue;
        }
#endif
        desc.id1 = it->first;
        desc.id2 = it->second;
#if DEBUG
        cout << trank(worker) << " added " << desc.id1 << "," << desc.id2 << endl;
#endif
        utcs[worker]->addTask(&desc, sizeof(desc));
    }
}


int main(int argc, char **argv)
{
    MPI_Comm comm = MPI_COMM_NULL;
    vector<string> all_argv;
    unsigned long nCk;
    time_t t1 = 0;                  /* start timer */
    time_t t2 = 0;                  /* stop timer */

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
    TascelConfig::initialize(NUM_WORKERS_DEFAULT, comm);
    utcs = new UniformTaskCollSplitHybrid*[NUM_WORKERS];
    stats = new AlignStats[NUM_WORKERS];
    treestats = new TreeStats[NUM_WORKERS];
#if OUTPUT_EDGES
    edge_results = new vector<EdgeResult>[NUM_WORKERS];
#endif
#if defined(GLOBAL_DUPLICATES)
    pairs = new set<pair<size_t,size_t> >[NUM_WORKERS];
#endif
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

    sequences = new SequenceDatabaseArmci(all_argv[1],
            parse_memory_budget(all_argv[2].c_str()), comm, NUM_WORKERS, DOLLAR);

    /* how many combinations of sequences are there? */
    nCk = binomial_coefficient(sequences->get_global_count(), 2);
    if (0 == trank(0)) {
        printf("brute force %lu C 2 has %lu combinations\n",
                sequences->get_global_count(), nCk);
    }

    double selectivity = 1.0;
    unsigned long nalignments = (long)(0.5+selectivity*nCk);
    unsigned long ntasks = nalignments;
    unsigned long global_num_workers = nprocs*NUM_WORKERS;
    unsigned long max_tasks_per_worker = ntasks / global_num_workers;
#if 1
    max_tasks_per_worker += ntasks % global_num_workers;
    max_tasks_per_worker *= 10;
    unsigned long GB = 1073741824;
    unsigned long GB_2 = 536870912;
    unsigned long GB_4 = 268435456;
    max_tasks_per_worker = std::min(max_tasks_per_worker, GB_2/sizeof(task_desc_align));
#else
    max_tasks_per_worker = max_tasks_per_worker * 0.001; /* approx. selectivity */
#endif
    if (0 == trank(0)) {
        printf("selectivity=%lf\n", selectivity);
        printf("nalignments=%lu\n", nalignments);
        printf("ntasks=%lu\n", ntasks);
        printf("global_num_workers=%lu\n", global_num_workers);
        printf("max_tasks_per_worker=%lu\n", max_tasks_per_worker);
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

    MPI_Barrier(comm);

    /* the tascel part */
    (void) time(&t1);
    suffix_buckets = new SuffixBuckets(sequences, parameters, comm);
    (void) time(&t2);
    if (0 == trank(0)) {
        printf("Bucketing finished in <%lld> secs\n", (long long)(t2-t1));
    }
    for (int worker=0; worker<NUM_WORKERS; ++worker)
    {
#if OUTPUT_EDGES
        edge_results[worker].reserve(max_tasks_per_worker);
#endif
        UniformTaskCollSplitHybrid*& utc = utcs[worker];
        TslFuncRegTbl *frt = new TslFuncRegTbl();
        TslFunc tf = frt->add(alignment_task);
        TslFunc tf2 = frt->add(tree_task);
        TaskCollProps props;
        props.functions(tf, frt, tf2)
            .taskSize(sizeof(task_desc_align))
            .maxTasks(max_tasks_per_worker)
            .taskSize2(sizeof(task_desc_tree))
            .maxTasks2(long(pow(26.0,parameters.window_size))/2); // TODO too big
        utc = new UniformTaskCollSplitHybrid(props, worker);
    }

    size_t even_split = suffix_buckets->bucket_size_total / NUM_WORKERS;
    int worker = 0;
    double poptimer = MPI_Wtime();
    for (size_t i=suffix_buckets->first_bucket;
            i < suffix_buckets->buckets_size; ++i) {
        if (suffix_buckets->bucket_owner[i] == rank) {
            if (NULL != suffix_buckets->buckets[i].suffixes) {
                task_desc_tree desc;
                desc.id1 = i;
#if DEBUG
                cout << trank(worker) << " added " << desc.id1 << endl;
#endif
                utcs[worker]->addTask2(&desc, sizeof(task_desc_tree));
                treestats[worker].trees += 1;
                treestats[worker].suffixes += suffix_buckets->buckets[i].size;
                if (treestats[worker].suffixes >= even_split) {
                    treestats[worker].times = MPI_Wtime() - poptimer;
                    poptimer = MPI_Wtime();
                    ++worker;
                    if (worker >= NUM_WORKERS) {
                        worker = NUM_WORKERS-1;
                    }
                }
            }
        }
        else {
            break;
        }
    }
    treestats[worker].times = MPI_Wtime() - poptimer;

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
    TreeStats * gtreestats = new TreeStats[NUM_WORKERS*nprocs];
    MPI_Gather(treestats, sizeof(TreeStats)*NUM_WORKERS, MPI_CHAR, 
	       gtreestats, sizeof(TreeStats)*NUM_WORKERS, MPI_CHAR, 
	       0, comm);

    /* synchronously print alignment stats all from process 0 */
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

    AlignStats * rstats = new AlignStats[NUM_WORKERS*nprocs];
    MPI_Gather(stats, sizeof(AlignStats)*NUM_WORKERS, MPI_CHAR, 
	       rstats, sizeof(AlignStats)*NUM_WORKERS, MPI_CHAR, 
	       0, comm);

    /* synchronously print alignment stats all from process 0 */
    if (0 == rank) {
        AlignStats totals;
        AlignStats mins;
        AlignStats maxs;
        cout << setw(4) << right << "pid" << rstats[0].getHeader() << endl;
        for(int i=0; i<nprocs*NUM_WORKERS; i++) {
            totals += rstats[i];
            mins < rstats[i];
            maxs > rstats[i];
            cout << std::setw(4) << right << i << " " << rstats[i] << endl;
        }
        cout << "==============================================" << endl;
        cout << setw(4) << right << "TOT " << totals << endl;
        cout << setw(4) << right << "MIN " << mins << endl;
        cout << setw(4) << right << "MAX " << maxs << endl;
        cout << setw(4) << right << "AVG " << (totals/(nprocs*NUM_WORKERS)) << endl;
        cout << setw(4) << right << "STD " << (totals.stddev(nprocs*NUM_WORKERS, rstats)) << endl;
        cout << "==============================================" << endl;
    }
    delete [] rstats;
    rstats=NULL;

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

    delete suffix_buckets;
    delete [] utcs;
#if OUTPUT_EDGES
    delete [] edge_results;
#endif
#if defined(GLOBAL_DUPLICATES)
    delete [] pairs;
#endif
    delete sequences;

    TascelConfig::finalize();
    MPI_Comm_free(&comm);
    MPI_Finalize();

    return 0;
}
