/**
 * @author jeff.daily@pnnl.gov
 *
 * A first attempt at a brute force alignment of an input dataset using work
 * stealing. Each MPI task reads the input file.
 */
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
#include <stdexcept>
#include <string>
#include <strstream>
#include <vector>

#include "AlignStats.h"
#include "combinations.h"
#include "dynamic.h"
#include "mpix.h"

using namespace std;
using namespace tascel;

#define DEBUG 0

#ifndef SORT_TASKS_LOCALLY
#define SORT_TASKS_LOCALLY 0
#endif

#ifndef SORT_SEQUENCES
#define SORT_SEQUENCES 0
#endif

#ifndef DUMP_COSTS
#define DUMP_COSTS 0
#endif

int rank = 0;
int nprocs = 0;
cell_t **tbl[NUM_WORKERS];
int **del[NUM_WORKERS];
int **ins[NUM_WORKERS];
vector<string> sequences;
size_t max_seq_len = 0;
ProcGroup* pgrp = NULL;
UniformTaskCollSplitHybrid* utcs[NUM_WORKERS];
#if DUMP_COSTS
ofstream out[NUM_WORKERS];
#endif
ofstream edges[NUM_WORKERS];
AlignStats stats[NUM_WORKERS];
// Synchronization for worker threads
pthread_barrier_t workersStart, workersEnd;
// Synchronization for server thread
pthread_barrier_t serverStart, serverEnd;
static pthread_t threadHandles[NUM_WORKERS + NUM_SERVERS];
static unsigned threadRanks[NUM_WORKERS + NUM_SERVERS];
volatile bool serverEnabled = true;
is_edge_param_t param;


bool longest_string_first(const string &i, const string &j)
{
    return i.size() > j.size();
}


/* each process indexes the file_buffer */
void parse_sequence_buffer(char *file_buffer, long file_size, size_t &max_seq_len)
{
    long seg_count = 0;
    istrstream input_stream(const_cast<const char*>(file_buffer), file_size);
    string line;
    string sequence;

    /* each process counts how many '>' characters are in the file_buffer */
    for (int i=0; i<file_size; ++i) {
        if (file_buffer[i] == '>') {
            ++seg_count;
        }
    }

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
}


void sort_sequences(MPI_Comm comm)
{
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

    int ret = pthread_setaffinity_np(thread, sizeof(cpu_set_t), &cpuset);
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


static string get_edges_filename(int thd)
{
    ostringstream str;
    str << "edges." << trank(thd) << ".txt";
    return str.str();
}


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


static void alignment_task(
        UniformTaskCollection *utc,
        void *_bigd, int bigd_len,
        void *pldata, int pldata_len,
        vector<void *> data_bufs, int thd) {
    task_description *desc = (task_description*)_bigd;
    unsigned long seq_id[2];
    cell_t result;
    int is_edge_answer = 0;
    double t = 0;
    int sscore;
    int maxLen;

    try {
    seq_id[0] = desc->id1;
    seq_id[1] = desc->id2;
    t = MPI_Wtime();
    affine_gap_align(
            sequences.at(seq_id[0]).c_str(),
            sequences.at(seq_id[0]).size(),
            sequences.at(seq_id[1]).c_str(),
            sequences.at(seq_id[1]).size(),
            &result, tbl[thd], del[thd], ins[thd]);
    is_edge_answer = is_edge(result,
            sequences.at(seq_id[0]).c_str(),
            sequences.at(seq_id[0]).size(),
            sequences.at(seq_id[1]).c_str(),
            sequences.at(seq_id[1]).size(),
            param, &sscore, &maxLen);
    ++stats[thd].align_counts;
    } catch (std::out_of_range &ex) {
        cerr << "[" << trank(thd) << "] seq_id[0] = " << seq_id[0] << endl;
        cerr << "[" << trank(thd) << "] seq_id[1] = " << seq_id[1] << endl;
        cerr << "[" << trank(thd) << "] sequences.size() = " << sequences.size() << endl;
        throw;
    }

    if (is_edge_answer) {
#if DEBUG
        cout << trank(thd)
            << ": aligned " << seq_id[0] << " " << seq_id[1]
            << ": (score,ndig,alen)=("
            << result.score << ","
            << result.ndig << ","
            << result.alen << ")"
            << ": (alen/max_len,nmatch/alen,score/sscore)=("
            << 1.0*result.alen/maxLen << ","
            << 1.0*result.ndig/result.alen << ","
            << 1.0*result.score/sscore << ")"
            << endl;
#endif
        edges[thd] << seq_id[0] << "\t" << seq_id[1] << "\t"
            << 1.0*result.alen/maxLen << "\t"
            << 1.0*result.ndig/result.alen << "\t"
            << 1.0*result.score/sscore << endl;
        ++stats[thd].edge_counts;
    }
    t = MPI_Wtime() - t;
#if DUMP_COSTS
    out[thd]
        << seq_id[0] << "\t"
        << seq_id[1] << "\t"
        << sequences[seq_id[0]].size() << "\t"
        << sequences[seq_id[1]].size() << "\t"
        << sequences[seq_id[0]].size()*sequences[seq_id[1]].size() << "\t"
        << t << endl;
#endif
    stats[thd].align_times_tot += t;
    stats[thd].calc_min(t);
    stats[thd].calc_max(t);
}


unsigned long populate_tasks(
        int ntasks, int tasks_per_worker, int worker, char *pair_file_buffer)
{
    task_description desc;
    int wrank = trank(worker);
    unsigned long count = 0;
    unsigned long i;
    unsigned long lower_limit = wrank*tasks_per_worker;
    unsigned long upper_limit = lower_limit + tasks_per_worker;
    unsigned long remainder = ntasks % (nprocs*NUM_WORKERS);
    unsigned long datum_size = sizeof(int)*2 + sizeof(char);

    for (i=lower_limit; i<upper_limit; ++i) {
        int id1 = *((int*)(pair_file_buffer+i*datum_size));
        int id2 = *((int*)(pair_file_buffer+i*datum_size+sizeof(int)));
        char newline = *(pair_file_buffer+i*datum_size+sizeof(int)*2);
        assert(id1 >= 0 && id1 < sequences.size());
        assert(id2 >= 0 && id2 < sequences.size());
        assert(newline == '\n');
        desc.id1 = id1;
        desc.id2 = id2;
#if DEBUG
        cout << wrank << " added " << desc.id1 << "," << desc.id2 << endl;
#endif
        utcs[worker]->addTask(&desc, sizeof(desc));
        count++;
    }
    /* if I'm the last worker, add the remainder of the tasks */
    if (wrank == nprocs*NUM_WORKERS-1) {
        for (/*ignore*/; i<upper_limit+remainder; ++i) {
            int id1 = *((int*)(pair_file_buffer+i*datum_size));
            int id2 = *((int*)(pair_file_buffer+i*datum_size+sizeof(int)));
            char newline = *(pair_file_buffer+i*datum_size+sizeof(int)*2);
            assert(id1 >= 0 && id1 < sequences.size());
            assert(id2 >= 0 && id2 < sequences.size());
            assert(newline == '\n');
            desc.id1 = id1;
            desc.id2 = id2;
#if DEBUG
            cout << wrank << " added " << desc.id1 << "," << desc.id2 << endl;
#endif
            utcs[worker]->addTask(&desc, sizeof(desc));
            count++;
        }
    }

    return count;
}


int main(int argc, char **argv)
{
    MPI_Comm comm = MPI_COMM_NULL;
    int provided;
    vector<string> all_argv;
    long file_size = -1;
    char *file_buffer = NULL;
    long pair_file_size = -1;
    char *pair_file_buffer = NULL;
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
    TascelConfig::initialize(NUM_WORKERS, comm);
    pgrp = ProcGroup::construct();
    for (int worker=0; worker<NUM_WORKERS; ++worker) {
        threadRanks[worker] = worker;
    }

    /* initialize dynamic code */
    init_map(SIGMA);
    param.AOL = 80;
    param.SIM = 40;
    param.OS = 30;

    /* MPI standard does not guarantee all procs receive argc and argv */
    mpix_bcast_argv(comm, argc, argv, all_argv);

#if DEBUG
    /* print the command line arguments */
    mpix_print_sync(comm, "all_argv", all_argv);
#endif

    /* sanity check that we got the correct number of arguments */
    if (all_argv.size() <= 1 || all_argv.size() >= 4) {
        if (0 == rank) {
            if (all_argv.size() <= 1) {
                printf("missing input file\n");
            }
            else if (all_argv.size() <= 2) {
                printf("missing pairs file\n");
            }
            else {
                printf("too many arguments\n");
            }
            printf("usage: brute_force sequence_file pairs_file\n");
        }
        TascelConfig::finalize();
        MPI_Comm_free(&comm);
        MPI_Finalize();
        return 1;
    }

    /* read sequence file on all procs */
    mpix_read_file(comm, all_argv[1], file_buffer, file_size);

    /* each process indexes the file_buffer */
    parse_sequence_buffer(file_buffer, file_size, max_seq_len);
    delete [] file_buffer;

#if SORT_SEQUENCES
    sort_sequences(comm);
#endif

#if DEBUG
    mpix_print_sync(comm, "sequences.size()", sequences.size());
    mpix_print_sync(comm, "sequences.capacity()", sequences.capacity());
    mpix_print_sync(comm, "max_seq_len", max_seq_len);
#endif

#if DEBUG
    /* print the first and last sequence on each process */
    mpix_print_sync(comm, "sequences[0]", sequences[0]);
    mpix_print_sync(comm, "sequences.back()", sequences.back());
#endif

    /* how many combinations of sequences are there? */
    nCk = binomial_coefficient(sequences.size(), 2);
    if (0 == trank(0)) {
        printf("brute force %lu C 2 has %lu combinations\n",
                sequences.size(), nCk);
    }

    /* read pair file on all procs */
    mpix_read_file(comm, all_argv[2], pair_file_buffer, pair_file_size);
    double size_check = 1.0 * pair_file_size / (sizeof(int)*2+sizeof(char));
    mpix_print_sync(comm, "size_check", size_check);
    assert(floor(size_check) == size_check);

    unsigned long ntasks = (unsigned long)size_check;
    unsigned long global_num_workers = nprocs*NUM_WORKERS;
    unsigned long tasks_per_worker = ntasks / global_num_workers;
    unsigned long max_tasks_per_worker = ntasks / global_num_workers;
    max_tasks_per_worker += ntasks % global_num_workers;
    max_tasks_per_worker *= 10;
    if (0 == trank(0)) {
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
    for (int worker=0; worker<NUM_WORKERS; ++worker)
    {
#if DUMP_COSTS
        out[worker].open(get_filename(worker).c_str());
#endif
        edges[worker].open(get_edges_filename(worker).c_str());
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
        count[worker] = populate_tasks(ntasks, tasks_per_worker, worker,
                pair_file_buffer);
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
//#if defined(SET_AFFINITY)
    cpu_set_t cpuset;
    pthread_t thread;
    CPU_ZERO(&cpuset);
    thread = pthread_self();
    CPU_SET(0, &cpuset);
    pthread_setaffinity_np(thread, sizeof(cpu_set_t), &cpuset);
//#endif

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
    for(unsigned i=0; i<NUM_WORKERS; i++) {
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

    /* clean up */
    for (int worker=0; worker<NUM_WORKERS; ++worker) {
        free_tbl(tbl[worker], NROW);
        free_int(del[worker], NROW);
        free_int(ins[worker], NROW);
        delete utcs[worker];
#if DUMP_COSTS
        out[worker].close();
#endif
        edges[worker].close();
    }
    delete pgrp;

    TascelConfig::finalize();
    MPI_Comm_free(&comm);
    MPI_Finalize();

    return 0;
}
