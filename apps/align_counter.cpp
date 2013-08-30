/**
 * @author jeff.daily@pnnl.gov
 *
 * Alignment of an input dataset using work stealing. Each MPI task reads the
 * input file.
 */
#include "config.h"

#include <mpi.h>
extern "C" {
#include <armci.h>
}

#include <fcntl.h>
#include <sys/errno.h>
#include <sys/mman.h>
#include <sys/time.h>
#include <unistd.h>
#include <semaphore.h>

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
#include "constants.h"
#include "combinations.h"
#include "alignment.hpp"
#include "mpix.hpp"

using namespace std;
using namespace pgraph;

#define DEBUG 1
#define SEP ","
#define ALL_RESULTS 1
#define GB (1073741824)

typedef struct {
    const char *str;
    size_t size;
} sequence_t;

static string _get_edges_filename(int rank);
static unsigned long _mytime();
static string _generate_shm_name();
static void* _shm_create(MPI_Comm comm, size_t size, string &name);
static sem_t* _sem_create(MPI_Comm comm, string &sem_name);
static void _pack_and_index_fasta(char *buffer, size_t size,
        char delimiter, vector<sequence_t> &sequences, size_t max_seq_size);
static void _index_fasta(const char *buffer, size_t size,
        char delimiter, vector<sequence_t> &sequences, size_t max_seq_size);


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

#if 0
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
        sequences->align(seq_id[0], seq_id[1], result.score, result.matches, result.length, open, gap, thd);
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
                << ": (score,ndig,alen)=("
                << result.score << ","
                << result.ndig << ","
                << result.alen << ")"
                << ": edge? " << is_edge_answer << endl;
#endif
            edge_results[thd].push_back(EdgeResult(
                        seq_id[0], seq_id[1],
#if 0
                        1.0*result.alen/max_len,
                        1.0*result.ndig/result.alen,
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
        t = MPI_Wtime() - t;
        stats[thd].align_times_tot += t;
        stats[thd].calc_min(t);
        stats[thd].calc_max(t);
    }
}
#endif


#if 0
unsigned long populate_tasks(
        unsigned long ntasks, unsigned long tasks_per_worker, int worker)
{
    task_description desc;
    int wrank = trank(worker);
    unsigned long count = 0;
    unsigned long i;
    unsigned long lower_limit = wrank*tasks_per_worker;
    unsigned long upper_limit = lower_limit + tasks_per_worker;
    unsigned long remainder = ntasks % (size_world*NUM_WORKERS);

    unsigned long seq_id[2];
    k_combination(lower_limit, 2, seq_id);
    for (i=lower_limit; i<upper_limit; ++i) {
        desc.id1 = seq_id[0];
        desc.id2 = seq_id[1];
#if DEBUG
        cout << wrank << " added " << seq_id[0] << "," << seq_id[1] << endl;
#endif
        next_combination(2, seq_id);
        utcs[worker]->addTask(&desc, sizeof(desc));
        count++;
    }
    /* if I'm the last worker, add the remainder of the tasks */
    if (wrank == size_world*NUM_WORKERS-1) {
        for (/*ignore*/; i<upper_limit+remainder; ++i) {
            count++;
            desc.id1 = seq_id[0];
            desc.id2 = seq_id[1];
            next_combination(2, seq_id);
            utcs[worker]->addTask(&desc, sizeof(desc));
        }
    }

    return count;
}
#endif


int main(int argc, char **argv)
{
    MPI_Comm comm_world = MPI_COMM_NULL;
    MPI_Comm comm_node = MPI_COMM_NULL;
    MPI_Comm comm_master = MPI_COMM_NULL;
    int rank_world = 0;
    int rank_node = 0;
    int rank_master = 0;
    int size_world = 0;
    int size_node = 0;
    int size_master = 0;
    AlignStats stats;
    vector<EdgeResult> edge_results;
    vector<string> all_argv;
    MPI_Offset file_size = -1;
    char *file_buffer = NULL;
    unsigned long nCk = 0;
    size_t sequence_count = 0;
    int retval = 0;
    vector<sequence_t> sequences;
    size_t max_seq_size = 0;
    cell_t **tbl = NULL;
    int **ins = NULL;
    int **del = NULL;
    sem_t *sem = NULL;
    string sem_name;
    string shm_name;

    /* initialize MPI */
    MPI_CHECK_C(MPI_Init(&argc, &argv));
    MPI_CHECK_C(MPI_Comm_dup(MPI_COMM_WORLD, &comm_world));
    MPI_CHECK_C(MPI_Comm_rank(comm_world, &rank_world));
    MPI_CHECK_C(MPI_Comm_size(comm_world, &size_world));
    
    /* split comm_world based on hostid */
    MPI_CHECK_C(MPI_Comm_split(comm_world, gethostid(), rank_world, &comm_node));
    MPI_CHECK_C(MPI_Comm_rank(comm_node, &rank_node));
    MPI_CHECK_C(MPI_Comm_size(comm_node, &size_node));

    /* split comm_world based on node master */
    if (0 == rank_node) {
        MPI_CHECK_C(MPI_Comm_split(comm_world, 0, rank_world, &comm_master));
        MPI_CHECK_C(MPI_Comm_rank(comm_master, &rank_master));
        MPI_CHECK_C(MPI_Comm_size(comm_master, &size_master));
    }
    else {
        MPI_CHECK_C(MPI_Comm_split(comm_world, 1, rank_world, &comm_master));
        MPI_CHECK_C(MPI_Comm_free(&comm_master));
        comm_master = MPI_COMM_NULL;
        rank_master = -1;
        size_master = -1;
    }

    mpix_print_sync("rank_world", rank_world, comm_world);
    mpix_print_sync("rank_node",  rank_node, comm_world);
    mpix_print_sync("rank_master",rank_master, comm_world);
    mpix_print_sync("size_world", size_world, comm_world);
    mpix_print_sync("size_node",  size_node, comm_world);
    mpix_print_sync("size_master",size_master, comm_world);


    /* MPI standard does not guarantee all procs receive argc and argv */
    if (0 == rank_world) {
        MPI_CHECK_C(MPI_Bcast(&argc, 1, MPI_INT, 0, comm_world));
        for (int i=0; i<argc; ++i) {
            int length = strlen(argv[i])+1;
            MPI_CHECK_C(MPI_Bcast(&length, 1, MPI_INT, 0, comm_world));
            MPI_CHECK_C(MPI_Bcast(argv[i], length, MPI_CHAR, 0, comm_world));
            all_argv.push_back(argv[i]);
        }
    } else {
        int all_argc;
        MPI_CHECK_C(MPI_Bcast(&all_argc, 1, MPI_INT, 0, comm_world));
        for (int i=0; i<all_argc; ++i) {
            int length;
            char buffer[ARG_LEN_MAX];
            MPI_CHECK_C(MPI_Bcast(&length, 1, MPI_INT, 0, comm_world));
            MPI_CHECK_C(MPI_Bcast(buffer, length, MPI_CHAR, 0, comm_world));
            all_argv.push_back(buffer);
        }
    }

#if DEBUG
    /* print the command line arguments */
    for (int i=0; i<size_world; ++i) {
        if (i == rank_world) {
            for (size_t j=0; j<all_argv.size(); ++j) {
                printf("[%d] argv[%zd]=%s\n", rank_world, j, all_argv[j].c_str());
            }
        }
        MPI_Barrier(comm_world);
    }
#endif

    /* sanity check that we got the correct number of arguments */
    if (all_argv.size() <= 1 || all_argv.size() >= 3) {
        if (0 == rank_world) {
            if (all_argv.size() <= 1) {
                printf("missing input file\n");
            }
            else if (all_argv.size() >= 3) {
                printf("too many arguments\n");
            }
            printf("usage: align sequence_file\n");
        }
        MPI_Comm_free(&comm_world);
        MPI_Finalize();
        return 1;
    }

    sem = _sem_create(comm_node, sem_name);

    /* rank_world 0 on each node will open the file into shared memory */
    if (0 == rank_node) {
        file_size = mpix_get_file_size(all_argv[1], comm_master);
    }
    mpix_bcast(file_size, 0, comm_node);
    file_buffer = (char*)_shm_create(comm_node, file_size, shm_name);

    /* rank_world 0 on each node will open the file into shared memory */
    if (0 == rank_node) {
        mpix_read_file(all_argv[1], file_buffer, file_size, GB, comm_master);
        _pack_and_index_fasta(file_buffer, file_size, '\0', sequences, max_seq_size);
        // unlock all other ranks on this node
        for (int i=0; i<size_node-1; ++i) {
            retval = sem_post(sem);
            if (-1 == retval) {
                perror("sem_post");
                MPI_Abort(MPI_COMM_WORLD, -1);
            }
        }
    }
    else {
        retval = sem_wait(sem);
        if (-1 == retval) {
            perror("sem_wait");
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
        while (file_buffer[0] != '>') {
            cout << "waiting for shm file" << endl;
            sleep(1);
            cout << string(file_buffer, file_size) << endl;
        }
        _index_fasta(file_buffer, file_size, '\0', sequences, max_seq_size);
    }
    MPI_Barrier(comm_world);

    mpix_print_sync("max_seq_size", max_seq_size, comm_world);

    {
        size_t max_seq_size_p1 = max_seq_size + 1;
        tbl = new cell_t*[2];
        tbl[0] = new cell_t[max_seq_size_p1];
        tbl[1] = new cell_t[max_seq_size_p1];
        del = new int*[2];
        del[0] = new int[max_seq_size_p1];
        del[1] = new int[max_seq_size_p1];
        ins = new int*[2];
        ins[0] = new int[max_seq_size_p1];
        ins[1] = new int[max_seq_size_p1];
#if USE_MEMSET
        (void)memset(tbl[0], 0, sizeof(cell_t)*max_seq_size_p1);
        (void)memset(tbl[1], 0, sizeof(cell_t)*max_seq_size_p1);
        (void)memset(del[0], 0, sizeof(int)*max_seq_size_p1);
        (void)memset(del[1], 0, sizeof(int)*max_seq_size_p1);
        (void)memset(ins[0], 0, sizeof(int)*max_seq_size_p1);
        (void)memset(ins[1], 0, sizeof(int)*max_seq_size_p1);
#endif
    }

    /* how many combinations of sequences are there? */
    nCk = binomial_coefficient(sequence_count, 2);
    if (0 == rank_world) {
        printf("brute force %lu C 2 has %lu combinations\n",
                sequence_count, nCk);
    }

    double selectivity = 1.0;
    if (all_argv.size() == 3) {
      selectivity  = fabs(min(1.0,atof(all_argv[2].c_str())));
    }
    unsigned long nalignments = (long)(0.5+selectivity*nCk);
    if (0 == rank_world) {
        printf("selectivity=%lf\n", selectivity);
        printf("nalignments=%lu\n", nalignments);
    }
    MPI_Barrier(comm_world);

    double mytimer = MPI_Wtime();
    // TODO TAKS COUNTER AND ALIGNMENT LOOP
    mytimer = MPI_Wtime() - mytimer;
    if (0 == rank_world) {
        cout << "mytimer=" << mytimer << endl;
    }

    MPI_Barrier(comm_world);

    AlignStats * rstats = new AlignStats[size_world];
    MPI_Gather(&stats, sizeof(AlignStats), MPI_CHAR, 
	       rstats, sizeof(AlignStats), MPI_CHAR, 
	       0, comm_world);

    /* synchronously print alignment stats all from process 0 */
    if (0 == rank_world) {
        AlignStats totals;
        cout << " pid" << rstats[0].getHeader() << endl;      
        for(int i=0; i<size_world; i++) {
            totals += rstats[i];
            cout << std::setw(4) << std::right << i << rstats[i] << endl;
        }
        cout << "==============================================" << endl;
        cout << setw(4) << right << "TOT" << totals << endl;
    }
    delete [] rstats;
    rstats=NULL;

    MPI_Barrier(comm_world);

    /* unmap the memory */
    retval = munmap(file_buffer, file_size);
    if (-1 == retval) {
        perror("munmap");
        //MPI_Abort(MPI_COMM_WORLD, retval);
    }

    /* remove the shared memory object */
    if (0 == rank_world) {
        retval = shm_unlink(shm_name.c_str());
        if (-1 == retval) {
            perror("shm_unlink");
            //MPI_Abort(MPI_COMM_WORLD, retval);
        }
    }

    retval = sem_close(sem);
    if (-1 == retval) {
        perror("sem_close");
    }
    if (0 == rank_node) {
        retval = sem_unlink(sem_name.c_str());
        if (-1 == retval) {
            perror("sem_unlink");
        }
    }

    if (MPI_COMM_NULL != comm_master) {
        MPI_Comm_free(&comm_master);
    }
    MPI_Comm_free(&comm_node);
    MPI_Comm_free(&comm_world);
    MPI_Finalize();

    return 0;
}


static string _get_edges_filename(int rank)
{
    ostringstream str;
    str << "edges." << rank << ".txt";
    return str.str();
}


#if 0
static size_t _parse_memory_budget(const string& value)
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
#endif


static void* _shm_create(MPI_Comm comm, size_t size, string &name)
{
    void *mapped = NULL;
    int fd = 0;
    int retval = 0;
    int rank = 0;

    MPI_CHECK_C(MPI_Comm_rank(comm, &rank));

    if (0 == rank) {
        name = _generate_shm_name();

        /* create shared memory segment */
        fd = shm_open(name.c_str(), O_CREAT|O_EXCL|O_RDWR, S_IRUSR|S_IWUSR);
        if (-1 == fd) {
            perror("_shm_create: shm_open");
            MPI_Abort(MPI_COMM_WORLD, -1);
        }

        /* set the size of my shared memory object */
        retval = ftruncate(fd, size);
        if (-1 == retval) {
            perror("_shm_create: ftruncate");
            MPI_Abort(MPI_COMM_WORLD, -1);
        }

        /* map into local address space */
        mapped = mmap(NULL, size, PROT_READ|PROT_WRITE, MAP_SHARED, fd, 0);
        if (MAP_FAILED == mapped) {
            perror("_shm_map: mmap");
            MPI_Abort(MPI_COMM_WORLD, -1);
        }

        /* close file descriptor */
        retval = close(fd);
        if (-1 == retval) {
            perror("_shm_create: close");
            MPI_Abort(MPI_COMM_WORLD, retval);
        }

        /* broadcast the name */
        mpix_bcast(name, 0, comm);
    }
    else {
        /* broadcast the name */
        mpix_bcast(name, 0, comm);

        /* attach to shared memory segment */
        fd = shm_open(name.c_str(), O_RDWR, S_IRUSR|S_IWUSR);
        if (-1 == fd) {
            perror("_shm_create (attach): shm_open");
            MPI_Abort(MPI_COMM_WORLD, -1);
        }

        /* map into local address space */
        mapped = mmap(NULL, size, PROT_READ|PROT_WRITE, MAP_SHARED, fd, 0);
        if (MAP_FAILED == mapped) {
            perror("_shm_map: mmap");
            MPI_Abort(MPI_COMM_WORLD, -1);
        }

        /* close file descriptor */
        retval = close(fd);
        if (-1 == retval) {
            perror("_shm_create (attach): close");
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
    }

    return mapped;
}


static sem_t* _sem_create(MPI_Comm comm, string &name)
{
    const char *cname = NULL;
    sem_t *sem = NULL;
    int retval = 0;
    int rank = 0;

    MPI_CHECK_C(MPI_Comm_rank(comm, &rank));

    if (0 == rank) {
        name = _generate_shm_name();
        cname = name.c_str();

        /* create shared semaphore */
        printf("_sem_create(%s, O_CREAT|O_EXCL, S_IRUSR|S_IWUSR, 0)\n", cname);
        sem = sem_open(cname, O_CREAT|O_EXCL, S_IRUSR|S_IWUSR, 0);
        if (SEM_FAILED == sem) {
            perror("_sem_create: sem_open");
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
        printf("_sem_create(%s, O_CREAT|O_EXCL, S_IRUSR|S_IWUSR, 0) OKAY\n", cname);

        /* broadcast the name */
        mpix_bcast(name, 0, comm);
    }
    else {
        /* broadcast the name */
        mpix_bcast(name, 0, comm);
        cname = name.c_str();

        /* open shared semaphore */
        printf("_sem_create(%s, 0)\n", cname);
        sem = sem_open(cname, 0);
        if (SEM_FAILED == sem) {
            perror("_sem_create: sem_open");
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
        printf("_sem_create(%s, 0) OKAY\n", cname);
    }

    return sem;
}


static unsigned long _mytime()
{
    struct timeval tv;
    int retval = 0;

    retval = gettimeofday(&tv, NULL);
    assert(0 == retval);

    return ((unsigned long)tv.tv_sec)*1000000UL + ((unsigned long)tv.tv_usec);
}


static string _generate_shm_name()
{
    ostringstream os;
    os << "/pgraph" << _mytime() << getpid();
    return os.str();
}


static void _pack_and_index_fasta(char *buffer, size_t size,
        char delimiter, vector<sequence_t> &sequences, size_t max_seq_size)
{
    size_t r = 0;           /* read index */
    size_t w = 0;           /* write index */
    size_t last_gt = 0;     /* index of last '>' seen */
    size_t last_hash = 0;   /* index of last '#' inserted */
    sequence_t sequence;

    /* We skip all newlines unless it is the newline that terminates the
     * sequence ID. We replace the ID-terminating newline with '#'. We replace
     * the sequence terminating newline with the given delimiter. */

    assert(buffer[0] == '>');

    while (r<size) {
        if (buffer[r] == '>') {
            last_gt = w;
            while (r<size && buffer[r] != '\n') {
                buffer[w++] = buffer[r++];
            }
            last_hash = w;
            buffer[w++] = '#';
            r++;
        }
        else {
            while (r<size && buffer[r] != '\n') {
                buffer[w++] = buffer[r++];
            }
            /* peek at next character, either EOF or '>' */
            if (r<size) {
                if (r+1 >= size || buffer[r+1] == '>') {
                    buffer[w++] = delimiter;
                    size_t sequence_offset = last_hash - last_gt + 1;
                    size_t sequence_length = w - last_hash - 1;
                    if (delimiter == '\0') {
                        sequence_length -= 1;
                    }
                    sequence.str = &buffer[sequence_offset];
                    sequence.size = sequence_length;
                    sequences.push_back(sequence);
                    max_seq_size = sequence.size > max_seq_size ?
                        sequence.size : max_seq_size;
                }
            }
            /* or file wasn't terminated with a newline */
            else {
                /* hopefully space enough to append the delimiter anyway */
                assert(w<r);
                {
                    buffer[w++] = delimiter;
                    size_t sequence_offset = last_hash - last_gt + 1;
                    size_t sequence_length = w - last_hash - 1;
                    if (delimiter == '\0') {
                        sequence_length -= 1;
                    }
                    sequence.str = &buffer[sequence_offset];
                    sequence.size = sequence_length;
                    sequences.push_back(sequence);
                    max_seq_size = sequence.size > max_seq_size ?
                        sequence.size : max_seq_size;
                }
                assert(0);
            }
            r++;
        }
    }

#if 0
    tbl = new cell_t**[num_threads];
    del = new int**[num_threads];
    ins = new int**[num_threads];
    for (size_t worker=0; worker<num_threads; ++worker) {
#if 1
        size_t max_seq_size_p1 = max_seq_size + 1;
        tbl[worker] = new cell_t*[2];
        tbl[worker][0] = new cell_t[max_seq_size_p1];
        tbl[worker][1] = new cell_t[max_seq_size_p1];
        del[worker] = new int*[2];
        del[worker][0] = new int[max_seq_size_p1];
        del[worker][1] = new int[max_seq_size_p1];
        ins[worker] = new int*[2];
        ins[worker][0] = new int[max_seq_size_p1];
        ins[worker][1] = new int[max_seq_size_p1];
#if USE_MEMSET
        (void)memset(tbl[worker][0], 0, sizeof(cell_t)*max_seq_size_p1);
        (void)memset(tbl[worker][1], 0, sizeof(cell_t)*max_seq_size_p1);
        (void)memset(del[worker][0], 0, sizeof(int)*max_seq_size_p1);
        (void)memset(del[worker][1], 0, sizeof(int)*max_seq_size_p1);
        (void)memset(ins[worker][0], 0, sizeof(int)*max_seq_size_p1);
        (void)memset(ins[worker][1], 0, sizeof(int)*max_seq_size_p1);
#endif
#else
        tbl[worker] = allocate_cell_table(2, max_seq_size);
        del[worker] = allocate_int_table(2, max_seq_size);
        ins[worker] = allocate_int_table(2, max_seq_size);
#endif
    }
#endif
}


static void _index_fasta(const char *buffer, size_t size,
        char delimiter, vector<sequence_t> &sequences, size_t max_seq_size)
{
    size_t r = 0;           /* read index */
    size_t w = 0;           /* write index */
    size_t last_gt = 0;     /* index of last '>' seen */
    size_t last_hash = 0;   /* index of last '#' inserted */
    sequence_t sequence;

    /* We skip all newlines unless it is the newline that terminates the
     * sequence ID. We replace the ID-terminating newline with '#'. We replace
     * the sequence terminating newline with the given delimiter. */

    assert(buffer[0] == '>');

    while (r<size) {
        if (buffer[r] == '>') {
            last_gt = w;
            while (r<size && buffer[r] != '\n') {
                w++;
                r++;
            }
            last_hash = w;
            w++;
            r++;
        }
        else {
            while (r<size && buffer[r] != '\n') {
                w++;
                r++;
            }
            /* peek at next character, either EOF or '>' */
            if (r<size) {
                if (r+1 >= size || buffer[r+1] == '>') {
                    w++;
                    size_t sequence_offset = last_hash - last_gt + 1;
                    size_t sequence_length = w - last_hash - 1;
                    if (delimiter == '\0') {
                        sequence_length -= 1;
                    }
                    sequence.str = &buffer[sequence_offset];
                    sequence.size = sequence_length;
                    sequences.push_back(sequence);
                    max_seq_size = sequence.size > max_seq_size ?
                        sequence.size : max_seq_size;
                }
            }
            /* or file wasn't terminated with a newline */
            else {
                /* hopefully space enough to append the delimiter anyway */
                assert(w<r);
                {
                    w++;
                    size_t sequence_offset = last_hash - last_gt + 1;
                    size_t sequence_length = w - last_hash - 1;
                    if (delimiter == '\0') {
                        sequence_length -= 1;
                    }
                    sequence.str = &buffer[sequence_offset];
                    sequence.size = sequence_length;
                    sequences.push_back(sequence);
                    max_seq_size = sequence.size > max_seq_size ?
                        sequence.size : max_seq_size;
                }
                assert(0);
            }
            r++;
        }
    }
#if 0
    tbl = new cell_t**[num_threads];
    del = new int**[num_threads];
    ins = new int**[num_threads];
    for (size_t worker=0; worker<num_threads; ++worker) {
#if 1
        size_t max_seq_size_p1 = max_seq_size + 1;
        tbl[worker] = new cell_t*[2];
        tbl[worker][0] = new cell_t[max_seq_size_p1];
        tbl[worker][1] = new cell_t[max_seq_size_p1];
        del[worker] = new int*[2];
        del[worker][0] = new int[max_seq_size_p1];
        del[worker][1] = new int[max_seq_size_p1];
        ins[worker] = new int*[2];
        ins[worker][0] = new int[max_seq_size_p1];
        ins[worker][1] = new int[max_seq_size_p1];
#if USE_MEMSET
        (void)memset(tbl[worker][0], 0, sizeof(cell_t)*max_seq_size_p1);
        (void)memset(tbl[worker][1], 0, sizeof(cell_t)*max_seq_size_p1);
        (void)memset(del[worker][0], 0, sizeof(int)*max_seq_size_p1);
        (void)memset(del[worker][1], 0, sizeof(int)*max_seq_size_p1);
        (void)memset(ins[worker][0], 0, sizeof(int)*max_seq_size_p1);
        (void)memset(ins[worker][1], 0, sizeof(int)*max_seq_size_p1);
#endif
#else
        tbl[worker] = allocate_cell_table(2, max_seq_size);
        del[worker] = allocate_int_table(2, max_seq_size);
        ins[worker] = allocate_int_table(2, max_seq_size);
#endif
    }
#endif
}

