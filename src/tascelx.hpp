/**
 * @file tascelx.hpp
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 *
 * Helper functions for routine tascel boilerplate code.
 */
#ifndef _TASCEL_X_H_
#define _TASCEL_X_H_

#include <pthread.h>

#include <tascel.h>
#include <tascel/UniformTaskCollection.h>

#include "Bootstrap.hpp"

#if defined(__xlc__) || defined(__xlC__) || defined(__IBMC__) || defined(__IBMCPP__)
#define MFENCE __asm__ __volatile__  ("sync" ::: "memory");
#else
#define MFENCE asm("mfence");
#endif

#if !HAVE_PTHREAD_BARRIER_T && defined(THREADED)
#include "pthread_fixes.h"
#endif

using namespace tascel;
/* For now, the pgraph namespace is only used in pthread_fixes.h header.
 * Once we properly namespace all of pgraph, this #if can be removed. */

static pthread_t *threadHandles = 0;
static unsigned *threadRanks = 0;
// Synchronization for worker threads
pthread_barrier_t workersStart;
pthread_barrier_t workersEnd;
// Synchronization for server thread
pthread_barrier_t serverStart;
pthread_barrier_t serverEnd;
volatile bool serverEnabled = true;

/**
 * Sets the current thread's affinity using pthread_setaffinity_np.
 *
 * @param[in] the core to bind pthread_self() to.
 */
static void set_affinity(const unsigned int rank)
{
#if defined(SET_AFFINITY)
    cpu_set_t cpuset;
    pthread_t thread;
    CPU_ZERO(&cpuset);
    thread = pthread_self();
    CPU_SET(rank, &cpuset);
    (void)pthread_setaffinity_np(thread, sizeof(cpu_set_t), &cpuset);
    //(void)sched_setaffinity(0, sizeof(cpu_set_t), &cpuset);
#endif
}


/**
 * The arguments passed to the server_thread when it is created.
 */
typedef struct server_thread_args {
    volatile bool *server_enabled;
    pthread_barrier_t *server_start;
    pthread_barrier_t *server_end;

    server_thread_args(
        volatile bool *server_enabled,
        pthread_barrier_t *server_start,
        pthread_barrier_t *server_end)
        :   server_enabled(server_enabled)
        ,   server_start(server_start)
        ,   server_end(server_end)
    { }
} server_thread_args;


/**
 * The server thread.
 *
 * Handles all active messages, data transfers, etc.
 *
 * @param[in] args a server_thread_args instance
 */
static void *server_thread(void *_args)
{
    server_thread_args *args = (server_thread_args *)_args;

    set_affinity(NUM_WORKERS);

    // When enabled execute any active messages that arrive
    pthread_barrier_wait(args->server_start);
    while (*(args->server_enabled) || !serverDispatcher.empty()) {
        AmListenObjCodelet<NullMutex> *lcodelet;
        if ((lcodelet = theAm().amListeners[0]->progress()) != NULL) {
            lcodelet->execute();
        }
        Codelet *codelet;
        if ((codelet = serverDispatcher.progress()) != NULL) {
            codelet->execute();
            // Assume that codelet was an AmRequest and needs to be freed
            delete reinterpret_cast<AmRequest *>(codelet);
        }
    }
    pthread_barrier_wait(args->server_end);

    delete args;

    return NULL;
}


/**
 * The arguments passed to the worker_thread when it is created.
 */
typedef struct worker_thread_args {
    unsigned int rank;
    UniformTaskCollection *utc;
    pthread_barrier_t *workers_start;
    pthread_barrier_t *workers_end;

    worker_thread_args(
        unsigned int rank,
        UniformTaskCollection *utc,
        pthread_barrier_t *workers_start,
        pthread_barrier_t *workers_end)
        :   rank(rank)
        ,   utc(utc)
        ,   workers_start(workers_start)
        ,   workers_end(workers_end)
    {}

} worker_thread_args;


/**
 * The worker thread.
 *
 * Simply begins processing the task collection.
 *
 * @param[in] args a worker_thread_args instance
 */
static void *worker_thread(void *_args)
{
    worker_thread_args *args = (worker_thread_args *)_args;

    set_affinity(args->rank);

    pthread_barrier_wait(args->workers_start);
    args->utc->process();
    pthread_barrier_wait(args->workers_end);

    delete args;

    return NULL;
}


/** barrier */
static void amBarrierThd()
{
    int epoch = pgrp->signalBarrier();
    while (!pgrp->testBarrier(epoch)) { }
}


static void allocate_threads()
{
    threadHandles = new pthread_t[NUM_WORKERS + NUM_SERVERS];
    threadRanks = new unsigned[NUM_WORKERS + NUM_SERVERS];
    for (int worker=0; worker<NUM_WORKERS; ++worker) {
        threadRanks[worker] = worker;
    }
    pthread_barrier_init(&workersStart, 0, NUM_WORKERS);
    pthread_barrier_init(&workersEnd, 0, NUM_WORKERS);
    pthread_barrier_init(&serverStart, 0, NUM_SERVERS + 1);
    pthread_barrier_init(&serverEnd, 0, NUM_SERVERS + 1);
    MFENCE
}


static void initialize_threads(UniformTaskCollection **utcs)
{
    set_affinity(0);
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
    MPI_Barrier(pgraph::comm);
    pthread_barrier_wait(&workersStart);
}


static void finalize_threads()
{
    pthread_barrier_wait(&workersEnd);
    amBarrierThd();
    for (int i = 1; i < NUM_WORKERS; ++i) {
        pthread_join(threadHandles[i], NULL);
    }

    serverEnabled = false;
    MFENCE
    pthread_barrier_wait(&serverEnd);
    amBarrier();
    pthread_join(threadHandles[NUM_WORKERS], NULL);

    delete [] threadHandles;
    delete [] threadRanks;
}


#endif /* _TASCEL_X_H_ */
