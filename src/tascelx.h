/**
 * @file tascelx.h
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
#include <tascel/UniformTaskCollSplitHybrid.h>

using namespace tascel;
/* For now, the pgraph namespace is only used in pthread_fixes.h header.
 * Once we properly namespace all of pgraph, this #if can be removed. */
#if !HAVE_PTHREAD_BARRIER_T
using namespace pgraph;
#endif

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
    int ret = pthread_setaffinity_np(thread, sizeof(cpu_set_t), &cpuset);
    //int ret = sched_setaffinity(0, sizeof(cpu_set_t), &cpuset);
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
    server_thread_args *args = (server_thread_args*)_args;

    set_affinity(NUM_WORKERS);

    // When enabled execute any active messages that arrive
    while (1) {
        pthread_barrier_wait(args->server_start);
        while (*(args->server_enabled)) {
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
        pthread_barrier_wait(args->server_end);
    }

    delete args;

    return NULL;
}


/**
 * The arguments passed to the worker_thread when it is created.
 */
typedef struct worker_thread_args {
    unsigned int rank;
    UniformTaskCollSplitHybrid *utcs;
    pthread_barrier_t *workers_start;
    pthread_barrier_t *workers_end;
    
    worker_thread_args(
            unsigned int rank,
            UniformTaskCollSplitHybrid *utcs,
            pthread_barrier_t *workers_start,
            pthread_barrier_t *workers_end)
        :   rank(rank)
        ,   utcs(utcs)
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
    worker_thread_args *args = (worker_thread_args*)_args;

    set_affinity(args->rank);

    while (1) {
        pthread_barrier_wait(args->workers_start);
        args->utcs->process(args->rank);
        pthread_barrier_wait(args->workers_end);
    }

    delete args;

    return NULL;
}


/** active message progress barrier */
static void amBarrier()
{
    int epoch = pgrp->signalBarrier();
    while(!pgrp->testBarrier(epoch)) {
        AmListenObjCodelet<NullMutex>* codelet;
        if((codelet=theAm().amListeners[0]->progress()) != NULL) {
            codelet->execute();
        }
    }
}


/** barrier */
static void amBarrierThd()
{
    int epoch = pgrp->signalBarrier();
    while(!pgrp->testBarrier(epoch)) { }
}

#endif /* _TASCEL_X_H_ */
