#include "config.h"

#include <unistd.h>

#include <cassert>

#include <mpi.h>

#include "Bootstrap.hpp"
#include "mpix.hpp"

#define PAUSE_ON_ERROR 0
#ifdef PAUSE_ON_ERROR
#include "SigSegvHandler.hpp"
#endif

#if HAVE_ARMCI
/* ANL's armci does not extern "C" inside the header */
extern "C" {
#include <armci.h>
}
#endif

namespace pgraph {

MPI_Comm comm = MPI_COMM_NULL;
int comm_rank = 0;
int comm_size = 0;
int *argc_ = NULL;
char ***argv_ = NULL;
bool is_armci_initialized = false;

void initialize(int &argc, char **&argv)
{
#ifdef PAUSE_ON_ERROR
    TrapSigSegv();
#endif

    argc_ = &argc;
    argv_ = &argv;
    is_armci_initialized = false;

    /* initialize MPI */
#if defined(THREADED)
    mpix::init_thread(argc, argv, MPI_THREAD_MULTIPLE);
#else
    mpix::init(argc, argv);
#endif

    /* get comm rank and size */
    comm = mpix::comm_dup(MPI_COMM_WORLD);
    comm_rank = mpix::comm_rank(comm);
    comm_size = mpix::comm_size(comm);
    ::std::cout << "Hello from " << comm_rank << " of " << comm_size << endl;
}


void finalize()
{
    finalize_armci();
    mpix::comm_free(comm);
    mpix::finalize();
}


void initialize_armci()
{
#if HAVE_ARMCI
    if (!is_armci_initialized) {
        ARMCI_Init_args(argc_,argv_);
        is_armci_initialized = true;
    }
#endif
}


void finalize_armci()
{
#if HAVE_ARMCI
    if (is_armci_initialized) {
        ARMCI_Finalize();
    }
#endif
}

};

