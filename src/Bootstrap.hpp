#ifndef _PGRAPH_BOOTSTRAP_H_
#define _PGRAPH_BOOTSTRAP_H_

#include <mpi.h>

namespace pgraph {
    extern MPI_Comm comm;
    extern int rank;
    extern int nprocs;
    extern int *argc_;
    extern char ***argv_;
    extern bool is_armci_initialized;
    void initialize(int *argc, char ***argv);
    void finalize();

    void initialize_armci();
    void finalize_armci();
};

#endif /* _PGRAPH_BOOTSTRAP_H_ */
