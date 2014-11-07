/**
 * @file mpix_types.cpp
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 *
 * MPI_Datatype creation functions for pgraph structures.
 */
#include <mpi.h>

#include <tascel/Counter.h>
#include <tascel/StealingStats.h>
#include <tascel/Timer.h>

#include "AlignStats.hpp"
#include "DbStats.hpp"
#include "DupStats.hpp"
#include "Stats.hpp"
#include "Suffix.hpp"
#include "TreeStats.hpp"

#include "mpix.hpp"
#include "mpix_helper.hpp"

using ::pgraph::AlignStats;
using ::pgraph::DbStats;
using ::pgraph::DupStats;
using ::pgraph::Stats;
using ::pgraph::Suffix;
using ::pgraph::TreeStats;

using ::tascel::Counter;
using ::tascel::StealingStats;
using ::tascel::Timer;

namespace mpix {

MPI_Datatype mpi_datatype_Counter;
MPI_Datatype mpi_datatype_Timer;
MPI_Datatype mpi_datatype_StealingStats;
MPI_Datatype mpi_datatype_Stats;
MPI_Datatype mpi_datatype_TreeStats;
MPI_Datatype mpi_datatype_DbStats;
MPI_Datatype mpi_datatype_DupStats;
MPI_Datatype mpi_datatype_AlignStats;
MPI_Datatype mpi_datatype_Suffix;

MPI_Datatype get_mpi_datatype(Counter object)       { return mpi_datatype_Counter; }
MPI_Datatype get_mpi_datatype(Timer object)         { return mpi_datatype_Timer; }
MPI_Datatype get_mpi_datatype(StealingStats object) { return mpi_datatype_StealingStats; }
MPI_Datatype get_mpi_datatype(Stats object)         { return mpi_datatype_Stats; }
MPI_Datatype get_mpi_datatype(TreeStats object)     { return mpi_datatype_TreeStats; }
MPI_Datatype get_mpi_datatype(DbStats object)       { return mpi_datatype_DbStats; }
MPI_Datatype get_mpi_datatype(DupStats object)      { return mpi_datatype_DupStats; }
MPI_Datatype get_mpi_datatype(AlignStats object)    { return mpi_datatype_AlignStats; }
MPI_Datatype get_mpi_datatype(Suffix object)        { return mpi_datatype_Suffix; }


static void build_mpi_datatype_Counter()
{
#if 1
    type_contiguous(1, MPI_LONG, mpi_datatype_Counter);
#else
    Counter object;
    MPI_Datatype type[1] = {
        get_mpi_datatype(object.num)
    };
    int blocklen[1] = {1};
    MPI_Aint disp[1] = {
        MPI_Aint(&object.num) - MPI_Aint(&object)
    };
    type_create_struct(1, blocklen, disp, type, mpi_datatype_Counter);
#endif
    type_commit(mpi_datatype_Counter);
}


static void build_mpi_datatype_Timer()
{
#if 1
    type_contiguous(2, MPI_DOUBLE, mpi_datatype_Timer);
#else
    Timer object;
    MPI_Datatype type[2] = {
        get_mpi_datatype(MPI_DOUBLE),
        get_mpi_datatype(MPI_DOUBLE)
    };
    int blocklen[2] = {1,1};
    MPI_Aint disp[2] = {
        MPI_Aint(&object.ttime) - MPI_Aint(&object),
        MPI_Aint(&object.stime) - MPI_Aint(&object)
    };
    type_create_struct(2, blocklen, disp, type, mpi_datatype_Timer);
#endif
    type_commit(mpi_datatype_Timer);
}


static void build_mpi_datatype_StealingStats()
{
    StealingStats object;
    MPI_Datatype type[9] = {
        get_mpi_datatype(object.numStealAttempts),
        get_mpi_datatype(object.numSteals),
        get_mpi_datatype(object.numTasks),
        get_mpi_datatype(object.numTasks2),
        get_mpi_datatype(object.stealTime),
        get_mpi_datatype(object.taskTime),
        get_mpi_datatype(object.taskTime2),
        get_mpi_datatype(object.tdTime),
        get_mpi_datatype(object.miscTime)
    };
    int blocklen[9] = {1,1,1,1,1,1,1,1,1};
    MPI_Aint disp[9] = {
        MPI_Aint(&object.numStealAttempts)   - MPI_Aint(&object),
        MPI_Aint(&object.numSteals)          - MPI_Aint(&object),
        MPI_Aint(&object.numTasks)           - MPI_Aint(&object),
        MPI_Aint(&object.numTasks2)          - MPI_Aint(&object),
        MPI_Aint(&object.stealTime)          - MPI_Aint(&object),
        MPI_Aint(&object.taskTime)           - MPI_Aint(&object),
        MPI_Aint(&object.taskTime2)          - MPI_Aint(&object),
        MPI_Aint(&object.tdTime)             - MPI_Aint(&object),
        MPI_Aint(&object.miscTime)           - MPI_Aint(&object)
    };
    type_create_struct(9, blocklen, disp, type, mpi_datatype_StealingStats);
    type_commit(mpi_datatype_StealingStats);
}


static void build_mpi_datatype_Stats()
{
    Stats object;
    MPI_Datatype type[6] = {
        get_mpi_datatype(object._n),
        get_mpi_datatype(object._mean),
        get_mpi_datatype(object._M2),
        get_mpi_datatype(object._sum),
        get_mpi_datatype(object._min),
        get_mpi_datatype(object._max)
    };
    int blocklen[6] = {1,1,1,1,1,1};
    MPI_Aint disp[6] = {
        MPI_Aint(&object._n)    - MPI_Aint(&object),
        MPI_Aint(&object._mean) - MPI_Aint(&object),
        MPI_Aint(&object._M2)   - MPI_Aint(&object),
        MPI_Aint(&object._sum)  - MPI_Aint(&object),
        MPI_Aint(&object._min)  - MPI_Aint(&object),
        MPI_Aint(&object._max)  - MPI_Aint(&object)
    };
    type_create_struct(6, blocklen, disp, type, mpi_datatype_Stats);
    type_commit(mpi_datatype_Stats);
}


static void build_mpi_datatype_TreeStats()
{
    TreeStats object;
    MPI_Datatype type[12] = {
        get_mpi_datatype(object.trees),
        get_mpi_datatype(object.suffixes),
        get_mpi_datatype(object.size),
        get_mpi_datatype(object.size_internal),
        get_mpi_datatype(object.fanout),
        get_mpi_datatype(object.depth),
        get_mpi_datatype(object.suffix_length),
        get_mpi_datatype(object.pairs),
        get_mpi_datatype(object.time_build),
        get_mpi_datatype(object.time_process),
        get_mpi_datatype(object.time_first),
        get_mpi_datatype(object.time_last)
    };
    int blocklen[12] = {1,1,1,1,1,1,1,1,1,1,1,1};
    MPI_Aint disp[12] = {
        MPI_Aint(&object.trees)         - MPI_Aint(&object),
        MPI_Aint(&object.suffixes)      - MPI_Aint(&object),
        MPI_Aint(&object.size)          - MPI_Aint(&object),
        MPI_Aint(&object.size_internal) - MPI_Aint(&object),
        MPI_Aint(&object.fanout)        - MPI_Aint(&object),
        MPI_Aint(&object.depth)         - MPI_Aint(&object),
        MPI_Aint(&object.suffix_length) - MPI_Aint(&object),
        MPI_Aint(&object.pairs)         - MPI_Aint(&object),
        MPI_Aint(&object.time_build)    - MPI_Aint(&object),
        MPI_Aint(&object.time_process)  - MPI_Aint(&object),
        MPI_Aint(&object.time_first)    - MPI_Aint(&object),
        MPI_Aint(&object.time_last)     - MPI_Aint(&object)
    };
    type_create_struct(12, blocklen, disp, type, mpi_datatype_TreeStats);
    type_commit(mpi_datatype_TreeStats);
}


static void build_mpi_datatype_DbStats()
{
    DbStats object;
#if 0
    mpi_datatype_DbStats = type_contiguous(3, get_mpi_datatype(object.time));
#else
    MPI_Datatype type[3] = {
        get_mpi_datatype(object.time),
        get_mpi_datatype(object.bytes),
        get_mpi_datatype(object.cum)
    };
    int blocklen[3] = {1,1,1};
    MPI_Aint disp[3] = {
        MPI_Aint(&object.time)  - MPI_Aint(&object),
        MPI_Aint(&object.bytes) - MPI_Aint(&object),
        MPI_Aint(&object.cum)   - MPI_Aint(&object)
    };
    type_create_struct(3, blocklen, disp, type, mpi_datatype_DbStats);
#endif
    type_commit(mpi_datatype_DbStats);
}


static void build_mpi_datatype_DupStats()
{
    DupStats object;
#if 0
    mpi_datatype_DupStats = type_contiguous(3, get_mpi_datatype(object.time));
#else
    MPI_Datatype type[3] = {
        get_mpi_datatype(object.time),
        get_mpi_datatype(object.checked),
        get_mpi_datatype(object.returned)
    };
    int blocklen[3] = {1,1,1};
    MPI_Aint disp[3] = {
        MPI_Aint(&object.time)     - MPI_Aint(&object),
        MPI_Aint(&object.checked)  - MPI_Aint(&object),
        MPI_Aint(&object.returned) - MPI_Aint(&object)
    };
    type_create_struct(3, blocklen, disp, type, mpi_datatype_DupStats);
#endif
    type_commit(mpi_datatype_DupStats);
}


static void build_mpi_datatype_AlignStats()
{
    AlignStats object;
    MPI_Datatype type[8] = {
        get_mpi_datatype(object.edge_counts),
        get_mpi_datatype(object.align_counts),
        get_mpi_datatype(object.align_skipped),
        get_mpi_datatype(object.time_align),
        get_mpi_datatype(object.time_kcomb),
        get_mpi_datatype(object.time_total),
        get_mpi_datatype(object.work),
        get_mpi_datatype(object.work_skipped)
    };
    int blocklen[8] = {1,1,1,1,1,1,1,1};
    MPI_Aint disp[8] = {
        MPI_Aint(&object.edge_counts)   - MPI_Aint(&object),
        MPI_Aint(&object.align_counts)  - MPI_Aint(&object),
        MPI_Aint(&object.align_skipped) - MPI_Aint(&object),
        MPI_Aint(&object.time_align)    - MPI_Aint(&object),
        MPI_Aint(&object.time_kcomb)    - MPI_Aint(&object),
        MPI_Aint(&object.time_total)    - MPI_Aint(&object),
        MPI_Aint(&object.work)          - MPI_Aint(&object),
        MPI_Aint(&object.work_skipped)  - MPI_Aint(&object)
    };
    type_create_struct(8, blocklen, disp, type, mpi_datatype_AlignStats);
    type_commit(mpi_datatype_AlignStats);
}


static void build_mpi_datatype_Suffix()
{
    Suffix object;
    MPI_Datatype type[5] = {
        get_mpi_datatype(object.sid),
        get_mpi_datatype(object.pid),
        get_mpi_datatype(object.bid),
        get_mpi_datatype(object.k),
        MPI_UNSIGNED_LONG /* void* */
    };
    int blocklen[5] = {1,1,1,1,1};
    MPI_Aint disp[5] = {
        MPI_Aint(&object.sid)  - MPI_Aint(&object),
        MPI_Aint(&object.pid)  - MPI_Aint(&object),
        MPI_Aint(&object.bid)  - MPI_Aint(&object),
        MPI_Aint(&object.k)    - MPI_Aint(&object),
        MPI_Aint(&object.next) - MPI_Aint(&object)
    };
    type_create_struct(5, blocklen, disp, type, mpi_datatype_Suffix);
    type_commit(mpi_datatype_Suffix);
}


void init_types()
{
    build_mpi_datatype_Counter();
    build_mpi_datatype_Timer();
    build_mpi_datatype_StealingStats();
    build_mpi_datatype_Stats();
    build_mpi_datatype_TreeStats();
    build_mpi_datatype_DbStats();
    build_mpi_datatype_DupStats();
    build_mpi_datatype_AlignStats();
    build_mpi_datatype_Suffix();
}


void free_types()
{
    type_free(mpi_datatype_Counter);
    type_free(mpi_datatype_Timer);
    type_free(mpi_datatype_StealingStats);
    type_free(mpi_datatype_Stats);
    type_free(mpi_datatype_TreeStats);
    type_free(mpi_datatype_DbStats);
    type_free(mpi_datatype_DupStats);
    type_free(mpi_datatype_AlignStats);
    type_free(mpi_datatype_Suffix);
}


MPIX_IMPL_ALL(Stats)
MPIX_IMPL_ALL(TreeStats)
MPIX_IMPL_ALL(DbStats)
MPIX_IMPL_ALL(DupStats)
MPIX_IMPL_ALL(AlignStats)
MPIX_IMPL_ALL(StealingStats)

}

