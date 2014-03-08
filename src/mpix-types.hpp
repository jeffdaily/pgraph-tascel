/**
 * @file mpix-types.hpp
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 *
 * MPI_Datatype creation functions for pgraph structures.
 */
#ifndef _MPIX_TYPES_H_
#define _MPIX_TYPES_H_

#include <tascel.h>

#include "AlignStats.hpp"
#include "Stats.hpp"
#include "TreeStats.hpp"

using ::pgraph::AlignStats;
using ::pgraph::Stats;
using ::pgraph::TreeStats;

using ::tascel::Counter;
using ::tascel::StealingStats;
using ::tascel::Timer;

namespace mpix {

template <>
inline MPI_Datatype get_mpi_datatype<Counter>(const Counter&)
{
    return MPI_LONG;
}

template <>
inline MPI_Datatype build_mpi_datatype<Timer>(const Timer&)
{
    return type_contiguous(2, MPI_DOUBLE);
}

template <>
inline MPI_Datatype build_mpi_datatype<StealingStats>(const StealingStats &object)
{
    MPI_Datatype result;
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
        MPI_Aint(&object.numStealAttempts) - MPI_Aint(&object),
        MPI_Aint(&object.numSteals) - MPI_Aint(&object),
        MPI_Aint(&object.numTasks) - MPI_Aint(&object),
        MPI_Aint(&object.numTasks2) - MPI_Aint(&object),
        MPI_Aint(&object.stealTime) - MPI_Aint(&object),
        MPI_Aint(&object.taskTime) - MPI_Aint(&object),
        MPI_Aint(&object.taskTime2) - MPI_Aint(&object),
        MPI_Aint(&object.tdTime) - MPI_Aint(&object),
        MPI_Aint(&object.miscTime) - MPI_Aint(&object)
    };
    result = type_create_struct(9, blocklen, disp, type);
    type_commit(result);
    return result;
}

template <>
inline MPI_Datatype build_mpi_datatype<Stats>(const Stats &object)
{
    MPI_Datatype result;
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
        MPI_Aint(&object._n) - MPI_Aint(&object),
        MPI_Aint(&object._mean) - MPI_Aint(&object),
        MPI_Aint(&object._M2) - MPI_Aint(&object),
        MPI_Aint(&object._sum) - MPI_Aint(&object),
        MPI_Aint(&object._min) - MPI_Aint(&object),
        MPI_Aint(&object._max) - MPI_Aint(&object)
    };
    result = type_create_struct(6, blocklen, disp, type);
    type_commit(result);
    return result;
}

template <>
inline MPI_Datatype build_mpi_datatype<TreeStats>(const TreeStats &object)
{
    MPI_Datatype result;
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
        MPI_Aint(&object.trees) - MPI_Aint(&object),
        MPI_Aint(&object.suffixes) - MPI_Aint(&object),
        MPI_Aint(&object.size) - MPI_Aint(&object),
        MPI_Aint(&object.size_internal) - MPI_Aint(&object),
        MPI_Aint(&object.fanout) - MPI_Aint(&object),
        MPI_Aint(&object.depth) - MPI_Aint(&object),
        MPI_Aint(&object.suffix_length) - MPI_Aint(&object),
        MPI_Aint(&object.pairs) - MPI_Aint(&object),
        MPI_Aint(&object.time_build) - MPI_Aint(&object),
        MPI_Aint(&object.time_process) - MPI_Aint(&object),
        MPI_Aint(&object.time_first) - MPI_Aint(&object),
        MPI_Aint(&object.time_last) - MPI_Aint(&object)
    };
    result = type_create_struct(12, blocklen, disp, type);
    type_commit(result);
    return result;
}

template <>
inline MPI_Datatype build_mpi_datatype<AlignStats>(const AlignStats &object)
{
    MPI_Datatype result;
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
        MPI_Aint(&object.edge_counts) - MPI_Aint(&object),
        MPI_Aint(&object.align_counts) - MPI_Aint(&object),
        MPI_Aint(&object.align_skipped) - MPI_Aint(&object),
        MPI_Aint(&object.time_align) - MPI_Aint(&object),
        MPI_Aint(&object.time_kcomb) - MPI_Aint(&object),
        MPI_Aint(&object.time_total) - MPI_Aint(&object),
        MPI_Aint(&object.work) - MPI_Aint(&object),
        MPI_Aint(&object.work_skipped) - MPI_Aint(&object)
    };
    result = type_create_struct(8, blocklen, disp, type);
    type_commit(result);
    return result;
}

}

#endif /* _MPIX_TYPES_H_ */
