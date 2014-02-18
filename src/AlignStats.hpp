/**
 * @file AlignStats.hpp
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 */
#ifndef ALIGNSTATS_H_
#define ALIGNSTATS_H_

#include <iomanip>
#include <iostream>
#include <string>

#include "Stats.hpp"

using std::fixed;
using std::ostream;
using std::setw;
using std::showpoint;
using std::string;

#define MIN(x, y) (((x)<(y))? (x) : (y))
#define MAX(x, y) (((x)>(y))? (x) : (y))

namespace pgraph {

class AlignStats
{
    public:
        double edge_counts;
        double align_counts;
        double align_skipped;
        Stats time_align;
        double time_kcomb;
        double time_total;
        double work;
        double work_skipped;

        AlignStats()
            : edge_counts(0)
            , align_counts(0)
            , align_skipped(0)
            , time_align()
            , time_kcomb(0.0)
            , time_total(0.0)
            , work(0)
            , work_skipped(0)
        { }

        static string header() {
            ostringstream os;
            os << right;
            os << setw(8) << "Edges";
            os << setw(12) << "Alignments";
            os << setw(12) << "Align_skip";
            os << Stats::header("TAlign");
            os << setw(14) << "TimeKcomb";
            os << setw(14) << "TimeTotal";
            os << setw(22) << "Cell_Updates";
            os << setw(22) << "Cell_Updates_Skipped";
            return os.str();
        }

        friend ostream &operator << (ostream &os, const AlignStats &stats) {
            int p = os.precision();
            os.precision(2);
            os << setw(8) << stats.edge_counts
               << setw(12) << stats.align_counts
               << setw(12) << stats.align_skipped
               << stats.time_align
               << setw(14) << fixed << showpoint << stats.time_kcomb
               << setw(14) << fixed << showpoint << stats.time_total
               << setw(22) << fixed << showpoint << stats.work
               << setw(22) << fixed << showpoint << stats.work_skipped;
            os.precision(p); // undo state change
            return os;
        }
};

};

#undef MIN
#undef MAX

#endif /* ALIGNSTATS_H_ */
