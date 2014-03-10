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
#include <sstream>
#include <string>

#include "Stats.hpp"

using ::std::fixed;
using ::std::ostream;
using ::std::ostringstream;
using ::std::right;
using ::std::setprecision;
using ::std::setw;
using ::std::showpoint;
using ::std::string;

namespace pgraph {

class AlignStats
{
    public:
        unsigned long edge_counts;
        unsigned long align_counts;
        unsigned long align_skipped;
        Stats time_align;
        double time_kcomb;
        double time_total;
        unsigned long work;
        unsigned long work_skipped;

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
            os << setw(10) << "TimeKcomb";
            os << setw(11) << "TimeTotal";
            os << setw(21) << "Cell_Updates";
            os << setw(21) << "Cell_Updates_Skipped";
            return os.str();
        }

        friend ostream &operator << (ostream &os, const AlignStats &stats) {
            os << setprecision(3) << fixed << showpoint;
            os << setw(8) << stats.edge_counts
               << setw(12) << stats.align_counts
               << setw(12) << stats.align_skipped
               << stats.time_align
               << setw(10) << stats.time_kcomb
               << setw(11) << stats.time_total
               << setw(21) << stats.work
               << setw(21) << stats.work_skipped;
            return os;
        }
};

};

#endif /* ALIGNSTATS_H_ */
