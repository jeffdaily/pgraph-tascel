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

using std::fixed;
using std::ostream;
using std::setw;
using std::showpoint;
using std::string;

#define MIN(x, y) (((x)<(y))? (x) : (y))
#define MAX(x, y) (((x)>(y))? (x) : (y))

class AlignStats
{
    public:
        unsigned long edge_counts;
        unsigned long align_counts;
        unsigned long align_skipped;
        double total_times_tot;
        double align_times_tot;
        double align_times_min;
        double align_times_max;
        double kcomb_times_tot;
        unsigned long work;
        unsigned long work_skipped;
        bool empty;

        AlignStats()
            : edge_counts(0)
            , align_counts(0)
            , align_skipped(0)
            , total_times_tot(0.0)
            , align_times_tot(0.0)
            , align_times_min(0.0)
            , align_times_max(0.0)
            , kcomb_times_tot(0.0)
            , work(0)
            , work_skipped(0)
            , empty(true)
        { }

        AlignStats &operator +=(const AlignStats &other) {
            if (!empty) {
                edge_counts += other.edge_counts;
                align_counts += other.align_counts;
                align_skipped += other.align_skipped;
                total_times_tot += other.total_times_tot;
                align_times_tot += other.align_times_tot;
                kcomb_times_tot += other.kcomb_times_tot;
                align_times_min = MIN(align_times_min, other.align_times_min);
                align_times_max = MAX(align_times_max, other.align_times_max);
                work += other.work;
                work_skipped += other.work_skipped;
            }
            else {
                edge_counts = other.edge_counts;
                align_counts = other.align_counts;
                align_skipped = other.align_skipped;
                total_times_tot = other.total_times_tot;
                align_times_tot = other.align_times_tot;
                kcomb_times_tot = other.kcomb_times_tot;
                align_times_min = other.align_times_min;
                align_times_max = other.align_times_max;
                work = other.work;
                work_skipped = other.work_skipped;
                empty = false;
            }
            return *this;
        }

        void calc_min(double value) {
            align_times_min = MIN(align_times_min, value);
        }

        void calc_max(double value) {
            align_times_max = MAX(align_times_max, value);
        }

        string getHeader() const {
            return 
                "  Edges"
                " Alignments"
                " Align_skip"
                "   Total_Time"
                "   Kcomb_Time"
                "   Align_Time"
                "  Min_Time"
                "  Max_Time"
                "  Avg_Time"
                "         Cell_Updates"
                " Cell_Updates_Skipped"
                ;
        }

        friend ostream &operator << (ostream &os, const AlignStats &stats) {
            int p = os.precision();
            os.precision(4);
            os << setw(7) << stats.edge_counts
               << setw(11) << stats.align_counts
               << setw(11) << stats.align_skipped
               << setw(13) << fixed << showpoint << stats.total_times_tot
               << setw(13) << fixed << showpoint << stats.kcomb_times_tot
               << setw(13) << fixed << showpoint << stats.align_times_tot
               << setw(10) << fixed << showpoint << stats.align_times_min
               << setw(10) << fixed << showpoint << stats.align_times_max
               << setw(10) << fixed << showpoint << (stats.align_times_tot / stats.align_counts)
               << setw(21) << fixed << showpoint << stats.work
               << setw(21) << fixed << showpoint << stats.work_skipped;
            os.precision(p); // undo state change
            return os;
        }
};

#undef MIN
#undef MAX

#endif /* ALIGNSTATS_H_ */
