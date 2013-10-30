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
        double edge_counts;
        double align_counts;
        double align_skipped;
        double total_times_tot;
        double align_times_tot;
        double align_times_min;
        double align_times_max;
        double kcomb_times_tot;
        double work;
        double work_skipped;
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

        AlignStats &operator >(const AlignStats &other) {
            if (!empty) {
                edge_counts = MAX(edge_counts, other.edge_counts);
                align_counts = MAX(align_counts, other.align_counts);
                align_skipped = MAX(align_skipped, other.align_skipped);
                total_times_tot = MAX(total_times_tot, other.total_times_tot);
                align_times_tot = MAX(align_times_tot, other.align_times_tot);
                kcomb_times_tot = MAX(kcomb_times_tot, other.kcomb_times_tot);
                align_times_min = MIN(align_times_min, other.align_times_min);
                align_times_max = MAX(align_times_max, other.align_times_max);
                work = MAX(work, other.work);
                work_skipped = MAX(work_skipped, other.work_skipped);
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

        AlignStats &operator <(const AlignStats &other) {
            if (!empty) {
                edge_counts = MIN(edge_counts, other.edge_counts);
                align_counts = MIN(align_counts, other.align_counts);
                align_skipped = MIN(align_skipped, other.align_skipped);
                total_times_tot = MIN(total_times_tot, other.total_times_tot);
                align_times_tot = MIN(align_times_tot, other.align_times_tot);
                kcomb_times_tot = MIN(kcomb_times_tot, other.kcomb_times_tot);
                align_times_min = MIN(align_times_min, other.align_times_min);
                align_times_max = MAX(align_times_max, other.align_times_max);
                work = MIN(work, other.work);
                work_skipped = MIN(work_skipped, other.work_skipped);
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

        AlignStats operator /(const int &other) const {
            AlignStats result;
            result.edge_counts = this->edge_counts / other;
            result.align_counts = this->align_counts / other;
            result.align_skipped = this->align_skipped / other;
            result.total_times_tot = this->total_times_tot / other;
            result.align_times_tot = this->align_times_tot / other;
            result.kcomb_times_tot = this->kcomb_times_tot / other;
            result.align_times_min = this->align_times_min;
            result.align_times_max = this->align_times_max;
            result.work = this->work / other;
            result.work_skipped = this->work_skipped / other;
            return result;
        }

        AlignStats stddev(const int &n, const AlignStats *trials) const {
            AlignStats avg = *this / n;
            AlignStats result;
            result.empty = false;
            for (int i=0; i<n; ++i) {
                result.edge_counts += pow(trials[i].edge_counts - avg.edge_counts, 2);
                result.align_counts += pow(trials[i].align_counts - avg.align_counts, 2);
                result.align_skipped += pow(trials[i].align_skipped - avg.align_skipped, 2);
                result.total_times_tot += pow(trials[i].total_times_tot - avg.total_times_tot, 2);
                result.align_times_tot += pow(trials[i].align_times_tot - avg.align_times_tot, 2);
                result.kcomb_times_tot += pow(trials[i].kcomb_times_tot - avg.kcomb_times_tot, 2);
                result.work += pow(trials[i].work - avg.work, 2);
                result.work_skipped += pow(trials[i].work_skipped - avg.work_skipped, 2);
            }
            result = result / n;
            result.edge_counts = pow(result.edge_counts, 0.5);
            result.align_counts = pow(result.align_counts, 0.5);
            result.align_skipped = pow(result.align_skipped, 0.5);
            result.total_times_tot = pow(result.total_times_tot, 0.5);
            result.align_times_tot = pow(result.align_times_tot, 0.5);
            result.kcomb_times_tot = pow(result.kcomb_times_tot, 0.5);
            result.work = pow(result.work, 0.5);
            result.work_skipped = pow(result.work_skipped, 0.5);
            return result;
        }

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
                "   Edges"
                "  Alignments"
                "  Align_skip"
                "    Total_Time"
                "    Kcomb_Time"
                "    Align_Time"
                "   Min_Time"
                "   Max_Time"
                "   Avg_Time"
                "          Cell_Updates"
                "  Cell_Updates_Skipped"
                ;
        }

        friend ostream &operator << (ostream &os, const AlignStats &stats) {
            int p = os.precision();
            os.precision(4);
            os << setw(8) << stats.edge_counts
               << setw(12) << stats.align_counts
               << setw(12) << stats.align_skipped
               << setw(14) << fixed << showpoint << stats.total_times_tot
               << setw(14) << fixed << showpoint << stats.kcomb_times_tot
               << setw(14) << fixed << showpoint << stats.align_times_tot
               << setw(11) << fixed << showpoint << stats.align_times_min
               << setw(11) << fixed << showpoint << stats.align_times_max
               << setw(11) << fixed << showpoint << (stats.align_times_tot / stats.align_counts)
               << setw(22) << fixed << showpoint << stats.work
               << setw(22) << fixed << showpoint << stats.work_skipped;
            os.precision(p); // undo state change
            return os;
        }
};

#undef MIN
#undef MAX

#endif /* ALIGNSTATS_H_ */
