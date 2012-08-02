#ifndef ALIGNSTATS_H_
#define ALIGNSTATS_H_

#include <iomanip>
#include <iostream>

using std::fixed;
using std::setw;
using std::showpoint;

#define MIN(x, y) (((x)<(y))? (x) : (y))
#define MAX(x, y) (((x)>(y))? (x) : (y))

class AlignStats {
    public:
        unsigned long edge_counts;
        unsigned long align_counts;
        double align_times_tot;
        double align_times_min;
        double align_times_max;

        AlignStats()
            : edge_counts(0)
            , align_counts(0)
            , align_times_tot(0.0)
            , align_times_min(DBL_MAX)
            , align_times_max(DBL_MIN)
        { }

        AlignStats& operator +=(const AlignStats &other) {
            edge_counts += other.edge_counts;
            align_counts += other.align_counts;
            align_times_tot += other.align_times_tot;
            align_times_min = MIN(align_times_min,other.align_times_min);
            align_times_max = MAX(align_times_max,other.align_times_max);
            return *this;
        }

        void calc_min(double value) {
            align_times_min = MIN(align_times_min, value);
        }

        void calc_max(double value) {
            align_times_max = MAX(align_times_max, value);
        }

        string getHeader() const {
            return "  Edges Alignments   Total_Time     Min_Time     Max_Time     Avg_Time";
        }

        friend ostream& operator << (ostream &os, const AlignStats &stats) {
            int p = os.precision();
            os.precision(4);
            os << setw(7) << stats.edge_counts
               << setw(11) << stats.align_counts
               << setw(13) << fixed << showpoint << stats.align_times_tot
               << setw(13) << fixed << showpoint << stats.align_times_min
               << setw(13) << fixed << showpoint << stats.align_times_max
               << setw(13) << fixed << showpoint << (stats.align_times_tot/stats.align_counts);
            os.precision(p); // undo state change
            return os;
        }
};

#undef MIN
#undef MAX

#endif /* ALIGNSTATS_H_ */
