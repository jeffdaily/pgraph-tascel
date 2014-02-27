/**
 * @file TreeStats.hpp
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 */
#ifndef TREESTATS_H_
#define TREESTATS_H_

#include <iomanip>
#include <iostream>
#include <string>

#include "Stats.hpp"

using ::std::fixed;
using ::std::ostream;
using ::std::setw;
using ::std::showpoint;
using ::std::string;

#define MIN(x, y) (((x)<(y))? (x) : (y))
#define MAX(x, y) (((x)>(y))? (x) : (y))

namespace pgraph {

class TreeStats
{
    public:
        double trees;
        Stats suffixes;
        Stats size;
        Stats size_internal;
        Stats fanout;
        Stats depth;
        Stats suffix_length;
        Stats time_build;
        Stats time_process;
        double time_first;
        double time_last;

        TreeStats()
            : trees(0.0)
            , suffixes()
            , size()
            , size_internal()
            , fanout()
            , depth()
            , suffix_length()
            , time_build()
            , time_process()
            , time_first(0.0)
            , time_last(0.0)
        { }

        static string header() {
            return 
                "         Trees"
                "      Suffixes"
                "         Nodes"
                " InternalNodes"
                "        Fanout"
                "         Depth"
                "  SuffixLength"
                "    Time_Build"
                "  Time_Process"
                "    Time_First"
                "    Time_Last"
                ;
        }

        friend ostream &operator << (ostream &os, const TreeStats &stats) {
            int p = os.precision();
            os.precision(1);
            os << setw(19) << right << "Suffixes" << stats.suffixes << endl;
            os << setw(19) << right << "Nodes" << stats.size << endl;
            os << setw(19) << right << "InternalNodes" << stats.size_internal << endl;
            os << setw(19) << right << "Fanout" << stats.fanout << endl;
            os << setw(19) << right << "Depth" << stats.depth << endl;
            os << setw(19) << right << "SuffixLength" << stats.suffix_length << endl;
            os << setw(19) << right << "TimeBuild" << stats.time_build << endl;
            os << setw(19) << right << "TimeProcess" << stats.time_process << endl;
            os << setw(19) << right << "Trees" << setw(15) << stats.trees << endl;
            os.precision(p); // undo state change
            return os;
        }

        TreeStats& operator += (const TreeStats &stats) {
            if (trees == 0.0) {
                *this = stats;
            }
            else {
                trees += stats.trees;
                suffixes.push_back(stats.suffixes);
                size.push_back(stats.size);
                size_internal.push_back(stats.size_internal);
                fanout.push_back(stats.fanout);
                depth.push_back(stats.depth);
                suffix_length.push_back(stats.suffix_length);
                time_build.push_back(stats.time_build);
                time_process.push_back(stats.time_process);
                time_first = time_first < stats.time_first ? time_first : stats.time_first;
                time_last = time_last > stats.time_last ? time_last : stats.time_last;
            }

            return *this;
        }
};

}; /* namespace pgraph */

#undef MIN
#undef MAX

#endif /* TREESTATS_H_ */
