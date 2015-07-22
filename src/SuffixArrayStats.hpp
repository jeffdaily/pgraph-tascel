/**
 * @file SuffixArrayStats.hpp
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 */
#ifndef SUFFIXARRAYSTATS_H_
#define SUFFIXARRAYSTATS_H_

#include <iomanip>
#include <iostream>
#include <string>

#include "Stats.hpp"

using ::std::ostream;
using ::std::setw;
using ::std::string;

namespace pgraph {

class SuffixArrayStats
{
    public:
        unsigned long arrays;
        Stats suffixes;
        Stats pairs;
        Stats time_build;
        Stats time_process;
        double time_first;
        double time_last;

        SuffixArrayStats()
            : arrays(0U)
            , suffixes()
            , pairs()
            , time_build()
            , time_process()
            , time_first(0.0)
            , time_last(0.0)
        { }

        static string header() {
            return 
                "        Arrays"
                "      Suffixes"
                "         Pairs"
                "    Time_Build"
                "  Time_Process"
                "    Time_First"
                "    Time_Last"
                ;
        }

        friend ostream &operator << (ostream &os, const SuffixArrayStats &stats) {
            os << setw(19) << right << "Suffixes" << stats.suffixes << endl;
            os << setw(19) << right << "Pairs" << stats.pairs << endl;
            os << setw(19) << right << "TimeBuild" << stats.time_build << endl;
            os << setw(19) << right << "TimeProcess" << stats.time_process << endl;
            os << setw(19) << right << "Arrays" << setw(Stats::width()) << stats.arrays << endl;
            return os;
        }

        SuffixArrayStats& operator += (const SuffixArrayStats &stats) {
            if (arrays == 0U) {
                *this = stats;
            }
            else {
                arrays += stats.arrays;
                suffixes.push_back(stats.suffixes);
                pairs.push_back(stats.pairs);
                time_build.push_back(stats.time_build);
                time_process.push_back(stats.time_process);
                time_first = time_first < stats.time_first ? time_first : stats.time_first;
                time_last = time_last > stats.time_last ? time_last : stats.time_last;
            }

            return *this;
        }
};

}; /* namespace pgraph */

#endif /* SUFFIXARRAYSTATS_H_ */
