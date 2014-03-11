/**
 * @file DupStats.hpp
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 */
#ifndef _DUPSTATS_H_
#define _DUPSTATS_H_

#include <iomanip>
#include <iostream>
#include <string>

#include "Stats.hpp"

using ::std::endl;
using ::std::ostream;
using ::std::setw;
using ::std::string;

namespace pgraph {

class DupStats
{
    public:
        Stats time;
        Stats checked;
        Stats returned;

        DupStats()
            : time()
            , checked()
            , returned()
        { }

        static string header() {
            return 
                "          Time"
                "       Checked"
                "      Returned"
                ;
        }

        friend ostream &operator << (ostream &os, const DupStats &stats) {
            os << setw(19) << right << "Time" << stats.time << endl;
            os << setw(19) << right << "Checked" << stats.checked << endl;
            os << setw(19) << right << "Returned" << stats.returned << endl;
            return os;
        }

        DupStats& operator += (const DupStats &stats) {
            time.push_back(stats.time);
            checked.push_back(stats.checked);
            returned.push_back(stats.returned);
            return *this;
        }
};

}; /* namespace pgraph */

#endif /* _DUPSTATS_H_ */
