/**
 * @file DbStats.hpp
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 */
#ifndef _DBSTATS_H_
#define _DBSTATS_H_

#include <iomanip>
#include <iostream>
#include <string>

#include "Stats.hpp"

using ::std::endl;
using ::std::ostream;
using ::std::setw;
using ::std::string;

namespace pgraph {

class DbStats
{
    public:
        Stats time;
        Stats bytes;
        bool cum;

        DbStats()
            : time()
            , bytes()
            , cum(false)
        { }

        static string header() {
            return 
                "     DbGetTime"
                "    DbGetBytes"
                ;
        }

        friend ostream &operator << (ostream &os, const DbStats &stats) {
            if (stats.cum) {
                os << setw(19) << right << " TotDbGetTime" << stats.time << endl;
                os << setw(19) << right << "TotDbGetBytes" << stats.bytes << endl;
                os << setw(19) << right << "TotDbGetCount" << " " << stats.bytes.n() << endl;
            }
            else {
                os << setw(19) << right << " DbGetTime" << stats.time << endl;
                os << setw(19) << right << "DbGetBytes" << stats.bytes << endl;
                os << setw(19) << right << "DbGetCount" << " " << stats.bytes.n() << endl;
            }
            return os;
        }

        DbStats& operator += (const DbStats &stats) {
            time.push_back(stats.time);
            bytes.push_back(stats.bytes);
            return *this;
        }
};

}; /* namespace pgraph */

#endif /* _DBSTATS_H_ */
