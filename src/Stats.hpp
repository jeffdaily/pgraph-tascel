/**
 * @file Stats.hpp
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 */
#ifndef _PGRAPH_STATS_H_
#define _PGRAPH_STATS_H_

#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

namespace pgraph {

class Stats
{
    public:
        Stats()
            : _n(0UL)
            , _mean(0.0)
            , _M2(0.0)
            , _sum(0.0)
            , _min(0.0)
            , _max(0.0)
        {}

        virtual ~Stats() {}

        void push_back(const double &x) {
            double delta = 0;

            /* extra stats */
            _sum = _sum + x;
            if (0UL == _n) {
                _min = x;
                _max = x;
            }
            else {
                _min = _min < x ? _min : x;
                _max = _max > x ? _max : x;
            }

            /* Knuth's online algorithm */
            _n = _n + 1UL;
            delta = x - _mean;
            _mean = _mean + delta/_n;
            _M2 = _M2 + delta * (x - _mean);
        }

        double size() const { return _n; }
        double mean() const { return _mean; }
        double variance() const { return _M2/(_n-1); }
        double stddev() const { return ::std::pow(variance(),0.5); }
        double sum() const { return _sum; }
        double min() const { return _min; }
        double max() const { return _max; }

        friend ::std::ostream& operator << (::std::ostream &os, const Stats &obj);

        static ::std::string header() {
            ostringstream os;

            os << ::std::right;
            os << ::std::setw(15);
            os << "Size";
            os << ::std::right;
            os << ::std::setw(15);
            os << "Mean";
            os << ::std::right;
            os << ::std::setw(15);
            os << "Variance";
            os << ::std::right;
            os << ::std::setw(15);
            os << "StdDev";
            os << ::std::right;
            os << ::std::setw(15);
            os << "Sum";
            os << ::std::right;
            os << ::std::setw(15);
            os << "Min";
            os << ::std::right;
            os << ::std::setw(15);
            os << "Max";

            return os.str();
        }
        
    protected:
        unsigned long _n;
        double _mean;
        double _M2;
        double _sum;
        double _min;
        double _max;
};

inline ::std::ostream& operator << (::std::ostream &os, const Stats &obj)
{
    ::std::streamsize width = os.width();
    ::std::streamsize precision = os.precision();
    ::std::ios_base::fmtflags flags = os.flags();

    os << ::std::fixed;
    os << ::std::showpoint;
    os << ::std::right;
    os << ::std::setprecision(2);

    os << ::std::setw(15);
    os << obj.size();
    os << ::std::setw(15);
    os << obj.mean();
    os << ::std::setw(15);
    os << obj.variance();
    os << ::std::setw(15);
    os << obj.stddev();
    os << ::std::setw(15);
    os << obj.sum();
    os << ::std::setw(15);
    os << obj.min();
    os << ::std::setw(15);
    os << obj.max();

    os.width(width);
    os.precision(precision);
    os.flags(flags);

    return os;
}

};

#endif /* _PGRAPH_STATS_H_ */
