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

#define WIDTH 15

namespace pgraph {

/** http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance */
class Stats
{
    public:
        unsigned long _n;
        double _mean;
        double _M2;
        double _sum;
        double _min;
        double _max;

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

            _n = _n + 1UL;
            delta = x - _mean;
            _mean = _mean + delta/_n;
            _M2 = _M2 + delta * (x - _mean);
        }

        void push_back(const Stats &B) {
            const Stats &A = *this;

            if (B.n() == 0) {
                return;
            }
            else if (A.n() == 0) {
                this->_n = B.n();
                this->_mean = B.mean();
                this->_M2 = B.M2();
                this->_sum = B.sum();
                this->_min = B.min();
                this->_max = B.max();
            }
            else {
                double delta = B.mean() - A.mean();
                double X_n = A.n() + B.n();
                //double X_mean = A.mean() + delta*(B.n()/X_n);
                double X_mean = (A.n()*A.mean() + B.n()*B.mean()) / X_n;
                double X_M2 = A.M2() + B.M2() + delta*delta*A.n()*B.n()/X_n;

                this->_sum += B.sum();
                this->_min = this->_min < B.min() ? this->_min : B.min();
                this->_max = this->_max > B.max() ? this->_max : B.max();
                this->_n = X_n;
                this->_mean = X_mean;
                this->_M2 = X_M2;
                this->_sum += B.sum();
            }
        }

        double n() const { return _n; }
        double mean() const { return _mean; }
        double M2() const { return _M2; }
        double variance() const { return _M2/(_n-1); }
        double stddev() const { return ::std::pow(variance(),0.5); }
        double sum() const { return _sum; }
        double min() const { return _min; }
        double max() const { return _max; }

        friend ::std::ostream& operator << (::std::ostream &os, const Stats &obj);

        static ::std::string header(const ::std::string &prefix="") {
            ::std::ostringstream os;

            //os << ::std::right;
            //os << ::std::setw(WIDTH);
            //os << prefix + "Size";
            os << ::std::right;
            os << ::std::setw(WIDTH);
            os << prefix + "Mean";
            //os << ::std::right;
            //os << ::std::setw(WIDTH);
            //os << prefix + "Variance";
            os << ::std::right;
            os << ::std::setw(WIDTH);
            os << prefix + "StdDev";
            os << ::std::right;
            os << ::std::setw(WIDTH);
            os << prefix + "Sum";
            os << ::std::right;
            os << ::std::setw(WIDTH);
            os << prefix + "Min";
            os << ::std::right;
            os << ::std::setw(WIDTH);
            os << prefix + "Max";

            return os.str();
        }
};

inline ::std::ostream& operator << (::std::ostream &os, const Stats &obj)
{
    ::std::streamsize width = os.width();
    ::std::ios_base::fmtflags flags = os.flags();

    os << ::std::fixed;
    os << ::std::showpoint;
    os << ::std::right;

    //os << ::std::setw(WIDTH);
    //os << obj.n();
    os << ::std::setw(WIDTH);
    os << obj.mean();
    //os << ::std::setw(WIDTH);
    //os << obj.variance();
    os << ::std::setw(WIDTH);
    os << obj.stddev();
    os << ::std::setw(WIDTH);
    os << obj.sum();
    os << ::std::setw(WIDTH);
    os << obj.min();
    os << ::std::setw(WIDTH);
    os << obj.max();

    os.width(width);
    os.flags(flags);

    return os;
}

#undef WIDTH

};

#endif /* _PGRAPH_STATS_H_ */
