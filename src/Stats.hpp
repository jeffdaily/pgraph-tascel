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
#include <ostream>
#include <sstream>
#include <string>

using ::std::ostream;
using ::std::pow;
using ::std::string;
using ::std::ostringstream;
using ::std::right;
using ::std::setw;
using ::std::streamsize;
using ::std::ios_base;
using ::std::fixed;
using ::std::showpoint;

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
                unsigned long X_n = A.n() + B.n();
                //double X_mean = A.mean() + delta*(B.n()/X_n);
                double X_mean = (A.n()*A.mean() + B.n()*B.mean()) / X_n;
                double X_M2 = A.M2() + B.M2() + delta*delta*A.n()*B.n()/X_n;

                this->_n = X_n;
                this->_mean = X_mean;
                this->_M2 = X_M2;
                this->_sum += B.sum();
                this->_min = this->_min < B.min() ? this->_min : B.min();
                this->_max = this->_max > B.max() ? this->_max : B.max();
            }
        }

        unsigned long n() const { return _n; }
        double mean() const { return _mean; }
        double M2() const { return _M2; }
        double variance() const { return _M2/(_n-1); }
        double stddev() const { return pow(variance(),0.5); }
        double sum() const { return _sum; }
        double min() const { return _min; }
        double max() const { return _max; }

        friend ostream& operator << (ostream &os, const Stats &obj);

        static string header(const string &prefix="") {
            ostringstream os;
            const int WIDTH = width();

            //os << right;
            //os << setw(WIDTH);
            //os << prefix + "Size";
            os << right;
            os << setw(WIDTH);
            os << prefix + "Mean";
            //os << right;
            //os << setw(WIDTH);
            //os << prefix + "Variance";
            os << right;
            os << setw(WIDTH);
            os << prefix + "SDev";
            os << right;
            os << setw(WIDTH);
            os << prefix + "Sum";
            os << right;
            os << setw(WIDTH);
            os << prefix + "Min";
            os << right;
            os << setw(WIDTH);
            os << prefix + "Max";

            return os.str();
        }

        static int width(const int &value=-1) {
            static int _width_ = 15;
            if (value >= 0) {
                _width_ = value;
            }
            return _width_;
        }
};


inline ostream& operator << (ostream &os, const Stats &obj)
{
    const int WIDTH = Stats::width();
    streamsize width = os.width();
    ios_base::fmtflags flags = os.flags();

    os << fixed;
    os << showpoint;
    os << right;

    //os << setw(WIDTH);
    //os << obj.n();
    os << setw(WIDTH);
    os << obj.mean();
    //os << setw(WIDTH);
    //os << obj.variance();
    os << setw(WIDTH);
    os << obj.stddev();
    os << setw(WIDTH);
    os << obj.sum();
    os << setw(WIDTH);
    os << obj.min();
    os << setw(WIDTH);
    os << obj.max();

    os.width(width);
    os.flags(flags);

    return os;
}

#undef WIDTH

};

#endif /* _PGRAPH_STATS_H_ */
