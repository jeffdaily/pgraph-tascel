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

using std::fixed;
using std::ostream;
using std::setw;
using std::showpoint;
using std::string;

#define MIN(x, y) (((x)<(y))? (x) : (y))
#define MAX(x, y) (((x)>(y))? (x) : (y))

class TreeStats
{
    public:
        double trees;
        double suffixes;
        double times;

        TreeStats()
            : trees(0.0)
            , suffixes(0.0)
            , times(0.0)
        { }

        TreeStats &operator >(const TreeStats &other) {
            trees = MAX(trees, other.trees);
            suffixes = MAX(suffixes, other.suffixes);
            times = MAX(times, other.times);
            return *this;
        }

        TreeStats &operator <(const TreeStats &other) {
            trees = MIN(trees, other.trees);
            suffixes = MIN(suffixes, other.suffixes);
            times = MIN(times, other.times);
            return *this;
        }

        TreeStats operator /(const int &other) const {
            TreeStats result;
            result.trees = this->trees / other;
            result.suffixes = this->suffixes / other;
            result.times = this->times / other;
            return result;
        }

        TreeStats stddev(const int &n, const TreeStats *trials) const {
            TreeStats avg = *this / n;
            TreeStats result;
            for (int i=0; i<n; ++i) {
                result.trees += pow(trials[i].trees - avg.trees, 2);
                result.suffixes += pow(trials[i].suffixes - avg.suffixes, 2);
                result.times += pow(trials[i].times - avg.times, 2);
            }
            result = result / n;
            result.trees = pow(result.trees, 0.5);
            result.suffixes = pow(result.suffixes, 0.5);
            result.times = pow(result.times, 0.5);
            return result;
        }

        TreeStats &operator +=(const TreeStats &other) {
            trees += other.trees;
            suffixes += other.suffixes;
            times += other.times;
            return *this;
        }

        string getHeader() const {
            return 
                "         Trees"
                "      Suffixes"
                "         Times"
                ;
        }

        friend ostream &operator << (ostream &os, const TreeStats &stats) {
            int p = os.precision();
            os.precision(4);
            os << setw(14) << fixed << showpoint << stats.trees
               << setw(14) << fixed << showpoint << stats.suffixes
               << setw(14) << fixed << showpoint << stats.times;
            os.precision(p); // undo state change
            return os;
        }
};

#undef MIN
#undef MAX

#endif /* TREESTATS_H_ */
