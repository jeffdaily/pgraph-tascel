#ifndef _PGRAPH_EDGERESULT_H_
#define _PGRAPH_EDGERESULT_H_

#include <ostream>

using std::ostream;

namespace pgraph {

class EdgeResult {
    public:
        static const char * SEP;

        unsigned long id1;
        unsigned long id2;
        double a;
        double b;
        double c;
        bool is_edge;

        EdgeResult(
                unsigned long id1, unsigned long id2,
                double a, double b, double c, bool is_edge)
            : id1(id1)
            , id2(id2)
            , a(a)
            , b(b)
            , c(c)
            , is_edge(is_edge)
        {}

        friend ostream& operator << (ostream &os, const EdgeResult &edge) {
            os << edge.id1
                << SEP << edge.id2
                << SEP << edge.a
                << SEP << edge.b
                << SEP << edge.c
                /*<< SEP << edge.is_edge*/
                ;
            return os;
        }
};

const char * EdgeResult::SEP = ",";

} /* namespace pgraph */

#endif /* _PGRAPH_EDGERESULT_H_ */
