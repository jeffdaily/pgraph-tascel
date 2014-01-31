#include "config.h"

#include <sys/time.h>

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <iterator>
#include <list>
#include <set>
#include <string>
#include <utility>
#include <vector>

#if HAVE_CXX_UNORDERED_SET
#include <unordered_set>
using std::unordered_set;
#elif HAVE_CXX_TR1_UNORDERED_SET
#include <tr1/unordered_set>
using std::tr1::unordered_set;
#endif
#if HAVE_CXX_UNORDERED_SET || HAVE_CXX_TR1_UNORDERED_SET
typedef std::pair<long,long> MyType;
namespace std { namespace tr1
{
    template <> struct hash<MyType> : public unary_function<MyType, size_t>
    {
        size_t operator()(const MyType& v) const
        {
            size_t seed = 0;
            seed ^= hash<long>()(v.first)  + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            seed ^= hash<long>()(v.second) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            return seed;
        }
    };
}}
#endif

using namespace std;

namespace {

    class StackPrinter
    {
        public:
            explicit StackPrinter(const char* msg) :
                msMsg(msg)
            {
                fprintf(stdout, "%s: --begin\n", msMsg.c_str());
                mfStartTime = getTime();
            }

            ~StackPrinter()
            {
                double fEndTime = getTime();
                fprintf(stdout, "%s: --end (duration: %g sec)\n", msMsg.c_str(), (fEndTime-mfStartTime));
            }

        private:
            double getTime() const
            {
                timeval tv;
                gettimeofday(&tv, NULL);
                return tv.tv_sec + tv.tv_usec / 1000000.0;
            }

            string msMsg;
            double mfStartTime;
    };

    class RandomLongIntPairGenerator
    {
        public:
            explicit RandomLongIntPairGenerator(size_t size) :
                index(0),
                pairs()
            {
                srandom(0);
                pairs.reserve(size);
                for (size_t i=0; i<size; ++i) {
                    pairs[i] = make_pair(random(), random());
                }
            }

            ~RandomLongIntPairGenerator()
            {
            }

            const MyType& getNext()
            {
                return pairs[index++];
            }

            void reset()
            {
                index = 0;
            }

        private:
            size_t index;
            vector<MyType> pairs;
    };

    template <typename T>
    class VectorInserter
    {
        public:
            explicit VectorInserter(vector<T> &v) :
                v(v)
            {
            }

            ~VectorInserter() {}

            void operator() (const T &item)
            {
                v.push_back(item);
            }

        private:
            vector<T> &v;
    };
}


int main(int argc, char **argv)
{
    size_t store_size = 50000000;
    RandomLongIntPairGenerator generator(store_size);

    if (argc == 2) {
        long param = atol(argv[1]);
        if (param > 0 && size_t(param) < store_size) {
            store_size = size_t(param);
        }
    }

    {
        StackPrinter __stack_printer__("vector non-reserved");
        string* ptr = 0x00000000;
        vector<void*> store;
        for (size_t i = 0; i < store_size; ++i)
            store.push_back(ptr++);
    }

    {
        StackPrinter __stack_printer__("vector non-reserved of random pairs");
        generator.reset();
        vector<MyType> store;
        for (size_t i = 0; i < store_size; ++i)
            store.push_back(generator.getNext());
    }

    {
        StackPrinter __stack_printer__("vector non-reserved of random pairs callback");
        generator.reset();
        vector<MyType> store;
        VectorInserter<MyType> new_store(store);
        for (size_t i = 0; i < store_size; ++i)
            new_store(generator.getNext());
    }

    {
        StackPrinter __stack_printer__("vector non-reserved of random pairs callback sorted");
        generator.reset();
        vector<MyType> store;
        VectorInserter<MyType> new_store(store);
        for (size_t i = 0; i < store_size; ++i)
            new_store(generator.getNext());
        sort(store.begin(), store.end());
        vector<MyType>::iterator it;
        it = unique(store.begin(), store.end());
        store.resize(distance(store.begin(), it));
    }

    {
        StackPrinter __stack_printer__("vector reserved");
        string* ptr = 0x00000000;
        vector<void*> store;
        store.reserve(store_size);
        for (size_t i = 0; i < store_size; ++i)
            store.push_back(ptr++);
    }

    {
        StackPrinter __stack_printer__("vector reserved of random pairs");
        generator.reset();
        vector<MyType> store;
        store.reserve(store_size);
        for (size_t i = 0; i < store_size; ++i)
            store.push_back(generator.getNext());
    }

    {
        StackPrinter __stack_printer__("vector reserved of random pairs callback");
        generator.reset();
        vector<MyType> store;
        VectorInserter<MyType> new_store(store);
        store.reserve(store_size);
        for (size_t i = 0; i < store_size; ++i)
            new_store(generator.getNext());
    }

    {
        StackPrinter __stack_printer__("vector reserved of random pairs callback sorted");
        generator.reset();
        vector<MyType> store;
        VectorInserter<MyType> new_store(store);
        store.reserve(store_size);
        for (size_t i = 0; i < store_size; ++i)
            new_store(generator.getNext());
        sort(store.begin(), store.end());
        vector<MyType>::iterator it;
        it = unique(store.begin(), store.end());
        store.resize(distance(store.begin(), it));
    }

    {
        StackPrinter __stack_printer__("list");
        string* ptr = 0x00000000;
        list<void*> store;
        for (size_t i = 0; i < store_size; ++i)
            store.push_back(ptr++);
    }

    {
        StackPrinter __stack_printer__("list of random pairs");
        generator.reset();
        list<MyType> store;
        for (size_t i = 0; i < store_size; ++i)
            store.push_back(generator.getNext());
    }

    {
        StackPrinter __stack_printer__("set");
        string* ptr = 0x00000000;
        set<void*> store;   
        for (size_t i = 0; i < store_size; ++i)
            store.insert(ptr++);
    }

    {
        StackPrinter __stack_printer__("set of random pairs");
        generator.reset();
        set<MyType> store;   
        for (size_t i = 0; i < store_size; ++i)
            store.insert(generator.getNext());
    }

#if HAVE_CXX_UNORDERED_SET || HAVE_CXX_TR1_UNORDERED_SET
    {
        StackPrinter __stack_printer__("unordered set");
        string* ptr = 0x00000000;
        unordered_set<void*> store;
        for (size_t i = 0; i < store_size; ++i)
            store.insert(ptr++);
    }

    {
        StackPrinter __stack_printer__("unordered set of random pairs");
        generator.reset();
        unordered_set<MyType> store;
        for (size_t i = 0; i < store_size; ++i)
            store.insert(generator.getNext());
    }
#endif
}
