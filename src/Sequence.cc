#include <cassert>
#include <cstring>
#include <iostream>
#include <string>

#include "Sequence.h"

using std::string;
using std::strlen;


Sequence::Sequence()
    :   is_owner(false)
    ,   data(NULL)
    ,   size(0)
{
}


Sequence::Sequence(const Sequence &seq)
    :   is_owner(false)
    ,   data(seq.data)
    ,   size(seq.size)
{
}


Sequence::Sequence(const char *seq)
    :   is_owner(false)
    ,   data(seq)
    ,   size(strlen(seq))
{
}


Sequence::Sequence(const char *seq, size_t len)
    :   is_owner(false)
    ,   data(seq)
    ,   size(len)
{
}


Sequence::~Sequence()
{
    if (is_owner) {
        delete [] data;
    }
    data = NULL;
}


ostream& operator << (ostream &os, const Sequence &s)
{
    for (size_t i=0; i<s.size; ++i) {
        os << s.data[i];
    }
    return os;
}
