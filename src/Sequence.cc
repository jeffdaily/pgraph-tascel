#include <cassert>
#include <string>

#include "Sequence.h"

using std::string;


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


Sequence::Sequence(const char *seq, size_t n)
    :   is_owner(false)
    ,   data(seq)
    ,   size(n)
{
}


Sequence::~Sequence()
{
    if (is_owner) {
        delete [] data;
    }
    data = NULL;
}

