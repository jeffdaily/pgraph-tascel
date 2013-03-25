#include <cassert>
#include <cstring>
#include <iostream>
#include <string>

#include "dynamic.h"

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


void Sequence::align(const Sequence &that, int &score, int &ndig, int &alen)
{
    cell_t result;
    cell_t **tbl = alloc_tbl(2, this->size > that.size ? this->size : that.size);
    int **del = alloc_int(2, this->size > that.size ? this->size : that.size);
    int **ins = alloc_int(2, this->size > that.size ? this->size : that.size);

    affine_gap_align(this->data, this->size, that.data, that.size, &result,
            tbl, del, ins);

    free_tbl(tbl, 2);
    free_int(del, 2);
    free_int(ins, 2);

    /* return */
    score = result.score;
    ndig = result.ndig;
    alen = result.alen;
}


Sequence::operator string() const
{
    return string(data, size);
}


ostream& operator << (ostream &os, const Sequence &s)
{
    for (size_t i=0; i<s.size; ++i) {
        os << s.data[i];
    }
    return os;
}
