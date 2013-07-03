/**
 * @file Sequence.cc
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 */
#include <cassert>
#include <cstring>
#include <iostream>
#include <string>

#include "dynamic.h"

#include "Sequence.hpp"

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
    cell_t result = {0,0,0};
    size_t bigger = 0;
    cell_t **tbl = NULL;
    int **del = NULL;
    int **ins = NULL;
    sequence_t s1 = { "", this->data, this->size };
    sequence_t s2 = { "", that.data, that.size };

    bigger = this->size > that.size ? this->size : that.size;
    tbl = pg_alloc_tbl(2, bigger);
    del = pg_alloc_int(2, bigger);
    ins = pg_alloc_int(2, bigger);

    pg_affine_gap_align_blosum(&s1, &s2, &result, tbl, del, ins, 10, 1);

    pg_free_tbl(tbl, 2);
    pg_free_int(del, 2);
    pg_free_int(ins, 2);

    /* return */
    score = result.score;
    ndig = result.ndig;
    alen = result.alen;
}


Sequence::operator string() const
{
    return string(data, size);
}


ostream &operator << (ostream &os, const Sequence &s)
{
    for (size_t i = 0; i < s.size; ++i) {
        os << s.data[i];
    }
    return os;
}
