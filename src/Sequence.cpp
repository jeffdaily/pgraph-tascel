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

using std::endl;
using std::string;
using std::strlen;


Sequence::Sequence()
    :   is_owner(false)
    ,   id(NULL)
    ,   id_size(0)
    ,   data(NULL)
    ,   data_size(0)
{
}


Sequence::Sequence(const Sequence &that)
    :   is_owner(false)
    ,   id(that.id)
    ,   id_size(that.id_size)
    ,   data(that.data)
    ,   data_size(that.data_size)
{
}


Sequence::Sequence(const char *data)
    :   is_owner(false)
    ,   id(NULL)
    ,   id_size(0)
    ,   data(data)
    ,   data_size(strlen(data))
{
}


Sequence::Sequence(const char *data, size_t data_size)
    :   is_owner(false)
    ,   id(NULL)
    ,   id_size(0)
    ,   data(data)
    ,   data_size(data_size)
{
}


Sequence::Sequence(const char *id, size_t id_size,
                   const char *data, size_t data_size)
    :   is_owner(false)
    ,   id(id)
    ,   id_size(id_size)
    ,   data(data)
    ,   data_size(data_size)
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
    sequence_t s1 = { "", this->data, this->data_size };
    sequence_t s2 = { "", that.data, that.data_size };

    bigger = this->data_size > that.data_size ? this->data_size : that.data_size;
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
    assert(NULL != data);
    assert(data_size > 0);

    return string(data, data_size);
}


ostream &operator << (ostream &os, const Sequence &s)
{
    for (size_t i = 0; i < s.id_size; ++i) {
        os << s.id[i];
    }
    os << endl;
    for (size_t i = 0; i < s.data_size; ++i) {
        os << s.data[i];
    }
    os << endl;
    return os;
}
