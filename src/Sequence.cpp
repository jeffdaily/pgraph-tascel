/**
 * @file Sequence.cc
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 */
#include "config.h"

#include <cassert>
#include <cstddef>
#include <cstring>
#include <iostream>
#include <string>

#include "alignment.hpp"

#include "Sequence.hpp"

using std::endl;
using std::size_t;
using std::string;
using std::strlen;

namespace pgraph {


Sequence::Sequence()
    :   is_owner(false)
    ,   data(NULL)
    ,   id_offset(0)
    ,   id_length(0)
    ,   sequence_offset(0)
    ,   sequence_length(0)
{
}


Sequence::Sequence(const Sequence &that)
    :   is_owner(false)
    ,   data(that.data)
    ,   id_offset(that.id_offset)
    ,   id_length(that.id_length)
    ,   sequence_offset(that.sequence_offset)
    ,   sequence_length(that.sequence_length)
{
}


Sequence::Sequence(const char *data, const bool &owns)
    :   is_owner(owns)
    ,   data(data)
    ,   id_offset(0)
    ,   id_length(0)
    ,   sequence_offset(0)
    ,   sequence_length(strlen(data))
{
}


Sequence::Sequence(const char *data, const size_t &length, const bool &owns)
    :   is_owner(owns)
    ,   data(data)
    ,   id_offset(0)
    ,   id_length(0)
    ,   sequence_offset(0)
    ,   sequence_length(length)
{
}


Sequence::Sequence(const char *data,
                   const size_t &id_offset,
                   const size_t &id_length,
                   const size_t &sequence_offset,
                   const size_t &sequence_length,
                   const bool &owns)
    :   is_owner(owns)
    ,   data(data)
    ,   id_offset(id_offset)
    ,   id_length(id_length)
    ,   sequence_offset(sequence_offset)
    ,   sequence_length(sequence_length)
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

    bigger = this->sequence_length > that.sequence_length ?
        this->sequence_length :
        that.sequence_length;
    tbl = allocate_cell_table(2, bigger);
    del = allocate_int_table(2, bigger);
    ins = allocate_int_table(2, bigger);

    result = affine_gap_align_blosum(
            &this->data[sequence_offset], this->sequence_length,
            &that.data[sequence_offset], that.sequence_length,
            tbl, del, ins, 10, 1);

    free_cell_table(tbl, 2);
    free_int_table(del, 2);
    free_int_table(ins, 2);

    /* return */
    score = result.score;
    ndig = result.matches;
    alen = result.length;
}


Sequence::operator string() const
{
    assert(NULL != data);
    assert(sequence_length > 0);

    return string(&data[sequence_offset], sequence_length);
}


ostream &operator << (ostream &os, const Sequence &s)
{
    for (size_t i = 0; i < s.id_length; ++i) {
        os << s.data[s.id_offset + i];
    }
    os << endl;
    for (size_t i = 0; i < s.sequence_length; ++i) {
        os << s.data[s.sequence_offset + i];
    }
    os << endl;
    return os;
}

}; /* namespace pgraph */

