/**
 * @file Sequence.cpp
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

#include "Sequence.hpp"

using std::endl;
using std::size_t;
using std::string;
using std::strlen;

namespace pgraph {


Sequence::Sequence()
    :   is_owner(false)
    ,   has_delimiter(false)
    ,   buffer(NULL)
    ,   id_offset(0)
    ,   id_length(0)
    ,   sequence_offset(0)
    ,   sequence_length(0)
{
}


Sequence::Sequence(const Sequence &that)
    :   is_owner(false)
    ,   has_delimiter(that.has_delimiter)
    ,   buffer(that.buffer)
    ,   id_offset(that.id_offset)
    ,   id_length(that.id_length)
    ,   sequence_offset(that.sequence_offset)
    ,   sequence_length(that.sequence_length)
{
}


Sequence::Sequence(const char *data, const bool &owns)
    :   is_owner(owns)
    ,   has_delimiter(false)
    ,   buffer(data)
    ,   id_offset(0)
    ,   id_length(0)
    ,   sequence_offset(0)
    ,   sequence_length(strlen(data))
{
}


Sequence::Sequence(const char *data, const size_t &length, const bool &owns)
    :   is_owner(owns)
    ,   has_delimiter(false)
    ,   buffer(data)
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
    ,   has_delimiter(false)
    ,   buffer(data)
    ,   id_offset(id_offset)
    ,   id_length(id_length)
    ,   sequence_offset(sequence_offset)
    ,   sequence_length(sequence_length)
{
}


ostream &operator << (ostream &os, const Sequence &s)
{
    for (size_t i = 0; i < s.id_length; ++i) {
        os << s.buffer[s.id_offset + i];
    }
    os << endl;
    for (size_t i = 0; i < s.sequence_length; ++i) {
        os << s.buffer[s.sequence_offset + i];
    }
    os << endl;
    return os;
}

}; /* namespace pgraph */

