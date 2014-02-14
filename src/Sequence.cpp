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

#include "alignment.hpp"

#include "Sequence.hpp"

using std::endl;
using std::size_t;
using std::string;
using std::strlen;

namespace pgraph {


Sequence::Sequence()
    :   is_owner(false)
    ,   buffer(NULL)
    ,   id_offset(0)
    ,   id_length(0)
    ,   sequence_offset(0)
    ,   sequence_length(0)
{
}


Sequence::Sequence(const Sequence &that)
    :   is_owner(false)
    ,   buffer(that.buffer)
    ,   id_offset(that.id_offset)
    ,   id_length(that.id_length)
    ,   sequence_offset(that.sequence_offset)
    ,   sequence_length(that.sequence_length)
{
}


Sequence::Sequence(const char *data, const bool &owns)
    :   is_owner(owns)
    ,   buffer(data)
    ,   id_offset(0)
    ,   id_length(0)
    ,   sequence_offset(0)
    ,   sequence_length(strlen(data))
{
}


Sequence::Sequence(const char *data, const size_t &length, const bool &owns)
    :   is_owner(owns)
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

cell_t align_global_affine(
        const Sequence &s1,
        const Sequence &s2,
        match_t match,
        int open, int gap,
        cell_t * const restrict * const restrict tbl,
        int * const restrict * const restrict del,
        int * const restrict * const restrict ins)
{
    return align_global_affine(s1.data(), s1.size(), s2.data(), s2.size(),
            match, open, gap, tbl, del, ins);
}

cell_t align_global_affine(
        const Sequence &s1,
        const Sequence &s2,
        const int * const restrict * const restrict sub,
        const int * const restrict map, char first,
        int open, int gap,
        cell_t * const restrict * const restrict tbl,
        int * const restrict * const restrict del,
        int * const restrict * const restrict ins)
{
    return align_global_affine(s1.data(), s1.size(), s2.data(), s2.size(),
            sub, map, first, open, gap, tbl, del, ins);
}

cell_t align_global_affine(
        const Sequence &s1,
        const Sequence &s2,
        int match, int mismatch,
        int open, int gap,
        cell_t * const restrict * const restrict tbl,
        int * const restrict * const restrict del,
        int * const restrict * const restrict ins)
{
    return align_global_affine(s1.data(), s1.size(), s2.data(), s2.size(),
            match, mismatch, open, gap, tbl, del, ins);
}

cell_t align_global_affine(
        const Sequence &s1,
        const Sequence &s2,
        int open, int gap,
        cell_t * const restrict * const restrict tbl,
        int * const restrict * const restrict del,
        int * const restrict * const restrict ins)
{
    return align_global_affine(s1.data(), s1.size(), s2.data(), s2.size(),
            open, gap, tbl, del, ins);
}

cell_t align_semi_affine(
        const Sequence &s1,
        const Sequence &s2,
        match_t match,
        int open, int gap,
        cell_t * const restrict * const restrict tbl,
        int * const restrict * const restrict del,
        int * const restrict * const restrict ins)
{
    return align_semi_affine(s1.data(), s1.size(), s2.data(), s2.size(),
            match, open, gap, tbl, del, ins);
}

cell_t align_semi_affine(
        const Sequence &s1,
        const Sequence &s2,
        const int * const restrict * const restrict sub,
        const int * const restrict map, char first,
        int open, int gap,
        cell_t * const restrict * const restrict tbl,
        int * const restrict * const restrict del,
        int * const restrict * const restrict ins)
{
    return align_semi_affine(s1.data(), s1.size(), s2.data(), s2.size(),
            sub, map, first, open, gap, tbl, del, ins);
}

cell_t align_semi_affine(
        const Sequence &s1,
        const Sequence &s2,
        int match, int mismatch,
        int open, int gap,
        cell_t * const restrict * const restrict tbl,
        int * const restrict * const restrict del,
        int * const restrict * const restrict ins)
{
    return align_semi_affine(s1.data(), s1.size(), s2.data(), s2.size(),
            match, mismatch, open, gap, tbl, del, ins);
}

cell_t align_semi_affine(
        const Sequence &s1,
        const Sequence &s2,
        int open, int gap,
        cell_t * const restrict * const restrict tbl,
        int * const restrict * const restrict del,
        int * const restrict * const restrict ins)
{
    return align_semi_affine(s1.data(), s1.size(), s2.data(), s2.size(),
            open, gap, tbl, del, ins);
}

cell_t align_local_affine(
        const Sequence &s1,
        const Sequence &s2,
        match_t match,
        int open, int gap,
        cell_t * const restrict * const restrict tbl,
        int * const restrict * const restrict del,
        int * const restrict * const restrict ins)
{
    return align_local_affine(s1.data(), s1.size(), s2.data(), s2.size(),
            match, open, gap, tbl, del, ins);
}

cell_t align_local_affine(
        const Sequence &s1,
        const Sequence &s2,
        const int * const restrict * const restrict sub,
        const int * const restrict map, char first,
        int open, int gap,
        cell_t * const restrict * const restrict tbl,
        int * const restrict * const restrict del,
        int * const restrict * const restrict ins)
{
    return align_local_affine(s1.data(), s1.size(), s2.data(), s2.size(),
            sub, map, first, open, gap, tbl, del, ins);
}

cell_t align_local_affine(
        const Sequence &s1,
        const Sequence &s2,
        int match, int mismatch,
        int open, int gap,
        cell_t * const restrict * const restrict tbl,
        int * const restrict * const restrict del,
        int * const restrict * const restrict ins)
{
    return align_local_affine(s1.data(), s1.size(), s2.data(), s2.size(),
            match, mismatch, open, gap, tbl, del, ins);
}

cell_t align_local_affine(
        const Sequence &s1,
        const Sequence &s2,
        int open, int gap,
        cell_t * const restrict * const restrict tbl,
        int * const restrict * const restrict del,
        int * const restrict * const restrict ins)
{
    return align_local_affine(s1.data(), s1.size(), s2.data(), s2.size(),
            open, gap, tbl, del, ins);
}

bool is_edge(
        const cell_t &result,
        const Sequence &s1,
        const Sequence &s2,
        int AOL,
        int SIM,
        int OS,
        int &self_score,
        size_t &max_len,
        const int ** const restrict sub,
        const int * const restrict map, char first)
{
    return is_edge(result, s1.data(), s1.size(), s2.data(), s2.size(),
            AOL, SIM, OS, self_score, max_len, sub, map, first);
}

bool is_edge(
        const cell_t &result,
        const Sequence &s1,
        const Sequence &s2,
        int AOL,
        int SIM,
        int OS,
        int &self_score,
        size_t &max_len,
        match_t match)
{
    return is_edge(result, s1.data(), s1.size(), s2.data(), s2.size(),
            AOL, SIM, OS, self_score, max_len, match);
}

bool is_edge(
        const cell_t &result,
        const Sequence &s1,
        const Sequence &s2,
        int AOL,
        int SIM,
        int OS,
        int &self_score,
        size_t &max_len,
        int match)
{
    return is_edge(result, s1.data(), s1.size(), s2.data(), s2.size(),
            AOL, SIM, OS, self_score, max_len, match);
}

bool is_edge(
        const cell_t &result,
        const Sequence &s1,
        const Sequence &s2,
        int AOL,
        int SIM,
        int OS,
        int &self_score,
        size_t &max_len)
{
    return is_edge(result, s1.data(), s1.size(), s2.data(), s2.size(),
            AOL, SIM, OS, self_score, max_len);
}

}; /* namespace pgraph */

