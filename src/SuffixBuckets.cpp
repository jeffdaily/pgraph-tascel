/**
 * @file SuffixBuckets.cpp
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2010 Washington State University. All rights reserved.
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 */
#include "config.h"

#include <mpi.h>

#include <sys/types.h>
#include <regex.h>

#include <cassert>
#include <cstddef>
#include <limits>
#include <iostream>
#include <vector>

#include <tascel.h>

#include "mpix.hpp"
#include "Parameters.hpp"
#include "SequenceDatabase.hpp"
#include "SuffixBuckets.hpp"
#include "Suffix.hpp"

using ::std::cerr;
using ::std::endl;
using ::std::numeric_limits;
using ::std::size_t;
using ::std::vector;

#ifndef SIZE_MAX
#define SIZE_MAX (size_t(-1))
#endif

namespace pgraph {

const size_t SuffixBuckets::npos(-1);

/** power function for size_t */
static size_t powz(size_t base, size_t n)
{
    size_t p = 1;

    for(/*empty*/; n > 0; --n) {
        assert(p < SIZE_MAX/base);
        p *= base;
    }

    return p;
}


size_t SuffixBuckets::bucket_index(const char *kmer)
{
    size_t value = 0;

    for (int i = 0; i < param.window_size; ++i) {
        size_t index = alphabet_table[(unsigned char)(kmer[i])];
        if (npos == index) {
            cerr << "kmer[" << i << "]='"
                << kmer[i] << "' not in alphabet" << endl;
            MPI_Abort(comm, -1);
        }
        value = value*SIGMA + index;
    }

    return value;
}


string SuffixBuckets::bucket_kmer(size_t bid)
{
    string ret;

    for (int i = param.window_size-1; i >= 0; --i) {
        size_t tmp = powz(SIGMA,i);
        size_t quo = bid / tmp;
        size_t rem = bid % tmp;
        assert(quo <= param.alphabet.size());
        ret += param.alphabet[quo];
        bid = rem;
    }

    return ret;
}


bool SuffixBuckets::filter_out(const char *kmer)
{
    bool result = false;

    string str(kmer, param.window_size);

    for (size_t i=0; i<param.re.size(); ++i) {
        regex_t *preq = param.re[i];
        int retval = regexec(preq, str.c_str(), 0, NULL, 0);
        result = result || (0 == retval);
        if (result) {
            //cout << "filter_out " << str << endl;
            break;
        }
    }

    return result;
}


SuffixBuckets::SuffixBuckets(SequenceDatabase *sequences,
                             const Parameters &param,
                             MPI_Comm comm)
    :   comm(comm)
    ,   comm_rank(mpix::comm_rank(comm))
    ,   comm_size(mpix::comm_size(comm))
    ,   n_buckets(powz(param.alphabet.size(), param.window_size))
    ,   sequences(sequences)
    ,   param(param)
    ,   SIGMA(param.alphabet.size())
    ,   alphabet_table(numeric_limits<unsigned char>::max(), npos)
{
    for (size_t i=0; i<SIGMA; ++i) {
        alphabet_table[(unsigned char)(param.alphabet[i])] = i;
    }
}


SuffixBuckets::~SuffixBuckets()
{
}

}; /* namespace pgraph */

