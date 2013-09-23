/**
 * @file Parameters.cpp
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2010 Washington State University. All rights reserved.
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 */
#include "config.h"

#include <cassert>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>

#include "Parameters.hpp"
#include "mpix.hpp"

#define COMMENT '#'

using std::cerr;
using std::endl;
using std::getline;
using std::ifstream;
using std::istringstream;
using std::map;
using std::string;

namespace pgraph {

Parameters::Parameters()
    : AOL(8)
    , SIM(4)
    , OS(3)
    , exact_match_len(4)
    , window_size(3)
    , open(-10)
    , gap(-1)
{
}


Parameters::Parameters(const char *parameters_file, MPI_Comm comm)
    : AOL(8)
    , SIM(4)
    , OS(3)
    , exact_match_len(4)
    , window_size(3)
    , open(-10)
    , gap(-1)
{
    parse(parameters_file, comm);
}


void Parameters::parse(const char *parameters_file, MPI_Comm comm)
{
    int comm_rank = 0;  /* rank 0 will open the file */
    int comm_size = 0;
    int status = 0;     /* for MPI routine return codes */
    map<string,int*> kv_int;
    map<string,string> kv_str;

    status = MPI_Comm_rank(comm, &comm_rank);
    MPI_CHECK_IERR(status, comm_rank, comm);
    status = MPI_Comm_size(comm, &comm_size);
    MPI_CHECK_IERR(status, comm_rank, comm);

    /* create mapping for parser */
    kv_int["AlignOverLongerSeq"] = &AOL;
    kv_int["MatchSimilarity"] = &SIM;
    kv_int["OptimalScoreOverSelfScore"] = &OS;
    kv_int["SlideWindowSize"] = &window_size;
    kv_int["ExactMatchLen"] = &exact_match_len;
    kv_str["CWD"] = "";

    if (0 == comm_rank) {
        ifstream is(parameters_file);
        string line;

        while (getline(is, line)) {
            string key;
            istringstream line_(line);

            if (line.empty() || *line.begin() == COMMENT) {
                continue;
            }
            line_ >> key;
            if (kv_int.find(key) != kv_int.end()) {
                int value;
                line_ >> value;
                *kv_int[key] = value;
            }
            else if (kv_str.find(key) != kv_str.end()) {
                string value;
                line_ >> value;
                kv_str[key] = value;
            }
            else {
                if (0 == comm_rank) {
                    cerr << "invalid config file param '" << key << "'" << endl;
                }
                MPI_Abort(comm, -1);
            }
        }
    }

    mpix_bcast(this, 0, comm);
}

}; /* namespace pgraph */

