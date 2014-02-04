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

#include <yaml-cpp/yaml.h>

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

#define KB (1024U)
#define MB (KB*KB)
#define GB (KB*KB*KB)

static size_t parse_memory_budget(const string& value)
{
    long budget_signed = 0;
    size_t budget = 0;
    char budget_multiplier = 0;
    istringstream iss(value);

    if (isdigit(*value.rbegin())) {
        iss >> budget_signed;
    }
    else {
        iss >> budget_signed >> budget_multiplier;
    }

    if (budget_signed <= 0) {
        cerr << "memory budget must be positive real number" << endl;
        assert(budget_signed > 0);
    }
    budget = size_t(budget_signed);

    if (budget_multiplier == 'b' || budget_multiplier == 'B') {
        budget *= 1; /* byte */
    }
    else if (budget_multiplier == 'k' || budget_multiplier == 'K') {
        budget *= KB; /* kilobyte */
    }
    else if (budget_multiplier == 'm' || budget_multiplier == 'M') {
        budget *= MB; /* megabyte */
    }
    else if (budget_multiplier == 'g' || budget_multiplier == 'G') {
        budget *= GB; /* gigabyte */
    }
    else if (budget_multiplier != 0) {
        cerr << "unrecognized size multiplier" << endl;
        assert(0);
    }

    assert(budget > 0);
    return size_t(budget);
}


Parameters::Parameters()
    : AOL(8)
    , SIM(4)
    , OS(3)
    , exact_match_len(4)
    , window_size(3)
    , open(-10)
    , gap(-1)
    , mem_worker(512U*MB)
    , mem_sequences(2U*GB)
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
    , mem_worker(512U*MB)
    , mem_sequences(2U*GB)
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
    map<string,size_t*> kv_zu;

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
    kv_zu["MemoryWorker"] = &mem_worker;
    kv_zu["MemorySequences"] = &mem_sequences;

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
            else if (kv_zu.find(key) != kv_zu.end()) {
                string value;
                line_ >> value;
                *kv_zu[key] = parse_memory_budget(value);
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


void Parameters::parse_yaml(const char *parameters_file, MPI_Comm comm)
{
    int comm_rank = 0;  /* rank 0 will open the file */
    int status = 0;     /* for MPI routine return codes */

    status = MPI_Comm_rank(comm, &comm_rank);
    MPI_CHECK_IERR(status, comm_rank, comm);

    if (0 == comm_rank) {
        const YAML::Node config = YAML::LoadFile(parameters_file);
        AOL = config["AlignOverLongerSeq"].as<int>(8);
        SIM = config["MatchSimilarity"].as<int>(4);
        OS = config["OptimalScoreOverSelfScore"].as<int>(3);
        exact_match_len = config["ExactMatchLen"].as<int>(8);
        window_size = config["SlideWindowSize"].as<int>(3);
        open = config["Open"].as<int>(-10);
        gap = config["Gap"].as<int>(-1);
        mem_worker = config["MemoryWorker"].as<size_t>(512U*MB);
        mem_sequences = config["MemorySequences"].as<size_t>(2U*GB);
    }

    //mpix_bcast(*this, 0, comm);
    status = MPI_Bcast(this, int(sizeof(Parameters)), MPI_CHAR, 0, comm);
    MPI_CHECK_IERR(status, comm_rank, comm);
}

ostream& operator<< (ostream &os, const Parameters &p)
{
    os << "AlignOverLongerSeq: " << p.AOL << endl;
    os << "MatchSimilarity: " << p.SIM << endl;
    os << "OptimalScoreOverSelfScore: " << p.OS << endl;
    os << "ExactMatchLen: " << p.exact_match_len << endl;
    os << "SlideWindowSize: " << p.window_size << endl;
    os << "Open: " << p.open << endl;
    os << "Gap: " << p.gap << endl;
    os << "MemoryWorker: " << p.mem_worker << endl;
    os << "MemorySequences: " << p.mem_sequences << endl;

    return os;
}

}; /* namespace pgraph */

