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

/* Keys */
const string Parameters::KEY_ALIGN_OVER_LONGER_SEQUENCE("AlignOverLongerSequence");
const string Parameters::KEY_MATCH_SIMILARITY("MatchSimilarity");
const string Parameters::KEY_OPTIMAL_SCORE_OVER_SELF_SCORE("OptimalScoreOverSelfScore");
const string Parameters::KEY_EXACT_MATCH_LENGTH("ExactMatchLength");
const string Parameters::KEY_SLIDE_WINDOW_SIZE("SlideWindowSize");
const string Parameters::KEY_OPEN("Open");
const string Parameters::KEY_GAP("Gap");
const string Parameters::KEY_MEMORY_WORKER("MemoryWorker");
const string Parameters::KEY_MEMORY_SEQUENCES("MemorySequences");
const string Parameters::KEY_SKIP_PREFIXES("SkipPrefixes");
const string Parameters::KEY_OUTPUT_ALL("OutputAll");
const string Parameters::KEY_DISTRIBUTE_SEQUENCES("DistributeSequences");
const string Parameters::KEY_USE_LENGTH_FILTER("UseLengthFilter");
const string Parameters::KEY_USE_ITERATOR("UseIterator");
/* Defaults */
const int Parameters::DEF_ALIGN_OVER_LONGER_SEQUENCE(8);
const int Parameters::DEF_MATCH_SIMILARITY(4);
const int Parameters::DEF_OPTIMAL_SCORE_OVER_SELF_SCORE(3);
const int Parameters::DEF_EXACT_MATCH_LENGTH(4);
const int Parameters::DEF_SLIDE_WINDOW_SIZE(3);
const int Parameters::DEF_OPEN(-10);
const int Parameters::DEF_GAP(-1);
const size_t Parameters::DEF_MEMORY_WORKER(512U*MB);
const size_t Parameters::DEF_MEMORY_SEQUENCES(2U*GB);
const vector<string> Parameters::DEF_SKIP_PREFIXES;
const bool Parameters::DEF_OUTPUT_ALL(false);
const bool Parameters::DEF_DISTRIBUTE_SEQUENCES(false);
const bool Parameters::DEF_USE_LENGTH_FILTER(true);
const bool Parameters::DEF_USE_ITERATOR(false);


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

    return budget;
}


static bool ends_with(string const &fullString, string const &ending)
{
    bool retval = false;

    if (fullString.length() >= ending.length()) {
        retval = (0 == fullString.compare(
                    fullString.length() - ending.length(),
                    ending.length(),
                    ending));
    }

    return retval;
}



Parameters::Parameters()
    : AOL(DEF_ALIGN_OVER_LONGER_SEQUENCE)
    , SIM(DEF_MATCH_SIMILARITY)
    , OS(DEF_OPTIMAL_SCORE_OVER_SELF_SCORE)
    , exact_match_length(DEF_EXACT_MATCH_LENGTH)
    , window_size(DEF_SLIDE_WINDOW_SIZE)
    , open(DEF_OPEN)
    , gap(DEF_GAP)
    , memory_worker(DEF_MEMORY_WORKER)
    , memory_sequences(DEF_MEMORY_SEQUENCES)
    , skip_prefixes()
    , output_all(DEF_OUTPUT_ALL)
    , distribute_sequences(DEF_DISTRIBUTE_SEQUENCES)
    , use_length_filter(DEF_USE_LENGTH_FILTER)
    , use_iterator(DEF_USE_ITERATOR)
{
}


Parameters::Parameters(const char *parameters_file, MPI_Comm comm)
    : AOL(DEF_ALIGN_OVER_LONGER_SEQUENCE)
    , SIM(DEF_MATCH_SIMILARITY)
    , OS(DEF_OPTIMAL_SCORE_OVER_SELF_SCORE)
    , exact_match_length(DEF_EXACT_MATCH_LENGTH)
    , window_size(DEF_SLIDE_WINDOW_SIZE)
    , open(DEF_OPEN)
    , gap(DEF_GAP)
    , memory_worker(DEF_MEMORY_WORKER)
    , memory_sequences(DEF_MEMORY_SEQUENCES)
    , skip_prefixes()
    , output_all(DEF_OUTPUT_ALL)
    , distribute_sequences(DEF_DISTRIBUTE_SEQUENCES)
    , use_length_filter(DEF_USE_LENGTH_FILTER)
    , use_iterator(DEF_USE_ITERATOR)
{
    parse(parameters_file, comm);
}


void Parameters::parse(const char *parameters_file, MPI_Comm comm)
{
    int comm_rank = 0;  /* rank 0 will open the file */
    int status = 0;     /* for MPI routine return codes */

    status = MPI_Comm_rank(comm, &comm_rank);
    MPI_CHECK_IERR(status, comm_rank, comm);

    if (!ends_with(parameters_file, ".yaml")) {
        cerr << "unrecognized config file format, expecting *.yaml" << endl;
        MPI_Abort(comm, -1);
    }

    if (0 == comm_rank) {
        const YAML::Node config = YAML::LoadFile(parameters_file);
        AOL = config[KEY_ALIGN_OVER_LONGER_SEQUENCE].as<int>(
                DEF_ALIGN_OVER_LONGER_SEQUENCE);
        SIM = config[KEY_MATCH_SIMILARITY].as<int>(
                DEF_MATCH_SIMILARITY);
        OS = config[KEY_OPTIMAL_SCORE_OVER_SELF_SCORE].as<int>(
                DEF_OPTIMAL_SCORE_OVER_SELF_SCORE);
        exact_match_length = config[KEY_EXACT_MATCH_LENGTH].as<int>(
                DEF_EXACT_MATCH_LENGTH);
        window_size = config[KEY_SLIDE_WINDOW_SIZE].as<int>(
                DEF_SLIDE_WINDOW_SIZE);
        open = config[KEY_OPEN].as<int>(
                DEF_OPEN);
        gap = config[KEY_GAP].as<int>(
                DEF_GAP);
        memory_worker = config[KEY_MEMORY_WORKER].as<size_t>(
                DEF_MEMORY_WORKER);
        memory_sequences = config[KEY_MEMORY_SEQUENCES].as<size_t>(
                DEF_MEMORY_SEQUENCES);
        output_all = config[KEY_OUTPUT_ALL].as<bool>(
                DEF_OUTPUT_ALL);
        distribute_sequences = config[KEY_DISTRIBUTE_SEQUENCES].as<bool>(
                DEF_DISTRIBUTE_SEQUENCES);
        use_length_filter = config[KEY_USE_LENGTH_FILTER].as<bool>(
                DEF_USE_LENGTH_FILTER);
        use_iterator = config[KEY_USE_ITERATOR].as<bool>(
                DEF_USE_ITERATOR);

        const YAML::Node filter_node = config[KEY_SKIP_PREFIXES];
        for (size_t i=0; i<filter_node.size(); ++i) {
            const string filter = filter_node[i].as<string>();
            if (filter.size() != ((unsigned)window_size)) {
                cerr << "skip prefix length must match slide window size" << endl;
                cerr << "'" << filter << "' len=" << filter.size() << " SlideWindowSize=" << window_size << endl;
                MPI_Abort(comm, -1);
            }
            skip_prefixes.push_back(filter);
        }
    }

    /* slow, but correct */
    mpix_bcast(AOL, 0, comm);
    mpix_bcast(SIM, 0, comm);
    mpix_bcast(OS, 0, comm);
    mpix_bcast(exact_match_length, 0, comm);
    mpix_bcast(window_size, 0, comm);
    mpix_bcast(open, 0, comm);
    mpix_bcast(gap, 0, comm);
    mpix_bcast(memory_worker, 0, comm);
    mpix_bcast(memory_sequences, 0, comm);
    mpix_bcast(output_all, 0, comm);
    mpix_bcast(distribute_sequences, 0, comm);
    mpix_bcast(use_length_filter, 0, comm);
    mpix_bcast(use_iterator, 0, comm);
    mpix_bcast(skip_prefixes, 0, comm);
}


ostream& operator<< (ostream &os, const Parameters &p)
{
    YAML::Emitter out;
    out << YAML::BeginMap;
    out << YAML::Key << Parameters::KEY_ALIGN_OVER_LONGER_SEQUENCE << YAML::Value << p.AOL;
    out << YAML::Key << Parameters::KEY_MATCH_SIMILARITY << YAML::Value << p.SIM;
    out << YAML::Key << Parameters::KEY_OPTIMAL_SCORE_OVER_SELF_SCORE << YAML::Value << p.OS;
    out << YAML::Key << Parameters::KEY_EXACT_MATCH_LENGTH << YAML::Value << p.exact_match_length;
    out << YAML::Key << Parameters::KEY_SLIDE_WINDOW_SIZE << YAML::Value << p.window_size;
    out << YAML::Key << Parameters::KEY_OPEN << YAML::Value << p.open;
    out << YAML::Key << Parameters::KEY_GAP << YAML::Value << p.gap;
    out << YAML::Key << Parameters::KEY_MEMORY_WORKER << YAML::Value << p.memory_worker;
    out << YAML::Key << Parameters::KEY_MEMORY_SEQUENCES << YAML::Value << p.memory_sequences;
    if (!p.skip_prefixes.empty()) {
        out << YAML::Key << Parameters::KEY_SKIP_PREFIXES << YAML::Value << p.skip_prefixes;
    }
    out << YAML::Key << Parameters::KEY_OUTPUT_ALL << YAML::Value << p.output_all;
    out << YAML::Key << Parameters::KEY_DISTRIBUTE_SEQUENCES << YAML::Value << p.distribute_sequences;
    out << YAML::Key << Parameters::KEY_USE_LENGTH_FILTER << YAML::Value << p.use_length_filter;
    out << YAML::Key << Parameters::KEY_USE_ITERATOR << YAML::Value << p.use_iterator;
    out << YAML::EndMap;
    os << out.c_str();
    return os;
}

}; /* namespace pgraph */

