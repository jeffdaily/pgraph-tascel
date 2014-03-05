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
const string Parameters::KEY_OUTPUT_TO_DISK("OutputToDisk");
const string Parameters::KEY_DISTRIBUTE_SEQUENCES("DistributeSequences");
const string Parameters::KEY_USE_LENGTH_FILTER("UseLengthFilter");
const string Parameters::KEY_USE_ITERATOR("UseIterator");
const string Parameters::KEY_USE_COUNTER("UseCounter");
const string Parameters::KEY_USE_TREE("UseTree");
const string Parameters::KEY_USE_TREE_DYNAMIC("UseTreeDynamic");
const string Parameters::KEY_USE_TREE_HYBRID("UseTreeHybrid");
const string Parameters::KEY_PRINT_STATS("PrintStats");
const string Parameters::KEY_ALPHABET("Alphabet");
const string Parameters::KEY_ALPHABET_BEGIN("AlphabetBegin");
const string Parameters::KEY_ALPHABET_DOLLAR("AlphabetDollar");
const string Parameters::KEY_DUPLICATES_LOCAL("DuplicatesLocal");
const string Parameters::KEY_DUPLICATES_SEMILOCAL("DuplicatesSemiLocal");
const string Parameters::KEY_DUPLICATES_SMP("DuplicatesSmp");
const string Parameters::KEY_DUPLICATES_GLOBAL("DuplicatesGlobal");
/* Defaults */
const int Parameters::DEF_ALIGN_OVER_LONGER_SEQUENCE(80);
const int Parameters::DEF_MATCH_SIMILARITY(40);
const int Parameters::DEF_OPTIMAL_SCORE_OVER_SELF_SCORE(30);
const int Parameters::DEF_EXACT_MATCH_LENGTH(4);
const int Parameters::DEF_SLIDE_WINDOW_SIZE(3);
const int Parameters::DEF_OPEN(-10);
const int Parameters::DEF_GAP(-1);
const size_t Parameters::DEF_MEMORY_WORKER(512U*MB);
const size_t Parameters::DEF_MEMORY_SEQUENCES(2U*GB);
const vector<string> Parameters::DEF_SKIP_PREFIXES;
const bool Parameters::DEF_OUTPUT_ALL(false);
const bool Parameters::DEF_OUTPUT_TO_DISK(true);
const bool Parameters::DEF_DISTRIBUTE_SEQUENCES(false);
const bool Parameters::DEF_USE_LENGTH_FILTER(true);
const bool Parameters::DEF_USE_ITERATOR(false);
const bool Parameters::DEF_USE_COUNTER(false);
const bool Parameters::DEF_USE_TREE(false);
const bool Parameters::DEF_USE_TREE_DYNAMIC(false);
const bool Parameters::DEF_USE_TREE_HYBRID(false);
const bool Parameters::DEF_PRINT_STATS(true);
const string Parameters::DEF_ALPHABET("ABCDEFGHIJKLMNOPQRSTUVWXYZ");
const char Parameters::DEF_ALPHABET_BEGIN('O');
const char Parameters::DEF_ALPHABET_DOLLAR('U');
const bool Parameters::DEF_DUPLICATES_LOCAL(true);
const bool Parameters::DEF_DUPLICATES_SEMILOCAL(false);
const bool Parameters::DEF_DUPLICATES_SMP(false);
const bool Parameters::DEF_DUPLICATES_GLOBAL(false);


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
    , output_to_disk(DEF_OUTPUT_TO_DISK)
    , distribute_sequences(DEF_DISTRIBUTE_SEQUENCES)
    , use_length_filter(DEF_USE_LENGTH_FILTER)
    , use_iterator(DEF_USE_ITERATOR)
    , use_counter(DEF_USE_COUNTER)
    , use_tree(DEF_USE_TREE)
    , use_tree_dynamic(DEF_USE_TREE_DYNAMIC)
    , use_tree_hybrid(DEF_USE_TREE_HYBRID)
    , print_stats(DEF_PRINT_STATS)
    , alphabet(DEF_ALPHABET)
    , alphabet_begin(DEF_ALPHABET_BEGIN)
    , alphabet_dollar(DEF_ALPHABET_DOLLAR)
    , dup_local(DEF_DUPLICATES_LOCAL)
    , dup_semilocal(DEF_DUPLICATES_SEMILOCAL)
    , dup_smp(DEF_DUPLICATES_SMP)
    , dup_global(DEF_DUPLICATES_GLOBAL)
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
    , output_to_disk(DEF_OUTPUT_TO_DISK)
    , distribute_sequences(DEF_DISTRIBUTE_SEQUENCES)
    , use_length_filter(DEF_USE_LENGTH_FILTER)
    , use_iterator(DEF_USE_ITERATOR)
    , use_counter(DEF_USE_COUNTER)
    , use_tree(DEF_USE_TREE)
    , use_tree_dynamic(DEF_USE_TREE_DYNAMIC)
    , use_tree_hybrid(DEF_USE_TREE_HYBRID)
    , print_stats(DEF_PRINT_STATS)
    , alphabet(DEF_ALPHABET)
    , alphabet_begin(DEF_ALPHABET_BEGIN)
    , alphabet_dollar(DEF_ALPHABET_DOLLAR)
    , dup_local(DEF_DUPLICATES_LOCAL)
    , dup_semilocal(DEF_DUPLICATES_SEMILOCAL)
    , dup_smp(DEF_DUPLICATES_SMP)
    , dup_global(DEF_DUPLICATES_GLOBAL)
{
    parse(parameters_file, comm);
}


void Parameters::parse(const char *parameters_file, MPI_Comm comm)
{
    int comm_rank = 0;  /* rank 0 will open the file */
    int status = 0;     /* for MPI routine return codes */
    char *file_buffer = NULL;
    MPI_Offset file_size = 0;
    long chunk_size = 1073741824;

    comm_rank = mpix_rank(comm);

    if (!ends_with(parameters_file, ".yaml")) {
        cerr << "unrecognized config file format, expecting *.yaml" << endl;
        MPI_Abort(comm, -1);
    }

    mpix_read_file_bcast(parameters_file, file_buffer, file_size, chunk_size, comm);

    try {
        const YAML::Node config = YAML::Load(file_buffer);
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
        output_all = config[KEY_OUTPUT_ALL].as<bool>(
                DEF_OUTPUT_ALL);
        output_to_disk = config[KEY_OUTPUT_TO_DISK].as<bool>(
                DEF_OUTPUT_TO_DISK);
        distribute_sequences = config[KEY_DISTRIBUTE_SEQUENCES].as<bool>(
                DEF_DISTRIBUTE_SEQUENCES);
        use_length_filter = config[KEY_USE_LENGTH_FILTER].as<bool>(
                DEF_USE_LENGTH_FILTER);
        use_iterator = config[KEY_USE_ITERATOR].as<bool>(
                DEF_USE_ITERATOR);
        use_counter = config[KEY_USE_COUNTER].as<bool>(
                DEF_USE_COUNTER);
        use_tree = config[KEY_USE_TREE].as<bool>(
                DEF_USE_TREE);
        use_tree_dynamic = config[KEY_USE_TREE_DYNAMIC].as<bool>(
                DEF_USE_TREE_DYNAMIC);
        use_tree_hybrid = config[KEY_USE_TREE_HYBRID].as<bool>(
                DEF_USE_TREE_HYBRID);
        print_stats = config[KEY_PRINT_STATS].as<bool>(
                DEF_PRINT_STATS);
        dup_local = config[KEY_DUPLICATES_LOCAL].as<bool>(
                DEF_DUPLICATES_LOCAL);
        dup_semilocal = config[KEY_DUPLICATES_SEMILOCAL].as<bool>(
                DEF_DUPLICATES_SEMILOCAL);
        dup_smp = config[KEY_DUPLICATES_SMP].as<bool>(
                DEF_DUPLICATES_SMP);
        dup_global = config[KEY_DUPLICATES_GLOBAL].as<bool>(
                DEF_DUPLICATES_GLOBAL);

        string val;
        val = config[KEY_MEMORY_WORKER].as<string>("");
        if (val.empty()) {
            memory_worker = DEF_MEMORY_WORKER;
        }
        else {
            memory_worker = parse_memory_budget(val);
        }

        val = config[KEY_MEMORY_SEQUENCES].as<string>("");
        if (val.empty()) {
            memory_sequences = DEF_MEMORY_SEQUENCES;
        }
        else {
            memory_sequences = parse_memory_budget(val);
        }

        const YAML::Node filter_node = config[KEY_SKIP_PREFIXES];
        for (size_t i=0; i<filter_node.size(); ++i) {
            const string filter = filter_node[i].as<string>();
            skip_prefixes.push_back(filter);
        }
    } catch (YAML::Exception &e) {
        cerr << "unable to parse YAML file" << endl;
        cerr << e.what() << endl;
        MPI_Abort(comm, -1);
    }

#if 0
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
    mpix_bcast(skip_prefixes, 0, comm);
    mpix_bcast(output_all, 0, comm);
    mpix_bcast(output_to_disk, 0, comm);
    mpix_bcast(distribute_sequences, 0, comm);
    mpix_bcast(use_length_filter, 0, comm);
    mpix_bcast(use_iterator, 0, comm);
    mpix_bcast(use_counter, 0, comm);
    mpix_bcast(use_tree, 0, comm);
    mpix_bcast(use_tree_dynamic, 0, comm);
    mpix_bcast(use_tree_hybrid, 0, comm);
    mpix_bcast(print_stats, 0, comm);
    mpix_bcast(alphabet, 0, comm);
    mpix_bcast(alphabet_begin, 0, comm);
    mpix_bcast(alphabet_dollar, 0, comm);
    mpix_bcast(dup_local, 0, comm);
    mpix_bcast(dup_semilocal, 0, comm);
    mpix_bcast(dup_smp, 0, comm);
    mpix_bcast(dup_global, 0, comm);
#endif
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
    out << YAML::Key << Parameters::KEY_OUTPUT_TO_DISK << YAML::Value << p.output_to_disk;
    out << YAML::Key << Parameters::KEY_DISTRIBUTE_SEQUENCES << YAML::Value << p.distribute_sequences;
    out << YAML::Key << Parameters::KEY_USE_LENGTH_FILTER << YAML::Value << p.use_length_filter;
    out << YAML::Key << Parameters::KEY_USE_ITERATOR << YAML::Value << p.use_iterator;
    out << YAML::Key << Parameters::KEY_USE_COUNTER << YAML::Value << p.use_counter;
    out << YAML::Key << Parameters::KEY_USE_TREE << YAML::Value << p.use_tree;
    out << YAML::Key << Parameters::KEY_USE_TREE_DYNAMIC << YAML::Value << p.use_tree_dynamic;
    out << YAML::Key << Parameters::KEY_USE_TREE_HYBRID << YAML::Value << p.use_tree_hybrid;
    out << YAML::Key << Parameters::KEY_PRINT_STATS << YAML::Value << p.print_stats;
    out << YAML::Key << Parameters::KEY_ALPHABET << YAML::Value << p.alphabet;
    out << YAML::Key << Parameters::KEY_ALPHABET_BEGIN << YAML::Value << p.alphabet_begin;
    out << YAML::Key << Parameters::KEY_ALPHABET_DOLLAR << YAML::Value << p.alphabet_dollar;
    out << YAML::Key << Parameters::KEY_DUPLICATES_LOCAL << YAML::Value << p.dup_local;
    out << YAML::Key << Parameters::KEY_DUPLICATES_SEMILOCAL << YAML::Value << p.dup_semilocal;
    out << YAML::Key << Parameters::KEY_DUPLICATES_SMP << YAML::Value << p.dup_smp;
    out << YAML::Key << Parameters::KEY_DUPLICATES_GLOBAL << YAML::Value << p.dup_global;
    out << YAML::EndMap;
    os << out.c_str();
    return os;
}

}; /* namespace pgraph */

