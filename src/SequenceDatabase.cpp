/**
 * @file SequenceDatabase.cc
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 */
#include <sstream>

#include "Sequence.hpp"
#include "SequenceDatabase.hpp"

using std::ostringstream;


SequenceDatabaseException::SequenceDatabaseException(
    const char *file,
    int line,
    const char *function,
    const char *message) throw()
{
    ostringstream oss;
    oss << file << ": " << function << ": " << line << ": " << message;
    this->message = oss.str();
}


SequenceDatabaseException::~SequenceDatabaseException() throw()
{

}


const char *SequenceDatabaseException::what() const throw()
{
    return this->message.c_str();
}



SequenceDatabase::SequenceDatabase()
{
}


SequenceDatabase::~SequenceDatabase()
{
}

