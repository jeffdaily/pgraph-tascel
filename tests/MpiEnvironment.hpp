#include <fstream>
#include <iostream>
#include <string>

#include <mpi.h>

#include "gtest/gtest.h"

using ::std::ofstream;
using ::std::ostream;
using ::std::streambuf;
using ::std::string;
using ::testing::Environment;

class MpiEnvironment : public Environment
{
    public:
        MpiEnvironment(int *argc, char ***argv);
        virtual ~MpiEnvironment();

        virtual void SetUp();
        virtual void TearDown();

    private:
        string get_program_path() const;
        string get_filename(const string &ext) const;
        streambuf* redirect_ostream(ostream &new_stream, ostream &old_stream);

        int *argc;
        char ***argv;
        MPI_Comm comm;
        int rank;
        int size;
        string program_name;
        ofstream file_out;
        ofstream file_err;
        ofstream file_ggl;
        streambuf *streambuf_out;
        streambuf *streambuf_err;
};
