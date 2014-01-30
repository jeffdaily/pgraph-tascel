#include "config.h"

#if HAVE_PROGNAME
extern const char * PROGNAME;
#endif

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include <mpi.h>

#include "gtest/gtest.h"

#include "tests/MpiEnvironment.hpp"
#include "tests/StreamPrettyPrinter.hpp"

using ::std::cerr;
using ::std::cout;
using ::std::ofstream;
using ::std::ostringstream;
using ::std::streambuf;
using ::std::string;
using ::testing::Environment;
using ::testing::UnitTest;
using ::testing::TestEventListeners;

#define PT_PATH_MAX 4096


MpiEnvironment::MpiEnvironment(int *argc, char ***argv)
    :   Environment()
    ,   argc(argc)
    ,   argv(argv)
    ,   comm(MPI_COMM_NULL)
    ,   rank(-1)
    ,   size(-1)
    ,   program_name((*argv)[0])
    ,   file_out()
    ,   file_err()
    ,   file_ggl()
    ,   streambuf_out(NULL)
    ,   streambuf_err(NULL)
{
    string path;
    ostringstream name;
    size_t last_slash = 0;
    int mpi_return_code = 0;

    mpi_return_code = MPI_Init(argc, argv);
    EXPECT_EQ(MPI_SUCCESS, mpi_return_code);

    mpi_return_code = MPI_Comm_dup(MPI_COMM_WORLD, &comm);
    EXPECT_EQ(MPI_SUCCESS, mpi_return_code);

    mpi_return_code = MPI_Comm_rank(comm, &rank);
    EXPECT_EQ(MPI_SUCCESS, mpi_return_code);

    mpi_return_code = MPI_Comm_size(comm, &size);
    EXPECT_EQ(MPI_SUCCESS, mpi_return_code);
    
    path = get_program_path();
    last_slash = path.find_last_of('/');
    if (last_slash != string::npos) {
        program_name = path.substr(last_slash+1);
    }
    /* else {
     *  program_name was already set in initializer list using argv[0]
    } */

    file_out.open(get_filename(".out").c_str());
    file_err.open(get_filename(".err").c_str());
    file_ggl.open(get_filename(".ggl").c_str());

    streambuf_out = redirect_ostream(file_out, cout);
    streambuf_err = redirect_ostream(file_err, cerr);

    UnitTest &unit_test = *UnitTest::GetInstance();
    TestEventListeners &listeners = unit_test.listeners();
    delete listeners.Release(listeners.default_result_printer());
    listeners.Append(new StreamPrettyPrinter(file_ggl));
}


MpiEnvironment::~MpiEnvironment()
{
    int mpi_return_code = 0;

    file_out.close();
    file_err.close();
    file_ggl.close();

    cout.rdbuf(streambuf_out);
    cerr.rdbuf(streambuf_err);

    mpi_return_code = MPI_Finalize();
    EXPECT_EQ(MPI_SUCCESS, mpi_return_code);
}


void MpiEnvironment::SetUp()
{
}


void MpiEnvironment::TearDown()
{
}


string MpiEnvironment::get_program_path() const
{
    string path_string;

    /* try /proc/self/exe first */
    if (path_string.empty()) {
        char *path = new char[PT_PATH_MAX];

        if (path != NULL) {
            ssize_t count = 0;

            count = readlink("/proc/self/exe", path, PT_PATH_MAX);
            if (count != -1) {
                path_string.assign(path, count);
            }

            delete [] path;
            path = NULL;
        }
    }

    /* try /proc/<pid>/exe next */
    if (path_string.empty()) {
        char *path = new char[PT_PATH_MAX];

        if (path != NULL) {
            ostringstream link_name;
            ssize_t count = 0;

            link_name << "/proc/" << getpid() << "/exe";
            count = readlink(link_name.str().c_str(), path, PT_PATH_MAX);
            if (count != -1) {
                path_string.assign(path, count);
            }

            delete [] path;
            path = NULL;
        }
    }

#if HAVE_PROGNAME
    if (path_string.empty()) {
        path_string = PROGNAME;
    }
#endif

    return path_string;
}


string MpiEnvironment::get_filename(const string &ext) const
{
    ostringstream name;
    name << program_name << '.' << rank << ext;
    return name.str();
}


streambuf* MpiEnvironment::redirect_ostream(ostream &new_stream, ostream &old_stream)
{
    // Save old buffer here
    streambuf* old_streambuf = old_stream.rdbuf();

    // Redirect cout to our stringstream buffer or any other ostream
    old_stream.rdbuf(new_stream.rdbuf());

    return old_streambuf;
}

