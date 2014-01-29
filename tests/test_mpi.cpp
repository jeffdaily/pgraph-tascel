#include <unistd.h>

#include <mpi.h>

#include "gtest/gtest.h"

#include "tests/StreamPrettyPrinter.hpp"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

using std::cerr;
using std::cout;
using std::endl;
using std::ofstream;
using std::ostream;
using std::ostringstream;
using std::streambuf;
using std::string;


string program_path()
{
    string path_string;
    char *path = new char[PATH_MAX];
    if (path != NULL) {
        ssize_t count = 0;
        count = readlink("/proc/self/exe", path, PATH_MAX);
        if (count != -1) {
            path_string.assign(path, count);
        }
        delete [] path;
        path = NULL;
    }
    return path_string;
}


string get_filename(const string &ext)
{
    int mpi_retval = 0;
    int mpi_initialized = 0;
    string path;
    int rank = -1;
    ostringstream name;
    size_t last_slash = 0;

    mpi_retval = MPI_Initialized(&mpi_initialized);
    assert(MPI_SUCCESS == mpi_retval);
    assert(mpi_initialized);

    mpi_retval = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    assert(MPI_SUCCESS == mpi_retval);
    assert(rank >= 0);

    path = program_path();
    last_slash = path.find_last_of('/');
    assert(last_slash != string::npos);

    name << path.substr(last_slash+1) << '.' << rank << ext;
    return name.str();
}


streambuf* redirect_ostream(ostream &new_stream, ostream &old_stream)
{
    // Save old buffer here
    streambuf* old_streambuf = old_stream.rdbuf();

    // Redirect cout to our stringstream buffer or any other ostream
    old_stream.rdbuf(new_stream.rdbuf());

    return old_streambuf;
}


#if 0
void restore_stdout()
{
    buffer.close();

    // When done redirect cout to its old self
    cout.rdbuf(sbuf);
}
#endif


TEST(MpiTest, Hello)
{
    int retval = 0;
    int flag = 0;
    string path = program_path();
    int rank = -1;

    retval = MPI_Initialized(&flag);
    EXPECT_EQ(MPI_SUCCESS, retval);
    ASSERT_TRUE(flag);

    retval = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    EXPECT_EQ(MPI_SUCCESS, retval);
    ASSERT_GE(rank, 0);

    cout << '[' << rank << "] Hello from " << path << endl;
}


int main(int argc, char* argv[])
{
    int result = 0;
    streambuf *streambuf_out = NULL;
    streambuf *streambuf_err = NULL;
    streambuf *streambuf_ggl = NULL;
    ofstream file_out;
    ofstream file_err;
    ofstream file_ggl;

    ::testing::InitGoogleTest(&argc, argv);

    MPI_Init(&argc, &argv);

    // Create file to buffer output
    file_out.open(get_filename(".out").c_str());
    file_err.open(get_filename(".err").c_str());
    file_ggl.open(get_filename(".ggl").c_str());

    streambuf_out = redirect_ostream(file_out, cout);
    streambuf_err = redirect_ostream(file_err, cerr);
    ::testing::UnitTest &unit_test = *UnitTest::GetInstance();
    ::testing::TestEventListeners &listeners = unit_test.listeners();
    delete listeners.Release(listeners.default_result_printer());
    listeners.Append(new StreamPrettyPrinter(file_ggl));

    result = RUN_ALL_TESTS();

    file_out.close();
    file_err.close();

    cout.rdbuf(streambuf_out);
    cerr.rdbuf(streambuf_err);

    MPI_Finalize();

    return result;
}

