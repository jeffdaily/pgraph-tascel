#include "config.h"

#include <mpi.h>

#include "gtest/gtest.h"

#include "tests/MpiEnvironment.hpp"
#include "tests/StreamPrettyPrinter.hpp"

using std::cerr;
using std::cout;
using std::endl;


TEST(MpiTest, Hello)
{
    int retval = 0;
    int flag = 0;
    int rank = -1;

    retval = MPI_Initialized(&flag);
    EXPECT_EQ(MPI_SUCCESS, retval);
    ASSERT_TRUE(flag);

    retval = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    EXPECT_EQ(MPI_SUCCESS, retval);
    ASSERT_GE(rank, 0);

    cout << '[' << rank << "] Hello stdout" << endl;
    cerr << '[' << rank << "] Hello stderr" << endl;
}


int main(int argc, char* argv[])
{
    int result = 0;

    ::testing::InitGoogleTest(&argc, argv);

    ::testing::Environment* const env = ::testing::AddGlobalTestEnvironment(
            new MpiEnvironment(&argc, &argv));

    result = RUN_ALL_TESTS();

    return result;
}

