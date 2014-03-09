#include "config.h"

#include <iostream>

#include <mpi.h>

#include "Bootstrap.hpp"
#include "Parameters.hpp"

using ::std::cout;
using ::std::endl;
using ::pgraph::Parameters;

int main(int argc, char **argv)
{
    Parameters parameters;

    pgraph::initialize(argc, argv);

    cout << parameters << endl;

    if (argc == 2) {
        parameters.parse(argv[1], MPI_COMM_WORLD);
        cout << endl;
        cout << "Parsing " << argv[1] << " from command line" << endl;
        cout << endl;
        cout << parameters << endl;
    }

    pgraph::finalize();

    return 0;
}
