#include "config.h"

#include <time.h>
#include <sys/time.h>

#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#endif

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <iomanip>

using namespace std;

#include "combinations.h"

static double wtime()
{
#if 0
    struct timeval tv;
    (void)gettimeofday(&tv, NULL);
    return double(tv.tv_sec) + double(tv.tv_usec) / 1000000.0;
#else
    struct timespec ts;
#   if defined(__MACH__) // OS X does not have clock_gettime, use clock_get_time
    clock_serv_t cclock;
    mach_timespec_t mts;
    host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
    clock_get_time(cclock, &mts);
    mach_port_deallocate(mach_task_self(), cclock);
    ts.tv_sec = mts.tv_sec;
    ts.tv_nsec = mts.tv_nsec;
#   else
    // clock_gettime(CLOCK_MONOTONIC, &ts); // Works on FreeBSD
    clock_gettime(CLOCK_REALTIME, &ts); // Works on Linux
#   endif
    return double(ts.tv_sec) + double(ts.tv_nsec) / 1000000000.0;
#endif
}


int main(int argc, char **argv)
{
    unsigned long result[2] = {0,0};
    unsigned long index = 0;
    unsigned long stop = 0;
    double timer = 0.0;

    assert(argc == 2);
    stop = (unsigned long)atol(argv[1]);

    cout << "stop after " << stop << " iterations" << endl;

    timer = wtime();
    for (index=0; index<stop; ++index) {
        k_combination2(index, result);
    }
    timer = wtime() - timer;

    cout << setprecision(6) << timer << endl;

    return EXIT_SUCCESS;
}
