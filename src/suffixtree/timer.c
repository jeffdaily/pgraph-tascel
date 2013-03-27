#include <sys/mta_task.h>

double timer(){
  return ((double)mta_get_clock(0) / (double)mta_clock_freq());
}
