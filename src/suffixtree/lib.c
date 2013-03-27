#include "lib.h"

#pragma mta parallel off
/* int power function */
int power(int base, int n){
    int p;
    
    for(p = 1; n > 0; n--){
        assert(p < INT_MAX/base);
        p *= base;
    }

    return p;
}

