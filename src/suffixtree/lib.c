/**
 * @file lib.c
 *
 * @author andy.cj.wu@gmail.com
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2010 Washington State University. All rights reserved.
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 */
#include <assert.h>
#include <limits.h>
#include <stdint.h>
#include <stddef.h>

#include "lib.h"

int power(int base, int n)
{
    int p;

    for (p = 1; n > 0; --n) {
        assert(p < INT_MAX / base);
        p *= base;
    }

    return p;
}

#ifndef SIZE_MAX
#define SIZE_MAX ((size_t)-1)
#endif

size_t zpower(size_t base, size_t n)
{
    size_t p;

    for (p = 1; n > 0; --n) {
        assert(p < SIZE_MAX / base);
        p *= base;
    }

    return p;
}
