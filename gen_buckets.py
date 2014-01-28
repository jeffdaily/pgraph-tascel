#!/usr/bin/env python

import itertools
import sys

if len(sys.argv) != 2:
    print "gen_buckets.py window_length"
    sys.exit()

k = int(sys.argv[1])
E_SIGMA = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
SIGMA = len(E_SIGMA)

def entry_index(kmer, k):
    value = 0;
    for i in range(k):
        tmp = ord(kmer[i]) - ord('A');
        assert(tmp >= 0);
        assert(tmp < SIGMA);
        value = value * SIGMA + tmp;
    return value;

iter = (''.join(i) for i in itertools.product(E_SIGMA,repeat=k))
for i in iter:
    print i, entry_index(i,k)
