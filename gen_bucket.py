#!/usr/bin/env python

import itertools
import sys

if len(sys.argv) != 2:
    print "gen_bucket.py str"
    sys.exit()

ID = sys.argv[1]
k = len(ID)
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

print ID,entry_index(ID,k)
