#!/usr/bin/env python

import itertools
import sys

if len(sys.argv) != 3:
    print "gen_bucket.py num k"
    sys.exit()

ID = int(sys.argv[1])
k = int(sys.argv[2])
E_SIGMA = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
SIGMA = len(E_SIGMA)
val = ""
for i in range(k-1,-1,-1):
    tmp = SIGMA**i
    quo = ID / tmp
    rem = ID % tmp
    val += chr(ord('A')+quo)
    ID = rem

print val
