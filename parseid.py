#!/usr/bin/env python

import sys

if len(sys.argv) != 2:
    print "parseid.py file"
    sys.exit()

E_SIGMA = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
SIGMA = len(E_SIGMA)

def entry_string(index, k):
    val = ""
    for i in range(k-1,-1,-1):
        tmp = SIGMA**i
        quo = index / tmp
        rem = index % tmp
        val += chr(ord('A')+quo)
        index = rem
    return val

for line in open(sys.argv[1]):
    line = line.strip()
    parts = line.split()
    for part in parts:
        if 'id=' in part:
            key = entry_string(int(part[3:]), 6)
            if 'X' not in key:
                print part, key, "<---------- UH OH"
            else:
                print part, key
    
