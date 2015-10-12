#! /usr/bin/env python

import sys

with open(sys.argv[1]) as fd:
    genome = fd.next().strip()

g = 0
c = 0
min_skew = 0
skew_list = list()
for i in range(len(genome)):
    if genome[i] == 'G':
        g += 1
    if genome[i] == 'C':
        c += 1
    skew = g - c
    if skew < min_skew:
        skew_list = [str(i+1)]
        min_skew = skew

    if skew == min_skew and str(i+1) not in skew_list:
        skew_list.append(str(i+1))

print " ".join(skew_list)
