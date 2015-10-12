#! /usr/bin/env python

import sys
import regex

with open(sys.argv[1]) as fd:
    pattern = fd.next().strip()
    genome = fd.next().strip()
    d = int(fd.next().strip())

output = list()
for i in range(len(genome) - len(pattern)):
    subseq = genome[i : i + len(pattern)]
    prsd_pttrn = regex.compile("(?:" + pattern + ")" + 
            "{i<1,d<1,s<=" + str(d) + "}", regex.BESTMATCH) #allow only substitution
    if prsd_pttrn.fullmatch(subseq):
#        print prsd_pttrn.fullmatch(subseq).fuzzy_counts
        output.append(str(i))

print " ".join(output)
