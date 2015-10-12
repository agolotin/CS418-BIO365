#! /usr/bin/env python

import sys

with open(sys.argv[1]) as fd:
    pattern = fd.next().strip()
    genome = fd.next().strip()

output = list()
k = len(pattern)
for i in range(len(genome)):
    ptrn = genome[i : i + k]
    if ptrn == pattern:
        output.append(str(i))

print " ".join(output)
