#! /usr/bin/env python

import sys

with open(sys.argv[1]) as fd:
    k = int(fd.next().strip())
    dna = fd.next().strip()

kmers = list()
for i in xrange(len(dna)- k+1):
    kmers.append(dna[i:i+k])
    

fd = open('output.txt', 'w')
fd.write("\n".join(sorted(kmers)))
fd.close()
