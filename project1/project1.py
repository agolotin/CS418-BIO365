#! /usr/bin/env python

import sys
from skew import minSkew
from freqStrWithMismatch import mainKmerFinder

with open(sys.argv[1]) as fd:
    info = fd.readline().strip()
    dna = "".join([line.strip() for line in fd.readlines()])

#cut down search space by looking for the mininum skew - where dnaa box may be
skew_list = minSkew(dna)
#look for freq words with mismatches
start_pos = int(skew_list[0])
window_size = 500
genome_piece = dna[start_pos : start_pos + window_size]
k = 9
d = 1

print skew_list
print skew_list[0]
print int(skew_list[0])+500
print "\n".join(mainKmerFinder(genome_piece, k, d))
