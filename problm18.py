#! /usr/bin/env python

import sys

with open(sys.argv[1]) as fd:
    dna = [line.strip() for line in fd.readlines()]

output = dict()
for seq1 in dna:
    for seq2 in dna:
        if seq1[1:len(seq1)] == seq2[0:len(seq2)-1]:
            output[seq1] = seq2

fd = open('output.txt', 'w')
for seqs in output.iteritems():
    fd.write(" -> ".join(seqs))
    fd.write("\n")
fd.close()
