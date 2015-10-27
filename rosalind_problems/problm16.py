#! /usr/bin/env python

import sys
from most_probable_kmer import profile_most_probable_kmer
from problm15 import score

with open(sys.argv[1]) as fd:
    k, t = map(int, fd.readline().strip().split(' '))
    dna = [line.strip() for line in fd.readlines()]

def profile(_motifs):
    columns = ["".join(seq) for seq in zip(*_motifs)]
    return [[float(col.count(nuc) + 1) / float(len(col) + 4) for nuc in "ACGT"] for col in columns]

def split(seq, length):
    return [seq[i:i+length] for i in xrange(len(seq) - length + 1)]

def greedyMotifSearch(_dna, _k, _t):
    bestScore = _k * _t
    bestMotifs = list()
    motifs = list()

    for motif in split(_dna[0], _k):
        motifs = [motif]
        for i in xrange(1, _t):
            current_profile = profile(motifs)
            motifs.append(profile_most_probable_kmer(_dna[i], k, current_profile))

        curr_score = score(motifs)
        if curr_score < bestScore:
            bestScore = curr_score
            bestMotifs = motifs

    return bestMotifs

print "\n".join(greedyMotifSearch(dna, k, t))
