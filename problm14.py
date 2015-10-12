#! /usr/bin/env python

import sys
from operator import mul

#with open(sys.argv[1]) as fd:
#    dna = fd.readline().strip()
#    k = int(fd.readline().strip())
    profile_matrix = [[float(val.strip()) for val in line.split(' ')] for line in fd.readlines()]


def probability(_kmer, profile_matrix):
    s = list()
    for i in range(len(_kmer)):
        if _kmer[i] == 'A':
            s.append(profile_matrix[0][i])
        elif _kmer[i] == 'C':
            s.append(profile_matrix[1][i])
        elif _kmer[i] == 'G':
            s.append(profile_matrix[2][i])
        elif _kmer[i] == 'T':
            s.append(profile_matrix[3][i])

    return reduce(mul, s, 1)

def profile_most_probable_kmer(dna, k, profile_matrix):
    best_pattern = ""
    best_probability = float(0)

    for i in range(len(dna) - k + 1):
        kmer = dna[i:i+k]
        probabil = probability(kmer, profile_matrix)
        if probabil > best_probability:
            best_pattern = kmer
            best_probability = probabil
    return best_pattern

