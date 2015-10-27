#! /usr/bin/env python

import sys
import itertools

with open(sys.argv[1]) ad fd:
    K, D = map(int, fd.readline().strip().split())
    dna_list = [dna.strip() for dna in fd.readlines()]

def mutations(word, hamming_distance, charset='ATCG'):
    # this enumerates all the positions in word
    for indices in itertools.combinations( range( len( word ) ), hamming_distance ):
        for replacements in itertools.product(charset, repeat=hamming_distance):
            mutation = list(word)
            for index, replacement in zip( indices, replacements ):
                mutation[ index ] = replacement
            yield "".join( mutation )

def motifEnumeration(_dna_list, k, d):
    motif_sets = list()
    for dna in _dna_list:
        patterns = list()
        for i in xrange(len(dna) - k+1):
            for kmer in mutations(dna[i:i+k], d):
                patterns.append(kmer)
        motif_sets.append(set(patterns))

    lambda_func = lambda a,b: a & b
    return sorted(list(reduce(lambda_func, motif_sets)))


print " ".join(motifEnumeration(dna_list, K, D))
