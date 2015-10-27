#! /usr/bin/env python

import sys
import itertools
from reverseComp import reverse
from collections import defaultdict

def mutations(word, hamming_distance, charset='ATCG'):
    # this enumerates all the positions in word
    for indices in itertools.combinations( range( len( word ) ), hamming_distance ):
        for replacements in itertools.product(charset, repeat=hamming_distance):
            mutation = list(word)
            for index, replacement in zip( indices, replacements ):
                mutation[ index ] = replacement
            yield "".join( mutation )

def kmer_finder_mismatches(seq, k, d):
    kmer_freq = defaultdict(int)
    for i in range(len(seq)-k+1):
        kmer_freq[seq[i:i+k]] += 1

    mismatch_count = defaultdict(int) 
    for kmer, freq in kmer_freq.iteritems():
        temp_mutations = set()
        for mismatch in mutations(kmer, d):
            temp_mutations.add(mismatch)
        for mismatch in temp_mutations:
            mismatch_count[mismatch] += freq

        temp_mutations = set()
        for mismatch in mutations(reverse(kmer), d):
            temp_mutations.add(mismatch)
        for mismatch in temp_mutations:
            mismatch_count[mismatch] += freq

    max_count = max(mismatch_count.values())
    return sorted([kmer for kmer, count in mismatch_count.iteritems() if count == max_count])

def mainKmerFinder(_seq, _k, _d):
    kmers = kmer_finder_mismatches(_seq, _k, _d)
    seq_kmers = list()
    for i in xrange(len(_seq) - _k + 1):
        seq_kmers.append(_seq[i : i + _k])
        seq_kmers.append(reverse(_seq[i : i + _k]))

    output = list()
    for kmer in kmers:
        if kmer in seq_kmers:
            output.append(kmer)

    return output
