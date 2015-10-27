#! /usr/bin/env python

import sys

with open(sys.argv[1]) as input:
	seq = input.next().strip()
	k = int(input.next().strip())

counts = dict()
highest_count  = 0
for i in range(len(seq)):
	kmer = seq[i : i + k]
	if kmer in counts:
		counts[kmer] += 1
	else:
		counts[kmer] = 1

	if counts[kmer] > highest_count:
		highest_count = counts[kmer]

highest_freq_kmers = []
for key, val in counts.iteritems():
	if val == highest_count:
		highest_freq_kmers.append(key)

print " ".join(sorted(highest_freq_kmers))
