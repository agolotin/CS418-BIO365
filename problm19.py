#! /usr/bin/env python

import sys
from collections import OrderedDict

with open(sys.argv[1]) as fd:
    k = int(fd.next().strip())
    dna = fd.next().strip()


graph = dict()#OrderedDict()
kmer_list = list()
for i in xrange(0, len(dna) - k + 1):
    kmer_list.append(dna[i:i+k])


for kmer in kmer_list:
    if kmer[0:k-1] not in graph:
        graph[kmer[0:k-1]] = list()
    graph[kmer[0:k-1]].append(kmer[1:k])

fd = open('output.txt', 'w')
for key, value in graph.iteritems():
    value = [",".join(value)]
    value.insert(0, key)
    fd.write(" -> ".join(value))
    fd.write("\n")
fd.close()
