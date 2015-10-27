#! /usr/bin/env python

import sys
from collections import OrderedDict

#with open(sys.argv[1]) as fd:
#    seqs = sorted([seq.strip() for seq in fd.readlines()])

def deBrujinGraph(seqs):
    graph = dict()#OrderedDict()

    for kmer in seqs:
            k = len(kmer)
            if kmer[0:k-1] not in graph:
                    graph[kmer[0:k-1]] = list()
            graph[kmer[0:k-1]].append(kmer[1:k])

            for k in seqs:
                    if kmer[1:] == k[1:]:
                            if kmer[1:] not in graph[kmer[:-1]]:
                                    graph[kmer[:-1]].append(kmer[1:])
    
    return graph


#fd = open('output.txt', 'w')
#for key, value in graph.iteritems():
#    value = [",".join(value)]
#    value.insert(0, key)
#    fd.write(" -> ".join(value))
#    fd.write("\n")
#fd.close()
