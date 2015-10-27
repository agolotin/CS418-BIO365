#! /usr/bin/env python

import sys
from problm22 import eulerian_path as eulerPath

with open(sys.argv[1]) as fd:
    k = int(fd.readline().strip())
    kmer_list = sorted([seq.strip() for seq in fd.readlines()])


graph = dict()
for kmer in kmer_list:
    if kmer[0:k-1] not in graph:
        graph[kmer[0:k-1]] = list()
    graph[kmer[0:k-1]].append(kmer[1:k])

path = eulerPath(graph)

output = ""
for i in xrange(len(path)):
    if i == 0:
        output += path[i]
    else:
        output += path[i][k-2:]

fd = open('output.txt', 'w')
fd.write(output)
fd.close()
