#! /usr/bin/env python

import sys
from collections import defaultdict
from math import ceil
import operator

try:
    with open(sys.argv[1]) as fd:
        sequences = filter(lambda a : a[0] != '>', [seq.strip() for seq in fd])
        k = int(sys.argv[2])
        errors = True if sys.argv[3][0].lower() == 'e' else False
except:
    print "USAGE: python graph_assembler.py <input_file_path> <kmer_size> <errors : noerrors>"
    sys.exit()

#  ============================== ACCOUNT FOR THE ERROR IN READS ================================ ||
kmers = defaultdict(int)
for seq in sequences:
    for i in xrange(len(seq) - k + 1):
        kmers[seq[i:i+k]] += 1

if errors:
    # we are going to ignore 1% of all least occuring kmers
    amount_to_ignore = int(round(len(kmers) * 0.01))
    kmers = sorted(kmers.items(), key=operator.itemgetter(1))[amount_to_ignore + 1:]

kmer_list = list(zip(*kmers)[0])

#  ============================================================================================== ||

graph = dict()
for kmer in kmer_list:
    if kmer[:-1] not in graph:
        graph[kmer[:-1]] = list()
    graph[kmer[:-1]].append(kmer[1:])


balanced, unbalanced = list(), list()
out_values = reduce(lambda a,b: a+b, graph.values())
for node in (out_values + graph.keys()):
    in_value = 0
    out_value = out_values.count(node)
    if node in graph:
        in_value = len(graph[node])
    
    if in_value == out_value == 1:
        balanced.append(node)
    else:
        unbalanced.append(node)

#==========================================================

def getContigs(node_first, node_second):
    contig = []
    for vertex in graph[node_first]:
        if vertex not in balanced:
            contig.append(node_second + vertex[-1])
        else:
            contig += getContigs(vertex, node_second+vertex[-1])
    return contig

#==========================================================

contig_output = list()
for node in (set(unbalanced) & set(graph.keys())):
    contig_output += sorted(getContigs(node, node))

print " ".join(sorted(contig_output))

