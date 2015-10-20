#! /usr/bin/env python

import sys
from collections import defaultdict
from math import floor

try:
    with open(sys.argv[1]) as fd:
        sequences = filter(lambda a : a[0] != '>', [seq.strip() for seq in fd])
        k = int(sys.argv[2])
except:
    print "USAGE: python graph_assembler.py <file_path> <kmer size>"
    sys.exit()

#  ============================== ACCOUNT FOR THE ERROR IN READS ================================ ||
kmers = defaultdict(int)
for seq in sequences:
    for i in xrange(len(seq) - k + 1):
        kmers[seq[i:i+k]] += 1

max_occurences = max(kmers.values()) # get the maximum number of occurences of a kmer

# we will ignore only 1% of all of the values
values_to_ignore = [int(floor(float(max_occurences) * 0.01))]
if values_to_ignore[0] > 1:
    values_to_ignore = range(1, values_to_ignore[0] + 1)
# filter the key-value pairs to match only the ones we need 
kmers = filter(lambda pair: pair[1] not in values_to_ignore, kmers.iteritems())
# now we are going to create a list of only kmers out of our kmer:num_occurences dictionary of kmers
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

