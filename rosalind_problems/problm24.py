#! /usr/bin/env python

import sys

with open(sys.argv[1]) as fd:
    kmer_list = [seq.strip() for seq in fd.readlines()]

# Get the graph from list of kmers
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
    contig_output += getContigs(node, node)

print " ".join(sorted(contig_output))

