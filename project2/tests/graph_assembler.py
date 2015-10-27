#! /usr/bin/env python

import sys
from collections import defaultdict
import operator


def buildKmers(sequences, k):
    kmers = defaultdict(int)
    for seq in sequences:
        for i in xrange(len(seq) - k + 1):
            kmers[seq[i:i+k]] += 1
    return kmers

#  ============================== ACCOUNT FOR THE ERROR IN READS ================================ ||

def checkForErrors(kmers, errors, percentline):
    ''' If there are errors in reads the function will get rid of 
        least occuring kmers based on the number of
		a single kmer occuring '''
    kmer_list = list()

    if errors:
        # get the maximum number of occurences of a kmer
        max_occurences = max(kmers.values()) 
        values_to_ignore = [int(float(max_occurences) * percentline)]
        if values_to_ignore[0] > 1:
            values_to_ignore = range(1, values_to_ignore[0] + 1)
        kmers = filter(lambda pair: (pair[1] not in values_to_ignore), kmers.iteritems())
        kmer_list = list(zip(*kmers)[0])
    else:
        kmer_list = list(kmers.keys())

    return kmer_list

#  ============================================================================================== ||

def buildDeBrujinGraph(kmer_list):
    graph = dict()
    for kmer in kmer_list:
        if kmer[:-1] not in graph:
            graph[kmer[:-1]] = list()
        graph[kmer[:-1]].append(kmer[1:])
    return graph


def getEdgeBalance(graph):
    balanced, unbalanced = list(), list()
    out_values = reduce(lambda a,b: a+b, graph.values())
    for node in set((out_values + graph.keys())):
        in_value = 0
        out_value = out_values.count(node)
        if node in graph:
            in_value = len(graph[node])
        
        if in_value == out_value == 1:
            balanced.append(node)
        else:
            unbalanced.append(node)

    return balanced, unbalanced


def getContigs(node_first, node_second, balanced):
    contig = []
    for vertex in graph[node_first]:
        if vertex not in balanced:
            contig.append(node_second + vertex[-1])
        else:
            contig += getContigs(vertex, node_second+vertex[-1], balanced)
    return contig


#==========================================================

def printOutput(contig_output):
    total = 0
    total_length = 0
    largest_size = -1
    for contig in contig_output:
        total += 1
        total_length += len(contig)
        if len(contig) > largest_size:
            largest_size = len(contig)

        print ">Contig" + str(total) + "_length" + str(len(contig))
        print contig
    print
    print "Average contig size: " + str(float(total_length) / float(total))
    print "Number of contigs returned: " + str(total)
    print "Largest contig size: " + str(largest_size)

if __name__ == "__main__":
    sys.setrecursionlimit(1000000)
    #load the parameters
    try:
        with open(sys.argv[1]) as fd:
            sequences = filter(lambda a : a[0] != '>', [seq.strip() for seq in fd])
            k = int(sys.argv[2])
            errors = True if sys.argv[3][0].lower() == 'e' else False
            percentline = float(sys.argv[4]) if len(sys.argv) > 4 else 0.1
    except:
        print "USAGE: python graph_assembler.py <input_file_path> <kmer_size> <errors : noerrors> optional:<percentline>"
        sys.exit()
    
    ''' modify recursion depth '''
    kmers = buildKmers(sequences, k)
    kmer_list = checkForErrors(kmers, errors, percentline)
    graph = buildDeBrujinGraph(kmer_list)
    balanced, unbalanced = getEdgeBalance(graph)

    contig_output = list()
    for node in (set(unbalanced) & set(graph.keys())):
        contig_output += sorted(getContigs(node, node, balanced))

    printOutput(contig_output)
