#!/usr/bin/env python

import sys

nodes = dict()
with open(sys.argv[1]) as fd:
    for key, value in [line.strip().split(" -> ") for line in fd.readlines()]:
        nodes[int(key)] = map(int, value.split(","))

def buildGraph(cnode, graph, cycle):
    cycle.append(cnode)
    if len(graph[cnode]) == 0:
        return cycle
    else:
        while len(graph[cnode]) != 0:
            n = graph[cnode][0]
            graph[cnode].remove(n)
            t = buildGraph(n, graph, [])
            cycle = cycle[:1] + t + cycle[1:]
        return cycle 

path = buildGraph(nodes.keys()[0], nodes, [])
fd = open('output.txt', 'w')
fd.write("->".join(map(str, path)))
fd.close()
