#! /usr/bin/env python

from eulerian_cycle import eulerian_cycle as Cycle
import sys

nodes = dict()
with open(sys.argv[1]) as fd:
    for key, value in [line.strip().split(" -> ") for line in fd.readlines()]:
        nodes[int(key)] = map(int, value.split(","))

def eulerian_path(edge_dict):
    out_values = reduce(lambda a,b: a+b, edge_dict.values())
    for node in set(out_values+edge_dict.keys()):
        out_value = out_values.count(node)
        if node in edge_dict:
            in_value = len(edge_dict[node])
        else:
            in_value = 0

        if in_value < out_value:
            unbalanced_from = node
        elif out_value < in_value:
            unbalanced_to = node

    if unbalanced_from in edge_dict:
        edge_dict[unbalanced_from].append(unbalanced_to)
    else:
        edge_dict[unbalanced_from] = [unbalanced_to]

    cycle = Cycle(edge_dict)
    divide_point = filter(lambda i: cycle[i:i+2] == [unbalanced_from, unbalanced_to], xrange(len(cycle)-1))[0]

    return cycle[divide_point+1:]+cycle[1:divide_point+1]


path = eulerian_path(nodes)
fd = open('output.txt', 'w')
fd.write("->".join(map(str, path)))
fd.close()
