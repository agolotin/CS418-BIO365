#! /usr/bin/env python

import sys
from Trie import Trie 
from operator import sub

with open(sys.argv[1]) as fd:
    seqs = [line.strip() for line in fd.readlines()]

trie = Trie(seqs)

output = list()
for item in trie.edges.items():
    output.append("->".join(map(str, map(sub, item[0], (1,1)))) + ":" + item[1])

def getKey(item):
    return item[3] "1->2:A"

print "\n".join(sorted(output, key=getKey))

