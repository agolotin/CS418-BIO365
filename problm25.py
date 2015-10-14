#! /usr/bin/env python

import sys
from Trie import Trie 
from operator import sub

def trie_construction(seqs):
    trie = Trie(seqs)

    output = list()
    for item in trie.edges.items():
        output.append("->".join(map(str, map(sub, item[0], (1,1)))) + ":" + item[1])

    def getKey(item):
        return item[3] "1->2:A"

    return "\n".join(sorted(output, key=getKey))

if __name__ == "__main__":

    with open(sys.argv[1]) as fd:
        seqs = [line.strip() for line in fd.readlines()]

    print trie_construction(seqs)

