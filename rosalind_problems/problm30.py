#! /usr/bin/env python

import sys
from Trie import Trie

def findMultiplePatternMatches(seq, patterns):
    trie = Trie(patterns)
    indicies = [i for i in xrange(len(seq) - max(map(len, patterns)) + 1) if trie.prefix_in_trie(seq[i:],1) is True]

#    for i in xrange(len(seq) - max(map(len, patterns)) + 1):
#        if trie.prefix_in_trie(seq[i:],1):
#            indicies.append(i)

    return indicies


if __name__ == "__main__":
    with open(sys.argv[1]) as fd:
        seq = fd.readline().strip()
        patterns = [line.strip() for line in fd]

    pattern_indicies = map(str, findMultiplePatternMatches(seq, patterns))
    print " ".join(pattern_indicies)
