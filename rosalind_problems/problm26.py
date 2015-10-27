#! /usr/bin/env python

import sys
from Trie import Trie 

def prefix_matching(main_seq, seqs):
    trie = Trie(seqs)
    min_seq_len = min(map(len, seqs))

    indicies = list()
    for i in xrange(len(main_seq) - min_seq_len + 1):
        if trie.prefix_in_trie(main_seq[i:], 1) == True:
            indicies.append(i)

    return " ".join(map(str, indicies))


if __name__ == "__main__":
    with open(sys.argv[1]) as fd:
        main_seq = fd.readline().strip()
        seqs = [line.strip() for line in fd]

    print prefix_matching(main_seq, seqs)
