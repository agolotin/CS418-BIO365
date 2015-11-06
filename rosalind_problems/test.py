#! /usr/bin/env python

from suffix_tree import SuffixTree
import sys



if __name__ == "__main__":
    with open(sys.argv[1]) as fd:
        text = fd.readline().strip()

    tree = SuffixTree(len(text))
    for char in text:
        tree.add_char(char)

    print tree.depthFirstSearch()
