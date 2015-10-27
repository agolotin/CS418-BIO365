#! /usr/bin/env python

import sys
import operator

if __name__ == "__main__":
    with open(sys.argv[1]) as fd:
        word = fd.readline().strip()
    suffix_map = {word[i:]:i for i in xrange(len(word))}
    suffix_array = ", ".join(map(str, list(zip(*sorted(suffix_map.items(), key=operator.itemgetter(0)))[1])))
    print suffix_array
