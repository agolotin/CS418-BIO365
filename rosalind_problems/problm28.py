#! /usr/bin/env python

import sys
from collections import defaultdict

def stringFromBWT(text):
    original = _enumerate(list(text))
    ordered = _enumerate(sorted(list(text)))

    _map = {original[i]:ordered[i] for i in xrange(len(text))}
    start = _map['$0']
    
    original_text = ""
    for i in xrange(len(text)):
        original_text += start[0]
        start = _map[start]

    return original_text

def _enumerate(_list):
    enumerated = list()
    char_count = defaultdict(int)

    for char in _list:
        char_count[char] += 1
        enumerated.append(char + str(char_count[char]-1))

    return enumerated


if __name__ == "__main__":
    with open(sys.argv[1]) as fd:
        text = fd.readline().strip()
    print stringFromBWT(text)
