#! /usr/bin/env python

import sys

def BWT(text):
    _len = len(text)
    cyclic_rotation_nuc = lambda i, n: text[(n-i) % _len]
    transforms = list()
    _text = text
    for i in xrange(1, _len+1):
        _text = cyclic_rotation_nuc(i, _len) + _text[:_len-1]
        transforms.append(_text)
    return "".join(zip(*sorted(transforms))[-1])

if __name__ == "__main__":
    with open(sys.argv[1]) as fd:
        text = fd.readline().strip()
    print BWT(text)
