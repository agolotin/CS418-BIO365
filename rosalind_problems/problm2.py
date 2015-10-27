#! /usr/bin/env python

import sys

with open(sys.argv[1]) as fh:
    ptrn = fh.next().strip()

def reverse(pattern):
    pattern = pattern[::-1]
    nucs = {'A':'T', 'G':'C', 'T':'A', 'C':'G'}

    reverse = ""
    for nuc in pattern:
        reverse += nucs[nuc]
    return reverse

#print reverse(ptrn)
