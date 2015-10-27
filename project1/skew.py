#! /usr/bin/env python

import sys

def minSkew(genome):
    g = 0
    c = 0
    min_skew = 0
    skew_list = list()
    for i in xrange(len(genome)):
        if genome[i] == 'G':
            g += 1
        if genome[i] == 'C':
            c += 1
        skew = g - c
        if skew < min_skew:
            skew_list = [str(i+1)]
            min_skew = skew

        if skew == min_skew and str(i+1) not in skew_list:
            skew_list.append(str(i+1))

    return skew_list
