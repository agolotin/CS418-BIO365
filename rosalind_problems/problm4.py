#! /usr/bin/env python

import sys

with open(sys.argv[1]) as fd:
    genome = fd.next().strip()
    k, L, t = fd.next().strip().split(" ")

Lclumps = list()
for i in range(len(genome)):
    clump = genome[i : i + int(L)]
    if len(clump) == int(L):
        Lclumps.append(clump)

######begin function#######
def findKmers(times, size, _clmp):
    kmrs = dict()
    highest = 0
    output = list()
    for i in range(len(_clmp)):
        kmr = _clmp[i : i + size]
        if kmr in kmrs:
            kmrs[kmr] += 1
        else:
            kmrs[kmr] = 1

        if kmrs[kmr] > highest:
            highest = kmrs[kmr]

    if highest < times:
        return False

    for key, val in kmrs.iteritems():
        if int(val) >= times:
            output.append(key)

    return output
######end function#######

output = list()
for clmp in Lclumps:
    kmer = findKmers(int(t), int(k), clmp)
    if kmer != False:
        output.extend(kmer)

print " ".join(set(output))

