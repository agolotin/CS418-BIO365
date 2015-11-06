#! /usr/bin/env python

from suffix_tree import SuffixTree
from collections import defaultdict
from datetime import datetime
import itertools
import operator
import sys

class overlapFinder(object):

    def __init__(self, sa, fr, bwt, ltf):
        self.suffix_array = sa
        self.first_rotation = fr
        self.bw_transform  = bwt
        self.last_to_first = ltf



    def findOverlaps(self, read):
        reverse_read = read[1][::-1]

        lft_indicies = list()
        previous_bw = self.first_rotation

        for i in xrange(len(read[1])-1):
            cur_char = reverse_read[i]
            next_char = reverse_read[i+1]

            first_indicies = [index for index, value in enumerate(self.first_rotation) 
                    if value[0] == cur_char and value in previous_bw]

            bw_indicies = [index for index in range(first_indicies[0],first_indicies[-1]+1) 
                    if self.bw_transform[index][0] == next_char]
            previous_bw = [self.bw_transform[bw_indicies[0]], self.bw_transform[bw_indicies[-1]]]

            lft_indicies = [self.last_to_first[bw_indicies[0]], self.last_to_first[bw_indicies[-1]]]

        print "[Logging " + getTime() + "] Read '" + read[0] + "' has been mapped"
        return list(set([self.suffix_array[i] for i in lft_indicies]))

'''
ACT
T: ltf 10, 11
C: ltf 5, 6
sa[5] and sa[6] - the positions where text occurs

eq:  A     C     T     G     A    C    A    T    G     A    C     T    G     A    A    C     A    $
pos: 0     1     2     3     4    5    6    7    8     9    10    11   12    13   14   15    16   17
sa:  17    16    13    14    4    9    0    6    15    5    10    1    12    3    8    11    2    7
1st: $_1   A_1   A_2   A_3   A_4  A_5  A_6  A_7  C_1   C_2  C_3   C_4  G_1   G_2  G_3  T_1   T_2  T_3
BWT: A_1   C_1   G_1   A_2   G_2  G_3  $_1  C_2  A_3   A_4  A_5   A_6  T_1   T_2  T_3  C_3   C_4  A_7
ltf: 1     8     12    2     13   14   0    9    3     4    5     6    15    16   17   10    11   7



main_sequence = "ACTGACATGACTGAACA$"
reads = ["AACA"]
main_sequence = "AATCGGGTTCAATCGGGGT$"
reads = ["ATCG", "GGGT"]
'''



def constructSuffixArray(main_sequence):

    tree = SuffixTree(len(main_sequence))

    for char in main_sequence:
        tree.add_char(char)
#    tree.print_graphviz_tree()

    return tree.depthFirstSearch()

def _enumerate(sequence):
    enumerated = list()
    char_count = defaultdict(int)

    for char in sequence:
        char_count[char] += 1
        enumerated.append(char + str(char_count[char]))

    return enumerated

def BurrowsWheelerTransform(seq, suffix_array):
    first_rotation = list()
    bw_transform = list()
    last_char_position = len(seq)-1

    for index in suffix_array:
        first_rotation.append(seq[index])
        bw_transform.append(seq[ index-1 if index != 0 else last_char_position ])

    return _enumerate(first_rotation), _enumerate(bw_transform)

def constructLTF(first_rotation, bw_transform):
    last_to_first = [first_rotation.index(key) for key in bw_transform]

    return last_to_first


def createSam(reference, read_overlaps):
    output = ""

    for read_tuple, overlap_indicies in read_overlaps:
        for overlap_index in overlap_indicies:

            output += str(read_tuple[0]) + "\t0\t" 
            output += reference + "\t" + str(overlap_index+1)
            output += "\t" + str(255) + "\t" + str(len(read_tuple[1])) + "M"
            output += "\t" + "*\t0\t0\t" + read_tuple[1] + "\t*\n"

    return output

def getTime():
    current_time = datetime.now()

    year = current_time.year
    month = current_time.month
    hour = current_time.hour
    minute = current_time.minute
    second = current_time.second

    return "{0}/{1} {2}:{3}:{4}".format(year, month, hour, minute, second)

if __name__ == "__main__":

    try: 

        with open(sys.argv[1]) as fd:
            reference = fd.readline().strip()[1:]
            main_sequence = "".join([seq.strip().lower() for seq in fd]) + "$"
        with open(sys.argv[2]) as fd:
            reads = [(header.strip()[1:], seq.strip().lower()) for header, seq in itertools.izip_longest(*[fd]*2)]
    except IndexError:
        print "USAGE: python dna_mapper.py <chromosome_file> <reads_file>"
        sys.exit()

    try:
        print "[Logging " + getTime() + "] Files loaded"
        suffix_array = constructSuffixArray(main_sequence)
        assert len(suffix_array) == len(main_sequence)

        print "[Logging " + getTime() + "] Suffix array constructed"
        first_rotation, bw_transform = BurrowsWheelerTransform(main_sequence, suffix_array)
        assert len(first_rotation) == len(bw_transform) == len(main_sequence)

        print "[Logging " + getTime() + "] The first rotation and the Burrows-Wheeler Transform constructed"
        last_to_first = constructLTF(first_rotation, bw_transform)
        assert len(last_to_first) == len(main_sequence)

        print "[Logging " + getTime() + "] Last-to-first array constructed"
        # Custom data structure that holds all of the other necessary data structures.
        # Should be written to a file later
        finder = overlapFinder(suffix_array, first_rotation, bw_transform, last_to_first)

        # Find there overlaps exist in reads
        print "[bogging " + getTime() + "] Searching for read overlaps"
        read_overlaps = [(single_read, finder.findOverlaps(single_read)) for single_read in reads]

        print "[Logging " + getTime() + "] Read overlap map has been constructed"
        sam_output = createSam(reference, read_overlaps)
        
        with open(sys.argv[1] + ".sam", "w+") as fd:
            fd.write(sam_output)
        print "[Logging " + getTime() + "] SAM file was written into the working directory as " + sys.argv[1] + ".sam"

    # Error occured when doing assertion
    except AssertionError:
        print "[ERROR " + getTime() + "] Check assert statements. Lengths of certain data structures do not match"
        sys.exit()


