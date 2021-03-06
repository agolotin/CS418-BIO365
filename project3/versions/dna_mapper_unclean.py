#! /usr/bin/env python
from suffix_tree import SuffixTree 
from collections import defaultdict
from datetime import datetime
import itertools
import operator 
import pickle
import sys
import os

''' Get current date and time for the logger '''
def getTime():
    current_time = datetime.now()

    year = current_time.year
    month = current_time.month
    hour = current_time.hour
    minute = current_time.minute
    second = current_time.second

    return "{0}/{1} {2}:{3}:{4}".format(year, month, hour, minute, second)

class overlapFinder(object):

    def __init__(self, sa, fr, bwt, ltf):
        self.suffix_array = sa
        self.first_rotation = fr
        self.bw_transform  = bwt
        self.last_to_first = ltf




    def findOverlaps(self, read):
        reverse_read = read[1][::-1]

        lft_indicies = list()
        bwt_values = self.first_rotation

        for i in xrange(len(read[1])-1):
            cur_char = reverse_read[i]
            next_char = reverse_read[i+1]

#            for index in xrange(self.first_rotation.index(bwt_values[0]), self.first_rotation.index(bwt_values[-1])):
#                if self.bw_transform[index][0] == next_char and self.first_rotation[index][0] == cur_char:
#                    print cur_char
#                    print next_char
#                    print self.first_rotation[index]
#                    print self.bw_transform[index]
#                    print 
#
#            first_indicies = [index for index, value in enumerate(self.first_rotation)  #change it to bwt_values
#                    if value[0] == cur_char and value in bwt_values]

            print "length of seq is " + str(len(read[1])) + " i = " + str(i)
            print "cur char " + cur_char
            print "next_char " + next_char

            ''' First, find the range that the current nucleotide belongs to in the 1st rotation of the bwt 
                In the same step, filter all of the indicies that do not correspond to the correct nucleotide
                in the bwt '''
            print range(self.first_rotation.index(bwt_values[0]),
                    self.first_rotation.index(bwt_values[-1]) + 1)
            bw_indicies = [index for index in xrange(self.first_rotation.index(bwt_values[0]), 
                            self.first_rotation.index(bwt_values[-1]) + 1) 
                            if self.first_rotation[index][0] == cur_char
                            and self.bw_transform[index][0] == next_char]
#            first_rot_indicies = [index for index in xrange(self.first_rotation.index(bwt_values[0]), 
#                            self.first_rotation.index(bwt_values[-1]) + 1) 
#                            if self.first_rotation[index][0] == cur_char]
#                            #and self.bw_transform[index][0] == next_char 

#            print "first rotation indicies " + str(first_rot_indicies)
            ''' Second, find all of the indicies in the bwt that correspond to the indices in the 1st rotation,
                but also match the next character in sequence '''
#            bw_indicies = [index for index in xrange(first_rot_indicies[0], first_rot_indicies[-1]+1) 
#                if self.bw_transform[index][0] == next_char]
#            print self.bw_transform[first_rot_indicies[0]]
            print "value of the first rot at index bwt_value " + str(self.first_rotation[self.first_rotation.index(bwt_values[0])])
            print "value of the bwt at index bwt_value " + str(self.bw_transform[self.first_rotation.index(bwt_values[0])])
            print "bw_indicies " + str(bw_indicies)
            ''' Third, get the actual strings that correspnd to the first and the last bwt indicies
                found above '''
            bwt_values = [self.bw_transform[bw_indicies[0]], self.bw_transform[bw_indicies[-1]]]

            print "bwt_values " + str(bwt_values)
            ''' Last, find what indicies in the first to last array correspond to the indicies of 
                the first and last occurence of the bwt'''
            lft_indicies = [self.last_to_first[bw_indicies[0]], self.last_to_first[bw_indicies[-1]]]
            print "lft indicies " + str(lft_indicies)
            print

        print "[Logging {0}] {1} has been mapped".format(getTime(), read[0])
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
reads = [('R1', 'ACT')]
main_sequence = "AATCGGGTTCAATCGGGGT$"
reads = ["ATCG", "GGGT"]
'''


''' Function constructs a suffix array from 
    first constructing a suffix tree '''
def constructSuffixArray(main_sequence):

    tree = SuffixTree(len(main_sequence))

    for char in main_sequence:
        tree.add_char(char)
#    tree.print_graphviz_tree()

    return tree.depthFirstSearch()

''' Enumerates characters in a sequence
    If input is "AGCTA" the output will be 
    ['A1', 'G1', 'C1', 'T1' 'A2'] '''
def _enumerate(sequence):
    enumerated = list()
    char_count = defaultdict(int)

    for char in sequence:
        char_count[char] += 1
        enumerated.append(char + str(char_count[char]))

    return enumerated

''' Constructs Burrows-Wheeler Transform as well as the
    first lexicographic rotation of a given string from
    suffix array '''
def BurrowsWheelerTransform(seq, suffix_array):
    first_rotation = list()
    bw_transform = list()
    last_char_position = len(seq)-1

    for index in suffix_array:
        first_rotation.append(seq[index])
        bw_transform.append(seq[ index-1 if index != 0 else last_char_position ])

    return _enumerate(first_rotation), _enumerate(bw_transform)

''' Constructs last to first array by searching for the 
    index of a string from Burrows-Wheeler Transform in
    the first lexicographic rotation '''
def constructLTF(first_rotation, bw_transform):
    last_to_first = [first_rotation.index(key) for key in bw_transform]

    return last_to_first

''' Creates a .sam file to be outputted '''
def createSam(reference, read_overlaps, output=""):

    for read_tuple, overlap_indicies in read_overlaps:
        for overlap_index in overlap_indicies:

            output += str(read_tuple[0]) + "\t0\t" 
            output += reference + "\t" + str(overlap_index+1)
            output += "\t" + str(255) + "\t" + str(len(read_tuple[1])) + "M"
            output += "\t" + "*\t0\t0\t" + read_tuple[1] + "\t*\n"

    return output

def loadOverlapFinder(filename, main_sequence, reads):

    if os.path.isfile("saved_objects/{0}.p".format(filename)):
        print "[Logging {0}] Loading previously constructed objects for this project".format(getTime())
        return pickle.load(open("saved_objects/{0}.p".format(filename), "r"))
    else:
        print "[Logging {0}] Creating a new finder object for input files".format(getTime())
        return createNewOverlapFinder(filename, main_sequence, reads)

def createNewOverlapFinder(filename, main_sequence, reads):

        suffix_array = constructSuffixArray(main_sequence)
        assert len(suffix_array) == len(main_sequence)
        print "[Logging {0}] Suffix array constructed".format(getTime())

        first_rotation, bw_transform = BurrowsWheelerTransform(main_sequence, suffix_array)
        assert len(first_rotation) == len(bw_transform) == len(main_sequence)
        print "[Logging {0}] The first rotation and the Burrows-Wheeler Transform constructed".format(getTime())

        last_to_first = constructLTF(first_rotation, bw_transform)
        assert len(last_to_first) == len(main_sequence)
        print "[Logging {0}] Last-to-first array constructed".format(getTime())

        # Custom data structure that holds all of the other necessary data structures.
        finder = overlapFinder(suffix_array, first_rotation, bw_transform, last_to_first)
        print "[Logging {0}] Overlap finder object constructed. Saving to file for future use".format(getTime())

        # Saving the data structure for later use
        output_file = "saved_objects/{0}.p".format(filename)
        pickle.dump(finder, open(output_file, "wb+"))
        print "[Logging {0}] Object saved successfully as {1}".format(getTime(), output_file)

        return finder


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
        main_filename = sys.argv[1].split("/")[-1]
        pickle_filename = sys.argv[2].split("/")[-1]
        
        print "[Logging {0}] Files loaded".format(getTime())
        finder = loadOverlapFinder(pickle_filename, main_sequence, reads)

#        with open("visual.txt", "w+") as fd:
#            fd.write("seq: ")
#            for nuc in main_sequence:
#                fd.write(str(nuc) +"        ")
#            fd.write("\nsa:  ")
#            for nuc in finder.suffix_array:
#                fd.write(str(nuc) + "    ")
#            fd.write("\n1st: ")
#            for nuc in finder.first_rotation:
#                fd.write(str(nuc) + "       ")
#            fd.write("\nbwt: ")
#            for nuc in finder.bw_transform:
#                fd.write(str(nuc) + "       ")
#            fd.write("\nltf: ")
#            for nuc in finder.last_to_first:
#                fd.write(str(nuc) + "    ")

        # Find where overlaps exist in reads
        print "[Logging {0}] Searching for read overlaps".format(getTime())
        read_overlaps = [(single_read, finder.findOverlaps(single_read)) for single_read in reads]

        print "[Logging {0}] Read overlap map has been constructed".format(getTime())
        sam_output = createSam(reference, read_overlaps)
        
        with open("output_sam/{0}.sam".format(main_filename), "w+") as fd:
            fd.write(sam_output)

        print "[Logging {0}] SAM file was written into output_sam/ directory as {1}.sam".format(getTime(), main_filename)

    except AssertionError:
        print "[ERROR {0}] Check assert statements. Lengths of certain data structures do not match".format(getTime())
        sys.exit()


