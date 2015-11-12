#! /usr/bin/env python

from collections import defaultdict
from suffix_tree import SuffixTree 
from datetime import datetime
import itertools
import pickle
import sys
import os


class overlapFinder(object):

    def __init__(self, sa, fr, bwt, ltf):
        self.suffix_array = sa
        self.first_rotation = fr
        self.bw_transform  = bwt
        self.last_to_first = ltf

    def generateKmers(self, seq, k):
        return [seq[i:i+k] for i in xrange(len(seq) - k + 1)]
    
    '''
    Example: find ACT
    T: ltf 10, 11
    C: ltf 5, 6
    sa[5] and sa[6] - the positions where the string ACT 
        starts; precisely pos #9 and pos #0
    eq:  A     C     T     G     A    C    A    T    G     A    C     T    G     A    A    C     A    $
    pos: 0     1     2     3     4    5    6    7    8     9    10    11   12    13   14   15    16   17
    sa:  17    16    13    14    4    9    0    6    15    5    10    1    12    3    8    11    2    7
    1st: $_1   A_1   A_2   A_3   A_4  A_5  A_6  A_7  C_1   C_2  C_3   C_4  G_1   G_2  G_3  T_1   T_2  T_3
    BWT: A_1   C_1   G_1   A_2   G_2  G_3  $_1  C_2  A_3   A_4  A_5   A_6  T_1   T_2  T_3  C_3   C_4  A_7
    ltf: 1     8     12    2     13   14   0    9    3     4    5     6    15    16   17   10    11   7

    main_sequence = "ACTGACATGACTGAACA$"
    reads = [('R1', 'ACT')] '''
    
    ''' Function is used to actually find overlaps
        in reads. The main part of the project '''
    def findOverlaps(self, read, kmer_len, error):
        kmer_map = defaultdict(int)
        read_kmers = self.generateKmers(read[1], kmer_len)

        for offset, kmer in enumerate(read_kmers):
            reverse_kmer = kmer[::-1]
            lft_indicies = list()
            bwt_values = self.first_rotation

            try:
                for i in xrange(len(kmer)-1):
                    cur_char = reverse_kmer[i]
                    next_char = reverse_kmer[i+1]
                    ''' First, find the range that the current nucleotide belongs to in the 1st rotation of the bwt '''
                    ''' Second, find all of the indicies in the bwt that correspond to the indices in the 1st rotation,
                        and also match the next character in sequence. Second step is met with the AND statement'''
                    bw_indicies = [index for index in xrange(self.first_rotation.index(bwt_values[0]), 
                                    self.first_rotation.index(bwt_values[-1]) + 1) 
                                    if self.first_rotation[index][0] == cur_char
                                    and self.bw_transform[index][0] == next_char]
                    ''' Third, get the actual strings that correspnd to the first and the last bwt indicies
                        found above '''
                    bwt_values = [self.bw_transform[bw_indicies[0]], self.bw_transform[bw_indicies[-1]]]
                    ''' Fourth, find the last to first indicies that correspot to the first and the last burrows-wheeler 
                        rotation indicies in case we are done mapping the read '''
                    lft_indicies = [self.last_to_first[bw_indicies[0]], self.last_to_first[bw_indicies[-1]]]
                
                ''' Map part of the read to a dictionary of frequent kmers '''
                for q in set(lft_indicies):
                    kmer_map[self.suffix_array[q] - offset] += 1

            except:
                if error: 
                    ''' The kmer did not map to anything in the genome, but 
                        we know the dataset has errors, so we continue '''
                    continue 
                else:
                    ''' Kmer did not map to anything and the dataset should not have had errors '''
                    print "[Logging {0}] {1} does not map to anything in the genome".format(getTime(), read[0])
                    return None
        
        ''' Get all of the positions in the genome where the read occurs '''
        max_occur = max(kmer_map.values())
        use_filter = lambda (k,v): v in xrange(max_occur-kmer_len, max_occur+kmer_len)
        kmer_map = filter(use_filter, kmer_map.iteritems())

        print "[Logging {0}] {1} has been mapped".format(getTime(), read[0])
        return [_kmer[0] for _kmer in kmer_map]



''' Get current date and time for the logger '''
def getTime():
    current_time = datetime.now()

    year = current_time.year
    month = current_time.month
    day = current_time.day
    hour = current_time.hour
    minute = current_time.minute
    second = current_time.second

    return "{0}/{1}/{2} {3}:{4}:{5}".format(year, month, day, hour, minute, second)

''' Function constructs a suffix array by 
    first constructing a suffix tree and
    performing depth first search on it '''
def constructSuffixArray(main_sequence):
    tree = SuffixTree(len(main_sequence))
    for char in main_sequence:
        tree.add_char(char)

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
def createSam(genome_header, read_overlaps, output=""):
    for read_tuple, overlap_indicies in filter(lambda x: x[1] is not None, read_overlaps):
        for overlap_index in overlap_indicies:
            output += str(read_tuple[0]) + "\t0\t" 
            output += genome_header + "\t" + str(overlap_index+1)
            output += "\t" + str(255) + "\t" + str(len(read_tuple[1])) + "M"
            output += "\t" + "*\t0\t0\t" + read_tuple[1] + "\t*\n"

    return output

''' If we have already tried to map a sequence we might have a file with all of the
    necessary precomputed objects. This function will determine if such a file exists
    and will load it. Otherwise it will have to create all those objects
    and save them for the future use '''
def loadOverlapFinder(filename, main_sequence, reads):
    input_file = "saved_objects/{0}.p".format(filename)

    if os.path.isfile(input_file):
        print "[Logging {0}] Loading previously constructed objects for this project".format(getTime())
        return pickle.load(open(input_file, "r"))

    else:
        print "[Logging {0}] Creating a new finder object for input files".format(getTime())
        return createNewOverlapFinder(input_file, main_sequence)

''' Creates a new finder object with suffix array, bwt, 1st rotation, and lft indicies
    and saves it for future use '''
def createNewOverlapFinder(output_file, main_sequence):

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
        pickle.dump(finder, open(output_file, "wb+"))
        print "[Logging {0}] Object saved successfully as {1}".format(getTime(), output_file)

        return finder


if __name__ == "__main__":

    print "[Logging {0}] The program started. Loading the input files".format(getTime())
    try: 
        with open(sys.argv[1]) as fd:
            genome_header = fd.readline().strip()[1:]
            main_sequence = "".join([seq.strip().lower() for seq in fd]) + "$"

        with open(sys.argv[2]) as fd:
            reads = [(header.strip()[1:], seq.strip().lower()) 
                    for header, seq in itertools.izip_longest(*[fd]*2)]
        kmer_len = int(sys.argv[3])
        errors = True if sys.argv[4][0].lower() == 'e' else False

    except IndexError:
        print "USAGE: python dna_mapper.py <chromosome_file> <reads_file> <kmer_length> <errors | noerrors>"
        sys.exit()

    try:
        print "[Logging {0}] Files loaded".format(getTime())
        main_filename = sys.argv[1].split("/")[-1]

        finder = loadOverlapFinder(main_filename, main_sequence, reads)
        print "[Logging {0}] Searching for read overlaps".format(getTime())

        # Find where overlaps exist in reads. Output format: [((<read_name>,<sequence>), <position_in_genome>)] 
        read_overlaps = [(single_read, finder.findOverlaps(single_read, kmer_len, errors)) for single_read in reads]
        print "[Logging {0}] The genome has been successfully indexed".format(getTime())

        sam_output = createSam(genome_header, read_overlaps)
        
        with open("output_sam/{0}.sam".format(main_filename), "w+") as fd:
            fd.write(sam_output)
        print "[Logging {0}] SAM file was written as output_sam/{1}.sam".format(getTime(), main_filename)

    except AssertionError:
        print "[ERROR {0}] Check assert statements. Lengths of certain data structures do not match".format(getTime())
        sys.exit()
