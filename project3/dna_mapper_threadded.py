#! /usr/bin/env python

from collections import defaultdict
from suffix_tree import SuffixTree 
from datetime import datetime
import itertools
import pickle
import sys
import os

# Multithreading libraries
import threading
import Queue

# ===== GLOBAL VARIABLES =====
global exitFlag
global kmer_len 
global finder 
global queueLock 
global readLock
global workQueue 
# ============================

def getTime():
	current_time = datetime.now()

	year = current_time.year
	month = current_time.month
	day = current_time.day
	hour = current_time.hour
	minute = current_time.minute
	second = current_time.second

	return "{0}/{1}/{2} {3}:{4}:{5}".format(year, month, day, hour, minute, second)

# ======================= CLASS DEFENITIONS ==========================
class overlapFinder(object):

	def __init__(self, sa, fr, bwt, ltf):
		self.suffix_array = sa
		self.first_rotation = fr
		self.bw_transform  = bwt
		self.last_to_first = ltf

class overlapThread(threading.Thread):

	def __init__(self, tid, q):
		super(overlapThread, self).__init__()
		self.tid = tid
		self.queue = q

	# Called on thread.start()
	def run(self):
		Noccurs = True if "n1" in finder.first_rotation else False
		nucleotides = ["$1", "a1", "c1", "g1", "t1"]
		if Noccurs:
			nucleotides = ["$1", "a1", "c1", "g1", "n1", "t1"]

		while not exitFlag:
			queueLock.acquire()
			if not workQueue.empty():
				read = self.queue.get()
				queueLock.release()
				# Find overlaps 
				overlaps = self.findOverlaps(read, len(read[1]), Noccurs, nucleotides)

				# Populate the read overlaps map
				readLock.acquire()
				read_overlaps.append((read, overlaps))
				readLock.release()
			else:
				queueLock.release()


	# ==================== OVERLAP FINDER METHODS  =======================
	def generateKmers(self, seq, k):
		return [seq[i:i+k] for i in xrange(len(seq) - k + 1)]

	def findOverlaps(self, read, k_len, Noccurs, nucleotides):
		
		read_len = len(read[1])
		kmer_map = defaultdict(int)
		# At first do not split a read into kmers, but see if the read maps by itself
		read_kmers = self.generateKmers(read[1], k_len) if read_len != k_len else [ read[1] ] 

		for offset, kmer in enumerate(read_kmers):
			reverse_kmer = kmer[::-1].lower()
			lft_indicies = list()

			try:
				# If "n" appears in our kmer and our main sequence 
				# does not contain n's, then we ignore that kmer
				if "n" in kmer and not Noccurs:
					raise ValueError 

				if reverse_kmer[0] == "t":
					lft_indicies = [finder.first_rotation.index(reverse_kmer[0]+"1"), len(finder.first_rotation)-1] 
				else:
					next_nucleotide = nucleotides[ nucleotides.index(reverse_kmer[0]+"1")+1 ] 
					lft_indicies = [finder.first_rotation.index(reverse_kmer[0]+"1"), finder.first_rotation.index(next_nucleotide)]

				for i in xrange(len(kmer)-1):
					cur_char = reverse_kmer[i]
					next_char = reverse_kmer[i+1]

					bwt_indicies = [index for index in xrange(lft_indicies[0], lft_indicies[-1]+1)
									if finder.bw_transform[index][0] == next_char]

					lft_indicies = [finder.last_to_first[bwt_indicies[0]], finder.last_to_first[bwt_indicies[-1]]]
				
				''' Map part of the read to a dictionary of frequent kmers '''
				for q in xrange(lft_indicies[0], lft_indicies[-1]+1):
					kmer_map[finder.suffix_array[q] - offset] += 1

			except:
				''' The kmer did not map to anything in the genome, but 
					we know the dataset has errors, so we continue '''
				if read_len != k_len:
					continue
				else:
					return self.findOverlaps(read, kmer_len, Noccurs, nucleotides)

		if read_len != k_len:
			# Check to see if kmer_map was populated at all, aka at least one read maps to something
			max_occur = max(kmer_map.values()) if len(kmer_map) > 0 else 0
			''' If the most frequent occuring kmer occurs less than half of the time we can ignore the read altogether '''
			if max_occur < len(read_kmers) / 4:
				print "[Logging {0}] Ignored read {1}. Most occuring kmer of the read: occurs {2} times out of {3} kmers.".format(getTime(), read[0], max_occur, len(read_kmers))
				return None

			''' Get all of the positions in the genome where the read occurs based on a range '''
			use_filter = lambda (k, v): v in xrange(max_occur-(k_len/2), max_occur+(k_len/2))
			kmer_map = filter(use_filter, kmer_map.iteritems())

		print "[Logging {0}] {1} has been mapped".format(getTime(), read[0])
		if type(kmer_map) is not list:
			kmer_map = kmer_map.items()

		return [_kmer[0] for _kmer in kmer_map]
	# ====================================================================

# ====================== REST OF THE PROJECT =============================
def constructSuffixArray(main_sequence):
	tree = SuffixTree(len(main_sequence))
	for char in main_sequence:
		tree.add_char(char)

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

# ==================== LOADING PICKELED FILES ========================
def loadOverlapFinder(filename, main_sequence):
	input_file = "saved_objects/{0}.p".format(filename)

	if os.path.isfile(input_file):
		print "[Logging {0}] Loading previously constructed objects for this project".format(getTime())
		return pickle.load(open(input_file, "r"))

	else:
		print "[Logging {0}] Creating a new finder object for input files".format(getTime())
		return createNewOverlapFinder(input_file, main_sequence)

def createNewOverlapFinder(output_file, main_sequence):
		print "[Logging {0}] Constructing suffix array".format(getTime())
		suffix_array = constructSuffixArray(main_sequence)

		print "[Logging {0}] Constructing first rotation and Burrows-Wheeler Transform".format(getTime())
		first_rotation, bw_transform = BurrowsWheelerTransform(main_sequence, suffix_array)

		print "[Logging {0}] Constructing last-to-first array".format(getTime())
		last_to_first = constructLTF(first_rotation, bw_transform)

		# Custom data structure that holds all of the other necessary data structures.
		finder = overlapFinder(suffix_array, first_rotation, bw_transform, last_to_first)
		print "[Logging {0}] Overlap finder object constructed. Saving to file for future use".format(getTime())

		# Saving the data structure for later use
		pickle.dump(finder, open(output_file, "wb+"))
		print "[Logging {0}] Object saved successfully as saved_objects/{1}".format(getTime(), output_file)

		return finder


# ==================== CREATING OUTPUT FILE ==========================
def createSam(genome_header, read_overlaps, output=""):
	for read_tuple, overlap_indicies in filter(lambda x: x[1] is not None, read_overlaps):
		for overlap_index in overlap_indicies:
			output += str(read_tuple[0]) + "\t0\t" 
			output += genome_header + "\t" + str(overlap_index+1)
			output += "\t" + str(255) + "\t" + str(len(read_tuple[1])) + "M"
			output += "\t" + "*\t0\t0\t" + read_tuple[1] + "\t*\n"

	return output

# ============================= MAIN =================================
if __name__ == "__main__":

	try: 
		with open(sys.argv[1]) as fd:
			genome_header = fd.readline().strip()[1:]
			main_sequence = "".join([seq.strip().lower() for seq in fd]) + "$"

		with open(sys.argv[2]) as fd:
			reads = [(header.strip()[1:], seq.strip().lower()) 
					for header, seq, indentifier, quality in itertools.izip_longest(*[fd]*4)] 
		kmer_len = int(sys.argv[3])
		nproc = int(sys.argv[4])

	except IndexError:
		print "USAGE: python dna_mapper_threadded.py <chromosome_file> <reads_file> <kmer_length> <num_threads>"
		sys.exit()

	print "[Logging {0}] The program started. Input files loaded. Using {1} threads. Kmer size is {2}".format(getTime(), nproc, kmer_len)
	main_filename = sys.argv[1].split("/")[-1]

	finder = loadOverlapFinder(main_filename, main_sequence)
	read_overlaps = list() # Global variable to hold all read overlaps
	# ==================== MULTITHREADING ========================
	exitFlag = 0
	queueLock = threading.Lock()
	readLock = threading.Lock()
	workQueue = Queue.Queue(len(reads))

	tid = 0
	worker_threads = list()
	for i in xrange(nproc):
		worker = overlapThread(tid, workQueue)
		worker.start()
		worker_threads.append(worker)
		tid += 1
		
	print "[Logging {0}] Populating work queue".format(getTime())
	# Populate work queue
	queueLock.acquire()
	for single_read in reads:
		workQueue.put(single_read)
	queueLock.release()
	print "[Logging {0}] Work queue is populated".format(getTime())
	print "[Logging {0}] Searching for read overlaps".format(getTime())

	# Wait for the work queue to empty
	while not workQueue.empty():
		pass
	# Let go of the threads 
	exitFlag = 1

	# Join threads
	for thread in worker_threads:
		thread.join()
	# ===========================================================
	print "[Logging {0}] The genome has been successfully indexed".format(getTime())

	# Create output SAM file
	sam_output = createSam(genome_header, read_overlaps)
	with open("output_sam/nproc{0}_kmerlen{1}_{2}.sam".format(nproc, kmer_len, main_filename), "w+") as fd:
		fd.write(sam_output)

	print "[Logging {0}] SAM file was written as output_sam/nproc{1}_kmerlen{2}_{3}.sam".format(getTime(), nproc, kmer_len, main_filename)
