#! /usr/bin/env python

import os

kmer_size = [12, 17, 21]
for kmer in kmer_size:
	for chunk in xrange(1, 11):
		os.system("sbatch submit.sh -c {0} -k {1}".format(chunk, kmer))
