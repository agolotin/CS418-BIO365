#! /usr/bin/env python

import os

noerror_datasets = ['../synthetic.example.noerror.small.fasta', '../synthetic.noerror.small.fasta', '../synthetic.noerror.large.fasta']
real_datasets = ['real.error.small.fasta', 'real.error.large.fasta']

kmer_list = [7,8,9,10,11,12,13,14,15,16,17,18,19,20]

for _file in noerror_datasets:
    for kmer_size in kmer_list: 
        print "python assembler.py " + _file + " " + str(kmer_size) + " noerror > " + _file + ".output.kmer_size"+str(kmer_size)
        os.system("python assembler.py " + _file + " " + str(kmer_size) + " noerror > " + _file + ".output.kmer_size"+str(kmer_size))

#kmer_size = 10
#percentages = [0.01, 0.05, 0.1, 0.15]
#for _file in real_datasets:
#    for percentage in percentages:
#        print "python assembler ../" + _file + " " + str(kmer_size) + " noerror > " + _file + ".output.kmer_size"+str(kmer_size)
#        os.system("python assembler.py ../" + _file + " " + kmer_size " error " + percentage + " > " + _file + ".output.assembler"+percentage)
#        os.system("python graph_assembler.py ../" + _file + " " + kmer_size " error " + percentage + " > " + _file + ".output.graph_assember"+percentage)
#    sys.exit()
