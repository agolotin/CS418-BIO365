#! /usr/bin/env python

import os

noerror_datasets = ['../datasets/synthetic.example.noerror.small.fasta', '../datasets/synthetic.noerror.small.fasta', '../datasets/synthetic.noerror.large.fasta']
outputs = ['output/synthetic.example.noerror.small.fasta', 'output/synthetic.noerror.small.fasta', 'output/synthetic.noerror.large.fasta', 'output/real.error.small.fasta', 'output/real.error.large.fasta']
real_datasets = ['../datasets/real.error.small.fasta', '../datasets/real.error.large.fasta']

kmer_list = [11,13,15,17,19,21,23,25,27,29,31]


i = 0
for _file in noerror_datasets:
    for kmer_size in kmer_list: 
        print "python assembler.py " + _file + " " + str(kmer_size) + " noerror > " + outputs[i] + ".output.kmer_size"+str(kmer_size)
        os.system("python assembler.py " + _file + " " + str(kmer_size) + " noerror > " + outputs[i] + ".output.kmer_size"+str(kmer_size))
    i += 1

k = 0
kmer_size = [15, 31] # best kmer size found using kmer genie
percentages = [0.1, 0.15, 0.20, 0.25]
for _file in real_datasets:
    for percentage in percentages:
#        os.system("python assembler.py " + _file + " " + str(kmer_size[k]) + " error " + str(percentage) + " > " + outputs[i] + ".output.assembler.percent"+str(percentage))
        os.system("python graph_assembler.py " + _file + " " + str(kmer_size[k]) + " error " + str(percentage) + " > " + outputs[i] + ".output.graph_assember.percent"+str(percentage))
    k += 1
    i += 1
