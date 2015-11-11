#! /usr/bin/env python

import sys

def parseBlosum(matrix={}):

    with open('blosum.txt') as fd:
        first_row = fd.readline().strip().split(" ")
        scoring_lines = { line.strip().split(" ")[0]:line.strip().split(" ")[1:] for line in fd }

    for amino_acid, scores in scoring_lines.iteritems():
        for i, score in enumerate(scores):
            matrix[(amino_acid, first_row[i])] = int(score)

    return matrix


def globalAlignment(seq1, seq2, blosum, INDEL):
    #initialize scoring matrix
    scoring_matrix = [[0 for j in xrange(len(seq2)+1)] for i in xrange(len(seq1)+1)]
    backtrack = [[0 for j in xrange(len(seq2)+1)] for i in xrange(len(seq1)+1)]

    #initialize edges 
    for i in xrange(1, len(seq1)+1):
        scoring_matrix[i][0] = scoring_matrix[i-1][0] + INDEL
    for j in xrange(1, len(seq2)+1):
        scoring_matrix[0][j] = scoring_matrix[0][j-1] + INDEL

    for i in xrange(1, len(seq1)+1):
        for j in xrange(1, len(seq2)+1):
            scores = [ scoring_matrix[i-1][j] + INDEL, 
                       scoring_matrix[i][j-1] + INDEL, 
                       scoring_matrix[i-1][j-1] + blosum[seq1[i-1],seq2[j-1]]]

            scoring_matrix[i][j] = max(scores)
            backtrack[i][j] = scores.index(scoring_matrix[i][j])

    i, j = len(seq1), len(seq2)
    max_score = str(scoring_matrix[i][j])

    insert_indel = lambda word, i: word[:i] + "-" + word[i:]

    seq1_alignment, seq2_alignment = seq1, seq2

    while i * j != 0:
        if backtrack[i][j] == 0:
            i -= 1
            seq2_alignment = insert_indel(seq2_alignment, j)
        elif backtrack[i][j] == 1:
            j -= 1
            seq1_alignment = insert_indel(seq1_alignment, i)
        else:
            i -= 1
            j -= 1

    for k in range(i):
        seq2_alignment = insert_indel(seq2_alignment, 0)
    for k in range(j):
        seq1_alignment = insert_indel(seq1_alignment, 0)

    return max_score, seq1_alignment, seq2_alignment



if __name__ == "__main__":
    blosum = parseBlosum()
    with open(sys.argv[1]) as fd:
        seq1, seq2 = [seq.strip() for seq in fd]

    print "\n".join(globalAlignment(seq1, seq2, blosum, -5))
