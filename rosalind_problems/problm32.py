#! /usr/bin/env python

import sys

def parsePAM(matrix={}):

    with open("pam250.txt") as fh:
            first_row = fh.readline().strip().split()
            lines = fh.readlines()

    value_dict = {}
    for line in lines:
            values = line.strip().split()
            line_dict = {}
            amino_acid = values[0]
            for i in range(1, len(values)):
                    line_dict[ first_row[i-1] ] = int(values[i])
            value_dict[amino_acid] = line_dict

    return value_dict

def getHighestScore(scoring_matrix):
    total_max = max( map(max, scoring_matrix) )
    max_i, max_j = -1, -1

    for i in range(len(scoring_matrix)):
        if total_max in scoring_matrix[i]:
            max_i, max_j = i, scoring_matrix[i].index(total_max)

    return total_max, max_i, max_j

def localAlignment(seq1, seq2, pam250, INDEL):
    #initialize scoring matrix
    scoring_matrix = [[0 for j in xrange(len(seq2)+1)] for i in xrange(len(seq1)+1)]
    backtrack = [[-1 for j in xrange(len(seq2)+1)] for i in xrange(len(seq1)+1)]

    for i in xrange(1, len(seq1)+1):
        for j in xrange(1, len(seq2)+1):
            scores = [scoring_matrix[i-1][j] - INDEL, 
                      scoring_matrix[i][j-1] - INDEL, 
                      scoring_matrix[i-1][j-1] + pam250[seq1[i-1]][seq2[j-1]], 0]

            scoring_matrix[i][j] = max(scores)
            if scoring_matrix[i][j] == 0:
                backtrack[i][j] = 3
            else:
                backtrack[i][j] = scores.index(scoring_matrix[i][j])

    max_score, i, j = getHighestScore(scoring_matrix)

    insert_indel = lambda word, k: word[:k] + "-" + word[k:]
    seq1_alignment = seq1[:i]
    seq2_alignment = seq2[:j]

    while backtrack[i][j] != 3 and i * j != 0:
        if backtrack[i][j] == 0:
            i -= 1
            seq2_alignment = insert_indel(seq2_alignment, j)
        elif backtrack[i][j] == 1:
            j -= 1
            seq1_alignment = insert_indel(seq1_alignment, i)
        elif backtrack[i][j] == 2:
            i -= 1
            j -= 1

    seq1_alignment = seq1_alignment[i:]
    seq2_alignment = seq2_alignment[j:]

    return str(max_score), seq1_alignment, seq2_alignment


if __name__ == "__main__":
    pam250 = parsePAM()
    with open(sys.argv[1]) as fd:
        seq1, seq2 = [seq.strip() for seq in fd]

    print "\n".join(localAlignment(seq1, seq2, pam250, 5))
