#! /usr/bin/env python

import sys 

def getHighestScore(scoring_matrix):
    total_max = max( map(max, scoring_matrix) )
    max_i, max_j = -1, -1

    for i in range(len(scoring_matrix)):
        if total_max in scoring_matrix[i]:
            max_i, max_j = i, scoring_matrix[i].index(total_max)

    return total_max, max_i, max_j

def fittingAlignment(seq1, seq2):
    #initialize scoring matrix
    scoring_matrix = [[0 for j in xrange(len(seq2)+1)] for i in xrange(len(seq1)+1)]
    backtrack = [[-100 for j in xrange(len(seq2)+1)] for i in xrange(len(seq1)+1)]

    max_score = -1 * (len(seq1) + len(seq2))
    max_indicies = ()

    for i in xrange(1, len(seq1)+1):
        for j in xrange(1, len(seq2)+1):
            match = 1 if seq1[i-1] == seq2[j-1] else -2
            scores = [scoring_matrix[i-1][j] - 2, 
                      scoring_matrix[i][j-1] - 2, 
                      scoring_matrix[i-1][j-1] + match]

            scoring_matrix[i][j] = max(scores)
            backtrack[i][j] = scores.index(scoring_matrix[i][j])

            # Check to see if we found a maximum along the last row and column
            if i == len(seq1) or j == len(seq2):
                if scoring_matrix[i][j] > max_score:
                    max_score = scoring_matrix[i][j]
                    max_indicies = (i, j)

    i, j = max_indicies

    insert_indel = lambda word, k: word[:k] + "-" + word[k:]
    seq1_alignment = seq1[:i]
    seq2_alignment = seq2[:j]

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

    seq1_alignment = seq1_alignment[i:]
    seq2_alignment = seq2_alignment[j:]

    return str(max_score), seq1_alignment, seq2_alignment


if __name__ == "__main__":
    with open(sys.argv[1]) as fd:
        seq1, seq2 = [seq.strip() for seq in fd]

    print "\n".join(fittingAlignment(seq1, seq2))
