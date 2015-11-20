#! /usr/bin/env python

import sys 

def fittingAlignment(seq1, seq2):
    #initialize scoring matrix
    scoring_matrix = [[0 for j in xrange(len(seq2)+1)] for i in xrange(len(seq1)+1)]
    backtrack = [[-100 for j in xrange(len(seq2)+1)] for i in xrange(len(seq1)+1)]

    for i in xrange(1, len(seq1)+1):
        for j in xrange(1, len(seq2)+1):
            match = 1 if seq1[i-1] == seq2[j-1] else -1
            scores = [scoring_matrix[i-1][j] - 1, 
                      scoring_matrix[i][j-1] - 1, 
                      scoring_matrix[i-1][j-1] + match]

            scoring_matrix[i][j] = max(scores)
            backtrack[i][j] = scores.index(scoring_matrix[i][j])

    j = len(seq2)
    i = max(enumerate([scoring_matrix[row][j] for row in xrange(len(seq2), len(seq1))]),key=lambda x: x[1])[0] + len(seq2)

    max_score = str(scoring_matrix[i][j])

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
        elif backtrack[i][j] == 2:
            i -= 1
            j -= 1

    seq1_alignment = seq1_alignment[i:]

    return str(max_score), seq1_alignment, seq2_alignment


if __name__ == "__main__":
    with open(sys.argv[1]) as fd:
        seq1, seq2 = [seq.strip() for seq in fd]

    print "\n".join(fittingAlignment(seq1, seq2))
