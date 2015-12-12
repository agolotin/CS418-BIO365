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


def globalAlignment(seq1, seq2, blosum, open_penalty, gap_ext):
    #initialize 3 scoring matrixes: bottom, middle and upper matrices
    scoring_matrix = [[[0 for j in xrange(len(seq2)+1)] for i in xrange(len(seq1)+1)] for k in xrange(3)]
    backtrack = [[[-1 for j in xrange(len(seq2)+1)] for i in xrange(len(seq1)+1)] for k in xrange(3)]

    #initialize edges 
    for i in xrange(1, len(seq1)+1):
        scoring_matrix[0][i][0] = -open_penalty - (i-1)*gap_ext
        scoring_matrix[1][i][0] = -open_penalty - (i-1)*gap_ext
        scoring_matrix[2][i][0] = -10 * open_penalty
    for j in xrange(1, len(seq2)+1):
        scoring_matrix[0][0][j] = -10 * open_penalty
        scoring_matrix[1][0][j] = -open_penalty - (j-1)*gap_ext
        scoring_matrix[2][0][j] = -open_penalty - (j-1)*gap_ext


    for i in xrange(1, len(seq1)+1):
        for j in xrange(1, len(seq2)+1):
            lower_scores = [scoring_matrix[0][i-1][j] - gap_ext, scoring_matrix[1][i-1][j] - open_penalty]
            scoring_matrix[0][i][j] = max(lower_scores)
            backtrack[0][i][j] = lower_scores.index(scoring_matrix[0][i][j])
            
            upper_scores = [scoring_matrix[2][i][j-1] - gap_ext, scoring_matrix[1][i][j-1] - open_penalty]
            scoring_matrix[2][i][j] = max(upper_scores)
            backtrack[2][i][j] = upper_scores.index(scoring_matrix[2][i][j])

            middle_scores = [scoring_matrix[0][i][j], scoring_matrix[1][i-1][j-1] + blosum[seq1[i-1],seq2[j-1]], scoring_matrix[2][i][j]]
            scoring_matrix[1][i][j] = max(middle_scores)
            backtrack[1][i][j] = middle_scores.index(scoring_matrix[1][i][j])

    i, j = len(seq1), len(seq2)
    seq1_alignment, seq2_alignment = seq1, seq2
    matrix_scores = [scoring_matrix[0][i][j], scoring_matrix[1][i][j], scoring_matrix[2][i][j]]
    max_score = str(max(matrix_scores))
    backtrack_matrix = matrix_scores.index(int(max_score))

#    for score in scoring_matrix:
#        for x in score:
#            print x
#        print
#    print
#    for score in backtrack:
#        for x in score:
#            print x
#        print
#    print backtrack_matrix
#    sys.exit()
    insert_indel = lambda word, i: word[:i] + "-" + word[i:]


    while i * j != 0:
        if backtrack_matrix == 0:
            if backtrack[0][i][j] == 1:
                backtrack_matrix = 1
            i -= 1
            seq2_alignment = insert_indel(seq2_alignment, j)
        elif backtrack_matrix == 1:
            if backtrack[1][i][j] == 0:
                backtrack_matrix = 0
            elif backtrack[1][i][j] == 2:
                backtrack_matrix = 2
            else:
                i -= 1
                j -= 1
        else:
            if backtrack[2][i][j] == 1:
                backtrack_matrix = 1
            j -= 1
            seq1_alignment = insert_indel(seq1_alignment, i)

    for k in range(i):
        seq2_alignment = insert_indel(seq2_alignment, 0)
    for k in range(j):
        seq1_alignment = insert_indel(seq1_alignment, 0)

    return max_score, seq1_alignment, seq2_alignment



if __name__ == "__main__":
    blosum = parseBlosum()
    with open(sys.argv[1]) as fd:
        seq1, seq2 = [seq.strip() for seq in fd]

    print "\n".join(globalAlignment(seq1, seq2, blosum, 11, 1))
