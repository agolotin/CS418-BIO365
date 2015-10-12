

def profile_most_probable_kmer(dna, k, profile):
    nuc_loc = {nucleotide:index for index,nucleotide in enumerate('ACGT')}
    max_probability = -1

    for i in xrange(len(dna)-k+1):
        current_probability = 1
        for j, nucleotide in enumerate(dna[i:i+k]):
            current_probability *= profile[j][nuc_loc[nucleotide]]

        if current_probability > max_probability:
            max_probability = current_probability
            most_probable = dna[i:i+k]

    return most_probable
