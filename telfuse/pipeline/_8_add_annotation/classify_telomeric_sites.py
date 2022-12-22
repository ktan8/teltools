#!/usr/bin/python

import sys


candidate_file = sys.argv[1]

candidates = open(candidate_file, "r")

telomere_seq = "TTAGGG"


def reverseComplement(seq):
        '''
        https://stackoverflow.com/questions/25188968/reverse-complement-of-dna-strand-using-python
        '''
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        reverse_complement = "".join(complement.get(base, base) for base in reversed(seq))

        return reverse_complement


def generateMotifCombinations(sequence):
        '''
        Generate all the possible sequence combinations
        for a motif of interest after shifting the frame.
        e.g. TTAGGG, TAGGGT, AGGGTT, GGGTTA, GGTTAG, GTTAGG
        '''
        seqLength = len(sequence)
        seqCombinations = []
        for i in range(seqLength):
                newString = sequence[i:seqLength] + sequence[0:i]

                seqCombinations.append(newString)

        return seqCombinations


telomere_seq = "TTAGGGTTAGGG"
telomere_seq_revcom = reverseComplement(telomere_seq)

telomere_seq_combinations = generateMotifCombinations(telomere_seq)
fusion_seq_combinations = generateMotifCombinations(telomere_seq_revcom)

#print(fusion_seq_combinations)

for line in candidates:
    lineList = line.strip().split("\t")
    softclipped_seq = lineList[7]
    softclipped_seq_len = len(softclipped_seq)
    softclipped_seq_terminal = ""
    if(softclipped_seq_len >= 12):
        softclipped_seq_terminal = softclipped_seq[0:12]
    else:
        softclipped_seq_terminal = softclipped_seq

    classification = "NA"
    if softclipped_seq_terminal in telomere_seq_combinations:
        classification = "new_telomere"
    elif softclipped_seq_terminal in fusion_seq_combinations:
        classification = "arm_fusion"
    
    lineList = lineList + [softclipped_seq_terminal, classification]

    print("\t".join(lineList))








