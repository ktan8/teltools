#!/usr/bin/python

import sys
import re

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


def check_seq_for_multiple_motifs(sequence, motif_set):
    pattern = "|".join(motif_set)
    result = re.search(pattern, sequence)
 
    if result is None:
        return 0
    else:
        return 1

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
    elif check_seq_for_multiple_motifs(softclipped_seq, telomere_seq_combinations):
        classification = "new_telomere_like"
    elif check_seq_for_multiple_motifs(softclipped_seq, fusion_seq_combinations):
        classification = "arm_fusion_like"

    
    lineList = lineList + [softclipped_seq_terminal, classification]

    print("\t".join(lineList))








