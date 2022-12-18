import numpy as np

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



def searchForMotif(sequence, motif):
	'''
	Search for a motif in a sequence of interest
	'''
	motifCombinations = generateMotifCombinations(motif)
	motifCounts = np.zeros(len(sequence))

	stringLength = len(sequence)
	motifLength = len(motif)

	# Search through the full string
	for i in range(stringLength - motifLength):
		subString = sequence[i:i+motifLength]
		# Loop through all the possible motifs
		for motifCurrent in motifCombinations:
			if subString == motifCurrent:
				motifCounts[i] += 1
			else:
				pass

	return(motifCounts)