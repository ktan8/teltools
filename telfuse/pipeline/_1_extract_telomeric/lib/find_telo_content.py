import numpy as np
import re
import sys
from decimal import *

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



def searchForMotifRe(sequence, motif):
	'''
	Search for a motif in a sequence of interest
	'''
	motifCombinations = generateMotifCombinations(motif)
	motifVector = np.zeros(len(sequence))

	stringLength = len(sequence)
	motifLength = len(motif)

	# Loop through all the possible motifs
	for motifCurrent in motifCombinations:
		matches = re.finditer(motifCurrent, sequence)
		matches_span = [i.span() for i in matches]
		posn = [i[0] for i in matches_span]
		for start in posn:
			motifVector[start:start+len(motif)] = 1

	return(motifVector)


def calc_telo_content(sequence, motif):
	'''
	Calculate telomeric content in a sequence using
	a provided motif.
	'''
	seq_len = len(sequence)
	motif_len_in_seq = np.sum(searchForMotifRe(sequence, motif))
	telo_content = motif_len_in_seq / seq_len
	
	#telo_content = round(Decimal(telo_content, 5))
	result = [seq_len, motif_len_in_seq, telo_content]

	return(result)


def calc_telo_content_both(sequence, motif):
	'''
	Calc telomeric content for both forward and reverse
	complement of motif in the sequence
	'''
	motif_complement = reverseComplement(motif)
	result_fwd	= calc_telo_content(sequence, motif)
	result_revCom 	= calc_telo_content(sequence, motif_complement)
	
	return([result_fwd, result_revCom])


def process_table(file, colnum, header=True, motif="TTAGGG"):
	colnum = int(colnum)
	filehandle = open(file, "r") 
	if header:
		header = filehandle.readline()
		headerList = header.strip().split("\t")
		headerList = headerList + ["seqLen", "TTAGGG_len", "TTAGGG_content", "CCCTAA_len", "CCCTAA_content"]
		#print(headerList)
		print("\t".join(headerList))		

	for line in filehandle:
		lineList = line.strip().split("\t")
		softclipped_seq = lineList[colnum]
		result_fwd, result_revCom = calc_telo_content_both(softclipped_seq, motif)
		result_all = result_fwd + result_revCom[1:3]
		lineList = lineList + result_all
		print("\t".join(map(str, lineList)))

#searchForMotifRe("TTAGGGTTAGGGTTAGGGTTA", "TAGGGT")
#print(calc_telo_content_both("TTAGGGTTAGGGTTAGGGTTAGGGTTTTT", "TAGGGT"))
#process_table("3_Supplementary_Tables_table_S2.txt", 8)

if __name__ == "__main__":
	file 	= sys.argv[1]
	colnum 	= sys.argv[2]
	process_table(file, colnum)


