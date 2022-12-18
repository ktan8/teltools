#!/usr/bin/python

import sys
#from lib.estimateTelomereLength import *
from lib.motif_lib import *


def read_fasta(fp):
	'''
	Generator code from Biopython. Parses a fasta file
	and yields a single fasta line on each call
	https://stackoverflow.com/questions/7654971/parsing-a-fasta-file-using-a-generator-python
	'''
	name, seq = None, []
	for line in fp:
		line = line.rstrip()
		if line.startswith(">"):
			if name: yield (name, ''.join(seq))
			name, seq = line, []
		else:
			seq.append(line)
	if name: yield (name, ''.join(seq))

	#print "good"



def check_motif_present(sequence, motif):
	'''
	Check if a sequence has a motif, or
	a rearrangement of it.
	'''

	motif_combinations = generateMotifCombinations(motif)

	for motif in motif_combinations:
		if motif in sequence:
			return 1



	# Return a false if no motif was found
	# across all the combinations
	return 0


def main():
	fasta_file = sys.argv[1]
	output_file = sys.argv[2]
	log_file = output_file + ".log"

	short_motif = "TTAGGG"
	long_motif = short_motif + short_motif
	long_motif_revCom = reverseComplement(long_motif)
	
	total = 0
	hits = 0

	fasta_file_handle = open(fasta_file, 'r')
	fasta_output = open(output_file, 'w')
	logfile_output = open(log_file, 'w')

	# Loop through the fasta file item generator
	for seq_name, seq in read_fasta(fasta_file_handle):
		total += 1

		# Heuristic to remove reads without TTAGGG (at least 2 TTAGGG present)
		if not (check_motif_present(seq, long_motif) or check_motif_present(seq, long_motif_revCom)):
			continue


		fasta_output.write(seq_name + "\n")
		fasta_output.write(seq + "\n")


		hits += 1



	fasta_output.close()


	logfile_output.write("Total reads: %d\n" %total)
	logfile_output.write("Total hits: %d\n" %hits)
	logfile_output.close()
	

if __name__ == "__main__":

	main()

