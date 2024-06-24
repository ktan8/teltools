import sys
import os
import subprocess
from os import sys, path
from sam_cigar import *
from lib.samReadClass import *
from lib.fastqLib import *

'''
Extract all the soft clipped sites from
a bam file to generate a list of candidate
sites
'''


def runProcess(exe):
    p = subprocess.Popen(exe, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    for line in p.stdout:
        yield line

    # This is necessary to tell OS to wait for exit status or we will
    # get a zombie process
    os.wait()

def extract_softclipped_sites(bamfile, fastqDict):
	'''
	Identify candidate sites from softclipped positions
	in a bam file
	'''

	cmd = "samtools view -F2304 %s" %bamfile

	samgenerator = runProcess(cmd.split())

	for samline in samgenerator:
		currSamRead = samRead(samline.decode('ascii'))



		currReadsamcigar = sam_cigar(currSamRead.Cigar, currSamRead.Seq, \
			currSamRead.Chr, currSamRead.Posn)

		softclipped_left_genome_posn, softclipped_right_genome_posn \
		= currReadsamcigar.get_softclipped_genome_site()
		
		softclipped_left_end, softclipped_right_start \
		= currReadsamcigar.get_softclipped_posn()

		softclipped_left_seq, softclipped_right_seq \
		= currReadsamcigar.get_softclipped_seq()

		softclipped_chr = currSamRead.Chr
		softclipped_mapQ = currSamRead.mapQ 

		mate_sequence = ""
		if currSamRead.readNumInPair == 1:
			mate_sequence = fastqDict[currSamRead.readName][1]
		elif currSamRead.readNumInPair == 2:
			mate_sequence = fastqDict[currSamRead.readName][0]
		    


		if softclipped_left_genome_posn > 0:
			mate_seq_type = ""
			final_mate_sequence = ""
			currSamReadAlignmentOrientation = currSamRead.getAlignmentOrientation()
			if currSamReadAlignmentOrientation == "forward":
				final_mate_sequence = mate_sequence
				mate_seq_type = "nonextending_event"
			elif currSamReadAlignmentOrientation == "reverse":
				final_mate_sequence = reverseComplement(mate_sequence)
				mate_seq_type = "extending_event"

			result_left = [currSamRead.Chr, softclipped_left_genome_posn, \
			"left", softclipped_left_end, softclipped_left_seq, softclipped_mapQ,
			currSamRead.readName, currSamRead.readNumInPair, currSamReadAlignmentOrientation,
			mate_seq_type, final_mate_sequence]
			print("\t".join(map(str, result_left)))


		if softclipped_right_genome_posn > 0:
			mate_seq_type = ""
			final_mate_sequence = ""
			currSamReadAlignmentOrientation = currSamRead.getAlignmentOrientation()
			if currSamReadAlignmentOrientation == "forward":
				final_mate_sequence = reverseComplement(mate_sequence)
				mate_seq_type = "extending_event"
			elif currSamReadAlignmentOrientation == "reverse":
				final_mate_sequence = mate_sequence
				mate_seq_type = "nonextending_event"

			result_right = [currSamRead.Chr, softclipped_right_genome_posn, \
			"right", softclipped_right_start, softclipped_right_seq, softclipped_mapQ,
			currSamRead.readName, currSamRead.readNumInPair, currSamReadAlignmentOrientation,
			mate_seq_type, final_mate_sequence]
			print("\t".join(map(str, result_right)))



if __name__ == "__main__":
	bamfile = sys.argv[1]
	fastq1file = sys.argv[2]
	fastq2file = sys.argv[3]
	fastqDict = fastqToSeqDict(fastq1file, fastq2file)
	#print(len(fastqDict.keys()))
	extract_softclipped_sites(bamfile, fastqDict)




