import sys
import os
import subprocess
from os import sys, path
from sam_cigar import *
from lib.samReadClass import *


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

def extract_softclipped_sites(bamfile):
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



		if softclipped_left_genome_posn > 0:
			result_left = [currSamRead.Chr, softclipped_left_genome_posn, \
			"left", softclipped_left_end, softclipped_left_seq, softclipped_mapQ]
			print("\t".join(map(str, result_left)))

		if softclipped_right_genome_posn > 0:
			result_right = [currSamRead.Chr, softclipped_right_genome_posn, \
			"right", softclipped_right_start, softclipped_right_seq, softclipped_mapQ]
			print("\t".join(map(str, result_right)))



if __name__ == "__main__":
	bamfile = sys.argv[1]
	extract_softclipped_sites(bamfile)




