#!/usr/bin/python

import sys
import os
import subprocess
from lib.identify_telomeric_sequence import *
from lib.extract_readpairs_readname import *

def runProcess(exe):
	p = subprocess.Popen(exe, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

	for line in p.stdout:
		yield line
	
	# This is necessary to tell OS to wait for exit status or we will
	# get a zombie process
	os.wait()



class sam_chunk():
	'''
	A sam chunk class that allows us to easily
	store and access information on a bam read
	'''

	def __init__(self, sam_iterator):
		self.samline = next(sam_iterator).strip()
		self.samline = self.samline.decode()
		self.samlist = self.samline.split("\t")
		self.readname = self.samlist[0]
		self.sequence = self.samlist[9]


	def generate_sam_string(self):
		'''
		Generate the sam info
		'''
		return self.samline



def write_sam_to_file(sam_chunk, file_handle):
	'''
	Write sam to file handle
	'''
	sam_string_output = sam_chunk.generate_sam_string() + "\n"
	file_handle.write(sam_string_output)


def write_samfile_with_reqr_reads(bamfile, samfile_withnames, file_handle):
	'''
	Write samfile with all required readnames
	'''
	reqr_readnames = get_readnames_from_sam(samfile_withnames)
	bam_generator = extract_reqr_reads_from_bam(bamfile, reqr_readnames)
	while True:
		try:
			outputline = next(bam_generator) + "\n"
			file_handle.write(outputline)
		except StopIteration:
			break



def extract_telomeric_bam(bamfile, label):
	'''
	Extract telomeric reads from bam file
	'''

	fastq1_out = label + ".fq1"
	fastq2_out = label + ".fq2"

	# File names of intermediate sam files
	telomere_tmp_candidates_sam = label + ".telomere.tmp.candidates.sam"
	telomere_tmp_allreads_sam = label + ".telomere.tmp.allreads.sam"


	# File names of intermediate bam files
	allreads_bam = label + ".telomere.tmp.allreads.bam"
	allreads_sortedname_bam = label + ".telomere.tmp.allreads.sortedname.bam"


	# Initialize some file handles we need
	telomere_tmp_candidates_sam_handle = open(telomere_tmp_candidates_sam, "w")
	telomere_tmp_allreads_sam_handle = open(telomere_tmp_allreads_sam, "w")

	# Define intermediate files that we will remove
	intermediate_files = [telomere_tmp_candidates_sam, telomere_tmp_allreads_sam, allreads_bam, allreads_sortedname_bam]


	bamCmd = "samtools view -F2304 %s" %(bamfile)
	bam_iter = runProcess(bamCmd.split())



	while True:
		try:
			sam_chunk_curr = sam_chunk(bam_iter)

			if check_telomeric_sequence(sam_chunk_curr.sequence):
				#print("apple")
				write_sam_to_file(sam_chunk_curr, telomere_tmp_candidates_sam_handle)

		except StopIteration:
			break

	# Close the file handle to ensure that file
	# that is compressed is the complete file
	telomere_tmp_candidates_sam_handle.close()



	# Extract the required reads from the bam file, using the names 
	# from the sam file, and write it to the file handle
	write_samfile_with_reqr_reads(bamfile, telomere_tmp_candidates_sam, telomere_tmp_allreads_sam_handle)
	telomere_tmp_allreads_sam_handle.close()



	os.system("samtools view -bS %s > %s" %(telomere_tmp_allreads_sam, allreads_bam))
	os.system("samtools sort -n %s > %s" %(allreads_bam, allreads_sortedname_bam))
	os.system("samtools fastq -1 %s -2 %s %s" %(fastq1_out, fastq2_out, allreads_sortedname_bam))


	# Compress the fastq file into gzipped files
	os.system("gzip %s" %fastq1_out)
	os.system("gzip %s" %fastq2_out)


	# Remove intermediate files
	for file in intermediate_files:
		os.system("rm %s" %file)


if __name__ == "__main__":
	bamfile = sys.argv[1]
	label 	= sys.argv[2]

	extract_telomeric_bam(bamfile, label)
