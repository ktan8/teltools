#!/usr/bin/python

import sys
import os
import subprocess
from lib.identify_telomeric_sequence import *

def runProcess(exe):
	p = subprocess.Popen(exe, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

	for line in p.stdout:
		yield line
	
	# This is necessary to tell OS to wait for exit status or we will
	# get a zombie process
	os.wait()



class fastq_chunk():
	'''
	A fastq chunk class that allows us to easily
	store and access information on a fastq read
	'''

	def __init__(self, fastq_iterator):
		self.readname = next(fastq_iterator)
		self.sequence = next(fastq_iterator)
		self.qualname = next(fastq_iterator)
		self.quality = next(fastq_iterator)

	def generate_fastq_string(self):
		'''
		Generate the fastq info from the
		fastq chunk info.
		'''
		fastq_string = self.readname + self.sequence + self.qualname + self.quality

		return fastq_string



def process_fastq_gz(fastq1_gz, fastq2_gz):
	pass


def write_fastq_to_file(fastq_chunk, file_handle):
	'''
	Write fastq to file handle
	'''
	fastq_string_output = fastq_chunk.generate_fastq_string()
	file_handle.write(fastq_string_output)


fastq1 = sys.argv[1]
fastq2 = sys.argv[2]
label = sys.argv[3]
telomere_repeat_count = 6


# fastq1_out = "test.fq1"
# fastq2_out = "test.fq2"

fastq1_out = label + ".fq1"
fastq2_out = label + ".fq2"

fastq1_out_handle = open(fastq1_out, "w")
fastq2_out_handle = open(fastq2_out, "w")


fastq1Cmd = "gunzip -c -d %s" %(fastq1)
fastq2Cmd = "gunzip -c -d %s" %(fastq2)

fastq1_iter = runProcess(fastq1Cmd.split())
fastq2_iter = runProcess(fastq2Cmd.split())



while True:
	try:
		fastq1_chunk = fastq_chunk(fastq1_iter)
		fastq2_chunk = fastq_chunk(fastq2_iter)
		#print fastq1_chunk.sequence

		if check_telomeric_sequence(fastq1_chunk.sequence, repeat_count=telomere_repeat_count) or \
		check_telomeric_sequence(fastq2_chunk.sequence, repeat_count=telomere_repeat_count):
			write_fastq_to_file(fastq1_chunk, fastq1_out_handle)
			write_fastq_to_file(fastq2_chunk, fastq2_out_handle )

	except StopIteration:
		break

# Close the file handle to ensure that file
# that is compressed is the complete file
fastq1_out_handle.close()
fastq2_out_handle.close()

# Compress the fastq file into gz
os.system("gzip %s" %fastq1_out)
os.system("gzip %s" %fastq2_out)

#print("gzip %s" %fastq1_out)
