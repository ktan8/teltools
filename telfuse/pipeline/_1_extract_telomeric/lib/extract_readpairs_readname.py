#!/usr/bin/env python

import sys
import os
import subprocess
from .samReadClass import samRead

'''
Extract all the read pairs from a .bam file based 
on a given list of readnames
'''



def runProcess(exe):
	#p = subprocess.Popen(exe, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	p = subprocess.Popen(exe, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

	for line in p.stdout:
		yield line
	
	# This is necessary to tell OS to wait for exit status or we will
	# get a zombie process
	os.wait()


def get_readnames_list(file):
	'''
	Read in a list of readnames within a file
	and return that as a list.
	'''
	f = open(file, "r")
	readnameList = []
	for line in f:
		readname = line.strip()
		readnameList.append(readname)
	readnamesUnique = unique(readnameList)
	f.close()

	return readnamesUnique



def unique(seq):
	'''
	Get unique values of a list using sets
	'''
	seq_set = set(seq)
	return seq_set


def get_readnames_from_sam(samfile):
	'''
	Extract all the readnames from a sam file and return
	a list of unique readnames
	'''
	f = open(samfile, "r")
	readnameList = []
	for line in f:
		# Skip header lines
		if line[0] == "@":
			next
		else:
			currSamRead = samRead(line)
			readnameList.append(currSamRead.readName)
	#print(readnameList)
	readnamesUnique = unique(readnameList)
	#print(readnamesUnique)
	f.close()

	return readnamesUnique



def extract_reqr_reads_from_bam(bamfile, readNamesToExtract):
	'''
	Extract all required reads from a bamfile.
	This is a generator.
	'''
	bamfile = sys.argv[1]

	#posnString = str(chrom) + ':' + str(posnStart) + '-' + str(posnEnd)
	pileupCmd = "samtools view -h -F 2304 %s" %(bamfile)


	# Iterate over all lines in the sam file
	for samline in runProcess(pileupCmd.split()):
		samline = samline.decode()
		# Header line
		if samline[0] == '@':
			yield samline.strip()

		else:
			currSamRead = samRead(samline)
			if currSamRead.readName in readNamesToExtract:
				yield currSamRead.getSamOutput()

def print_samfile_with_reqr_reads(bamfile, samfile_withnames):
	reqr_readnames = get_readnames_from_sam(samfile_withnames)
	bam_generator = extract_reqr_reads_from_bam(bamfile, reqr_readnames)
	while True:
		try:
			print(next(bam_generator))
		except StopIteration:
			break


if __name__ == "__main__":
	bamfile = sys.argv[1]
	samfile_withnames = sys.argv[2]
	print_samfile_with_reqr_reads(bamfile, samfile_withnames)
