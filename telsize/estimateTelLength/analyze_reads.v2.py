#!/usr/bin/env python

import sys
import os
import argparse

def classify_site(telo_start, telo_end, read_length, bufferlen=100):
	'''
	Decide how to classify the telomeric reads
	'''
	classification = "NA"
	if telo_start == -1 and telo_end == -1:
		classification  = "non_telomeric"
	elif telo_start <= bufferlen and telo_end >= read_length - bufferlen:
		classification = "full_telomeric"
	elif telo_start <= bufferlen and telo_end < read_length - bufferlen:
		classification = "left_telomeric"
	elif telo_end >= read_length - bufferlen and telo_start > bufferlen:
		classification = "right_telomeric"
	elif telo_start > bufferlen and telo_end < read_length - bufferlen:
		classification = "intra_telomeric"

	return classification


def extract_nontelomeric_region(sequence, telo_start, telo_end, read_length):
	classification = classify_site(telo_start, telo_end, read_length)
	
	if classification == "left_telomeric":
		extracted_seq = [sequence[telo_end:read_length]]
	elif classification == "right_telomeric":
		extracted_seq = [sequence[0:telo_start]]
	elif classification == "intra_telomeric":
		extracted_seq_left = sequence[0:telo_start]
		extracted_seq_right = sequence[telo_end:read_length]
		extracted_seq = [extracted_seq_left, extracted_seq_right]
	else:
		extracted_seq = ["NNNNN"]

	return classification, extracted_seq

def generate_nontelomeric_fasta(telomere_length_file, outputfile):
	f = open(telomere_length_file, "r")
	#outputfile = telomere_length_file + ".nontelomeric.fasta"
	out = open(outputfile, "w")
	
	header = f.readline()
	for line in f.readlines():
		lineArr = line.strip().split("\t")
		#print(lineArr)
		readname   = ">" + lineArr[0]
		teloStart  = int(lineArr[1])
		teloEnd	   = int(lineArr[2])
		readLength = int(lineArr[6])
		sequence   = lineArr[9]
		
		classification, nontelomeric_sequence = extract_nontelomeric_region(sequence, teloStart, teloEnd, readLength)
		#print(nontelomeric_sequence)
		if len(nontelomeric_sequence) == 1:
			readname_new = readname + ":nontelomeric:" + classification
			out.write(readname_new + "\n")
			out.write(nontelomeric_sequence[0] + "\n")
		elif len(nontelomeric_sequence) == 2:
			readname_new_left = readname + ":nontelomeric:" + classification + "_left"
			out.write(readname_new_left + "\n")
			out.write(nontelomeric_sequence[0] + "\n")
			readname_new_right = readname + ":nontelomeric:" + classification + "_right"
			out.write(readname_new_right + "\n")
			out.write(nontelomeric_sequence[1] + "\n")
			
	f.close()
	out.close()

def map_nontemeric_fasta(non_telomere_fasta, genome, paf_file, platform="ont"):
	#paf_file = non_telomere_fasta + ".map_genome.paf"
	if platform == "ont":
		os.system("minimap2 -c -x map-ont -t 16 %s %s > %s" %(genome, non_telomere_fasta, paf_file)) # -c argument for base level alignment
	elif platform == "pb":
		os.system("minimap2 -c -x map-pb -t 16 %s %s > %s" %(genome, non_telomere_fasta, paf_file)) # -c argument for base level alignment
	

def map_full_fasta(full_fasta, genome):
	pass
 

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Extract non-telomeric sequence and map it to reference genome',
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('telLenFile', metavar='teloLengthFile', type=str, nargs=1,
		help='Telomere length file from teltools')
	parser.add_argument('outputlabel', metavar='outputlabel', type=str, nargs=1,
		help='Label for output file')
	parser.add_argument('genome', metavar='genome', type=str, nargs="?",
		help='Reference genome to map to (mapping will be skipped if genome is not specified)')
	parser.add_argument('-x', metavar='platform', type=str, default="ont", choices=['ont', 'pb'],
		help='Sequencing platform for long-reads (ont/pb)')

	args = parser.parse_args()


	#teloLengthFile  = sys.argv[1]
	#genome 		= sys.argv[2]
	teloLengthFile  = args.telLenFile[0]
	#genome 		= args.genome[0]
	#generate_nontelomeric_fasta(teloLengthFile)

	#genome_hg38 = "~/kartong/genome/Homo_sapiens_assembly38.fasta"
	#genome_chm13 = "~/kartong/genome/chm13.draft_v1.0.fasta"
	#non_telomere_fasta = teloLengthFile + ".nontelomeric.fasta"
	non_telomere_fasta = args.outputlabel[0] + ".nontelomeric.fasta"
	generate_nontelomeric_fasta(teloLengthFile, non_telomere_fasta)

	#map_hg38_paf = non_telomere_fasta + ".map_hg38.paf"
	#map_chm13_paf = non_telomere_fasta + ".map_chm13.paf"
	map_genome_paf = non_telomere_fasta + ".map_genome.paf"
	
	#map_nontemeric_fasta(non_telomere_fasta, genome_hg38, map_hg38_paf)
	#map_nontemeric_fasta(non_telomere_fasta, genome_chm13, map_chm13_paf)

	if args.genome is not None:
		map_nontemeric_fasta(non_telomere_fasta, args.genome[0], map_genome_paf, platform=args.platform)
	
# To be modified to:
# 1) introduce argparse
# 2) introduce genome as an argument
# 3) introduce parameters for different types of data


