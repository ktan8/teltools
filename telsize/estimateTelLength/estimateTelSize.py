#!/usr/bin/python

import sys
import os
import time
import argparse
import multiprocessing
from estimateTelomereLength import *
from analyze_reads import classify_site

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
			if name: yield (name[1:], ''.join(seq))
			name, seq = line, []
		else:
			seq.append(line)
	if name: yield (name[1:], ''.join(seq))





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


def analyze_read_parallel(data, env_params):
	(seq_name, seq) = data
	(short_motif, long_motif, long_motif_revCom, pic_format, \
        AverageTelomereSignalCutoff, movingAveWindow, medianKernelSize, folder, plot_fig) = env_params

	# Skip sequences without the motif	
	if not (check_motif_present(seq, long_motif) or check_motif_present(seq, long_motif_revCom)):
		TelomereStart = -1
		TelomereEnd = -1
		TelomereLength = 0
		TelomereRegions = 0
		AverageTelomereSignal = "NA"
		stringLength = len(seq)
		strand = "NA"
		classification = "non_telomeric"
		string = seq

		result = (seq_name, TelomereStart, TelomereEnd, TelomereLength, TelomereRegions,
		AverageTelomereSignal, stringLength, strand, classification, string)

		return result

	(TelomereStart, TelomereEnd, TelomereLength, TelomereRegions,
	AverageTelomereSignal, stringLength, strand, string) \
	= analyze_sequence_for_motif(seq, short_motif, seq_name, pic_format=pic_format, \
	threshold=AverageTelomereSignalCutoff, movingAveWindow=movingAveWindow, \
	medianKernelSize=medianKernelSize, folder=folder, plot_fig=plot_fig)


	classification = classify_site(TelomereStart, TelomereEnd, stringLength, bufferlen=50)
	result = (seq_name, TelomereStart, TelomereEnd, TelomereLength, TelomereRegions,
	AverageTelomereSignal, stringLength, strand, classification, string)

	return result

		


def main(fasta_file, short_motif = "TTAGGG", pic_format="png", AverageTelomereSignalCutoff=0.2,
	movingAveWindow=50, medianKernelSize=501, threads=12, folder="./", noseq=False,
	penalty=False, penaltyval=0.5, plot_fig=True):

	long_motif = short_motif + short_motif
	long_motif_revCom = reverseComplement(long_motif)
	
	total = 0
	hits = 0
	#fasta_file_handle = open(fasta_file, 'r')
	fasta_file_handle = os.popen("zcat -f %s" %fasta_file) # to support gzip file


	# Header
	header = ()
	if noseq:
		header = ("readname", "TelomereStart", "TelomereEnd", "TelomereLength", "TelomereRegions", 
		"AverageTelomereSignal", "stringLength", "strand", "classification")
	else:
		header = ("readname", "TelomereStart", "TelomereEnd", "TelomereLength", "TelomereRegions",
		"AverageTelomereSignal", "stringLength", "strand", "classification", "sequence")
	print("\t".join(header))

	# Create folder for output
	if not os.path.exists(folder):
		os.mkdir(folder)


	##########################
	### Replaced by parallel code
	##########################
#	for seq_name, seq in read_fasta(fasta_file_handle):
#		total += 1
#
#		# Heuristic to remove reads without TTAGGG (at least 2 TTAGGG present)
#		if not (check_motif_present(seq, long_motif) or check_motif_present(seq, long_motif_revCom)):
#			continue
#
#		(TelomereStart, TelomereEnd, TelomereLength, TelomereRegions, 
#			AverageTelomereSignal, stringLength, strand, string) \
#		= analyze_sequence_for_motif(seq, short_motif, seq_name, pic_format=pic_format, threshold=AverageTelomereSignalCutoff, movingAveWindow=movingAveWindow, medianKernelSize=medianKernelSize, folder=folder)		
#
#
#		# Print results
#		classification = classify_site(TelomereStart, TelomereEnd, stringLength, bufferlen=50)
#		result = (seq_name, TelomereStart, TelomereEnd, TelomereLength, TelomereRegions, 
#			AverageTelomereSignal, stringLength, strand, classification, string)
#		print("\t".join(map(str, result)))
#
#		if AverageTelomereSignal != "NA" and AverageTelomereSignal > AverageTelomereSignalCutoff:
#			hits += 1
	#########################
	#########################
	reads = [(seq_name, seq) for seq_name, seq in read_fasta(fasta_file_handle)]
	env_params = [(short_motif, long_motif, long_motif_revCom, pic_format, \
	AverageTelomereSignalCutoff, movingAveWindow, medianKernelSize, folder, plot_fig)] * len(reads)

	pool = multiprocessing.Pool(threads)
	results = pool.starmap(analyze_read_parallel, zip(reads, env_params))
	pool.close()
	pool.join()
	#print(results)
	#print("\t".join(map(str, [line_list for line_list in results])))
	#"\t".join(map(str, line_list)) for line_list in results
	for line_list in results:
		if noseq:
			result_list = line_list[:-1]
		else:
			result_list = line_list
		print("\t".join(map(str, result_list)))
		
	

	# Print results



	# print("Total sequences: " + str(total))
	# print("Total hits: " + str(hits))



if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Get length for repeats in fasta',
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	
	currtime = int(time.time())
	image_dir = "./img_" + str(currtime) + "/"


	parser.add_argument('fasta', metavar='fasta', type=str, nargs=1,
                    help='fasta file to evaluate (unzipped|gzipped)')
	parser.add_argument('--cutoff', metavar='cutoff', type=float, default=0.35,
					help='cutoff of repeat signal')
	parser.add_argument('--movave', metavar='movave', type=int, default=50,
                                        help='Window size for moving average window')
	parser.add_argument('--movmed', metavar='movmed', type=int, default=501,
                                        help='Window size for moving median window')
	parser.add_argument('--motif', metavar='motif', type=str, default='TTAGGG',
					help='motif to assess for repeat signal (e.g. TTAGGG)')
	parser.add_argument('--format', metavar='format', type=str, default='pdf',
					help='output file type (e.g. pdf, png)')
	parser.add_argument('--folder', metavar='folder', type=str, default=image_dir,
                                        help='folder for images')
	parser.add_argument('--noseq', action='store_true',
					help='do not print sequence of read in output')
	parser.add_argument('--nofig', action='store_false',
					help='do not plot figures for sequences')
	#parser.add_argument('--fastq', metavar='fastq', type=str, default='./',
        #                                help='folder for images')
	parser.add_argument('--threads', metavar='threads', type=int, default=12,
                                        help='number of threads to use')
	parser.add_argument('--penalty', action='store_true',
					help='calculate telomere length using penalty based method')
	parser.add_argument('--penaltyval', metavar='movave', type=int, default=0.5,
					help='penalty value to use for penalty method')
	args = parser.parse_args()


	main(args.fasta[0], pic_format=args.format, AverageTelomereSignalCutoff=args.cutoff,
		short_motif = args.motif, movingAveWindow=args.movave, medianKernelSize=args.movmed,
		folder=args.folder, threads = args.threads, noseq=args.noseq,
		penalty=args.penalty, penaltyval=args.penaltyval, plot_fig=args.nofig)

