#!/usr/bin/env python
import sys
from collections import defaultdict
from lib.find_telo_content import *


def Average(lst):
	return sum(lst) / len(lst)


class softclipped_reads:
	'''
	A class to store information on all
	the softclipped reads with the same
	softclipped position and orientation.
	'''

	def __init__(self, list_of_softclippedseqs, list_of_readanmes):
		#self.softclippedseqs = ["CCCAACCCCAACCCCAACA", "CCCAACCCCAACCCCACCCCC", "CCCAACCCCGACCCCCAACCCCACCCCAACCCCCCGCCCACCCCCGACCCCCAACTACCCCCA"]
		self.softclippedseqs = list_of_softclippedseqs
		self.readnames	     = list_of_readanmes
		#print(list_of_softclippedseqs)
		self.num_seqs = len(self.softclippedseqs)
		#self.num_seqs_dedup = self.count_numseq_dedup()
		self.num_seqs_dedup = self.count_numseq_dedup_readpairs()
		self.longest_seq_len = self.longest_seq_len(self.softclippedseqs)
		self.consensus_seq = self.get_consensus_seq()

		# self.
		# pass

	def most_frequent(self, List):
		'''
		Identify the most frequent item in a list.
		'''
		return max(set(List), key = List.count)


	def longest_seq_len(self, seqList):
		'''
		Identify the longest sequence in a list
		of sequences.
		'''
		max_seqlen = -1
		for seq in seqList:
			if len(seq) > max_seqlen:
				max_seqlen = len(seq)
		return max_seqlen

	def count_numseq_dedup(self):
		'''
		Count the number of non-duplicated
		reads. This does not consider if the 
		reads are pointing to the left or to
		the right. This is done by assuming that 
		if the softclipped length is the same,
		then the reads are duplicates
		'''
		seq_len_dict = defaultdict(int)
		for seq in self.softclippedseqs:
			seq_len_dict[len(seq)] += 1
		unique_seq_count = len(seq_len_dict.keys())

		return unique_seq_count

	def count_numseq_dedup_readpairs(self):
		'''
		Count the number of non-duplicated read pairs.
		Look through all the readnames and then if
		there are two reads with same readname, 
		'''

		seq_len_dict = defaultdict(int)
		readname_dict = dict()

		# Take the longest len for softclipped seq
		# if the same readname occurs twice
		for i in range(len(self.softclippedseqs)):
			readname = self.readnames[i]
			seq = self.softclippedseqs[i]
			seqLen = len(seq)
			if readname in readname_dict:
				newSeqLen = max(readname_dict[readname], seqLen)
				readname_dict[readname] = newSeqLen
			else:
				readname_dict[readname] = seqLen
		
		# Then check these reads for duplicates
		for seqLen in readname_dict.values():
			seq_len_dict[seqLen] += 1
		unique_seq_count = len(seq_len_dict.keys())

		return unique_seq_count


	def get_consensus_seq(self):
		'''
		Get the consensus position based on
		all the soft clipped reads. In this case,
		we assume that all the softclipped sequences
		are properly left aligned, and that there are
		no indels which might mess up the alignment.
		Even if there are indels in a few rare reads, 
		'''
		consensus_seq = ""

		for seq_posn in range(self.longest_seq_len):
			chars_posns = []
			for seq_num in range(self.num_seqs):
				try:
					curr_char = self.softclippedseqs[seq_num][seq_posn]
					chars_posns.append(curr_char)
				except:
					pass

			# Get consensus character for each position and append
			# it to the overall consensus sequence.
			consensus_character = self.most_frequent(chars_posns)
			consensus_seq = consensus_seq + consensus_character

		return consensus_seq


	def calc_sequence_identity(self):
		'''
		Calculate the sequence identify of each
		of the sequences with respect to the consensus
		sequence we derived from the sequences we have
		'''
		total_bases = 0
		total_matched_bases = 0
		seq_identity_all = []
		for seq_posn in range(self.longest_seq_len):
			# chars_posns = []
			matched_base_at_posn_count = 0
			total_base_at_posn_count = 0 

			for seq_num in range(self.num_seqs):
				try:
					curr_char = self.softclippedseqs[seq_num][seq_posn]
					# chars_posns.append(curr_char)

					# curr_char
					if curr_char == self.consensus_seq[seq_posn]:
						matched_base_at_posn_count += 1
						total_matched_bases += 1

					total_base_at_posn_count += 1
					total_bases += 1
				except:
					pass

			seq_identity_curr_posn = matched_base_at_posn_count / float(total_base_at_posn_count)
			# print(seq_identity_curr_posn)
			seq_identity_all.append(seq_identity_curr_posn)


		ave_seq_identity = sum(seq_identity_all) / self.longest_seq_len
		ave_weighted_seq_identity = float(total_matched_bases) / total_bases

		return ave_seq_identity, ave_weighted_seq_identity




class softclipped_aggregated_sites:
	def __init__(self):
		self.softclipped_dict = dict()
		self.querySite_dict = dict()
		self.mapQ_dict = dict()
		self.readinfo_dict = dict()


	def append(self, softclipline):
		'''
		Add the current soft clipped site infor to 
		this object
		'''
		lineList 	= softclipline.strip().split("\t")
		chrom 		= lineList[0]
		posn 		= lineList[1]
		orientation	= lineList[2]
		clippedQuerySite = lineList[3]
		clippedQuerySeq  = lineList[4]
		mapQ		= lineList[5]
		readname	= lineList[6]
		readNuminPair	= lineList[7]
		alignmentOrientation	= lineList[8]
		mate_seq_type	= lineList[9]
		final_mate_sequence = lineList[10]

		dict_key = str(chrom) + "|" + str(posn) + "|" + str(orientation)
		readinfo = [readname, readNuminPair, alignmentOrientation, mate_seq_type, final_mate_sequence]

		if dict_key in self.softclipped_dict:
			self.softclipped_dict[dict_key].append(clippedQuerySeq)
			self.querySite_dict[dict_key].append(clippedQuerySite)
			self.mapQ_dict[dict_key].append(mapQ)
			self.readinfo_dict[dict_key].append(readinfo)
		else:
			#print(list(clippedQuerySeq))
			initial_clippedQuerySeq_list = []
			initial_clippedQuerySeq_list.append(clippedQuerySeq)
			#self.softclipped_dict[dict_key] = list(clippedQuerySeq)
			self.softclipped_dict[dict_key] = initial_clippedQuerySeq_list

			initial_clippedQuerySite_list = []
			initial_clippedQuerySite_list.append(clippedQuerySite)
			self.querySite_dict[dict_key] = initial_clippedQuerySite_list

			initial_mapQ_list = []
			initial_mapQ_list.append(mapQ)
			self.mapQ_dict[dict_key] = initial_mapQ_list

			initial_readinfo_list = []
			initial_readinfo_list.append(readinfo)
			self.readinfo_dict[dict_key] = initial_readinfo_list
			
	
	def average_softclip_read_posn(self, dict_site_key):
		'''
		Calculate the average position of the softclipped
		site on the read.
		'''
		site_posn_list = self.querySite_dict[dict_site_key]
		average_read_posn = sum(map(int, site_posn_list)) / float(len(site_posn_list))

		return average_read_posn


	def mapQ_stats(self, dict_site_key):
		'''
		Calculate the statistics for the mapQ
		i.e. average, min, max
		'''
		site_mapQ_list = self.mapQ_dict[dict_site_key]
		average_mapQ = sum(map(int, site_mapQ_list)) / float(len(site_mapQ_list))
		min_maqQ = min(site_mapQ_list)
		max_mapQ = max(site_mapQ_list)

		return average_mapQ, min_maqQ, max_mapQ


	def calc_telomeric_content(self, dict_site_key, motif="TTAGGG"):
		'''
		Calculate the telomeric content of the other read (mate)
		in the readpair that falls in the softclipped region.
	        '''
		site_readinfo_list = self.readinfo_dict[dict_site_key]
		required_seq = []
		for read in site_readinfo_list:
			mate_seq_type = read[3]
			mate_sequence = read[4]
			
			if mate_seq_type == "extending_event":
				required_seq.append(mate_sequence)
		
		teloContent_fwd_list = []
		teloContent_rev_list = []
		for sequence in required_seq:
			(telmetric_fwd, telmetric_rev) = calc_telo_content_both(sequence, motif)
			[seq_len, motif_len_in_seq, telo_content_fwd] = telmetric_fwd
			[seq_len, motif_len_in_seq, telo_content_rev] = telmetric_rev
			teloContent_fwd_list.append(telo_content_fwd)
			teloContent_rev_list.append(telo_content_rev)


		extending_reads = len(teloContent_fwd_list)
		teloContent_fwd_mean = "NA"
		teloContent_rev_mean = "NA"
		if teloContent_fwd_list:
			teloContent_fwd_mean = Average(teloContent_fwd_list)
			teloContent_rev_mean = Average(teloContent_rev_list)

		
		return extending_reads, teloContent_fwd_mean, teloContent_rev_mean


	def get_all_readnames(self, dict_site_key):
		'''
		Get the readnames of all reads
		'''
		readnameList = []
		readinfoList = self.readinfo_dict[dict_site_key]
		for read in readinfoList:
			readname = read[0]
			readnameList.append(readname)

		return readnameList
		

	def softclipped_analysis_generator(self, motif="TTAGGG"):
		'''
		Generator for the softclipped read
		sequences. Will take in 
		'''

		for key in self.softclipped_dict:
			softclippedseqsforPosn = self.softclipped_dict[key]
			readnamesList = self.get_all_readnames(key)
			softclippedseqsforPosnObj = softclipped_reads(softclippedseqsforPosn, readnamesList)


			num_seqs_dedup = softclippedseqsforPosnObj.num_seqs_dedup
			consensus_seq = softclippedseqsforPosnObj.consensus_seq
			ave_seq_identity, ave_weighted_seq_identity = \
			softclippedseqsforPosnObj.calc_sequence_identity()
			softclip_ave_posn_on_read = self.average_softclip_read_posn(key)
			[ave_mapQ, min_mapQ, max_mapQ] = self.mapQ_stats(key)

			extendingReads, teloContentFwd, teloContentRev = self.calc_telomeric_content(key, motif=motif)


			# print(key)
			chrom, posn, orientation = key.split("|")
			result = [chrom, posn, orientation, num_seqs_dedup, ave_seq_identity, \
			ave_weighted_seq_identity, consensus_seq, softclip_ave_posn_on_read, \
			ave_mapQ, min_mapQ, max_mapQ, extendingReads, teloContentFwd, teloContentRev]

			yield result




# num_reads
# seq_identity
# weighted_seq_identity
# seq_length
# number_of_bases


# # Generate a dict of softclipped class.
# # Each soft clipped
# dict


# a = softclipped_reads()
# print(a.get_consensus_seq())
# print(a.calc_sequence_identity())


site_file = sys.argv[1]
#site_file = "HCC1954_illumina.alltelomeric.readname.hg38.softclipped.sites"

site_file_handle = open(site_file, "r")

softclipped_group = softclipped_aggregated_sites()

# Aggregate the sites by unique coordinates
for line in site_file_handle:
	softclipped_group.append(line)
site_file_handle.close()

# Convert the aggregated sites into a generator
# which can be readily analyzed.
softclipped_generator = softclipped_group.softclipped_analysis_generator()

for result in softclipped_generator:
	print("\t".join(map(str, result)))
