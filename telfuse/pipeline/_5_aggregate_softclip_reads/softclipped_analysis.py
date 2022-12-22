#!/usr/bin/env python
import sys

class softclipped_reads:
	'''
	A class to store information on all
	the softclipped reads with the same
	softclipped position and orientation.
	'''

	def __init__(self, list_of_softclippedseqs):
		#self.softclippedseqs = ["CCCAACCCCAACCCCAACA", "CCCAACCCCAACCCCACCCCC", "CCCAACCCCGACCCCCAACCCCACCCCAACCCCCCGCCCACCCCCGACCCCCAACTACCCCCA"]
		self.softclippedseqs = list_of_softclippedseqs
		#print(list_of_softclippedseqs)
		self.num_seqs = len(self.softclippedseqs)
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

	def append(self, softclipline):
		'''
		Add the current soft clipped site infor to 
		this object
		'''
		lineList 	= softclipline.strip().split("\t")
		chrom 		= lineList[0]
		posn 		= lineList[1]
		orientation = lineList[2]
		clippedQuerySite = lineList[3]
		clippedQuerySeq  = lineList[4]
		mapQ		= lineList[5]

		dict_key = str(chrom) + "|" + str(posn) + "|" + str(orientation)

		if dict_key in self.softclipped_dict:
			self.softclipped_dict[dict_key].append(clippedQuerySeq)
			self.querySite_dict[dict_key].append(clippedQuerySite)
			self.mapQ_dict[dict_key].append(mapQ)
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


	def softclipped_analysis_generator(self):
		'''
		Generator for the softclipped read
		sequences. Will take in 
		'''

		for key in self.softclipped_dict:
			softclippedseqsforPosn = self.softclipped_dict[key]
			softclippedseqsforPosnObj = softclipped_reads(softclippedseqsforPosn)


			num_seqs = softclippedseqsforPosnObj.num_seqs
			consensus_seq = softclippedseqsforPosnObj.consensus_seq
			ave_seq_identity, ave_weighted_seq_identity = \
			softclippedseqsforPosnObj.calc_sequence_identity()
			softclip_ave_posn_on_read = self.average_softclip_read_posn(key)
			[ave_mapQ, min_mapQ, max_mapQ] = self.mapQ_stats(key)

			# print(key)
			chrom, posn, orientation = key.split("|")
			result = [chrom, posn, orientation, num_seqs, ave_seq_identity, \
			ave_weighted_seq_identity, consensus_seq, softclip_ave_posn_on_read, \
			ave_mapQ, min_mapQ, max_mapQ]

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
