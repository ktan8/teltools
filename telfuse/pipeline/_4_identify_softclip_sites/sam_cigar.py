#!/usr/bin/env python

# def calc_discriminat(a, b, c):

import re
from lib.samReadClass import *

consume_query = ['M', 'I', 'S', '=', 'X']
consume_ref = ['M', 'D', 'N', '=', 'X']

class sam_cigar():
	'''
	
	'''
	def __init__(self, cigar_str, readseq, ref_chr, ref_posn):
		self.cigar_str = cigar_str
		self.readseq = readseq
		self.ref_chr = ref_chr
		self.ref_posn = int(ref_posn)


		self.cigar_list = self.analyze_cigar()
		# self.has_softclipped = 
		# self.left_softclipped = 
		# self.right_softclipped =
		self.readlen = len(self.readseq)

		# Indicates length, and position of mapping on 
		# reference genome. Start position is 1-based 
		# inclusive, end position is 1-based.
		# self.ref_len = 
		# self.ref_start = self.position
		# self.ref_end = self.position + self.ref_len


		# Indicates length, and position of query (read) 
		# which has aligned on the reference genome. 
		# Query start position is 0-based inclusive,
		# end position is 0-based exclusive.
		# self.query_len = 
		# self.query_start =
		# self.query_end = 
		

	def analyze_cigar(self):
		match = re.findall(r'(\d+)(\w)', self.cigar_str)

		# Clean up the numeric character and make sure
		# that it is a int.
		match_clean = []
		for pos, char in match:
			match_clean.append([int(pos), char])
		
		return match_clean

	def get_ref_info(self):
		pass

	# def get_query_info(self):
	# 	query_len = 0
	# 	query_start = 0
	# 	query_end = 0

	# 	query_first_char_size = self.cigar_list[0][0]
	# 	query_first_char_char = self.cigar_list[0][1]
	# 	query_last_char_size = self.cigar_list[-1][0]
	# 	query_last_char_char = self.cigar_list[-1][1]

	# 	if query_first_char_char == 'S':
	# 		query_start = query_first_char_size # because it is zerobased
	# 	else:
	# 		query_start = 0

	# 	if quert_last_char_char == 'S':
	# 		query_end
			
		
	# 	for (cigar_size, cigar_char) in self.cigar_list:
	# 		if cigar_char 
			
	# 	pass


	def get_mapping_length(self):
		'''
		Get the span of mapping of a 
		particular read.
		'''
		mapping_len = 0
		for (cigar_size, cigar_char) in self.cigar_list:
			if cigar_char in consume_ref:
				mapping_len += cigar_size

		return mapping_len


	def get_mapping_posn(self):
		'''
		Get the fullspan of mapping of the read on
		the reference genome. Provides the start and
		stop genomic coordinates in 1-based inclusive, 
		and 1-based exclusive coordinates
		'''
		map_start = self.ref_posn
		map_end = self.ref_posn + self.get_mapping_length()


		return(map_start, map_end)




	def get_softclipped_posn(self):
		'''
		Get the position of the softclipped
		sites in the read. Start and end positions
		are provided in 0-based inclusive, and 
		0-based exclusive coordiates respectively.
		(Mainly to ease sequence extraction)

		If an end of a read has no softclipping, the
		coordinates representing this end of the
		read will be indicated as -1.
		'''
		#print(self.cigar_list)
		query_first_char_size = self.cigar_list[0][0]
		query_first_char_char = self.cigar_list[0][1]
		query_last_char_size = self.cigar_list[-1][0]
		query_last_char_char = self.cigar_list[-1][1]

		if query_first_char_char == 'S':
			# softclipped_left_start = 0
			softclipped_left_end = query_first_char_size
		else:
			# Set left_end as zero to essentially
			# say that there is no softclip on the left
			# softclipped_left_start = 0
			softclipped_left_end = -1
		
		if query_last_char_char == 'S':
			softclipped_right_start = self.readlen - query_last_char_size
			# softclipped_right_end = self.readlen
		else:
			# Set right_end as zero to essentially
			# say that there is no softclip on the right
			softclipped_right_start = -1
			# softclipped_right_end = 0

		# return softclipped_left_start, softclipped_left_end, \
		# softclipped_right_start, softclipped_right_end

		return softclipped_left_end, softclipped_right_start


	def get_softclipped_seq(self):
		'''
		Get softclipped sequence.
		For left softclipped site, we are going to report
		the reverse complement sequence for ease of alignment
		and comparison. For right softclipped site, we will
		report the soft clipped sequence as it is.
		'''
		softclipped_left_start = 0
		softclipped_right_end = self.readlen

		softclipped_left_end, softclipped_right_start \
		= self.get_softclipped_posn()
		# print [softclipped_left_start,softclipped_left_end]
		# print self.readseq[softclipped_left_start:softclipped_left_end]

		softclipped_left_seq = self.readseq[softclipped_left_start:softclipped_left_end]
		softclipped_right_seq = self.readseq[softclipped_right_start:softclipped_right_end]


		softclipped_left_seq_revcom = reverseComplement(softclipped_left_seq)

		return softclipped_left_seq_revcom, softclipped_right_seq


	def get_softclipped_genome_site(self):
		'''
		Determine which site in the genome is 
		softclipped.

		An example:

		Read: ATTAGGTCCCCCATTTAGGGA
                  ||||||||||||||
		Ref : ATTAGGTCCCCCATTTAGGGA
		          |             |
				20001(L)      20015(R)

		In the example above, we will report
		the genomic coordinate as 20001 as
		softclipped left, and 20015 as soft
		clipped right position. Thus, if there
		is softclipping on both the left, and right,
		softclipped coords == mapping coords

		If there is no softclipping, we will report
		the softclipping coordinates as -1
		'''

		softclipped_left_genome_posn = -1
		softclipped_right_genome_posn = -1


		# softclipped_left_start, softclipped_left_end, \
		# softclipped_right_start, softclipped_right_end = \
		# self.get_softclipped_posn()

		softclipped_left_end, softclipped_right_start \
		= self.get_softclipped_posn()

		map_start, map_end = self.get_mapping_posn()

		# if softclipped_left_end > softclipped_left_start:
		# 	softclipped_left_genome_posn = map_start
		# if softclipped_right_end > softclipped_right_start:
		# 	softclipped_right_genome_posn = map_end


		if softclipped_left_end >= 0:
			softclipped_left_genome_posn = map_start
		if softclipped_right_start >= 0:
			softclipped_right_genome_posn = map_end

		return softclipped_left_genome_posn, softclipped_right_genome_posn

	
if __name__ == "__main__":

	# cigar_str = '2S40M25N5M3S'
	#cigar_str = '2S40M25N8M'
	cigar_str = '40M25N8M2S'
	read_seq = "TAGAGAGAGAGDHEGDGDFDFDFDEHHHHSSSSGEEBCSDGEYQLGQKNV"
	ref_chr = "chr2"
	ref_posn = 20001
	a = sam_cigar(cigar_str, read_seq, ref_chr, ref_posn)

	print(a.analyze_cigar()[0][1])

	print(a.get_softclipped_posn())
	print(a.get_softclipped_seq())
	print(a.get_mapping_posn())
	print(a.get_softclipped_genome_site())

