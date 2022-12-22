import sys
import re
from samReadClass import *



class samReadClassFilter(samRead):
	'''
	Inherit the samRead class. Additional
	functions to help filter relevant sam
	reads
	'''
	def __init__(self, samLine):
		samRead.__init__(self, samLine)

		self.cigar_char_list = re.findall(r'[0-9]+([MIDNSHPX=])', self.Cigar)
		self.cigar_value_list = re.findall(r'([0-9]+)[MIDNSHPX=]', self.Cigar)
		self.readmapped = self.check_read_mapped()


	def check_read_mapped(self):
		'''
		Check the flag to tell if the read is mapped
		'''
		# Read unmapped
		if int(self.Flag) & 4:
			return 0
		# Read mapped
		else:
			return 1


	def check_perfect_mapping(self):
		'''
		Check if the read is perfectly mapped
		onto the location. 
		'''
		allowedMappedChar = ['M', 'D', 'I']

		if check_read_mapped():
			if len(self.cigar_char_list) == 1 and self.cigar_char_list[0] == 'M':
				return 1
			elif len(self.cigar_char_list) > 1:
				'''
				Check that there is only mapped, deleted or inserted bases
				in the alignment, which indicates perfect mapping
				'''
				for cigar_char in self.cigar_char_list:
					if cigar_char not in allowedMappedChar:
						return 0
				return 1

		return 0



	def has_softclip(self):
		'''
		Check if part of the read contains softclipped
		regions.
		'''
		has_softlip = 0
		for cigarChar in self.cigar_char_list:
			if cigarChar == "S":
				return 1

		return 0


	def get_softclip_region(self):
		'''
		Identify and return where the soft clip
		region is.
		- start position is 0-based inclusive
		- end position is 0-based exclusive
		'''
		currStart = 0
		cigar_with_read_span = ["M", "S", "I"]
		softclipped_regions = []
		for i in range(0,len(self.cigar_char_list)):
			cigar_value = int(self.cigar_value_list[i])
			cigar_char = self.cigar_char_list[i]

			if cigar_char == "S":
				softclipped_start = currStart
				softclipped_end = softclipped_start + cigar_value

				region_string = str(softclipped_start) + "-" + str(softclipped_end)
				softclipped_regions.append(region_string)

			if cigar_char in cigar_with_read_span:
				currStart += cigar_value

		return softclipped_regions



	def check_softclip_lowqual(self, qualCutoff = 20):
		'''
		Check if soft clipped region is of low 
		quality. Note that we can have two possible
		regions corresponding to left and right ends
		of the read

		- What if the read is a secondary alignment and has 
		  no qual values?
		'''
		softclippedRegions = self.get_softclip_region()
		qualCutoff = float(qualCutoff)
		#print softclippedRegions

		lowQualResults = []
		for region in softclippedRegions:
			[start, end] = region.split("-")
			start 	= int(start)
			end 	= int(end)
			region_qual = self.QualVals[start:end]

			total_qualVals = 0
			for qual in region_qual:
				total_qualVals += ord(qual)
			#print len(region_qual)
			ave_qualVal = (float(total_qualVals) / len(region_qual)) - 33
			#print ave_qualVal

			if ave_qualVal < qualCutoff:
				lowQualResults.append(1)
			else:
				lowQualResults.append(0)

		#print lowQualResults
		for result in lowQualResults:
			# If at least one of the regions is not
			# of low quality.
			if result == 0:
				return 0

		return 1

	def is_candidate_read(self):
		'''
		Identify reads that are candidates for a 
		fusion event.
		'''

		if self.readmapped == 0:
			'''
			Unmapped reads
			'''
			return 1
		else:
			'''
			Mapped reads
			'''
			if self.has_softclip():
				if not self.check_softclip_lowqual():
					return 1

		return 0











