#!/usr/bin/python


###########################
# Object to store the counts of an
# object.
###########################

from sites import *


class SiteCount:
	def __init__(self, chrom, posn, orientation):
		#self.key =
		self.chrom = chrom
		self.posn = posn
		self.orientation = orientation
		self.samples = 0
		self.readcount = 0	

	def add_sample(self, read_count):
		'''
		Add information for a sample
		to a site_count object
		'''
		self.samples += 1
		self.readcount += read_count


	def generate_unique_key(self):
		'''
		Generate the key to unique identify
		the site of interest.
		'''
		data = [self.chrom, self.posn, self.orientation]
		unique_key = "_".join(map(str, data))

		return unique_key
		
	def generate_output(self):
		'''
		Generate output describing the site
		as a list
		'''
		output = [self.chrom, self.posn, self.orientation, self.samples, self.readcount]
		
		return output


class SiteCountCollection:
	def __init__(self):
		self.sites = dict()
		

	def add_site(self, chrom, posn, orientation, reads):
		data = [chrom, posn, orientation]
		unique_key = "_".join(map(str, data))
		
		# Check if site exists already
		if unique_key in self.sites:
			self.sites[unique_key].add_sample(reads)
		else:
			curr_site = SiteCount(chrom, posn, orientation)
			self.sites[unique_key] = curr_site
			self.sites[unique_key].add_sample(reads)
		
	def generate_output_line(self):
		'''
		Generate output for each
		site line by line
		'''
		for site in self.sites.values():
			output = site.generate_output()
			output_line = "\t".join(map(str, output))
			
			yield output_line
			
