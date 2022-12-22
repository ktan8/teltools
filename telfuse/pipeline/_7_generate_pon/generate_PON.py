#!/usr/bin/python


########
# Aggregate and generate PON
########

import sys
from sites import *
from aggregated_sites import *

import re



def filter_sequence(seq, pattern):
	'''
	Check if the sequence matches
	the required sequence
	'''
	result = re.search(pattern, seq)

	return result


#
#class SiteCount:
#	def __init__(self, chrom, posn, orientation):
#		#self.key =
#		self.chrom = chrom
#		self.posn = posn
#		self.orientation = orientation
#		self.samples = 0
#		self.readcount = 0	
#
#	def add_sample(self, read_count):
#		'''
#		Add information for a sample
#		to a site_count object
#		'''
#		self.samples += 1
#		self.readcount += read_count
#
#
#	def generate_unique_key(self):
#		'''
#		Generate the key to unique identify
#		the site of interest.
#		'''
#		data = [self.chrom, self.posn, self.orientation]
#		unique_key = "_".join(map(str, data))
#
#		return unique_key
#		
#	def generate_output(self):
#		'''
#		Generate output describing the site
#		as a list
#		'''
#		output = [self.chrom. self.posn, self.orientation, self.samples, self.readcount]
#		
#		return output
#
#
#class SiteCountCollection:
#	def __init__(self):
#		self.sites = dict()
#		
#
#	def add_site(self, chrom, posn, orientation, reads):
#		data = [chrom, posn, orientation]
#		unique_key = "_".join(map(str, data))
#		
#		# Check if site exists already
#		if unique_key in self.sites:
#			self.sites[unique_key].add_sample(reads)
#		else:
#			curr_site = SiteCount(chrom, posn, orientation)
#			self.sites[unique_key] = curr_site
#			self.sites[unique_key].add_sample(reads)
#		
#	def generate_output_line(self):
#		'''
#		Generate output for each
#		site line by line
#		'''
#		for site in self.sites.values():
#			output = site.generate_output()
#			output_line = "\t".join(map(str, output))
#			
#			yield output_line
#			



if __name__ == "__main__":
	aggregated_file = sys.argv[1]
	agg_site_file = AggregateSiteFile(aggregated_file)
	
	#site_collection = SiteCountCollection()

	
	#for site in agg_site_file:
	#	site_collection.add_site(site.chrom, site.posn, site.orientation, site.reads)
		
	
	#for agg_output in site_collection.generate_output_line():
	#	print agg_output

	for counts in agg_site_file.count_sites().generate_output_line():
		print(counts)

