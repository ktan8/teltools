#!/usr/bin/python


#rom lib.sites import *
from site_count_collection import *


class AggregateSiteLine:
	'''
	'''
	def __init__(self, line):
		try:
			self.line 	= line
			self.lineArr 	= line.strip().split("\t")
			self.sample 	= self.lineArr[0]
			self.chrom 	= self.lineArr[1]
			self.position 	= self.lineArr[2]
			self.orientation = self.lineArr[3]
			self.reads 		= int(self.lineArr[4])
			self.identity_ave 	= self.lineArr[5]
			self.identity_weighted 	= self.lineArr[6]
			self.sequence 		= self.lineArr[7]
		except:
			print(self.lineArr)


class AggregateSiteFile:
	def __init__(self, aggregated_file):
		self.aggregated_file = aggregated_file


	def line_generator(self):
		'''
		Generator to loop through line by line of
		the aggreated file. Returns an object
		representing each 
		'''
		agg_file = open(self.aggregated_file, "r")
		
		for line in agg_file:
			agg_site_line = AggregateSiteLine(line)
	
			yield agg_site_line

	def count_sites(self):
		'''
		Count number of samples and reads
		supporting a particular site
		'''
		site_count_collection = SiteCountCollection()
	
		for site in self.line_generator():
			#unique_key = chrom + "_" + posn + "_" + orientation
			#site_count[unique_key] += 1
			site_count_collection.add_site(site.chrom, site.position, site.orientation, site.reads)
	
		return site_count_collection


