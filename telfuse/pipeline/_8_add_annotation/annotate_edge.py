#!/usr/bin/env python

import sys

input_file = sys.argv[1]
fai_file = sys.argv[2]


class Fai():
	'''
	Create a fai file class
	'''
	def __init__(self, fai_file):
		'''
		
		'''
		self.fai_file = fai_file
		
		self.fai_dict = self.parse_fai_file()


	def parse_fai_file(self):
		'''
		Parses the fai file
		'''
		fai_dict = dict()
	
		f = open(self.fai_file, "r")
		for line in f:
			line_arr = line.strip().split("\t")
			chrom = line_arr[0]
			position = line_arr[1]
			fai_dict[chrom] = int(position)
			
		return(fai_dict)

class InputListLine():
	'''
	Class to reprsent what the input
	list line looks like
	'''
	def __init__(self, line):
		self.line = line
		self.lineArr = line.strip().split("\t")
		self.sample = self.lineArr[0]
		self.chrom = self.lineArr[1]
		self.position = int(self.lineArr[2])
		self.dist_left = None
		self.dist_right = None
		self.near_edge = None	
	
	def update_dist_variables(self, dist_left, dist_right):
		'''
		Update the variables for the 
		distance from left and right
		'''
		self.dist_left = dist_left
		self.dist_right = dist_right

	def update_nearedge(self, near_edge):
		'''
		Update whether the site is near
		the edge or not.
		'''
		self.near_edge = near_edge

	def get_output(self):
		'''
		Generate and return the output as a list
		'''
		output = self.lineArr + [self.dist_left, self.dist_right, self.near_edge]
		
		return output
		


class InputList():
	'''
	Input list of sites that can be extracted
	and analyzed
	'''
	def __init__(self, sitefile):
		self.sitefile = sitefile

	def calc_edge_distance(self, fai_dict):
		'''
		Calculate how far the 
		'''
		for line_class in self.site_file_generator():
			dist_left_edge = line_class.position
			dist_right_edge = fai_dict[line_class.chrom] - line_class.position

			line_class.update_dist_variables(dist_left_edge, dist_right_edge)
			
			yield line_class
			

	def site_file_generator(self):
		'''
		Generator for the site file
		'''
		f = open(self.sitefile, "r")
		
		for line in f:
			line_class = InputListLine(line)
			yield line_class


	def check_near_edge(self, edge_cutoff, fai_dict):
		'''
		Check if the site is near the edge
		or not
		'''
		for line_class in self.calc_edge_distance(fai_dict):
			near_edge = None
			if line_class.dist_left < edge_cutoff or \
			line_class.dist_right < edge_cutoff:
				near_edge = 1
			else:
				near_edge = 0

			line_class.update_nearedge(near_edge)
			
			yield line_class
	
	def check_near_edge_output(self, edge_cutoff, fai_dict):
		for line_class in self.check_near_edge(edge_cutoff, fai_dict):
			yield line_class.get_output()	


def annotate_filter_edges(fai_dict, site_list):
	'''
	Annotate sites that are at the edge
	to allow for them to be filtered out.
	'''
	edgecutoff = 500000

	for result in site_list.check_near_edge_output(edgecutoff, fai_dict):
		print("\t".join(map(str, result)))
	
	
	
if __name__ == "__main__":
	input_file = sys.argv[1]
	fai_file = sys.argv[2]
	#print Fai(fai_file).fai_dict
	fai_dict = Fai(fai_file).parse_fai_file()
	sites = InputList(input_file)
	
	annotate_filter_edges(fai_dict, sites)

