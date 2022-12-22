#!/usr/bin/env python

import sys




class Pon:
	def __init__(self, PON_file):
		self.PON_file = PON_file		


	def load_PON_dict(self):
		'''
		Load each entry in the PON
		into a dict
		'''
		f = open(self.PON_file, "r")
		PON_dict = dict()
		for line in f:
			lineArr = line.strip().split("\t")
			chrom  = lineArr[0]
			posn = lineArr[1]
			orientation = lineArr[2]
			#left_count = lineArr[3]
			#right_count = lineArr[4]
			samples = lineArr[3]
			reads = lineArr[4]
			#near_edge = lineArr[5]
			key = "_".join([chrom, posn, orientation])
			
			#result = [left_count, right_count, near_edge]
			result = [samples, reads]

			PON_dict[key] = result

		return PON_dict







# Add PON annotation
#def add_PON_anno




#def add_edge_anno():
	



input_file = sys.argv[1]
PON_file = sys.argv[2]


PON_dict = Pon(PON_file).load_PON_dict()

f = open(input_file, "r")
for line in f:
	result = ["NA", "NA"]
	lineArr = line.strip().split("\t")
	key = "_".join(lineArr[1:4])
	#print key
	if key in PON_dict:
		result = PON_dict[key]

	finalLine = lineArr + result
	print("\t".join(finalLine))





