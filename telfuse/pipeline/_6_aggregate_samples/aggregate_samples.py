#!/usr/bin/env python

################################
# Script used to aggregate the 
# different samples together
################################


import sys
import os
from sites import *

folder = sys.argv[1]



filelist = os.listdir(folder)
#print filelist


for sitefile in filelist:
	sitefile_full = folder + "/" + sitefile
	sitefile_list = sitefile.split(".")
	
	# Get sample label
	label = sitefile_list[0]

	sitefile = SiteFile(sitefile_full)

	for line_list in sitefile.line_generator(outputList = True):
		#print label
		output_list = [label] + line_list

		print("\t".join(output_list))

