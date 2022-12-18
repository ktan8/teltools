#!/usr/bin/python

import numpy as np
from lib.motif_lib import *
from scipy.signal import medfilt




def moving_average(a, n=3) :
	'''
	https://stackoverflow.com/questions/14313510/how-to-calculate-moving-average-using-numpy
	'''
	ret = np.cumsum(a, dtype=float)
	ret[n:] = ret[n:] - ret[:-n]

	return ret[n - 1:] / n


def denoiseData(motifCounts, movingAveWindow=50, medianKernelSize=501):
	'''
	Smoothern the data by getting the moving
	average and then the median filter
	'''
	motifCountsAve = moving_average(motifCounts, n=movingAveWindow)
	motifCountsAveMedian = medfilt(motifCountsAve, kernel_size=medianKernelSize)

	return(motifCountsAve, motifCountsAveMedian)


def calc_telomere_length(countsAcrossRead, threshold=0.05):
	'''
	Identify region above threshold which correspond
	to telomeres.
	'''

	regionAboveThreshold = countsAcrossRead > threshold


	falseToTrue = []
	trueToFalse = []
	for i in range(len(regionAboveThreshold) - 1):
		if regionAboveThreshold[i] == True and regionAboveThreshold[i+1] == False:
			trueToFalse.append(i)
		if regionAboveThreshold[i] == False and regionAboveThreshold[i+1] == True:
			falseToTrue.append(i)

	# Deal with case where the whole sequence is telomeric
	if regionAboveThreshold[0] == True:
		falseToTrue.append(0)
	if regionAboveThreshold[-1] == True:
			trueToFalse.append(len(regionAboveThreshold))


	if len(falseToTrue) > 0 or len(trueToFalse) > 0:
		# If there's no false to true
		if not falseToTrue:
			falseToTrue.append(1)
		# If there's no true to false
		if not trueToFalse:
			trueToFalse.append(len(regionAboveThreshold))
	else:
		falseToTrue.append(-1)
		trueToFalse.append(-1)


	# Final results
	telomereStart 	= falseToTrue[0]
	telomereEnd 	= trueToFalse[0]
	telomereLength 	= trueToFalse[0] - falseToTrue[0]
	telomereRegions = len(falseToTrue)


	averageTelomereSignal = np.average(countsAcrossRead[telomereStart:telomereEnd])

	return(telomereStart, telomereEnd, telomereLength, telomereRegions, averageTelomereSignal)




class TelomereSignalDetails():
	def __init__(self, sequence, motif):
		self.sequence 	= sequence
		self.motif 		= motif
		self.strand		= "NA"

		self.MotifCounts = searchForMotif(self.sequence, self.motif)
		self.MotifCountsAve, self.MotifCountsAveMedian \
		= denoiseData(self.MotifCounts)

		(self.TelomereStart, self.TelomereEnd, self.TelomereLength, 
		self.TelomereRegions, self.AverageTelomereSignal) \
		= calc_telomere_length(self.MotifCountsAveMedian)



class AnalysisInfo():
	def __init__(self, sequence, motif, sequence_name):
		self.sequence 		= sequence
		self.motif 			= motif
		self.revComMotif	= reverseComplement(motif)

		self.sequence_name 	= sequence_name
		self.sequenceLength = len(self.sequence)


		# Count for forward strand
		self.fwdSeqDetails = TelomereSignalDetails(self.sequence, self.motif)

		# Count for reverse strand
		self.revSeqDetails = TelomereSignalDetails(self.sequence, self.revComMotif)





	# self.fwdMotifCounts = searchForMotif(self.sequence, self.motif)
	# self.fwdMotifCountsAve, self.fwdMotifCountsAveMedian = denoiseData(self.fwdMotifCounts)
	# (self.fwdTelomereStart, self.fwdTelomereEnd, self.fwdTelomereLength, 
	# 	self.fwdTelomereRegions, self.fwdAverageTelomereSignal) 
	# = calc_telomere_length(fwdMotifCountsAveMedian)


	# # Count for reverse strand
	# self.revMotifCounts = searchForMotif(self.sequence, self.revComMotif)
	# self.revMotifCountsAve, self.revMotifCountsAveMedian 
	# = denoiseData(self.revMotifCounts)
	# (self.revTelomereStart, self.revTelomereEnd, self.revTelomereLength, 
	# 	self.revTelomereRegions, self.revAverageTelomereSignal) 
	# = calc_telomere_length(revMotifCountsAveMedian)



		self.finalResult = []

		if self.fwdSeqDetails.TelomereLength > self.revSeqDetails.TelomereLength:
			self.finalResult = self.fwdSeqDetails
			self.finalResult.strand = "forward"
		else:
			self.finalResult = self.revSeqDetails
			self.finalResult.strand = "reverse"


		# self.finalResult =  [self.revTelomereStart, self.revTelomereEnd, self.revTelomereLength, 
		# self.revTelomereRegions, self.revAverageTelomereSignal, self.stringLength, self."reverse"]
		# [self.fwdTelomereStart, self.fwdTelomereEnd, self.fwdTelomereLength, 
		# self.fwdTelomereRegions, self.fwdAverageTelomereSignal, self.stringLength, "forward",]



	def generate_output(self):
		result = [self.sequence_name, self.finalResult.TelomereStart, self.finalResult.TelomereEnd, 
				  self.finalResult.TelomereLength, self.finalResult.TelomereRegions, 
				  self.finalResult.AverageTelomereSignal, self.sequenceLength, 
				  self.finalResult.strand]

		return(result)




