#!/usr/bin/python


import os
import numpy as np
import matplotlib.pyplot as plt
import re
from scipy.signal import medfilt



def reverseComplement(seq):
	'''
	https://stackoverflow.com/questions/25188968/reverse-complement-of-dna-strand-using-python
	'''
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
	reverse_complement = "".join(complement.get(base, base) for base in reversed(seq))

	return reverse_complement


def generateMotifCombinations(sequence):
	'''
	Generate all the possible sequence combinations
	for a motif of interest after shifting the frame.
	e.g. TTAGGG, TAGGGT, AGGGTT, GGGTTA, GGTTAG, GTTAGG
	'''
	seqLength = len(sequence)
	seqCombinations = []
	for i in range(seqLength):
		newString = sequence[i:seqLength] + sequence[0:i]

		seqCombinations.append(newString)

	return seqCombinations



def searchForMotif(sequence, motif):
	'''
	Search for a motif in a sequence of interest
	'''
	motifCombinations = generateMotifCombinations(motif)
	motifCounts = np.zeros(len(sequence))

	stringLength = len(sequence)
	motifLength = len(motif)

	# Search through the full string
	for i in range(stringLength - motifLength):
		subString = sequence[i:i+motifLength]
		# Loop through all the possible motifs
		for motifCurrent in motifCombinations:
			if subString == motifCurrent:
				motifCounts[i] += 1
			else:
				pass

	return(motifCounts)



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


	# countsAcrossRead.append(0)
	# countsAcrossRead.append(0, 0)

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
		#falseToTrue.append(0)
		falseToTrue.insert(0,1)
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


	# Calculate average telomere signal in identified region
	#print("=====")
	#print(len(countsAcrossRead[telomereStart:telomereEnd]))
	#print(countsAcrossRead[telomereStart:telomereEnd])
	#print(countsAcrossRead)
	#print(telomereStart, telomereEnd)
	

	if countsAcrossRead[telomereStart:telomereEnd].size != 0: 
		averageTelomereSignal = np.average(countsAcrossRead[telomereStart:telomereEnd])
	else:
		averageTelomereSignal = "NA" # because start and end for arrays are -1 and will give an error
	#averageTelomereSignal = 10

	return(telomereStart, telomereEnd, telomereLength, telomereRegions, averageTelomereSignal)




def analyze_sequence_for_motif(string, motif, sequence_name, pic_format="png", 
	threshold=0.05, movingAveWindow=50, medianKernelSize=501, folder="./", plot_fig=True):
	'''
	Analyze an input sequence for a motif of interest.
	For instance, one can input the 'TTAGGG' sequence to look
	for telomeres
	'''
	if not os.path.exists(folder):	
		os.mkdir(folder)
	sequence_name_clean = re.sub("/", "_", sequence_name) # in case slash is in the readname
	output_image = folder + sequence_name_clean + "." + pic_format
	#print(output_image)	
	stringLength = len(string)
	motifRevComplement = reverseComplement(motif)


	# Count for forward strand
	fwdMotifCounts = searchForMotif(string, motif)
	fwdMotifCountsAve, fwdMotifCountsAveMedian = denoiseData(fwdMotifCounts, 
		movingAveWindow=movingAveWindow, medianKernelSize=medianKernelSize)
	(fwdTelomereStart, fwdTelomereEnd, fwdTelomereLength, 
		fwdTelomereRegions, fwdAverageTelomereSignal) = calc_telomere_length(fwdMotifCountsAveMedian, threshold=threshold)


	# Count for reverse strand
	revMotifCounts = searchForMotif(string, motifRevComplement)
	revMotifCountsAve, revMotifCountsAveMedian = denoiseData(revMotifCounts,
		movingAveWindow=movingAveWindow, medianKernelSize=medianKernelSize)
	(revTelomereStart, revTelomereEnd, revTelomereLength, 
		revTelomereRegions, revAverageTelomereSignal) = calc_telomere_length(revMotifCountsAveMedian, threshold=threshold)


	# print [fwdTelomereStart, fwdTelomereEnd, fwdTelomereLength, fwdTelomereRegions]
	# print [revTelomereStart, revTelomereEnd, revTelomereLength, revTelomereRegions]
	finalResult = []

	if fwdTelomereLength > revTelomereLength:
		finalResult =  [fwdTelomereStart, fwdTelomereEnd, fwdTelomereLength, 
		fwdTelomereRegions, fwdAverageTelomereSignal, stringLength, "forward",]
	elif fwdTelomereLength < revTelomereLength:
		finalResult =  [revTelomereStart, revTelomereEnd, revTelomereLength, 
		revTelomereRegions, revAverageTelomereSignal, stringLength, "reverse"]
	# Where both equal and zero.
	elif fwdTelomereLength == 0 & revTelomereLength == 0:
		finalResult =  [revTelomereStart, revTelomereEnd, revTelomereLength,
		revTelomereRegions, revAverageTelomereSignal, stringLength, "NA"]
	# Where both equal. Happens especially when both are non-telomeric
	else:
		finalResult =  [revTelomereStart, revTelomereEnd, revTelomereLength,
		revTelomereRegions, revAverageTelomereSignal, stringLength, "equal"]
	

	if plot_fig == True:	
    	    # For aesthetics. If it is forward, we plot forward
            # on top of reverse. If it is reverse, we plot reverse
	    # on top of forward.
	    plt.figure(figsize=(9,4.8))	

	    if finalResult[6] == "forward":
	            # Plot reverse motif call
		    plt.plot(revMotifCountsAve, color='mistyrose', ls='-')
		    plt.plot(revMotifCountsAveMedian, color='maroon', ls='-')

		    # Plot forward motif call
		    plt.plot(fwdMotifCountsAve, color='lightsteelblue', ls='-')
		    plt.plot(fwdMotifCountsAveMedian, color='royalblue', ls='-')

	    else:
		    # Plot forward motif call
		    plt.plot(fwdMotifCountsAve, color='lightsteelblue', ls='-')
		    plt.plot(fwdMotifCountsAveMedian, color='royalblue', ls='-')

		    # Plot reverse motif call
		    plt.plot(revMotifCountsAve, color='mistyrose', ls='-')
		    plt.plot(revMotifCountsAveMedian, color='maroon', ls='-')


	#plt.figure(figsize=(80,50))



	    # Plot vertical line
	    plt.axvline(x=finalResult[0], color='grey', ls='--')
	    plt.axvline(x=finalResult[1], color='grey', ls='--')


	    plt.xlabel("Position (bp)")
	    plt.ylabel("Telomeric repeat signal")
	    plt.ylim(-0.05,1.05)

	    #plt.show()

	    # build a rectangle in axes coords
	    left, width = 0.01*len(string) , len(string)
	    bottom, height = 0, 1
	    right = left + 0.55 * width
	    top = bottom + 0.98 * height


	    telomereLengthString = 'Telomere length (nt) = ' + str(finalResult[2])
	
	    # Decide where to place the text depending on whether
	    # the telomeric region is in the forward or reverse direction
	    if finalResult[6] == "forward":
		    plt.text(left, top, telomereLengthString,
	            horizontalalignment='left',
	            verticalalignment='bottom')
	    else:
		    plt.text(right, top, telomereLengthString,
			horizontalalignment='left',
			verticalalignment='bottom')

	    plt.savefig(output_image, format=pic_format)
	    plt.close()
	
	finalResult.append(string)	
	
	return finalResult





if __name__ == "__main__":
	string 	= "AGGGATTATTTGGGTGAAGAGTGGATCGTATTATTGATCACCCACCATGTACCTAAACCCACAGTGAGGATCTCTCTGTTCCCACAAGCCTTAAGTGGGGCATCCCAGTCGGGGGTAAGGCAGTAGCCATTGGAACAACAGAACAAAGGAAAAAGGTGGAGACTTGGAGCCAGACATCTCACGCGCAGGAAGTGAAGAAGTCCCAATAAAATTCCTGACAGGGACTCTTAGGCCTGTTTTAATGCACCGCTCAGCCACTCAATCCCATTTTTTCTACAAAAGCTATTTCACAACTTGGGTTGCTTTTGCAAATGAGTTATATGCCATTTGGTAATGCCTATTGGTGAAAACTTTTACTTCTTCAAAGTTAAACCAAGAAACTGGGACAATCGCCCTCCCTCGCTAATAGCTCGTTCAAGTTGATTTACAGAACTGATGGGGCTAATAAACGCGCTCTCTCTGGACTTTAGGGGCGCAGTGAGGCCGTAACACCAGGACCCAAGGGCCCTGCCTAGTCCCAATCTGCTCCGCAGGTGGCGTGCAGCCACGCGACACTGACAGCAATAAGGCCGGCAGTGTCATCATCGATGCAGGACAGGCGGCGTTACGGGCACCACACCATAATGCAAGATGACCAGCAGTGCCATGTCGTCGCTGCCACACACGGGAGCAAGAGGATCCTGGAGTGCTCCCAGAGGATAGAGCGGCGTGCCCGACTATACTGCAGGCAAGAGGGTCCTGGGCATTGTCCCAGCTGCAGGTAGGCGGGCGTGATGCCACTACACGTGAGCACAAGAGGGTCGGAAAGTGTCACAAGCTGACAGCAGGCGGGCGTGCTGCCACTACACTGTGTAGCAAGAGGGCCCCTGAAACCGTCCCGTAGCTGCAGCAGGCCGGCGGCCGCCACTATTTAGCGAGCAAGAGATCCTGAAGCGACCCGGAGACCAGCAGGGGACGCGTGCCACCAGGTAAGAAAGAGGGACTGCGGTGTCCTAGCCGCTAGCAGGGGGCGCAATGTGGAAAGCACCGCGGGTAGCGGGTCCTGTAGCTTGCACGGCTTGCAAGCAGGGGGCCCAGAGGACGGCTTTTCGGATTACGAGGTTTCAACCCGTCTCGCTGCCGCGCCCCGGGGAAGTGAGTCTCTGCCTGACAAAACGCTCCAACCCCCGCGCCTGTCCCGGTTGGCGTGTTGGACTGTGGATGGCGCGTTCGCGCGCTGCCCCCAGCCACCCAGGGGAACAGCGCAAAGGTGGATGCTAGACGAGCACTCCCCCGCCCCCACAGGGGAACGAGATCTCTGAGCCTGCGCCGGCGAGCCGAAGCCTCTCTGCGCCTGACAGAGGCGGCGAGGCCGAGCCTCTCTCTGCGCCTGCGCGGCGCGCCGGCCTATCTGCGCCTGCGAGGCGCAGCGGACGCCCTAGGTGCCGCCTGGCGGAGGCCGCCCGCGCCTCTGCGCCCGCGCCGGCGCGACGCCTTGGCCGCCTGCGCACGGCAGTCGGATTTATCAGCATCCTATTGTAGTGGTGTCACCGGCGGCGCGGCCGGCGAGCTCTCCTTTGTCTGTTCGTCTATGCCCACTTTCTACCTTGCCTGTGTCTCCCTTTCCATTGCTCCCGCTCTCTTCTCTCCCCCTTACTTGACAGCGCCATTCTGGCGTCGAGAGATTGCGCATCGGCGTCTGCTCTTGCTTGCGGCACAGGAATGCGAAACGGACGAGAGCACACGTGGGCAGCCTGCGCCGTCACTCGCGCCCTGCGGCGTGGCGCGCCGCGCAGGCTCTATGTCGCGTACGAGCCAGGCCGCGCCGAGCCGGTCTGCGCCCTGACGCCGGAGCGCCGACGCCTCTTGCGCCCTGCGCCGGGCCGCGCCGCGCCTGCTCTGCGCCCTGCGCCGCGCGCGCGGCCTCTCTGCGCCTGCGCCGGCGCGCCGCGCCTCTCATGCGCCGCGCCGCGCGGCCGCGGCGCGCGCTTCGCTGTTTACGCCTAGCGCGGGTGTGCCAGGCCGGCCGACGCCTGCTCGAGCCTGCGCCGTGGGCACGCCAAGCGCCTTCTCTGCGCCTGCGCGGGGCGGGACTGAGACCACTCTCATGGACTGCGCCGGCGCGCCGTGAGTTGCCTCTTCTGCGCCCCCTGCGCCGGGCGCGCCGAGCCTTATGCGCCATGCTGGCGGCGCCGCACGCGACTATCTGCGCACAGAGACGTTGCTGCGCCGACCTGGTGCGAGGGATCTTTGTGGAGTTGCAGTTCTCCTTAAGCACAGTAACCCGGAGAGCATCGCGAGGGGGAGCTGCGTTTGCTGGCACAGCCGGAGGATGCTGTGAAAGGGGATCAGCAGCGTGCTCCTCAGCACAGACCCCGCCGGGCGTGGTTATCCGGCCGGCACCGCGAAGGCGGGAATGCGTTCTGCCTCAGCAGCCACCCGGTGGTAAACAAAAGGTGGGAAGAAGCTTCTCGCAGCATCGACCCTGGTGTTGGGTTAGGGTAGGGGTTGGCGCTTGTAGGGTTCGGGTTAGGGGTTTAGGGGGTAAGGGTTCGGGTTAGGGTTAGGGTTAGGGTTCGGGTTAGGGTTCGGGTTAGGGTAGGGTTAGGGGTGTCGGGTTAGGGTTCGGGTAGGTGTTTGTGAGGGGTTAGGGGGTGCATAGGGGAGTAAGGTAGTTAGGGGGTTTAGGTTAGGGTTAGGGTTAGAGGTTGGGTTGTTAGGTGTTTAGGAGCGGTTACAGGGTTAGGGTTAAGCGGGTGTAGGGTAGTAGGGTAGGGTGAGGGTTAAGGGTTAGTGGTAGGGTTGGGGTTGGGGTGGGGGGTGAGGCAGTGAGGGTAGAGTGGTGAGGGTAGAGGGTGAGGGTGAGGAGTGAGGGTTAGGGTTAAGCGCTCCACAGGGTTTAGGGGTGCAAGAGTGACGGGTGAGGGTGAGGGATAGGGTGAGGGGTATAGTGGTGTAGGGTTTTTGTAATGGTCAGAGGGTCTAGGGGGTAGGGTTGCTAGGGTGAGTGGTAGGGGTAGGGCGTTATAGTTAGGGGTTTAGGGTTGGTTCGGGATTAGGGTTGGGTTAGGGTTCAGGGGTTAGTGGTTAGGGTTGGGTTAGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTGCTAGGGTTTTAGGGTCGTTAGGGTTAGGGTCTGAGGGGTTAGGGGCTGCAGGGGTGAGGGGATGAGGGTGGGGTTAGGGTGTAGGGTTAGCGGTGTTGGGTTAGGTGTTAGGGTTAGGGTGATAAGGGTTAGGGTTAGGGTATCTAGGGTTAGGGTTAGGGTCTTCGGGGTAGGGTTAAGGGGTTGGGTTAGGGTTAGGGAGGGTTAGGGTTTTTAGGGCTTAGGGTAGGGTTGGGGGTTGGGGTTGGCGGCTGTTGAGGGGCCACGGTGGGTGGGGGGTTTGTAGGCGGTGTAGGGTCGGGTAGGGTTGGGCGTTAGGGTTAGGGTAGGGTCTTAGGGTTAGGGTTGAGGGTTAGGGTGATGGTGGGGGTCTGCGGTGTACCTGACGAGGGTGGAGGTGGGGAGGGTGAAGGACCGTGAGGGTGGAAGGGTGAGGGAAGACGAGGTGCAGGGTAGGGTTAAGGGGGGTTAGGGTTCTGGGTTGGTGGGTGAGGTTGAGGGTGAGGGTGGGGTGCGGTGTGGGGTGAGGGTTGACGGGTTGAGGGTGAGGGATAGGGTTAGGGGTTAGGGTATGGGTTAGGGTATAGGGTTAGGGGTAGGGGGTAGGGTTAGGTAGGGTTAGGCGTTAGGGTGAGGGTTAGGGTAGGGTTAGAGTATAGGGTTTGGGTGTAGGGTTAGGGTTAGGGTAGGGTAGGGTAGGGTTTAGGGTAGGGTTAGGGTTGGGTTAGGGTTAGTTAGGGTAGGGTTAGGGTTAGGGTTAGGGTTTAGGGGTAGGGTAGGTTTTAGGGTTAGGGTGTATAGGTAGGGTTTAGGGTTAGGGTTAGGGTTAGGGTTAGGTGGGTATAGGGTTGTGGTTCGGTTAGGGTTCAGGGTAGGGTTAGGGTTTAGGGTTAGGGTAGGGTTAGGGTCTAGGGGTAGGGTGTCGGGTAGGGTTGGGTTAGGGTTAGGGTTAGGGATATTAGGTTAGGTTACGGGTTCAGGGTTGGGGTTAGGGTCTTAGGGTTTGGGTTAGGGGTAGGGTTAGGGTTAGGGTGGTTAGGGTTAGGGTTAGGGTTAGGGTTCGGGTTAGAAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGGTTAGGGTTAGGGGGGATGGTTAGGGTTAAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGGTAGGGCTTAGGTTAGCGGTGGATCAAGGGTTAGGTTTAGGTGTAAGGGTGTTACGCGGGGTTAGAGGGTTAGGGGTTAGGTGGGGTGTAGGTTGAGGTATAGGGTTAGGTGTTTAGGGGTTGTAGGGTTAGGTTAGCGGTTGTAGGGTTGTGGAGGTTGCGGTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTGGGTTAGGGGTTAGGGTTAGGGTTAGGGGCGAGGTTACTAGGCGTTTAGGTGTAGGGTTACTTGAGGGTTAGGGTTAGGGTGGGTGAGGGGGCTGTAGGGCTAGGGTATAGGGTTCGGGTAGGTTAGGGTTAGGGTTAGGGGGGTAGGGTATAGGGTAGGGGTTAGGCGTTCGGGTTAGGGTAGGAGTTAGGGTTAGGGTTAGGGTTTCGGGTTGGGTTAGGTTAGGTAGGGTTAGGGTGTAGGGTTAGGGTTAGGGTTAGGGTGTAGGGGTTGGGGAGTTAGGGTACGGGTTCGGGTTAGGAGTTAGGGTTAGGGTATAGGGTTCGGGTTCGGGTTAGGGTGTGTAGGGTTTGGTTTAGAGGTTAGGGTTAGGGTGAGGGTAGGGTTAGGGTTGGGATTAGGGGTTGGGTAATCGGGTTAGGGTTAGGTTAGGGAGCGGTTAGGTTAGGGTTAGGGGTTCGGGTTGGGTTAGGGTTAGGGTTGGGGGGTTGTAGGGTTAGGGTTAGGGTTAGGTTAGTGTTAGGGTAGGGTTCGGGTTAGGATAGGGGTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTATGGGTTGGGTTAGGGTTAGGGTTATGAGGCGTTAGGGTTAGGGTGTTAGGGTTAGGGTAGGGTTAGGGTTTAGGGTTAGGGTTAGGGATTGAGGGGTTAGGGTTTGGTTAGGGTTAGGGTGGTATTTAGGGTTAGGGTTAGGGTTCAGGTTCGGGGTTAGGGTTTTAGGGTTAGGGTAGGGTTAGGGTTAGGGGCGTAGGGTTAGGGGTAGGGTTAGGGTTAGGTTTAGGAGTTAGGTCTAAGGGGTTAGGGTGATGGGTTAGGGTTGGGTTAGGGTTAGGGTTAGGGTAGGCTTAGGTTAGGGTTGAGGGGTATGGTAGGGTTCAGGGTAGGGTTAAGCGGTTAGGTGTTAGGGTTGGTGGTGAGCGGGGTAGGGTAGGGTTAGGGTTTAGGGTTAGGGTTGGTTAGGGTAGGGTTAGGGTTAGGCTTAGGGTTAGGGTTGGGTTAGGGTGTAGGTTAGGGTCGCGGTTAGGGTTAGGGTTCAGGTCTTTAGGTTAGCATTAGGGTTAGCGTTAGGGTTAGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTTAGGTAGGGTAGGGGTTAGGGAGGGGTTAGGGTTTAGGGTTAGGGGTAGGGTTAGGCGAGGGTTAGGGTTAGGGTTAGGGTAGGGTTAGGGTGGGTTATGGGTAGGGTTAGGGTTAGGGTTAGGGTTAGTGTAGGGTTGGGTTAGGGTTAGGGTTAGGGCTTAGGGTTAGGGTTCGCGGTAGGGTTAGGGTTAGGGTATAAGGGTTAGGGTTAGTTAAAGGGTTCAGGTGTTAGGGTTAGGTTAGGTTTAGGGTAGGGTAGGGGTTAGGGGTTAGGGTTGGGTTAGGTTTTAGGTTAGGGGTAGGGTCTTAGGGGTAGGGTTGAGGGTTAAGGGTTAGGGTTAGGGTTAGGGTTAGGAGTCGGGTTGGGTAGGGTTAGGGTGTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTTAGGGTTAGGTGTGTATGGGGTGGGGGGTGAGGGTTAGGGTTAGGTTTAGGGTTAGGTTAGGCTTAGGGAGGGGTTAGGTTAGGGTTTAGGGTTAGGGTTAGGGTTAGGTTAGGGTTAGGGTTAGGGTAGGGTTGAGGGTAGGGTTAGGGTTAGGTAGGGGTAGCCGGTAGGGTTAGGGTTAGGTATAGGGTTAGGGTAGGTTAGGGTTAGGTTAGTTAGGGTTAGGGTTGGGTTAGGGGTAGGGTTAGGGTTAGGCGTTAGTGTGGGTTAGGGTACGGGTTAGGTTATGTTTAGGGTGGGTTAGGGGTGGGTTAGGCGTTAGGGTCTAGGGTTTTAGGTTCGGGTTAGGGTTAGGGTTAGGGTTAAGGGTAGGGTTGGGTGGGGGGGGTTAGGGTTAGGGTTGCTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTATGGTTTCGGGGTGTGAGGTAGCGTTCG"
	motif 	= "TTAGGG"
	string = str(string)
	motif = str(motif)
	threshold = 0.35
	
	print(analyze_sequence_for_motif(string, motif, "apple", threshold=threshold))





