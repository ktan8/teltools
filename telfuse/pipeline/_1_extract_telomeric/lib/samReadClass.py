


class samRead:
	'''
	An object type to represent a SAM format read
	'''
	
	def __init__(self, samLine):
		samLine = samLine.strip()
		samLineArr = samLine.split("\t")
		self.samLine = samLine
		self.readName	= samLineArr[0]
		self.Flag	= samLineArr[1]
		self.Chr	= samLineArr[2]
		self.Posn 	= samLineArr[3]
		self.mapQ	= samLineArr[4]
		self.Cigar	= samLineArr[5]
		self.rnext	= samLineArr[6]
		self.pnext	= samLineArr[7]
		self.tlen	= samLineArr[8]
		self.Seq 	= samLineArr[9]
		self.QualVals	= samLineArr[10]
		self.readNumInPair = self.getReadNumFromFlag()
		

	def getNucleotideAtPosn(self, reqrPositionRelativeChr):
		'''
		Get the nucleotide of the base at the position
		'''
		posnRelativeRead = int(float(reqrPositionRelativeChr)) - int(float(self.Posn))
		if posnRelativeRead < 0:
			return '-'
		elif posnRelativeRead > len(self.Seq) - 1:
			return '-'
		else:
			return self.Seq[posnRelativeRead:posnRelativeRead+1]
		
	def getReadSequence(self):
		'''
		Get the actual read sequence after fixing the reverse complement
		problem.
		'''
		# Check if reverse strand
		if(int(self.Flag) & 16 == 16):
			return reverseComplement(self.Seq)
		else:
			return self.Seq
		
		
	def getReadQual(self):
		if(int(self.Flag) & 16 == 16):
			return reverseString(self.QualVals)
		else:
			return self.QualVals


	def getSamOutput(self):
		'''
		Returns the normal SAM output as a string
		'''
		return self.samLine


	def getReadNumFromFlag(self):
		'''
		Deduce if it is the first or second read in
		a pair from its flag value
		'''
		if(int(self.Flag) & 64 == 64):
			return 1
		elif(int(self.Flag) & 128 == 128):
			return 2
		else:
			return 0


	def setReadNum(self, val):
		'''
		Set the read number within the pair
		'''
		self.readNumInPair = val 


	def getReadNum(self):
		'''
		Get the read number within the pair
		'''
		return self.ReadNumInPair
	
	def getAlignmentOrientation(self):
		'''
		Get the orientation of the alignment
		from the cigar flag
		'''
		if(int(self.Flag) & 16 == 16):
			return "reverse"
		else:
			return "forward"

	
def reverseString(string):
	'''
	Returns input string in reverse order
	'''
	return string[::-1]


def complementString(nucString):
	'''
	Get the complement of a nucleotide string
	'''
	complementTable = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'N':'N' }
	stringComplement = ''
	for char in nucString:
		stringComplement = stringComplement + complementTable[char]

	return stringComplement
		

def reverseComplement(nucString):
	'''
	Get the reverse complement of a nucleotide string
	'''
	nucStringRev = reverseString(nucString)
	nucStringRevCom = complementString(nucStringRev)
	
	return nucStringRevCom



#string = "ATGCN"
#print string
#print reverseComplement(string)

#apple = samRead(3, 'abcdefghij')
#print apple.getNucleotideAtPosn(130)
