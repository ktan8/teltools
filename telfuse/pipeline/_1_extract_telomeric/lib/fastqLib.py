import subprocess
import os

class fastq_chunk():
	'''
	A fastq chunk class that allows us to easily
	store and access information on a fastq read
	'''

	def __init__(self, fastq_iterator):
		self.readname = next(fastq_iterator)
		self.sequence = next(fastq_iterator)
		self.qualname = next(fastq_iterator)
		self.quality = next(fastq_iterator)

	def generate_fastq_string(self):
		'''
		Generate the fastq info from the
		fastq chunk info.
		'''
		fastq_string = self.readname + self.sequence + self.qualname + self.quality

		return fastq_string


def runProcess(exe):
	p = subprocess.Popen(exe, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

	for line in p.stdout:
		yield line.decode().strip()

	# This is necessary to tell OS to wait for exit status or we will
	# get a zombie process
	os.wait()


def fastqToSeqDict(fastq1, fastq2):
	'''
	Parse all sequences in the fastq files into a dict
	'''
	fastq1Cmd = "gunzip -c -d %s" %(fastq1)
	fastq2Cmd = "gunzip -c -d %s" %(fastq2)

	fastq1_iter = runProcess(fastq1Cmd.split())
	fastq2_iter = runProcess(fastq2Cmd.split())

	fastqDict = dict()

	while True:
		try:
			fastq1_chunk = fastq_chunk(fastq1_iter)
			fastq2_chunk = fastq_chunk(fastq2_iter)
			readSequences = [fastq1_chunk.sequence, fastq2_chunk.sequence]

			fastq_seq_name = fastq1_chunk.readname[1:]
			fastq_seq_name_clean = fastq_seq_name.split(" ")[0]
			#print(fastq_seq_name)
			fastqDict[fastq_seq_name_clean] = readSequences

		except StopIteration:
			break

	return fastqDict

