def read_fasta(fp):
	'''
	Generator code stolen from Biopython. Parses a fasta file
	and yields a single fasta line on each call
	https://stackoverflow.com/questions/7654971/parsing-a-fasta-file-using-a-generator-python
	'''
	name, seq = None, []
	for line in fp:
		line = line.rstrip()
		if line.startswith(">"):
			if name: yield (name, ''.join(seq))
			name, seq = line, []
		else:
			seq.append(line)
	if name: yield (name, ''.join(seq))