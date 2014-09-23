#!/usr/bin/python


def parse(filepath):
	"""
	Function for parsing a file containing fastv sequences.
	A fastv file is based of the fastq file format.
	The purpose of a fastv file is to assign a conservation score to each of the amino acids (or nucleotides) in the reference file to a set of homologs.
	
	First line starts with '@' and holds the reference sequence name.
	Second line holds the DNA or protein sequence.
	Third line starts with '+' and then gives how many sequences were used in the alignment (and therefore the calculation of the score)
	Fourth line holds uses ascii characters to assign conservation (in %) of that residue on a set of aligned sequences.
	After that follows the alignment of all sequences in fasta format.
	
	The input is the absolute filepath including the filename (as a string).
	The parsed data is returned as id, seq, qual_val tuples in a generator.
	"""
	infile = list(open(filepath))
	id, seq, num, qual_val = None, None, None, None
	alignments = []
	for i in range(len(infile)):
		line = infile[i].strip()
		if line: #line cannot be empty
			if line.startswith("@") and i==0: #first line, holds ID
				id = line[1:]
			elif id and not seq and not num and not qual_val and i == 1: #second line, holds the seq
				seq = line
			elif line.startswith("+") and i == 2: #third line, holds n used for the alignment
				num = line[1:]
			elif id and seq and num and i == 3:
				qual_val = line	
			elif i > 3:
				alignments.append(line)
	alignments = '\n'.join(alignments)
	#print(id, seq, num, qual_val)
	return id, seq, num, qual_val, alignments
					