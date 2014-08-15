#!/usr/bin/python

def parse(filepath):
	"""
	Function for parsing a file containing fastq sequences.
	The input is the absolute filepath including the filename (as a string).
	The parsed data is returned as id, dna, qual_val tuples in a generator.
	"""
	infile = open(filepath)
	id, dna, id2, qual_val = None, None, None, None
	for line in infile:
		line = line.strip()
		if line: #line cannot be empty
			if not id and line.startswith("@"): #first id line
				id = line[1:]
			elif id and line.startswith("@"): #first id line of the next entry
				assert len(dna) == len(qual_val), "Error, the length of DNA and the associated qual_val in %s are not the same." % id
				yield id, dna, id2, qual_val #return a generator
				id = line[1:]
				dna, id2, qual_val = None, None, None
			elif line.startswith("+"): #id2 line
				if len(line) == 1:
					id2 = id
				else:
					id2 = line[1:]
			elif not id2: #if id2 has not yet been assigned
				dna = line
			elif id and id2 and dna: #only qualifier values line remains
				qual_val = line		