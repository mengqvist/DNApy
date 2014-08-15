#!/usr/bin/python

from itertools import groupby

def parse(filepath):
	"""
	Function for parsing a file containing fasta sequences.
	The input is the absolute filepath including the filename (as a string).
	The parsed data is returned as id, dna tuples in a generator.
	"""
	infile = open(filepath)
	faiter = (x[1] for x in groupby(infile, lambda line: line[0] == ">")) #make generator
	for header in faiter:
		id = header.next()[1:].strip() #skip the >
		seq = "".join(s.strip() for s in faiter.next()) # join all sequence lines to one.
		yield id, seq #return a generator