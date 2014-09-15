#!/usr/bin/python

from itertools import groupby



def parseString(string):
	"""
	Function for parsing a file containing fasta sequences.
	The input is a string containing all the fasta records.
	The parsed data is returned as id, dna tuples in a generator.
	"""
#	faiter = (x[1] for x in groupby(string, lambda line: line[0] == ">")) #make generator
#	for header in faiter:
#		print('header', header)
		#id = header.next()[1:].strip() #skip the >
		#seq = "".join(s.strip() for s in faiter.next()) # join all sequence lines to one.
		#yield id, seq #return a generator
		
	records = string.split('>')
	records = [s for s in records if s != '']
	for record in records:
		temp = record.split('\n')
		temp = [s for s in temp if s != '']
		id = temp[0]
		seq = ''.join(temp[1:])
		yield id, seq
		
def parseFile(filepath):
	"""
	Function for parsing a file containing fasta sequences.
	The input is the absolute filepath including the filename (as a string).
	The parsed data is returned as id, dna tuples in a generator.
	"""
	f = open(filepath)
	string = f.read()
	return parseString(string)
