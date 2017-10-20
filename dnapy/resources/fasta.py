#!/usr/bin/env python3


#DNApy is a DNA editor written purely in python.
#The program is intended to be an intuitive and fully featured
#editor for molecular and synthetic biology.
#Enjoy!
#
#copyright (C) 2014-2015  Martin Engqvist |
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#LICENSE:
#
#DNApy is free software; you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation; either version 3 of the License, or
#(at your option) any later version.
#
#DNApy is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU Library General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program; if not, write to the Free Software Foundation,
#Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#Get source code at: https://github.com/mengqvist/DNApy
#


from dnapy.resources import bioseq
import re


def read_fasta(iterator):
	"""
	Function for iterating over lines to identify what constitutes the header
	and the entire sequence that goes with it.
	"""
	header, seq = None, []
	for line in iterator:
		line = line.rstrip()
		if line.startswith(">"):
			if header is not None:
				yield (header, ''.join(seq))
			header, seq = line, []
		else:
			seq.append(line)
	if header is not None:
		yield (header, ''.join(seq))


def parse_string(string):
	"""
	Function for parsing a string containing fasta sequences.
	The input is a string containing all the fasta records.
	The parsed data is returned as id, sequence tuples in a generator.
	"""
	f = (x for x in re.split('\n', string))
	return read_fasta(f)


def parse_file(filepath):
	"""
	Function for parsing a file containing fasta sequences.
	The input is the absolute filepath including the filename (as a string).
	The parsed data is returned as id, sequence tuples in a generator.
	"""
	if type(filepath) == str:
		with open(filepath) as f:
			for record in read_fasta(f):
				yield record
	elif type(filepath) == file:
		for record in read_fasta(filepath):
			yield record
	else:
		raise ValueError




class FreqTable(object):

	def __init__(self, fasta_file):
		self.make_codon_freq_table(fasta_file)


	def make_codon_freq_table(self, fasta_file):
		'''
		Input is a file path.
		Counts the usage of each codon in a FASTA file of DNA sequences.
		Then converts that as codon usage per 1000 codons.
		Good for generating codon tables.
		Output is a dictionary of codon frequencies per 1000 codons and the total number in brackets.
		'''

		records = fasta.parse_file(fasta_file)
		for record in records:
			cds = bioseq.DNA(record[1])
			codons = cds.count_codons()

		#sum codons
		codon_sum = sum(codons.values())

		#divide each by the sum and multiply by 1000
		freq_table = {}
		for key in list(codons.keys()):
			freq_table[key] = '%s(%s)' % (1000*(codons[key]/codon_sum), codons[key]) #ouput is following format: freq/thousand(number)
		self.table = freq_table
