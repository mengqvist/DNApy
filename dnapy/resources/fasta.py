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


def read_fasta(iterator, seqtype=None):
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
		if seqtype == 'dna':
			yield (header, bioseq.DNA(''.join(seq)))
		elif seqtype == 'rna':
			yield (header, bioseq.RNA(''.join(seq)))
		elif seqtype == 'protein':
			yield (header, bioseq.Protein(''.join(seq)))
		elif seqtype == 'None':
			#should guess type here
			print('Cannot guess sequence type, specify "dna", "rna", or "protein"')
			raise ValueError


def parse_string(string, seqtype=None):
	"""
	Function for parsing a string containing fasta sequences.
	The input is a string containing all the fasta records.
	The parsed data is returned as id, sequence tuples in a generator.
	"""
	f = (x for x in re.split('\n', string))
	return read_fasta(f)


def parse_file(filepath, seqtype=None):
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
