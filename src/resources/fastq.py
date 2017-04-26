#!/usr/bin/env python


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
