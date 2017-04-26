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


# this file is intended for sequence alignments of protein and DNA.
# The alignments are performed using MUSCLE which is an algorithm for global sequence alignment.
# You will need to have MUSCLE installed on your Linux system for this to work.
# For local sequence aligments use Water (TODO)



from io import StringIO
from muscle_wrapper import MuscleCommandline
import fasta


def align_fasta(fasta_string):
	'''
	Aligns all records from a string containing fasta entries.
	'''
	assert type(fasta_string) is str, 'Error, input must be a string.'
	#add more asserts

	muscle_cline = MuscleCommandline() #instantiate
	stdout, stderr = muscle_cline(stdin=fasta_string) #perform alignment
	aln = fasta.parseString(stdout) #parse output

	#add some assert to make sure the alignment went well

	#prepare for export
	output = []
	for entry in aln:
		output.append((entry[0], entry[1]))
	return output


def align_pairwise(id1, id2, seq1, seq2):
	'''
	Perform pairwise alignment.
	'''
	#add asserts for input...
	
	records = '>%s\n%s\n>%s\n%s' % (id1, seq1, id2, seq2) #prepare 'virtual' FASTA file
	return align_fasta(records)













