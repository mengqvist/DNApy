#!/usr/bin/python


#This file is part of DNApy. DNApy is a DNA editor written purely in python. 
#The program is intended to be an intuitive, fully featured, 
#extendible, editor for molecular and synthetic biology.  
#Enjoy!
#
#Copyright (C) 2014  Martin Engqvist | 
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#LICENSE:
#This file is part of DNApy.
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
#Get source code at: https://github.com/0b0bby0/DNApy
#

import dna
import random

def one_to_three(one_letter):
	'''
	Convert a one letter code amino acid to a three letter code.
	'''
	assert one_letter.upper() in 'FLSYCWPHERIMTNKVADQG*', 'Error, %s is not a valid amino acid' % one_letter
	
	AA = {'I':'Ile',
	'V':'Val',
	'L':'Leu',
	'F':'Phe',
	'C':'Cys',
	'M':'Met',
	'A':'Ala',
	'G':'Gly',
	'T':'Thr',
	'W':'Trp',
	'S':'Ser',
	'Y':'Tyr',
	'P':'Pro',
	'H':'His',
	'E':'Glu',
	'Q':'Gln',
	'D':'Asp',
	'N':'Asn',
	'K':'Lys',
	'R':'Arg',
	'*':'***'}	
	return AA[one_letter.upper()]
	
	
def three_to_one(three_letter):
	'''
	Convert a three letter code amino acid to a one letter code.
	'''
	assert three_letter.upper() in ['ILE','VAL','LEU','PHE','CYS','MET','ALA','GLY','THR','TRP','SER','TYR','PRO','HIS','GLU','GLN','ASP','ASN','LYS','ARG','STOP'], 'Error, %s is not a valid amino acid' % three_letter

	AA = {'ILE':'I',
	'VAL':'V',
	'LEU':'L',
	'PHE':'F',
	'CYS':'C',
	'MET':'M',
	'ALA':'A',
	'GLY':'G',
	'THR':'T',
	'TRP':'W',
	'SER':'S',
	'TYR':'Y',
	'PRO':'P',
	'HIS':'H',
	'GLU':'E',
	'GLN':'Q',
	'ASP':'D',
	'ASN':'N',
	'LYS':'K',
	'ARG':'R',
	'***':'*'}
	return AA[three_letter.upper()]
	
	
def one_to_full(one_letter):
	'''
	Convert one-letter amino acid code to full amino acid name.
	'''
	assert one_letter.upper() in 'FLSYCWPHERIMTNKVADQG*', 'Error, %s is not a valid amino acid' % one_letter
	AA = {'F':'Phenylalanine', 
	'L':'Leucine', 
	'S':'Serine', 
	'Y':'Tyrosine', 
	'*':'Stop', 
	'C':'Cysteine', 
	'W':'Tryptophan', 
	'P':'Proline', 
	'H':'Histidine', 
	'Q':'Glutamine', 
	'R':'Arginine', 
	'I':'Isoleucine', 
	'M':'Methionine', 
	'T':'Threonine', 
	'N':'Asparagine', 
	'K':'Lysine', 
	'V':'Valine', 
	'A':'Alanine', 
	'D':'Aspartic acid', 
	'E':'Glutamic acid', 
	'G':'Glycine'}
	return AA[one_letter.upper()]

	
def full_to_one(full):
	'''
	Convert full amino acid name to one-letter amino acid code.
	'''
	assert full.lower() in ['phenylalanine', 
							'leucine', 
							'serine', 
							'tyrosine', 
							'stop', 
							'cysteine', 
							'tryptophan', 
							'proline', 
							'histidine', 
							'glutamine', 
							'arginine', 
							'isoleucine', 
							'methionine', 
							'threonine', 
							'asparagine', 
							'lysine', 
							'valine', 
							'alanine', 
							'aspartic acid', 
							'glutamic acid', 
							'glycine'], 'Error, %s is not a valid amino acid' % full

	AA = {'phenylalanine':'F', 
			'leucine':'L', 
			'serine':'S', 
			'tyrosine':'Y', 
			'stop':'*', 
			'cysteine':'C', 
			'tryptophan':'W', 
			'proline':'P', 
			'histidine':'H', 
			'glutamine':'Q', 
			'arginine':'R', 
			'isoleucine':'I', 
			'methionine':'M', 
			'threonine':'T', 
			'asparagine':'N', 
			'lysine':'K', 
			'valine':'V', 
			'alanine':'A', 
			'aspartic acid':'D', 
			'glutamic acid':'E', 
			'glycine':'G'}
	return AA[full.lower()]
	
	
def three_to_full(three_letter):
	'''
	Convert amino acid three letter code to full amino acid names.
	'''
	return one_to_full(three_to_one(three_letter))

	
def full_to_three(full):
	'''
	Convert full amino acid names to three letter code.
	'''
	return one_to_three(full_to_one(full))
	
	
def count_aa(seq):
	'''
	Count occurrences of all amino acids in sequence. 
	The X character for unknown amino acid is allowed.
	Return as dictionary.
	'''
	seq = seq.upper()
	assert all([s in 'XIVLFCMAGTWSYPHEQDNKR*' for s in seq]) is True, 'Error, unknown amino acids %s in sequence: %s' % (str([s for s in seq if s not in 'XIVLFCMAGTWSYPHEQDNKR*']), seq)
	
	AA = {'I':seq.count('I'),
	'V':seq.count('V'),
	'L':seq.count('L'),
	'F':seq.count('F'),
	'C':seq.count('C'),
	'M':seq.count('M'),
	'A':seq.count('A'),
	'G':seq.count('G'),
	'T':seq.count('T'),
	'W':seq.count('W'),
	'S':seq.count('S'),
	'Y':seq.count('Y'),
	'P':seq.count('P'),
	'H':seq.count('H'),
	'E':seq.count('E'),
	'Q':seq.count('Q'),
	'D':seq.count('D'),
	'N':seq.count('N'),
	'K':seq.count('K'),
	'R':seq.count('R'),
	'*':seq.count('*'),
	'X':seq.count('X')}	
	return AA
	

	
def reverse_translate(prot_seq, table=1):
	'''
	Reverse translates protein sequence to DNA sequence using the specified codon table.
	For each amino acid the DNA codon is chosen randomly from those that encode the specified amino acid.
	table defaults to the standard codon table 1
	
	Input is a protein sequence in one-letter code as a string.
	Output is a DNA sequence as a string.
	table defaults to the standard codon table 1	
	'''
	assert type(prot_seq) is str, 'Error, the input must be a string containing amino acids in single letter code.'
	dna_seq = []
	for AA in prot_seq:
		possible = dna.GetCodons(AA, table=table, separate=False)
		dna_seq.append(random.choice(possible))
		
	return ''.join(dna_seq)
	