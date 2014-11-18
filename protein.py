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
	three_letter = AA[one_letter.upper()]
	return three_letter
	
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
	one_letter = AA[three_letter.upper()]
	return one_letter
	
def count_aa(seq):
	'''
	Count occurrences of all amino acids in sequence. Return as dictionary.
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
	'X':seq.count('R'),
	'*':seq.count('*')}	
	return AA