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
