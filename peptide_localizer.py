#!/usr/bin/env python


#This file is part of DNApy. DNApy is a DNA editor written purely in python. 
#The program is intended to be an intuitive, fully featured, 
#extendable, editor for molecular and synthetic biology.  
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


#This script lends ideas and code snippets from the oligo finder written by Andrea Cabibbo andrea.cabibbo@uniroma2.it   
#http://www.cellbiol.com/scripts/oligo/oligo_motif_sequence_finder.php    


import re
import string

def re_const(peptide):
    '''This is the actual search function'''
    peptide=peptide.upper()
    re_out=[]
    for char in peptide:
        if char in 'ACGILMPSTVFRYKW':
            re_out.append('[%sX]' % char)
        elif char in 'ND':
            re_out.append('[%sBX]' % char)
        elif char in 'QE':
            re_out.append('[%sZX]' % char)
        elif char == 'B':
            re_out.append('[%sNDX]' % char)
        elif char == 'Z':
            re_out.append('[%sQEX]' % char)
        elif char == 'X':
            re_out.append('[ACGILMPSTVFNRYDEKQWBZX]')
    re_out=string.join(re_out,'')
    re_out="(?=(%s))" % re_out
    re_out_comp=re.compile(re_out,re.IGNORECASE)
    return re_out_comp


def match_peptide(seq,peptide,mismatches=0):
	'''Function for searching for a certain oligo'''
	print('seq', seq)
	re_peptide=re_const(peptide)
	location_out=[]
	for match in re_peptide.finditer(seq):
		location_out.append([match.start()+1, match.end()+len(match.group(1))])
	return location_out
