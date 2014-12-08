#!/usr/bin/env python



###
# name        : seqfiles.py
# description : Script for analyzing sequencing results
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LICENSE:
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#Copyright (C) 2014  Martin K. M. Engqvist | 
###


import os
import fnmatch
import dna


from StringIO import StringIO
import pprint
import copy
import re

import string


import ABIreader #parser
#import scf #parser
#import abi #parser
#import ztr #parser
import fastq #parser
import fasta #parser
import NeedlemanWunsch as NW #alignment
import pyalign




#TODO
#I need to make my own bindings to muscle so that I do not have to depend on biopython


##### Json style data structure used for this script #####
# [{
#       "reference": 
#		{
#		"name": "reference_dna",
#		"dna": "CACCGG",
#		"aln_dna": "C-ACCG-G",
#		},
#       "samples":
#		[
#		{
#		"name": "rrnb_primer1.seq",
#		"dna": "CTACGTG",
#		"orientation": "forward",
#		"aln_dna": "CTA-CGTG",	
#		"sign_missmatches": " * *  * ",
#		},
#		{
#		"name": "rrnb_primer2.seq",
#		"dna": "CACGTAG",
#		"aln_dna": "CTA-CGTG",	
#		"sign_missmatches": " * *  * ",
#		}
#		],
#		"contig": 
#		{
#		"name": "rrnb_contig"
#		"dna": "CTACGTG",
#		"aln_dna": "CTA-CGTG",
#		"sign_missmatches": " * *  * ",
#		"dna_missmatches": "list of dna missmatches",
#		"aa_missmatches": "list of aa missmatches",	
#		}
#		
#	
#     },
#     {
#       "reference": 
#		{
#		"name": "referene_dna",
#		"dna": "CACCGG",
#		"aln_dna": "C-ACCG-G",
#		},
#       "samples":
#		[
#		{
#		"name": "ccdb_primer1.seq",
#		"dna": "CTACGTG",
#		"orientation": "forward",
#		"aln_dna": "CTA-CGTG",
#		"sign_missmatches": " * *  * ",
#		},
#		{
#		"name": "ccdb_primer2.seq",
#		"dna": "CACGTAG",
#		"orientation": "reverse",
#		"aln_dna": "CTA-CGTG",	
#		"sign_missmatches": " * *  * ",
#		}
#		],
#	"contig": 
#		{
#		"name": "ccdb_contig"
#		"dna": "CTACGTG",
#		"aln_dna": "CTA-CGTG",
#		"sign_missmatches": " * *  * ",
#		"dna_missmatches": "list of dna missmatches"
#		"aa_missmatches": "list of aa missmatches",	
#		}
#		
#	
#     }
# ]
#######################################################





class SeqObj:
	'''
	Sequencing object for storing data related to a single Sanger sequencing run.
	filepath is the complete path to the file, including the file name.
	
	The file is parsed based on its file name ending.
	'''
	def __init__(self, filepath, name, primer, extension):
		#self.filepath = string.replace(filepath, "\\", "/") #replace slashes and hold the input filepath

		self.filepath = filepath #path to file
		self.filename = self.filepath.split('/').pop() #entire filename
		self.name = name #name of sequence
		self.primer = primer  #name of primer used to generate sequence
		self.input_type = extension.upper() #type of input file
		

		self.orientation = False #fw or rv orientation of the dna
		self.seq = False	#complete DNA sequence
		self.qual_val = False #contains the qualifier values for the sequence (if derived from an .ab1, .fastq or .scf file)
		self.trace = False #sequencing trace

		self.RC = False #keep track of whether the sequence has been reverse-complemented
		self.seq_clipped = False #clipped DNA sequence (removal of poor sequence (based on qual_val))

				
#		x =  '!"#$%%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~' #without escape characters: '!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~'

		self.getInput() #open the file
		return
	
	def setName(self, name):
		self.name = name #string
	def getName(self):
		return self.name #string
		
	def setPrimer(self, primer):
		self.primer = primer
	def getPrimer(self):
		return self.primer
		
	def setOrientation(self, orientation):
		self.orientation = orientation #string
	def getOrientation(self):
		return self.orientation #string
		
	def setDNA(self, seq):
		self.seq = seq #string
	def getDNA(self):
		return self.seq #string

	def setQualVal(self, qual_val):
		self.qual_val = qual_val #string
	def getQualVal(self):
		return self.qual_val #string	

	def setTrace(self, trace):
		self.trace = trace #list [G, A, T, C]
	def getTrace(self):
		return self.trace #list	[G, A, T, C]
		
	def setRC(self, bool):
		self.RC = bool #true or false
	def getRC(self):
		return self.RC #true or false
		
	def setDNAClipped(self, seq):
		self.seq_clipped = seq #string
	def getDNAClipped(self):
		return self.seq_clipped #string		

	
	def getInput(self):
		'''
		Open a single .seq, .fasta, .fastq, .ztr, .scf, .ab1 file (or even a text file with a DNA sequence) and set variables accordingly.
		'''

		#read the input
		if self.input_type in ['TXT', 'SEQ', 'SEQ.CLIPPED', None] and self.filename not in ['allseqs.txt']:
			f = open(self.filepath, 'r') 
			input = f.read() 
			f.close()
			self.setDNA(input.replace('\n', ''))
			
		elif self.input_type in ['AB1', 'ABI', 'ABIF']:
			ab1 = ABIreader.Trace(self.filepath, trimming=True) #optionally ', trimming=True'
			self.setDNA(ab1.seq)
			self.setQualVal(ab1.qual_val)
			self.setTrace([ab1.data['raw1'], ab1.data['raw2'], ab1.data['raw3'], ab1.data['raw4']]) #need to RC this too

#		elif self.input_type == 'ZTR':
#			print('Support for .ztr files has not yet been implemented')
			
#		elif self.input_type == 'SCF':
#			print('Support for .scf files has not yet been implemented')
			
		elif self.input_type is 'FASTA':
			id, seq = fasta.parseFile(self.filepath) #parse the fasta file. File should contain ONE entry
			self.setDNA(seq)
			
		elif self.input_type is 'FASTQ':
			id, seq, id2, qual_val = fastq.parse(self.filepath) #parse the fastq file. File should contain ONE entry
			self.setDNA(seq)
			self.setQualVal(qual_val)
		
		else:
			print('"%s" is not a .txt, .seq, .scf, .fasta, .fastq, .abif, .ab1, .abi or .ztr file' % self.filename)

			
			


class SeqAnalysis:
	'''
	Sequencing analysis object for assembling and analysing data from many Sanger sequencing runs.
	There should be one instance of this object for every physical real-world construct sequenced.
	
	For example: If you sequences two clones of the same construct with three primers each, 
	then there should be one SeqAnalysis object for each of these clones. 
	Each of the SeqAnalysis objects will hold data from the three sequencing runs done on that clone.
	
	'''
	def __init__(self, name):
		self.seqdata = {} #this is where all the sequence info gets stored
		self.reference = None #a SeqObj containing the reference DNA
		self.name = name #name of the construct that this SeqOverview object belongs to
		
		self.reference = None #reference to a SeqObj that holds the reference sequence
		self.aligned = False #keep track of whether sequences have been aligned or not
		
		self.consensus = None #reference to a SeqObj that represents the consensus sequence of which this sequence is a part
		self.sign_missmatches = None #sign (* *  *   *) mismatches vs the reference or consensus
		self.seq_missmatches = None #nucleotide mismatches vs the reference or consensus
		self.aa_missmatches	= None #amino acid mismatches vs the reference or consensus
		
		self.next = None #reference to a SeqObj that is sequentially before this one
		self.previous = None #reference to a SeqObj that is sequentially after this one

	
	def addSeq(self, seqobject):
		'''
		Adds a single SeqObj sequencing object to the data structure.
		'''

		#if a reference sequence has been added to the SeqAnalysis object, perpetuate that down to the single sequence
#		if self.getReference() != None:
#			seqobject.setReference(self.getReference())
			
		#add the object to the data structure.
		primer = seqobject.primer
		if self.seqdata.get(primer) == None:
			self.seqdata[primer] = {}
			self.seqdata[primer][seqobject.input_type] = seqobject
		else:
			self.seqdata[primer][seqobject.input_type] = seqobject


		## I should enforce that _ is an illegal character except for separating primer and construct names ##

		#this workaround works for construct names with _ in it. But it fails when the primer name has _ in it....
		### this is overly complicated to deal with filenames that has more than one '_' in the filename
#		nameparts = dictionary['name'].split('_')
#		nameparts.pop() #remove last
#		name = ''
#		for item in nameparts: #stitch them together
#			name += '%s_' % item
		###	
		
		#construct name and primer name MUST be separated by _ . Construct name CANNOT contain _ .
	#	name = seqobject.getName().split('_')[0] #get first part of name




	
	def setReference(self, seq):
		self.reference = seq #string
	def getReference(self):
		return self.reference #string		

	def setSignMissmatch(self, string):
		self.sign_missmatches = string #string
	def getSignMissmatch(self):
		return self.sign_missmatches #string	

	def setDNAMissmatch(self, string):
		self.seq_missmatches = string #string
	def getDNAMissmatch(self):
		return self.seq_missmatches #string

	def setAAMissmatch(self, string):
		self.aa_missmatches = string #string
	def getAAMissmatch(self):
		return self.aa_missmatches #string

	def setContig(self, seq):
		self.contig = seq #string
	def getContig(self):
		return self.contig #string		
		
	def setNext(self, seqobj):
		self.next = seqobj #SeqObj object
	def getNext(self):
		return self.next #SeqObj object

	def setPrevious(self, seqobj):
		self.previous = seqobj #SeqObj object
	def getPrevious(self):
		return self.previous #SeqObj object


		
	def assemble(self):
		if self.getReference() is None: #if no reference sequence is present, do de Novo assembly.
			self.assemble_de_novo()
		else: 
			self.assemble_to_ref() #if a reference sequence is present, align each sequence to it.
	
		#check missmatches and set the associated variables

		self.aligned = True

	def assemble_de_novo(self):
		'''
		
		'''
		
		
		
		#somewhere I need to pick which type of file to use if there are multiple. I'm gonna go with the AB1 data for now.
		
		#i'm starting with the simple case of two sequences
		seq1 = self.seqdata['p423-GPD-RV']['AB1'].getDNA()
		seq2 = self.seqdata['p423-GPD-FW']['AB1'].getDNA()
		print(self.getOverlap(seq1, seq2))
		print(self.getOverlap(seq1, dna.RC(seq2)))
	
			
		#### I may adapt the structure below #########
		#take one seq, find the best match, combine them. Start over.
#		for i in range(len(seqdata)): #for all of the sequence groups
#			fw_list = []
#			rv_list = []
#			seqdata[i]['contig']['dna'] = ''
#			samples = copy.deepcopy(seqdata[i]['samples']) #copy samples so that I can delete them as I go
#			overlap_found = False #has overlap been found?
#			whole_list = False #did we traverse the whole list?
#			while overlap_found is False and whole_list is False:
#				overlap_found = False #has overlap been found?
#				whole_list = False #did we traverse the whole list?
#				for n in range(len(samples)):
#					for o in range(n, len(samples)):
#						fw_dict = {}
#						seq1 = seqdata[i]['samples'][0]['dna']
#						seq2 = seqdata[i]['samples'][n]['dna']
#						forward = findoverlap(seq1, seq2)
#						if forward is not False:
#							overlap_found = True
#						fw_list.append(forward[2])
#						
#						seq1 = seqdata[i]['samples'][0]['dna']
#						seq2 = seqdata[i]['samples'][n]['dna']
#						reverse = findoverlap(seq1, dna.reversecomplement(seq2))
#						if reverse is not False:
#							overlap_found = True
#						rv_list.append(forward[2])
#					if overlap_found is True:
#						if max(fw_list) >= max(rv_list):
#							index, value = max(enumerate(fw_list))
#						elif max(fw_list) < max(rv_list):
#							index, value = max(enumerate(fw_list))
						
							
#				whole_list = True
		
	def assemble_to_ref(self):
		'''
		This method aligns many sequences to a reference, one sequence at a time, then makes sure that all gaps match.
		'''
		
		############# need to fix this method so that it works with the new data structure #############
		
		#################################################################################################
		
		
		temp_dictlist = []
		for n in range(len(seqdata[index]['samples'])): #make sure each dna has an empty aln_dna entry
			seqdata[index]['samples'][n]['aln_dna'] = '' 

		#make virtual FASTA file of two files and align. The first entry should be the reference.
		for n in range(len(seqdata[index]['samples'])): 
			records = ''
			if seqdata[index]['reference']['name'][0] != '>':
				records += '>%s\n%s\n' % (seqdata[index]['reference']['name'], seqdata[index]['reference']['dna']) #reference
			elif seqdata[index]['reference']['name'][0] == '>':
				records += '%s\n%s\n' % (seqdata[index]['reference']['name'], seqdata[index]['reference']['dna'])
			else:
				print('Muscle name error')	
				
			if seqdata[index]['samples'][n]['name'][0] != '>':
				records += '>%s\n%s\n' % (seqdata[index]['samples'][n]['name'], seqdata[index]['samples'][n]['dna']) #sample
			elif seqdata[index]['samples'][n]['name'][0] == '>':
				records += '%s\n%s\n' % (seqdata[index]['samples'][n]['name'], seqdata[index]['samples'][n]['dna'])
			else:
				print('Muscle name error')

			records_handle = StringIO(records) #turn string into a handle
			tempdata = records_handle.getvalue()
		
			#for seperate fasta entries
			muscle_cline = MuscleCommandline()
			stdout, stderr = muscle_cline(stdin=tempdata)
			stdout = parse_fasta(stdout)	


			#sort so that ref is first
			#is that needed? seems to be working fine
			

			if n == 0:
				seqdata[index]['reference']['aln_dna'] = stdout[0]['dna']
				seqdata[index]['samples'][n]['aln_dna'] = stdout[1]['dna']

			
			else:
			
				#compare sequences with temp dictlist, check for - in one and not the other. Make changes to all.
				temp_dictlist.append(stdout[0])
				temp_dictlist.append(stdout[1])
				i = 0
				while i <= len(seqdata[index]['reference']['aln_dna']):

					
					#check if reference sequences have the same spaces
					if len(seqdata[index]['reference']['aln_dna']) == i and len(temp_dictlist[0]['dna']) == i:
						i += 1 
					
					elif len(seqdata[index]['reference']['aln_dna']) == i and len(temp_dictlist[0]['dna']) > i: #for dealing with end of sequence
						seqdata[index]['reference']['aln_dna'] = seqdata[index]['reference']['aln_dna'] + '-'
						for x in range(len(seqdata[index]['samples'])): #make change in all sequences present
							if seqdata[index]['samples'][x]['aln_dna'] != '':
								seqdata[index]['samples'][x]['aln_dna'] += '-'
						i = 0

					elif len(temp_dictlist[0]['dna']) == i and len(seqdata[index]['reference']['aln_dna']) > i: #for dealing with end of sequence
						for entry in temp_dictlist:
							entry['dna'] = entry['dna'] + '-'
						i = 0				

					elif seqdata[index]['reference']['aln_dna'][i] == '-' and temp_dictlist[0]['dna'][i] != '-': # if gap in old alignment, but not in new
						for entry in temp_dictlist:
							entry['dna'] = entry['dna'][:i] + '-' + entry['dna'][i:]
						i = 0

					elif seqdata[index]['reference']['aln_dna'][i] != '-' and temp_dictlist[0]['dna'][i] == '-': #if gap in new alignment, but not in old
						seqdata[index]['reference']['aln_dna'] = seqdata[index]['reference']['aln_dna'][:i]	+ '-' + seqdata[index]['reference']['aln_dna'][i:]			
						for x in range(len(seqdata[index]['samples'])): #make change in all sequences present
							if seqdata[index]['samples'][x]['aln_dna'] != '':

								seqdata[index]['samples'][x]['aln_dna'] = seqdata[index]['samples'][x]['aln_dna'][:i] + '-' + seqdata[index]['samples'][x]['aln_dna'][i:]
						i = 0

					else:
						i += 1
				
				seqdata[index]['samples'][n]['aln_dna'] = temp_dictlist[1]['dna']
				
		#I should probably add something to go through and check that I don't have '-' for all sequences at some position

			
	def findoverlap(self, Seq1, Seq2, min_overlap=20):
		"""Function for finding overlaps of two sequences.
			Returns start of overlap on Seq1, start of overlap on Seq2, and length of overlap"""
		Seq1 = Seq1.upper()
		Seq1 = Seq1.replace('\n','')

		Seq2 = Seq2.upper()
		Seq2 = Seq2.replace('\n','')

		seq_matcher = difflib.SequenceMatcher(None, Seq1, Seq2)
		seq1_loc, seq2_loc, match_len = seq_matcher.find_longest_match(0, len(Seq1), 0, len(Seq2))
		if match_len < min_overlap: #the match is shorter than the minimum specified
			return False
		else:
			return seq1_loc, seq2_loc, match_len #return overlap DNA, overlap start on first seq, overlap start on second seq
		
##############################################################################################	
	def getOverlap(self, seq1, seq2):
		'''
		Find and return the maximum overlap length of two sequences.
		This only works when using the NeedlemanWunsch algorithm to get the overlap. (The algorithm is much too promiscous otherwise.
		If there is no overlap, return empty string.
		'''
		alignment = NW.PairwiseAlignment(seq1, seq2)
		score = alignment.score
		
		seq1aln = alignment.seq1aligned.upper()
		seq2aln = alignment.seq2aligned.upper()
		assert len(seq1aln) == len(seq2aln), 'Error, the sequences are not of the same length'
				
		overlap = False
		start = False
		end = len(seq1aln)
		length = False
		first = False
	
		#get start, end and length of overlap. Any number of N is tolerated and a double -- is tolerated.
		for i in range(len(seq1aln)):
			if overlap is False and seq1aln[i] in 'ATCGN' and seq2aln[i] in 'ATCGN': #first nucleotide of overlap
				start = copy.copy(i)
				overlap = True
				
			elif overlap is True and (seq1aln[i] == '-' or seq2aln[i] == '-'):
				if seq1aln[i+1:i+3] == '--' or seq2aln[i+1:i+3] == '--': #allow for two missing bases, but not more
					end = copy.copy(i)
					break
#		a = open('test.txt', 'a') #open it for writing	
#		a.write(seq1aln)
#		a.write(seq2aln)
#		a.close()
		#Find which sequence is the first (the leftmost)
		#possible topologies:

		#AAAAAAAAAAAAAATTTTT-----------
		#--------------AAAAACCCCCCCCCCC
		if (seq1aln[0] in 'ATCGN' and seq2aln[0] == '-') and (seq1aln[end+1] == '-' and seq2aln[end+1] in 'ATCGN'):
			first = 1 #seq1 first
			
		#--------------AAAAACCCCCCCCCCC
		#AAAAAAAAAAAAAATTTTT-----------
		elif (seq1aln[0] == '-' and seq2aln[0] in 'ATCGN') and (seq1aln[end+1] in 'ATCGN' and seq2aln[end+1] == '-'):
			first = 2 #seq2 first
		
		#CCCCCCAAAAAAAAAAAAAATTTTTCCCCCCCCCCC		
		#------AAAAAAAAAAAAAATTTTT-----------
		elif (seq1aln[0] in 'ATCGN' and seq2aln[0] == '-') and (seq1aln[end+1] in 'ATCGN' and seq2aln[end+1] == '-'):
			overlap = False #mark this as not an overlap

		#------AAAAAAAAAAAAAATTTTT-----------
		#CCCCCCAAAAAAAAAAAAAATTTTTCCCCCCCCCCC			
		elif (seq1aln[0] == '-' and seq2aln[0] in 'ATCGN') and (seq1aln[end+1] == '-' and seq2aln[end+1] in 'ATCGN'):
			overlap = False #mark this as not an overlap

		#CCCCCCAAAAAAAAAAAAAATTTTTCC---------
		#CCCCCCAAAAAAAAAAAAAATTTTTCCCCCCCCCCC		
		elif (seq1aln[0] in 'ATCGN' and seq2aln[0] in 'ATCGN') and (seq1aln[end+1] == '-' and seq2aln[end+1] in 'ATCGN'):
			first = 1 #seq1 first

		#CCCCCCAAAAAAAAAAAAAATTTTTCCCCCCCCCCC				
		#CCCCCCAAAAAAAAAAAAAATTTTTCC---------
		elif (seq1aln[0] in 'ATCGN' and seq2aln[0] in 'ATCGN') and (seq1aln[end+1] in 'ATCGN' and seq2aln[end+1] == '-'):
			first = 2 #seq2 first
			
		#-----CAAAAAAAAAAAAAATTTTTCCCCCCCCCCC
		#CCCCCCAAAAAAAAAAAAAATTTTTCCCCCCCCCCC
		elif (seq1aln[0] == '-' and seq2aln[0] in 'ATCGN') and (seq1aln[end+1] in 'ATCGN' and seq2aln[end+1] in 'ATCGN'):
			first = 2 #seq2 first

		#CCCCCCAAAAAAAAAAAAAATTTTTCCCCCCCCCCC
		#-----CAAAAAAAAAAAAAATTTTTCCCCCCCCCCC
		elif (seq1aln[0] in 'ATCGN' and seq2aln[0] == '-') and (seq1aln[end+1] in 'ATCGN' and seq2aln[end+1] in 'ATCGN'):
			first = 1 #seq1 first
			
		#CCCCCCAAAAAAAAAAAAAATTTTTCCCCCCCCCCC
		#CCCCCCAAAAAAAAAAAAAATTTTTCCCCCCCCCCC
		elif (seq1aln[0] in 'ATCGN' and seq2aln[0] in 'ATCGN') and (seq1aln[len(seq1aln)] in 'ATCGN' and seq2aln[len(seq2aln)] in 'ATCGN'):
			first = 1 #it does not matter which one is first, but let's pick seq1		

			
		if overlap is False:
			return False, False, False, False
		else:
			return seq1aln, seq2aln, seq1aln[start:end], first, score #return aligned seq1 (str), seq2 (str), the overlap (str), and an integer (1 or 2 ) that indicates which sequence is first in the alignment.


			
	def findFirst(self, alnscores):
		'''
		Find the leftmost sequence of alignment pairs.
		'''
		#find the leftmost sequence
		keys = alnscores.keys() #seqs with overlaps to the right
		leftmost = []
		for i in keys:
			present = False
			for j in keys:
				if i in alnscores[j]: #they cannot be present as the 'right' sequence of any other sequence
					present = True
					break
			if present is False and i not in leftmost:
				leftmost.append(i)
		print('leftmost', leftmost)
		return leftmost # a list


	def sortSeqs(self):
		'''
		Align each sequence with each other sequence. 
		Save the alignment scores in a matrix where the index (in the list) for each SeqObj is used as the identifier.
		'''

		alnscores = {} #store the alignment overlaps
		alnseqs = {} #store the alignment sequences
		for group in self.seqdata:
			for i in range(len(group['samples'])-1):
				for j in range(i+1, len(group['samples'])):
					seqobj1 = group['samples'][i]
					seqobj2 = group['samples'][j]
					alnseq1, alnseq2, overlap, first = self.getOverlap(seqobj1.getDNA(), seqobj2.getDNA()) #get the overlap and which sequence is first
	#					print('%s and %s' % (seqobj1.getName(), seqobj2.getName()), overlap)
					if overlap is not False and first == 1:
						if i in alnscores:
							alnscores[i].update({j:len(overlap)})
							alnseqs[i].update({j:(alnseq1, alnseq2)})
						else:
							alnscores[i] = {j:len(overlap)}
							alnseqs[i] = {j:(alnseq1, alnseq2)}
					elif overlap is not False and first == 2:
						if j in alnscores:
							alnscores[j].update({i:len(overlap)})
							alnseqs[j].update({i:(alnseq2, alnseq1)})
						else:
							alnscores[j] = {i:len(overlap)}
							alnseqs[j] = {i:(alnseq2, alnseq1)}
			print(alnscores)
			
			#get the leftmost (the first) sequence
			leftmost = self.findFirst(alnscores)[0] 
			
			#determine the order of sequences
			sequence = [leftmost]
			keys = alnscores.keys() 
			while sequence[-1] in keys:
				sequence.append(max(alnscores[sequence[-1]], key=alnscores[sequence[-1]].get)) #get the max value for the dictionary under the key
			print('sequence', sequence)

		

	def setReference(self, seq):
		'''
		Set a reference sequence.
		'''
		self.reference = seq
		
	def getReference(self):
		'''
		Return the reference sequence.
		'''
		return self.reference
		
		
	def printAlnScores(self):
		pass

		
	def setConsensus(self):
		'''
		Build a consensus sequence from sorted sequence reads.
		'''
		pass
			
	


	
			
################ Generating outputs ######################

	def writeseqfiles(self, path):
		'''Write all sequence entries from a list of dictionaries to textfile'''
		a = open(path + 'allseqs.txt', 'w') #open it for writing
		for i in range(len(self.seqdata)):
			for n in range(len(self.seqdata[i]['samples'])):
				a.write(self.seqdata[i]['samples'][n].getName())	
				a.write('\n')
				a.write(self.seqdata[i]['samples'][n].getDNA())
				a.write('\n')
		a.close()


##########################################################

			

class SeqOverview:
	'''
	The SeqOverview object manages one or several instances of SeqAnalysis and provides means to interact with these.
	It should be possible to load a whole folder of files into this class and have it organize the sequences by name and make 
	SeqAnalysis instances as appropriate.
	The aim is to have a GUI interact with this class.
	'''

	def __init__(self):
		self.data = {} #stores all the data in a dictionary
		
	def addFile(self, filepath):
		'''
		Add a single file.
		'''

		#get the filename. The filename may only have one _ character,
		#and that character MUST be used to separate the construct name from the primer name.
		filename = re.split('[\\\/]', filepath)[-1] #the entire filename
		extension = '.'.join(filename.split('.')[1:]) #get the extension
	
		if extension.upper() in ['TXT', 'AB1', 'SCF', 'FASTQ', 'SEQ', 'SEQ.CLIPPED']: #here I determine which files are allowed
			print('Adding: %s' % filename)
			name, primer = filename.replace('.' + extension, '').split('_') #get construct and primer names
			
			#check whether any other SeqAnalysis object with that name is already present, if not, make one.
			instance = self.data.get(name)
			if instance == None: #not present, so add new SeqAnalysis instance.
				self.data[name] = SeqAnalysis(name)
				self.data[name].addSeq(SeqObj(filepath, name, primer, extension)) #make a new SeqObj with the info and add it to the SeqAnalysis object.
				
			elif instance != None: #is present, so add sequence to the existing instance.
				instance.addSeq(SeqObj(filepath, name, primer, extension)) #make a new SeqObj with the info and add it to the SeqAnalysis object.
		else:
			print('Skipping: %s' % filename)
		
		
	def addFolder(self, path):
		'''
		Add an entire folder worth of files.
		'''
		file_list = sorted(os.listdir(path))
		for filename in file_list: 
			filepath = path+filename
			self.addFile(filepath)
	
	
	def addList(self, path, file_list):
		'''
		Add a specified files from a path.
		'''
		for filename in file_list: 
			filepath = path+'/'+filename
			self.addFile(filepath)		
	
	
	
	
def RemoveN(Seq):
	Seq = Seq.upper()
	BeginningCounter = 0
	EndCounter = len(Seq)
	for i in range(int(len(Seq)/2)):
		if 'N' in Seq[i]:
			BeginningCounter = i+1
	for i in range(int(len(Seq)/2)+1, len(Seq)):
		if 'N' in Seq[i]:
			EndCounter = i
			break
	
	Seq = Seq[BeginningCounter:EndCounter]
	#print(Seq)
	return Seq
#End of RemoveN function



		




def join(Seq1, Seq2):
	"""joins two sequences that have an overlap.
		Seq1 is assumed to come before Seq2.
		Both are assumed to be in the same orientation."""
	vars = findoverlap(Seq1, Seq2) #vars contains the start of overlap on seq1, start of overlap on seq2, and length of overlap
	#print(vars)
	if vars == False:
		print('No overlap')
		return False
	elif type(var[0]) == int and type(var[1]) == int and type(var[2]) == int:
		print('Overlap')
		JointSeq = Seq1[0:vars[0]] + Seq2[vars[1]:len(Seq2)] 
		return JointSeq
	else:
		print('error while joining')


		
	

#align function
def align(Seq1, Seq2): #aligns two sequences that have an overlap
	Seq1 = Seq1.replace('\n','')
	Seq2 = Seq2.replace('\n','')
	Seq1 = Seq1.upper()
	Seq2 = Seq2.upper()
	vars = findoverlap(Seq1, Seq2)
	#print(vars)
	if vars == False:
		#print('No overlap')
		return False
	elif type(vars[0]) == int and type(vars[1]) == int and type(vars[2]) == int:
		dash = '-'
		dash = dash * (var[1]-var[2]) #putting dashes in front of seq2 to fill up until the alignment
		Seq2 = dash + Seq2 
		alnvar = ''

		if len(Seq1) < len(Seq2):
			Seq1 = Seq1 + (len(Seq2) - len(Seq1)) *  '-' #fill up the back
			#print(len(Seq1))
			#print(len(Seq2))
		elif len(Seq1) > len(Seq2):
			Seq2 = Seq2 + (len(Seq1) - len(Seq2)) *  '-' #fill up the back
			#print(len(Seq1))
			#print(len(Seq2))

		for i in range(len(Seq1)):
			if Seq1[i] == Seq2[i]:
				alnvar = alnvar + '|'
			elif Seq1[i] != Seq2[i]:
				alnvar = alnvar + 'x'
			else:
				print('error making alnvar')
		
		return (Seq1, alnvar, Seq2)

	else:
		print('error while aligning')
#end of align function



	





#Function for finding a certain sequence in a string and delete everything in front. Used for finding beginning of gene in seqfile.
#Currently uses 10 nucleotides for the search
#Also deletes everything after the stop codon of the RefSeq
def findstart(RefSeq, Seq):
	RefSeq = RefSeq.upper()
	Seq = Seq.upper()
	
	for i in range(len(RefSeq)):
		if RefSeq[i:(i+3)] == 'ATG':
			RefSeq = RefSeq[i:len(RefSeq)]	
			#print(RefSeq[0:10])
			break
	for i in range(len(Seq)):
		if RefSeq[0:10] == Seq[i:(i+10)]: #Find the start here
			Seq = Seq[i:len(Seq)]
		if RefSeq[len(RefSeq)-11:len(RefSeq)-1] == Seq[i:(i+10)]: #find the end. Had to add an extra -1 since last character was \n
			Seq = Seq[0:(i+10)]
	#print(Seq[0:10])
	return(RefSeq, Seq)	
#End of findstart function





#Function for comparing the sequencing reaction to a reference sequence, one triplet at a time
def tripletalign(RefSeq, Seq):
	RefSeq = RefSeq.upper()
	Seq = Seq.upper()
	counter = 0
	for i in range(len(Seq)):
		if i%3==0 or i == 0:
			#print(RefSeq[i:i+3])
			#print(Seq[i:i+3])
			if RefSeq[i:(i+3)] == Seq[i:(i+3)]:
				counter = 0
				pass			
			elif i+3>len(Seq): #to make sure it's a complete codon... i.e. leave out last codon if it is shorter than 3
				pass	
			else:					
				if counter <= 5:
					print('DNA pos %d, %s mutated to %s --- %s%d%s' % (i+3-1, RefSeq[i:(i+3)], Seq[i:(i+3)], dna.translate(RefSeq[i:(i+3)]), int((i+3)/3), dna.translate(Seq[i:(i+3)])))
					counter += 1
				else:
					print('over %d consecutive mismatches, rest of construct is likely out of frame' % (counter-1))
					break
	print('\n')
#End of triplet align function


def check(seqdataentry):
	#for entry in Sequences:
	print(seqdataentry['contig']['name'])
	dna = findstart(seqdataentry['reference']['dna'], seqdataentry['contig']['dna'])
	tripletalign(dna[0], dna[1]) #(RefSeq, dna)


### Add a Seqalign function here that analyzes the current sequence vs refseq and spits out missmatches independent of codons...
#########






def getreference(path, index):
	'''Get the reference sequence for each group of sequences'''
	#open the RefSeq file
	#RefSeq file should contain the gene of interest, starting with ATG
	filepath = path + 'RefSeq.txt' #Open reference sequence file
	f = open(filepath, 'r') #open content of current file
	seqdata[index]['reference'] = dict(name = 'reference_sequence', dna = f.read()) 
	f.close()		#close current file	

#######################################################









#I don't think this is needed any more
def group_sequences():
	'''groups sequences based on name'''
	Sequences = []
	Currentsequence = []
	filepath = path + '/allseqs.txt' #Open reference sequence file
	f = open(filepath, 'r') #open content of current file
	for line in f:
		line = line.replace('\n', '')
		if line[0] == '>':
			if Currentsequence == []:
				name, rest = line.split('_')
				Currentsequence.append(name)
			
			elif Currentsequence[0] in line:
				Sequences.append(Currentsequence)
				
			else:
				Currentsequence = []
				name, rest = line.split('_')
				Currentsequence.append(name)
				
		elif line[0] == 'A' or line[0] == 'T' or line[0] == 'C' or line[0] == 'G':
			Currentsequence.append(line)
		elif line[0] == 'N':
			print('Sequence starts with N')
		else:
			print('Something is wrong')
	f.close()
	return Sequences




			
			
	
def analyze():
	###### here I call the functions #########

	alnresults = path + '/' + 'alnresults.txt'  #get the path for the results file (output file)
	f = open(alnresults, 'w') #open it for writing

	results = path + '/' + 'results.txt'  #get the path for the results file (output file)
	a = open(results, 'w') #open it for writing


	for i in range(len(seqdata)):
		f.write('>%str\n' % seqdata[i]['reference']['name'])
		f.write('%s\n' % seqdata[i]['reference']['aln_dna'])
		n = 1	
		while n < len(seqdata[i]['samples']):
			f.write('%s\n' % seqdata[i]['samples'][n]['name'])		
			f.write('%s\n' % seqdata[i]['samples'][n]['aln_dna'])		

####### fix the assembly method
####### save it in 'contig'

			assembled = join(entry[n], entry[i])
			if assembled==False:
				n = i
			else:		
				entry[i] = ''
				entry[n] = assembled
#		for i in range(len(entry)):
#			if i == 0:
#				a.write('>' + entry[i])
#				a.write('\n')
#			else:
#				a.write(entry[i])
#				a.write('\n')
#						

			#check whether any mutations are present

		#this is a shortcut for now
		seqdata[i]['contig'] = {'name': seqdata[i]['samples'][0]['name']}
		seqdata[i]['contig']['dna'] = seqdata[i]['samples'][0]['dna']
		check(seqdata[i])

	f.close()
	a.close()

	
	
if __name__ == '__main__': #if script is run by itself and not loaded	
	import sys
	assert len(sys.argv) == 2, 'Error, this script requires a path to a folder containing the sequencing files as an argument.'
	print('Opening %s' % str(sys.argv[1]))
	path = str(sys.argv[1]) #Path to folder that contains the sequences

	
	# There are two modes of running the software. Either by aligning sequences to a reference,
	# or by assembling the sequences by finding overlaps.
	
		
	#One should choose what type of sequence should be loaded. ab1>scf>fastq>fasta
	#In general, the sequences should be loaded and then clipped to remove bad sequence stretches. 
	#Otherwise these would really mess up the process of finding overlaps.
	#Then align to reference sequence OR assemble based on overlaps.
	
	#Then missmatches and N's need to be found. The process should depend on which mode is used.
	#If the sequences are aligned to reference, find missmatches and N's throughout the sequence.
	#If the sequences are assembled, then find missmatches in the overlap regions.
	
#	FindSingleGeneMissmatch = 'Y' #Test for one gene whether missmatches occur or not. It is expected that RefSeq starts with ATG
	x = SeqOverview()
	x.addFolder(path)
	
	#assemble the sequences
#	for key in x.keys():
#		x[key].assemble()
	
	
	

