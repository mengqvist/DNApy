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




import os
import fnmatch
import dna as DNA

#from Bio.Align.Applications import MuscleCommandline
from StringIO import StringIO
#from Bio import AlignIO
import pprint
import copy

import string

#import difflib

import ABIreader #parser
#import scf #parser
#import abi #parser
#import ztr #parser
import fastq #parser
import fasta #parser
import NeedlemanWunsch as NW #alignment




class SeqObj:
	'''
	Sequencing object for storing data related to a single Sanger sequencing run.
	filepath is the complete path to the file, including the file name
	'''
	def __init__(self, filepath):
		self.filepath = string.replace(filepath, "\\", "/") #replace slashes and hold the input filepath
		self.input_type = False #type of input file

		self.name = False #name of sequence
		self.orientation = False #fw or rv orientation of the dna
		self.dna = False	#complete DNA sequence
		self.qual_val = False #contains the qualifier values for the sequence (if derived from an .ab1 or .scf file)
		self.trace = False #sequencing trace
		
		self.dna_clipped = False #clipped DNA sequence (removal of poor sequence (based on qual_val))
		self.reference = False #reference to a SeqObj that holds the reference sequence
		self.sign_missmatches = False #sign (* *  *   *) mismatches vs the reference
		self.dna_missmatches = False #nucleotide mismatches vs the reference
		self.aa_missmatches	= False #amino acid mismatches vs the reference

		self.consensus = False #reference to a SeqObj that represents the consensus sequence of which this sequence is a part
		self.next = False #reference to a SeqObj that is sequentially before this one
		self.previous = False #reference to a SeqObj that is sequentially after this one
	
#		x =  '!"#$%%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~' #without escape characters: '!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~'

		self.getInput() #open the file
		return
	
	def setName(self, name):
		self.name = name #string
	def getName(self):
		return self.name #string
		
	def setOrientation(self, orientation):
		self.Orientation = orientation #string
	def getOrientation(self):
		return self.orientation #string
		
	def setDNA(self, dna):
		self.dna = dna #string
	def getDNA(self):
		return self.dna #string

	def setQualVal(self, qual_val):
		self.qual_val = qual_val #string
	def getQualVal(self):
		return self.qual_val #string	

	def setTrace(self, trace):
		self.trace = trace #list [G, A, T, C]
	def getTrace(self):
		return self.trace #list	[G, A, T, C]

	def setDNAClipped(self, dna):
		self.dna_clipped = dna #string
	def getDNAClipped(self):
		return self.dna_clipped #string		
	
	def setReference(self, dna):
		self.reference = dna #string
	def getReference(self):
		return self.refrence #string		

	def setSignMissmatch(self, string):
		self.sign_missmatches = string #string
	def getSignMissmatch(self):
		return self.sign_missmatches #string	

	def setDNAMissmatch(self, string):
		self.dna_missmatches = string #string
	def getDNAMissmatch(self):
		return self.dna_missmatches #string

	def setAAMissmatch(self, string):
		self.aa_missmatches = string #string
	def getAAMissmatch(self):
		return self.aa_missmatches #string

	def setContig(self, dna):
		self.contig = dna #string
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

	
	def getInput(self):
		'''Open a single .seq, .fasta, .fastq, .ztr, .scf, .ab1 file (or even a text file with a DNA sequence) and set variables accordingly.'''
		parts = self.filepath.split('/')
		filename = parts.pop() #get filename
		path = '/'.join(parts)+'/' #path to file
		
		#establish type of input file
		if '.' in filename: 
			self.input_type = filename.split('.')[-1].upper() 
			print('type', self.input_type)
		else:
			self.input_type = None
		
	
		#read the input
		if self.input_type in ['TXT', 'SEQ', None] and filename not in ['allseqs.txt']:
			self.setName(filename)
			f = open(self.filepath, 'r') 
			dna = f.read() 
			self.setDNA(dna.replace('\n', ''))
			#add an assert that there are only dna bases here
			f.close()
			#establish orientation of DNA
			if filename.split('.')[0][-2:].upper() == 'FW':
				self.setOrientation('fw')
			elif filename.split('.')[0][-2:].upper() == 'RV':
				self.setOrientation('rv')
			else:
				raise TypeError, 'The last two characters of the filename (before the .) must specify whether the sequence is fw or rv. Pleace rename file %s accordingly' % filename

			
		elif self.input_type in ['AB1', 'ABI', 'ABIF']:
			self.setName(filename)
			ab1 = ABIreader.Trace(self.filepath, trimming=False) #optionally ', trimming=True'
			self.setDNA(ab1.seq)
			self.setQualVal(ab1.qual_val)
			self.setTrace([ab1.data['raw1'], ab1.data['raw2'], ab1.data['raw3'], ab1.data['raw4']])
			#abi=dict(baseorder=ab1.data['baseorder'], qual_val=ab1.qual_val, G=str(AB1Trace.data['raw1']), A=str(AB1Trace.data['raw2']), T=str(AB1Trace.data['raw3']), C=str(AB1Trace.data['raw4']))
			#establish orientation of DNA
			if filename.split('.')[0][-2:].upper() == 'FW':
				self.setOrientation('fw')
			elif filename.split('.')[0][-2:].upper() == 'RV':
				self.setOrientation('rv')
			else:
				raise TypeError, 'The last two characters of the filename (before the .) must specify whether the sequence is fw or rv. Pleace rename file %s accordingly' % filename

#		elif self.input_type == 'ZTR':
#			print('Support for .ztr files has not yet been implemented')
			
#		elif self.input_type == 'SCF':
#			print('Support for .scf files has not yet been implemented')
			
		elif fnmatch.fnmatch(filename, '*.fasta'):
			self.setName(filename)
			id, dna = fasta.parse(self.filepath) #parse the fasta file. File should contain ONE entry
			self.setDNA(dna)
			#establish orientation of DNA
			if filename.split('.')[0][-2:].upper() == 'FW':
				self.setOrientation('fw')
			elif filename.split('.')[0][-2:].upper() == 'RV':
				self.setOrientation('rv')
			else:
				raise TypeError, 'The last two characters of the filename (before the .) must specify whether the sequence is fw or rv. Pleace rename file %s accordingly' % filename
			
		elif fnmatch.fnmatch(filename, '*.fastq'):
			self.setName(filename)
			id, dna, id2, qual_val = fastq.parse(self.filepath) #parse the fastq file. File should contain ONE entry
			self.setDNA(dna)
			self.setQualVal(qual_val)
			#establish orientation of DNA
			if filename.split('.')[0][-2:].upper() == 'FW':
				self.setOrientation('fw')
			elif filename.split('.')[0][-2:].upper() == 'RV':
				self.setOrientation('rv')
			else:
				raise TypeError, 'The last two characters of the filename (before the .) must specify whether the sequence is fw or rv. Pleace rename file %s accordingly' % filename
			
		else:
			pass
			print('"%s" is not a .txt, .seq, .scf, .fasta, .fastq, .abif, .ab1, .abi or .ztr file' % filename)

			
			


class SeqAnalysis:
	'''
	Sequencing analysis object for assembling and analysing data from many Sanger sequencing runs.
	'''
	def __init__(self):
		self.seqdata = [] #this is where all the sequence info gets stored

		
	def addSeq(self, seqobject):
		'''
		Adds a single SeqObj sequencing object to the data structure.
		'''
		if len(self.seqdata) == 0:
			self.seqdata.append(dict(samples = [seqobject]))

		elif len(self.seqdata) > 0:

			#this workaround works for construct names with _ in it. But it fails when the primer name has _ in it....
			### this is overly complicated to deal with filenames that has more than one '_' in the filename
	#		nameparts = dictionary['name'].split('_')
	#		nameparts.pop() #remove last
	#		name = ''
	#		for item in nameparts: #stitch them together
	#			name += '%s_' % item
			###	
			
			#construct name and primer name MUST be separated by _ . Construct name CANNOT contain _ .
			name = seqobject.getName().split('_')[0] #get first part of name
			
			for i in range(len(self.seqdata)):
				if name in self.seqdata[i]['samples'][0].getName():
					self.seqdata[i]['samples'].append(seqobject)
					break
				elif i == len(self.seqdata)-1:
					self.seqdata.append(dict(samples = [seqobject]))

		else:
			print('Error adding sequence to self.seqdata')	

				
	def addFolder(self, path):
		'''Opens many .Seq or .ab1 file and return id and sequence into a list of dictionaries'''
		#Here I move all the sequence data from all .Seq or .ab1 files in one folder into a list of dictionaries
		file_list = sorted(os.listdir(path))
		for filename in file_list: 
			print(filename)

			#check if there is an ab1 file for this name in the folder. I will also add support for .fasta, .scf and .fastq files
			#I need to change this to an internal ranking. For each name it should just pick the 'best' file
			if filename.split('.')[-1].upper() == 'SEQ':
				if [filename.split('.')[0]+'.ab1' in s for s in file_list] is False: #if it's not present, add the .seq file
					seqobject = SeqObj(path+'/'+filename) #make a SeqObj object
					print('if if', seqobject)
			else:
				seqobject = SeqObj(path+'/'+filename)	#make a SeqObj object

			if seqobject.getDNA() is not False and seqobject.getName() is not False: #a seqobject is only made if the input file is the appropriate file type
				self.addSeq(seqobject)
		

	def getOverlap(self, seq1, seq2):
		'''
		Find and return the maximum overlap length of two sequences.
		This only works when using the NeedlemanWunsch algorithm to get the overlap. (The algorithm is much too promiscous otherwise.
		If there is no overlap, return empty string.
		'''
		alignment = NW.PairwiseAlignment(seq1, seq2)
		seq1aln = alignment.seq1aligned.upper()
		seq2aln = alignment.seq2aligned.upper()
		assert len(seq1aln) == len(seq2aln), 'Error, the sequences are not of the same length'
		
				
		overlap = False
		start = False
		end = False
		length = False
	
		#get start, end and length of overlap. N is tolerated and double -- are tolerated.
		for i in range(len(seq1aln)):
			if overlap is False and seq1aln[i] in 'ATCGN' and seq2aln[i] in 'ATCGN': #first nucleotide of overlap
				start = copy.copy(i)
				overlap = True
				
			elif overlap is True and (seq1aln[i] == '-' or seq2aln[i] == '-'):
				if seq1aln[i+1:i+3] == '--' or seq2aln[i+1:i+3] == '--': #allow for two missing bases, but not more
					end = copy.copy(i)
					length = end-start
					break
					
		if overlap is False:
			return 0
		else:
			return length #make this to look like the difflib

			
	def sortSeqs(self):
		'''
		Align each sequence with each other sequence. Sort them in the order in which they should physically occur.
		'''
		for group in seqdata:
			alnmatrix = []
			for seqobj1 in group['samples']:
				for seqobj2 in group['samples']:
					if seqobj1 == seqobj2:
						pass
					else:
						temp = NW.PairwiseAlignment(seqobj1.getDNA(), seqobj2.getDNA()) #align the two
						#somehow score the alignment


		# sortlist = [0]

		# for i in range(1,len(entry)):
			# counter = 0
			

			# for i in range(len(temp)):
				# if temp[i] == '-':
					# counter += 1
				# elif temp[i] == 'x':
					# counter = 'x'	
				# elif temp[i] != '-':
					# break
			# sortlist.append(counter)
		
		# i = len(sortlist)
		# n = 0
		# while n < i:
			# if sortlist[n] == 'x':
				# del entry[n]
				# del sortlist[n]
				# i -= 1
				# n -= 1
			# else:
				# n += 1
		# sortedseqs = [x for (y,x) in sorted(zip(sortlist, entry))]
		# return sortedseqs

				
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


def findoverlap(Seq1, Seq2, min_overlap=20):
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


		
		
###function that aligns everything in fasta file
def M_align(dictlist):
	'''Muscle alignment, takes a list of dictionaries as input and aligns all at once'''
	#prepare input as a 'virtual' FASTA file
		
	records = ''
	for entry in dictlist:
		if entry['name'][0] != '>':
			records += '>%s\n%s\n' % (entry['name'], entry['dna'])
		elif entry['name'][0] == '>':
			records += '%s\n%s\n' % (entry['name'], entry['dna'])
		else:
			print('Muscle name error')

	records_handle = StringIO(records) #turn string into a handle
	tempdata = records_handle.getvalue()
	
	#for seperate fasta entries
	muscle_cline = MuscleCommandline()
	stdout, stderr = muscle_cline(stdin=tempdata)
	stdout = parse_fasta(stdout)

#	#for aligned fasta entries
#	muscle_cline = MuscleCommandline(clw=True)
#	stdout, stderr = muscle_cline(stdin=tempdata)
#	print(stdout)
	
#	#the clustalw-type output can be further formated
#	align = AlignIO.read(StringIO(stdout), "clustal")
#	print(align)

	return stdout

	
	
def M_align_iterative(index):
	'''This function aligns many sequences to a reference, one sequence at a time, then makes sure that all gaps match'''
	#prepare input as a 'virtual' FASTA file
	#first entry should be the reference

	temp_dictlist = []
	for n in range(len(seqdata[index]['samples'])): #make sure each dna has an empty aln_dna entry
		seqdata[index]['samples'][n]['aln_dna'] = '' 

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

				
			
	#I should probably add someting to go through and check that I don't have '-' for all sequences at some position


def parse_fasta(string):
	'''Convert fasta sting to a list of dictionaries'''
	sequence = ''
	name = ''
	output = []
	list = string.split('\n')
	for line in list:
		if line == '':
			pass
		elif line == '\n':
			pass
		elif line[0] == '>':
			if sequence != '':
				output.append(dict(name=str(name), dna=sequence))
			sequence = ''
			name = line
		else:
			sequence += line
	output.append(dict(name=str(name), dna=sequence))
	return output		

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
					print('DNA pos %d, %s mutated to %s --- %s%d%s' % (i+3-1, RefSeq[i:(i+3)], Seq[i:(i+3)], DNA.translate(RefSeq[i:(i+3)]), int((i+3)/3), DNA.translate(Seq[i:(i+3)])))
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
	'''grups sequences based on name'''
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



#write a function for de-novo assembly without refseq
def assemble():
	'''Assembles all of the sequences in every group (if possible)'''
	
	#take one seq, find the best match, combine them. Start over.
	for i in range(len(seqdata)): #for all of the sequence groups
		fw_list = []
		rv_list = []
		seqdata[i]['contig']['dna'] = ''
		samples = copy.deepcopy(seqdata[i]['samples']) #copy samples so that I can delete them as I go
		overlap_found = False #has overlap been found?
		whole_list = False #did we traverse the whole list?
		while overlap_found is False and whole_list is False:
			overlap_found = False #has overlap been found?
			whole_list = False #did we traverse the whole list?
			for n in range(len(samples)):
				for o in range(n, len(samples)):
					fw_dict = {}
					seq1 = seqdata[i]['samples'][0]['dna']
					seq2 = seqdata[i]['samples'][n]['dna']
					forward = findoverlap(seq1, seq2)
					if forward is not False:
						overlap_found = True
					fw_list.append(forward[2])
					
					seq1 = seqdata[i]['samples'][0]['dna']
					seq2 = seqdata[i]['samples'][n]['dna']
					reverse = findoverlap(seq1, DNA.reversecomplement(seq2))
					if reverse is not False:
						overlap_found = True
					rv_list.append(forward[2])
				if overlap_found is True:
					if max(fw_list) >= max(rv_list):
						index, value = max(enumerate(fw_list))
					elif max(fw_list) < max(rv_list):
						index, value = max(enumerate(fw_list))
					
						
			whole_list = True
			
			
	
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
	FindSingleGeneMissmatch = 'Y' #Test for one gene whether missmatches occur or not. It is expected that RefSeq starts with ATG
	
	#this grabs all sequences from a folder, removes N, joins them, and checks for mutations
	all_files_from_path(path) #to get all files from a folder
	writeseqfiles(path) #writes sequences to allseqs
	# pprint.pprint(seqdata)	
	# #for making alignment
	# for i in range(len(seqdata)): #for all entries (groups of sequences)
		# getreference(path, i)
		# M_align_iterative(i)
# #	print(seqdata[0]['samples'][0]['aln_dna'])
	
	# #call sortseqs here for when I have more sequences

	# analyze()


