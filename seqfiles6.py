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
import itertools
import string


from Bio import AlignIO
from Bio.Align.Applications import MuscleCommandline

import ABIreader #parser
#import scf #parser
#import abi #parser
#import ztr #parser
import fastq #parser
import fasta #parser
import NeedlemanWunsch as NW #alignment





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
		
		self.aln_dna = None #store dna sequence aligned to reference

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
		
	def setAlignment(self, seq):
		self.aln_dna = seq
	def getAlignment(self):
		return self.aln_dna
		
	def setRC(self, bool):
		self.RC = bool #true or false
	def getRC(self):
		return self.RC #true or false
		
	def setDNAClipped(self, seq):
		self.seq_clipped = seq #string
	def getDNAClipped(self):
		return self.seq_clipped #string		

	def RCObject(self):
		'''
		Reverse-complement the entire object.
		'''
		#this will need a lot of work to make it work right
		
		#update variable to keep track of whether it has been reverse-complemented or not
		if self.getRC() is True:
			self.setRC(False)
		elif self.getRC() is False:
			self.setRC(True)
			
		#update trace
		#how?
		
		#update qual val
		#how?
		
		#update sequence
		self.setDNA(dna.RC(self.getDNA()))
	
	
	def getInput(self):
		'''
		Open a single .seq, .fasta, .fastq, .ztr, .scf, .ab1 file (or even a text file with a DNA sequence) and set variables accordingly.
		'''
	
		#read the input
		if self.input_type in ['TXT', None] and self.filename not in ['allseqs.txt']:
			f = open(self.filepath, 'r') 
			input = f.read() 
			f.close()
			self.setDNA(input.replace('\n', ''))
			
		elif self.input_type in ['SEQ', 'SEQ.CLIPPED']:
			output = fasta.parseFile(self.filepath) #parse the fasta file. File should contain ONE entry
			for item in output:
				id, seq = item
			self.setDNA(seq)
			
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
		
		self.reference = None #reference DNA
		self.aligned = False #keep track of whether sequences have been aligned or not
		
		self.consensus = None #consensus sequence
		
		self.aln_dna = None #sequence of reference DNA aligned with all other sequenes
		self.sign_missmatches = None #sign (* *  *   *) mismatches of sequences vs the reference or consensus
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
		'''
		Set a reference sequence.
		'''
		self.reference = seq #string
		
	def getReference(self):
		'''
		Return the reference sequence.
		'''
		return self.reference #string		

	def setName(name):
		self.name = name
	def getName(self):
		return self.name
	
	def setAlignment(self, seq):
		self.aln_dna = seq
	def getAlignment(self):
		return self.aln_dna
	
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
	
	def setConsensus(self):
		'''
		Build a consensus sequence from sorted sequence reads.
		'''
		pass
					
		
		
		
	######################################################
	def parse_fasta(self, string):
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

	
	def assemble(self):
		print('Assembling %s' % self.getName())
		if self.getReference() is None: #if no reference sequence is present, do de Novo assembly.
			self.assemble_de_novo()
		else: 
			self.assemble_to_ref() #if a reference sequence is present, align each sequence to it.
	
		#check missmatches and set the associated variables

		self.aligned = True

	def assemble_de_novo(self, seq_type='SEQ.CLIPPED'):
		'''
		This method aims to achieve the de-novo assembly of sequence reads without using a reference sequence.
		The sequence reads are arranged from longest to shortest.
		The longest sequence is used as a starting point of the contig.
		The method then takes the contig and find which other sequence generates the largest overlap.
		This sequence is added to the growing contig.
		This is repeated until no other overlaps are found or untill all sequences are aligned.
		
		#it works pretty well but needs some further improvement
		#mainly the entire sequence objects need to be reverse-complemented for the reverse, not just the sequence.
		#the objects should also be stored in a list for later use.
		
		seq_type defines what type of sequence data should be used (AB1, SFC, SEQ, FASTQ...). There may be several sources for a single sequence.
		
		#I should improve the performance of this algorithm and then convert it to cython...
		'''
		#sort the sequence objects based on sequence length, longest first, shortest last
		seq_list = [self.seqdata[key] for key in self.seqdata.keys()]
		seq_list = sorted(seq_list, key=lambda k: k[seq_type].getDNA(), reverse=True) 
		
		#set up contig
		order = [seq_list[0]] #stores the sequence objects in order, not currently used.....
		contig = seq_list[0][seq_type].getDNA() #add the longest sequence
		used_index = [0]
		
		#take the longest seq, find the best match, combine them. Start over. Do until all sequences have been assembled.
		for n in range(len(seq_list)-1): #do it as many times as the list is long -1
			best_score = 0 #keep track of best score
			best = None #keep the best (with the best alignment score) sequence object
			index_of_best = None
			for o in range(len(seq_list)):
				if o in used_index:
					pass
				else:
					#get sequences
					seq1 = contig
					seq2 = seq_list[o][seq_type].getDNA()
				
					#align sequences in FW
					forward = self.getOverlap(seq1, seq2)
					
					#align sequences with one in RV
					reverse = self.getOverlap(seq1, dna.RC(seq2))
					
					#was the fw or rv alignment best?
					threshold_score = 10
#					print('FW score: %s' % forward['score'])
#					print('RV score: %s' % reverse['score'])
					if forward['score'] >= threshold_score > reverse['score']: #only the fw is above threshold
						if forward['score'] > best_score:
							best_score = forward['score']
							best = forward
							index_of_best = o
							
					elif forward['score'] < threshold_score <= reverse['score']: #only the rv is above threshold
						if reverse['score'] > best_score:
							best_score = reverse['score']
							best = reverse
							index_of_best = o
							
					elif forward['score'] > reverse['score'] >= threshold_score: #both are above, but fw is largest
						if forward['score'] > best_score:
							best_score = forward['score']
							best = forward
							index_of_best = o
						
					elif threshold_score <= forward['score'] < reverse['score']: #both are above, but rv is largest
						if reverse['score'] > best_score:
							best_score = reverse['score']
							best = reverse
							index_of_best = o
						
					elif forward['score'] < threshold_score and reverse['score'] < threshold_score: #both are below the threshold
						pass						
						
					elif forward['score'] == reverse['score'] >= threshold_score: #they can't both be above the threshold and be equal, because only one alignment is biologically correct
						print('Both FW and RV have the same score.')
						raise ValueError
						
					else: 
						print('Unanticipated scores for FW and RV.')
						raise ValueError
					
					#I should add the object to the "order" list. Reverse-complemented if necessary.
					#RCObject
				
			#make sure that something actually aligned
			if best == None:
				print('No further alignments possible even though sequences are left.')
				self.setContig(contig)
				break
				
			#now take the "best" sequence and add it to the contig
			bases = []
			for i,j in itertools.izip_longest(best['seq1'],best['seq2']):
				if i == j and i.upper() in 'GATCN':
					bases.append(i)
				elif i.upper() in 'GATC' and j.upper() in '-N':
					bases.append(i)
				elif i.upper() in '-N' and j.upper() in 'GATC':
					bases.append(j)
				elif i != j and i.upper() in 'GATC' and j.upper() in 'GATC': #both are GATC but they are not the same
					bases.append(i) #just pick one. I will change it later so that basecalls are used to pick the best one.
				else:
					print('unanticipated base value: %s, %s' % (i,j))
					raise ValueError
			contig = ''.join(bases)
			
			#add that sequence to the used index
			used_index.append(index_of_best)
	
		#add the finished contig
		self.setContig(contig)

		
	def assemble_to_ref(self, seq_type='SEQ.CLIPPED'):
		'''
		This method aligns many sequences to a reference, one sequence at a time, then makes sure that all gaps match.
		
		The method assumes that all sequences are in the FW direction.
		'''
		
		#need to sort the sequences first!!! They all need to be in the FW direction. Fix!
		


		#set dna as aligned sequence. Gaps will be introduced later.
		self.setAlignment(self.getReference())
		
		#make virtual FASTA file of reference and a sequence file and align. The first entry is always the reference.
		for key in self.seqdata.keys(): 
			records = '>%s_ref\n%s\n' % (self.getName().lstrip('>'), self.getReference()) #reference
			records += '>%s\n%s\n' % (self.seqdata[key][seq_type].getName().lstrip('>'), self.seqdata[key][seq_type].getDNA()) #sample

			#turn string into a handle
			records_handle = StringIO(records) 
			tempdata = records_handle.getvalue()
		
			#do the alignment using Muscle
			muscle_cline = MuscleCommandline()
			stdout, stderr = muscle_cline(stdin=tempdata)
			stdout = self.parse_fasta(stdout)	

			#set the aligned sequence
			self.seqdata[key][seq_type].setAlignment(stdout[1]['dna'])
			
			#set the aligned reference
			current_ref_aln = stdout[0]['dna']
			
			#Now compare the current aligned reference with the one stored in the data structure.
			#If there are gap differences, then update all stored alignments. 
			#It is enough to compare the reference sequences from each pairwise alignment to get the info needed to match all gaps.
			current_pos = 0
			while current_ref_aln.upper() != self.getAlignment().upper(): #as long as the entire sequences don't match
				if len(current_ref_aln) < current_pos-1: #if new alignment is too short, add - to the end
					print('adding to new')
					#update the new reference
					current_ref_aln += '-'
					#update the last aligned sequence
					new_seq = self.seqdata[key][seq_type].getAlignment() + '-'
					self.seqdata[key][seq_type].setAlignment(new_seq)
				
				elif len(self.getAlignment()) < current_pos-1: #if old alignment is too short, add - to the end
					print('adding to old')
					#update the old reference
					self.setAlignment(self.getAlignment() + '-' )
					#update all aligned sequences
					for key in self.seqdata.keys():
						if self.seqdata[key][seq_type].getAlignment() == None:
							pass
						elif len(self.seqdata[key][seq_type].getAlignment()) < current_pos-1:
							new_seq = self.seqdata[key][seq_type].getAlignment() + '-'
							self.seqdata[key][seq_type].setAlignment(new_seq)
					
				elif current_ref_aln[current_pos].upper() == self.getAlignment()[current_pos].upper(): #they match at this position
					current_pos += 1

				else: #they do not match at this position, so fix it.
					if self.getAlignment()[current_pos].upper() == '-' and current_ref_aln[current_pos].upper() != '-': # if gap in old alignment, but not in new
						#update the new reference
						new_seq = current_ref_aln[:current_pos] + '-' + current_ref_aln[current_pos:]
						current_ref_aln = new_seq
						#update the last aligned sequence
						new_seq = self.seqdata[key][seq_type].getAlignment()[:current_pos] + '-' + self.seqdata[key][seq_type].getAlignment()[current_pos:]
						self.seqdata[key][seq_type].setAlignment(new_seq)

					elif self.getAlignment()[current_pos].upper() != '-' and current_ref_aln[current_pos].upper() == '-': #if gap in new alignment, but not in old
						#update the old reference
						new_seq = self.getAlignment()[:current_pos] + '-' + self.getAlignment()[current_pos:]
						self.setAlignment(new_seq)
						#update all aligned sequences
						for key in self.seqdata.keys():
							if self.seqdata[key][seq_type].getAlignment() == None:
								pass
							elif self.seqdata[key][seq_type].getAlignment()[current_pos] != '-':
								new_seq = self.seqdata[key][seq_type].getAlignment()[:current_pos] + '-' + self.seqdata[key][seq_type].getAlignment()[current_pos:]
								self.seqdata[key][seq_type].setAlignment(new_seq)

					else:
						raise ValueError
				
			#make contig here....	


	def getOverlap(self, seq1, seq2):
		'''
		Align two sequences using the NeedlemanWunsch algorithm. 
		Return a dictionary of the aligned sequences and the alignment score.
		'''
		alignment = NW.PairwiseAlignment(seq1, seq2)
		score = alignment.score
		seq1aln = alignment.seq1aligned.upper()
		seq2aln = alignment.seq2aligned.upper()
		assert len(seq1aln) == len(seq2aln), 'Error, the sequences are not of the same length'

		return {'seq1':seq1aln, 'seq2':seq2aln, 'score':score} 

	def generate_output(self, path):
		with open('%s/%s_output.txt' % (path, self.getName()), 'w') as f:
			f.write('>%s\n' % self.getName())
			f.write(self.getContig())


	
			
################ Generating outputs ######################

#	def writeseqfiles(self, path):
#		'''Write all sequence entries from a list of dictionaries to textfile'''
#		a = open(path + 'allseqs.txt', 'w') #open it for writing
#		for i in range(len(self.seqdata)):
#			for n in range(len(self.seqdata[i]['samples'])):
#				a.write(self.seqdata[i]['samples'][n].getName())	
#				a.write('\n')
#				a.write(self.seqdata[i]['samples'][n].getDNA())
#				a.write('\n')
#		a.close()


##########################################################

			

class SeqOverview:
	'''
	The SeqOverview object manages one or several instances of SeqAnalysis and provides means to interact with these.
	It should be possible to load a whole folder of files into this class and have it organize the sequences by name and make SeqAnalysis instances as appropriate.
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

	
	def assembleAll(self):
		'''
		Assemble all 
		'''
		for key in x.data.keys():
			x.data[key].assemble()
			
	def outputAll(self, path):
		'''
		Make text output of the contigs.
		'''
		for key in x.data.keys():
			print(key)
			x.data[key].generate_output(path)

			
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
	x.assembleAll()
	x.outputAll(path)
	


	
	
	

