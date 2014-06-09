#!/usr/bin/python


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


def R(dna):
	"""Returns the reverse of a DNA string"""
	dna = dna.replace('\n','')
	return dna[::-1]  #makes the reverse of the input string
		
def C(dna):
	"""Returns the complement of a DNA string"""
	dna = dna.replace('\n','')
	compdna = ''
	for i in range(len(dna)):    #this loop makes the complement of a string
		if dna[i] == 'a':
			base = 't'
			
		elif dna[i] == 't':
			base = 'a'

		elif dna[i] == 'c':
			base = 'g'

		elif dna[i] == 'g':
			base = 'c'

		elif dna[i] == 'n':
			base = 'n'
			
		elif dna[i] == 'A':
			base = 'T'

		elif dna[i] == 'T':
			base = 'A'

		elif dna[i] == 'C':
			base = 'G'

		elif dna[i] == 'G':
			base = 'C'

		elif dna[i] == 'N':
			base = 'N'
			
		elif dna[i] == '\n': #\n is the end of line character found in many files
			base = ''	

		else:
			base = 'x'
		compdna = compdna + base
	return compdna
	
	
def RC(dna):
	"""Returns the reverse complement of a DNA string"""
	return R(C(dna))


def Translate(dna):
	"""Returns protein sequence from DNA string input"""
	F = ['TTT', 'TTC'];
	L = ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'];
	S = ['TCT', 'TCC', 'TCA', 'TCG', 'AGC', 'AGT'];
	Y = ['TAT', 'TAC'];
	stop = ['TAA', 'TAG', 'TGA'];
	C = ['TGT', 'TGC'];
	W = ['TGG'];
	P = ['CCT', 'CCA', 'CCG', 'CCC'];
	H = ['CAT', 'CAC'];
	E = ['GAA', 'GAG'];
	R = ['CGT', 'CGA', 'CGG', 'CGC', 'AGG', 'AGA'];
	I = ['ATT', 'ATC', 'ATA'];
	M = ['ATG'];
	T = ['ACG', 'ACC', 'ACA', 'ACT'];
	N = ['AAT', 'AAC'];
	K = ['AAA', 'AAG'];
	V = ['GTT', 'GTA', 'GTG', 'GTC'];
	A = ['GCG', 'GCC', 'GCA', 'GCT'];
	D = ['GAT', 'GAC'];
	Q = ['CAG', 'CAA'];
	G = ['GGA', 'GGT', 'GGC', 'GGG']

	protein = ''
	dna = dna.upper()
	dna = dna.replace('\n', '') #get rid of row breaks
	dna = dna.replace('U', 'T')

	for i in range(len(dna)):
		if i%3==0:
			if i+3>len(dna):
				pass
			elif any(dna[i:(i+3)] in s for s in F):
				protein = protein + 'F'
			elif any(dna[i:(i+3)] in s for s in L):
				protein = protein + 'L'
			elif any(dna[i:(i+3)] in s for s in S):
				protein = protein + 'S'
			elif any(dna[i:(i+3)] in s for s in Y):
				protein = protein + 'Y'
			elif any(dna[i:(i+3)] in s for s in stop):
				protein = protein + '*'
			elif any(dna[i:(i+3)] in s for s in C):
				protein = protein + 'C'
			elif any(dna[i:(i+3)] in s for s in W):
				protein = protein + 'W'
			elif any(dna[i:(i+3)] in s for s in P):
				protein = protein + 'P'
			elif any(dna[i:(i+3)] in s for s in H):
				protein = protein + 'H'
			elif any(dna[i:(i+3)] in s for s in E):
				protein = protein + 'E'
			elif any(dna[i:(i+3)] in s for s in R):
				protein = protein + 'R'
			elif any(dna[i:(i+3)] in s for s in I):
				protein = protein + 'I'
			elif any(dna[i:(i+3)] in s for s in M):
				protein = protein + 'M'
			elif any(dna[i:(i+3)] in s for s in T):
				protein = protein + 'T'
			elif any(dna[i:(i+3)] in s for s in N):
				protein = protein + 'N'
			elif any(dna[i:(i+3)] in s for s in K):
				protein = protein + 'K'
			elif any(dna[i:(i+3)] in s for s in V):
				protein = protein + 'V'
			elif any(dna[i:(i+3)] in s for s in A):
				protein = protein + 'A'
			elif any(dna[i:(i+3)] in s for s in D):
				protein = protein + 'D'
			elif any(dna[i:(i+3)] in s for s in Q):
				protein = protein + 'Q'
			elif any(dna[i:(i+3)] in s for s in G):
				protein = protein + 'G'
			else:
				protein = protein + '?'
	return protein	

def TranslateRC(dna):
	'''Translate the reverse complement of DNA'''
	dna = RC(dna)
	return Translate(dna)

#add randomize dna


#add abireader function

#add function for fetching dna from uniprot, ncbi...




def pair_ident(Seq1, Seq2, single_gaps):
	'''Takes two aligned sequences and returns their percent identity.
Assumes that Seq1 and Seq2 are sequence strings'''

	l=0.0 # counts alignment length, excluding identical gaps, but including single gaps
	n=0.0 # count number of single gaps
	i=0.0 # counts identity hits


	for j in range(len(Seq2)): 
		if Seq1[j] == '-' and Seq2[j] == '-': #DON'T count identical gaps towards alignment length
			pass
		else:
			if Seq2[j] == Seq1[j]:
				i += 1 #count matches
			elif Seq2[j] == '-' or Seq1[j] == '-':
				n += 1 #count number of single gaps
			l += 1 #count total length with single gaps

	if single_gaps == 'y': #include single gaps 
		percent = round(100*(i/l),1) #calculate identity
	elif single_gaps == 'n': #exclude single gaps
		if n == l: #if gaps same as total length identity is 0
			percent = 0.0
		else:	
			percent = round(100*(i/(l-n)),1) #calculate identity
	return percent




#### Analyze alignments ####
def single_ident(align_file, single_gaps):
	'''Get identities for all protines compared to one reference sequence in aligned FASTA file.
Assumes that align_file is a fasta-formatted alignment.
Assumes that reference sequence is the first one.
Returns a touple of the results.'''

	from Bio import AlignIO
	results = ()

	input_handle = open(align_file, "rU")
	alignment = AlignIO.read(input_handle, "fasta")
	input_handle.close()
	for record2 in range(1, len(alignment)):
		percent = pair_ident(alignment[0].seq, alignment[record2].seq, single_gaps)
		results += (alignment[0].id + ' vs ' + alignment[record2].id+' '+str(percent),)
	return results

def all_ident(align_file, single_gaps):
	'''Get identities for all combinations of proteins in aligned FASTA file.
Assumes that align_file is a fasta-formatted alignment.
Returns a touple of the results.'''

	from Bio import AlignIO
	results = ()

	input_handle = open(align_file, "rU")
	alignment = AlignIO.read(input_handle, "fasta")
	input_handle.close()
	for record in range(len(alignment)):
		for record2 in range(record+1, len(alignment)):
			percent = pair_ident(alignment[record].seq, alignment[record2].seq, single_gaps)
			results += (alignment[record].id + ' vs ' + alignment[record2].id+' '+str(percent),)
	return results

########


def blast(blast_type, database, seq): #function for blasting 
	'''Blasts a sequence against the NCBI database'''
	#"tblastn", "nt", seq
	#add asserts...
	from Bio.Blast import NCBIWWW
	from Bio.Blast import NCBIXML
	from copy import deepcopy

	result_handle = NCBIWWW.qblast(blast_type, database, seq) #perform the BLAST
	blast_record = NCBIXML.read(result_handle)
	string = ''

	for alignment in blast_record.alignments: #print all seq IDs and sequences in the fasta file
		for hsp in alignment.hsps:
			subject = hsp.sbjct[0:]
			nodash = subject.replace("-", "")
			string += ">"+alignment.title+"\n"+nodash+"\n"
	return string
	result_handle.close()
#	save_file = open(protname+"_my_blast_nt.xml", "w")
#	save_file.write(result_handle.read())
#	save_file.close()
	


def smithwaterman(dna1, dna2):
	"""Module aligns two DNA sequences and returns aligned seqA, alignment "bars (|)" and aligned seqB"""
	##got script from Kevin Kwok and modified...
	##http://kevinakwok.tumblr.com/post/18372160357/smith-waterman-algorithm
	seqA = dna1
	seqB = dna2
	seqB = seqB.upper()
	seqA = seqA.upper()
	seqB = seqB[:-1]
	seqA = "-"+seqA
	seqB = "-"+seqB
	
	#Sets the values for match, mismatch, and gap.
	match = 1
	mismatch = 0
	gap = -1

	row = len(seqA)
	col = len(seqB)

	#Function to print out matrix. Used for testing.
	def print_matrix():
		for i in range(row):
			print("{}".format(A[i]))

	#Creates blank matrix
	def create_matrix(row,col):
		A = [0] * row
		for i in range(row):
			A[i] = [0] * col
		return A
	#isMatch
	def isMatch(i,j):
		if seqA[i] == seqB[j]:
			matchVal = match
		else:
			matchVal = mismatch
		return matchVal

	#Returns the new value if diagonal is used
	def diag(i,j):
		return A[i-1][j-1] + isMatch(i,j)

	#Returns the new value if up is used
	def up(i,j):
		return A[i-1][j] + isMatch(i,j) + gap

	#Returns the new value if left is used
	def left(i,j):
		return A[i][j-1] + isMatch(i,j) + gap

	#Fills matrix with correct scores.
	def complete_matrix(row,col):
		for i in range(1,row):
			for j in range(1,col):
				A[i][j] = max(0,diag(i,j),up(i,j),left(i,j)) #here the scoring matrix is generated
		    				#it seems that maybe using the different scoring methods is a bad idea
		    				#in the matrix there are two identical max scores that represent
		    				#the same alignment
		return A

	#FInd the highest scoring cell.
	def get_max(A): #A is a matrix of scores, this function goes through to find the max score
		local_max = [[0,0]]
		for i in range(row):
			for j in range(col):
				if A[i][j] == A[local_max[0][0]][local_max[0][1]]:
					local_max.append([i,j])
				elif A[i][j] > A[local_max[0][0]][local_max[0][1]]:
					local_max = [[i,j]]
		return local_max

	#Gives you the next location.
	def get_next(A,location):
		i = location[0]
		j = location[1]
		maxVal = max(A[i-1][j-1],A[i-1][j]+gap,A[i][j-1]+gap)
		if A[i-1][j-1] == maxVal:
			return [i-1,j-1]
		#Is this the right ordering of the three?
		elif A[i][j-1]+gap == maxVal:
			return [i,j-1]
		else:
			return [i-1,j]

	#Traces the path back given starting location
	def trace_back(A,tracer):
		if A[tracer[len(tracer)-1][0]][tracer[len(tracer)-1][1]] == 0:
			return tracer
		next_cell = get_next(A,tracer[len(tracer)-1])
		#tracer.insert(0,next_cell)
		tracer.append(next_cell)
		return trace_back(A,tracer)

	#Uses tracer to return final sequence
	def get_seq(A,tracer,k,seq):
		if k == 0:
			original_sequence = seqA
		else:
			original_sequence = seqB
		N = len(tracer)
		for i in range(0,N-1):
			if tracer[i][k] == tracer[i+1][k]+1:
				seq = seq + original_sequence[tracer[i][k]]
			elif tracer[i][k] == tracer[i+1][k]:
				seq = seq + "-"
		return seq

	#Shows the relevant lines for matching pairs
	def get_middle(finalA,finalB):
		middle = ""
		for k in range(len(finalA)):
			mid = " "
			if finalA[k] == finalB[k]:
				mid = "|"
			middle = middle + mid
		return middle

	A = create_matrix(row,col)
	A = complete_matrix(row,col)
	num_answers = len(get_max(A))

	for i in range(num_answers):
		tracer = trace_back(A,[get_max(A)[i]])
		finalA = get_seq(A,tracer,0,"")
		for n in range(0 , i): #I had to add this because it was giving me two answers with the same alignment
			tracerOld = trace_back(A,[get_max(A)[n]])
			finalAOld = get_seq(A,tracerOld,0,"")
			if n == i:
				pass
			elif finalAOld == finalA:
				pass
	    		
			else:
				finalA = finalA.replace('\n', '')
				finalB = get_seq(A,tracer,1,"")
				finalA = finalA[::-1]
				finalB = finalB[::-1]
				middle = get_middle(finalA,finalB)
				return (finalA,middle,finalB)
				break

