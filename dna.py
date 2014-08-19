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

import fasta
import string
import random
from Bio import AlignIO #biopython package
from Bio.Align.Applications import MuscleCommandline #biopython package
from StringIO import StringIO

############# Basic DNA functions #####################

def CleanDNA(DNA, ambiguous=False, silent=True):
	'''
	Function for cleaning DNA of non-DNA characters.
	The variable ambiguous (bool) determines whether the full set of ambiguous codons are used or not.
	the variable silent determines whether an output is printed every time a non-DNA character is omitted.
	'''
	assert type(DNA) == str or type(DNA) == unicode, 'Error, input sequence must be a string or unicode'
	if ambiguous == False:
		bases = 'GATC'
	elif ambiguous == True:
		bases = 'GATCRYWSMKHBVDN'
	cleaned_seq = ''
	for char in DNA:
		if (char in string.digits) or (char in string.whitespace) or char==('/'or'\\'or ' ' or '  ') or (char.upper() not in bases):
			if silent == False:
				print('Character "%s" is not a valid DNA character and was deleted' % char)
		else:
			cleaned_seq += char
	return cleaned_seq


def R(DNA):
	"""Returns the reverse of a DNA string"""
	assert type(DNA) == str or type(DNA) == unicode, 'Error, input sequence must be a string or unicode'
	return DNA[::-1]  #makes the reverse of the input string

		
def C(DNA):
	"""Returns the complement of a DNA string"""
	assert type(DNA) == str or type(DNA) == unicode, 'Error, input sequence must be a string or unicode'
	complement = {'a':'t', 't':'a', 'c':'g', 'g':'c', 'y':'r', 'r':'y', 'w':'w', 's':'s', 'k':'m', 'm':'k', 'd':'h', 'v':'b', 'h':'d', 'b':'v', 'n':'n', 
'A':'T', 'T':'A', 'C':'G', 'G':'C', 'Y':'R', 'R':'Y', 'W':'W', 'S':'S', 'K':'M', 'M':'K', 'D':'H', 'V':'B', 'H':'D', 'B':'V', 'N':'N'}
	bases = list(DNA) 
	bases = [complement[base] for base in bases] 
	return ''.join(bases)

	
def RC(DNA):
	"""Returns the reverse complement of a DNA string"""
	assert type(DNA) == str or type(DNA) == unicode, 'Error, input sequence must be a string or unicode'
	return R(C(DNA))


def Translate(DNA):
	"""Returns protein sequence from DNA string input"""
	assert type(DNA) == str or type(DNA) == unicode, 'Error, input sequence must be a string or unicode'
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
	DNA = DNA.upper()
	DNA = DNA.replace('\n', '') #get rid of row breaks
	DNA = DNA.replace('U', 'T')

	for i in range(len(DNA)):
		if i%3==0:
			if i+3>len(DNA):
				pass
			elif any(DNA[i:(i+3)] in s for s in F):
				protein = protein + 'F'
			elif any(DNA[i:(i+3)] in s for s in L):
				protein = protein + 'L'
			elif any(DNA[i:(i+3)] in s for s in S):
				protein = protein + 'S'
			elif any(DNA[i:(i+3)] in s for s in Y):
				protein = protein + 'Y'
			elif any(DNA[i:(i+3)] in s for s in stop):
				protein = protein + '*'
			elif any(DNA[i:(i+3)] in s for s in C):
				protein = protein + 'C'
			elif any(DNA[i:(i+3)] in s for s in W):
				protein = protein + 'W'
			elif any(DNA[i:(i+3)] in s for s in P):
				protein = protein + 'P'
			elif any(DNA[i:(i+3)] in s for s in H):
				protein = protein + 'H'
			elif any(DNA[i:(i+3)] in s for s in E):
				protein = protein + 'E'
			elif any(DNA[i:(i+3)] in s for s in R):
				protein = protein + 'R'
			elif any(DNA[i:(i+3)] in s for s in I):
				protein = protein + 'I'
			elif any(DNA[i:(i+3)] in s for s in M):
				protein = protein + 'M'
			elif any(DNA[i:(i+3)] in s for s in T):
				protein = protein + 'T'
			elif any(DNA[i:(i+3)] in s for s in N):
				protein = protein + 'N'
			elif any(DNA[i:(i+3)] in s for s in K):
				protein = protein + 'K'
			elif any(DNA[i:(i+3)] in s for s in V):
				protein = protein + 'V'
			elif any(DNA[i:(i+3)] in s for s in A):
				protein = protein + 'A'
			elif any(DNA[i:(i+3)] in s for s in D):
				protein = protein + 'D'
			elif any(DNA[i:(i+3)] in s for s in Q):
				protein = protein + 'Q'
			elif any(DNA[i:(i+3)] in s for s in G):
				protein = protein + 'G'
			else:
				protein = protein + '?'
	return protein	


def TranslateRC(DNA):
	'''Translate the reverse complement of DNA'''
	assert type(DNA) == str or type(DNA) == unicode, 'Error, input sequence must be a string or unicode'
	DNA = RC(DNA)
	return Translate(DNA)

def GetCodons(AA):
	'''Get the codons for a specified AA. Returns a list of strings.'''
	assert AA.upper() in 'FLSYCWPHERIMTNKVADQG', 'Error, %s is not a valid amino acid' % str(AA)
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
	codons = eval(AA)
	return codons


############################################################

def combine(input_list, max_total=50000):
	'''
	Takes a list of lists (nucleotides or amino acids) and makes every combination of these while retaining the internal order.
	For example [['A'], ['T','A','G'], ['C', 'A']] will result in six sequences: ['ATC', 'AAC', 'AGC', 'ATA', 'AAA', 'AGA']
	max_total puts a limit on the maximum number of sequences that can be computed. The number tends to increase explosively....
	The output is a generator with a list of the resulting sequences.
	'''

	#make sure the input will not result in too many sequences
	list_of_nums = [len(x) for x in input_list]
	total = reduce(lambda x, y: x*y, list_of_nums)
	assert total < max_total, 'The sequence "%s" would result in %s total sequences, this is above the set maximum of %s sequences.' % (string, total, max_total)

	#now combine them
	#first copy the output list the number of times that there are nucleotides in the current position
	output = []
	for pos in input_list:
		output_len = len(output)
		if output_len == 0:
			output.extend(pos)
			output_len = len(output)
		else:
			output.extend(output*(len(pos)-1)) #duplicate output the number of times that there are new nucleotides to add
			for i in range(0, len(pos)): #for every nucleotide to be added
				for j in range(0, output_len): #add that nucleotide the number of times that the prevous output was long (before the last duplication)
					output[j+i*output_len] += pos[i]

	yield output #return a generator



def UnAmb(string):
	'''Converts an ambigous nucletotide sequence to a list of sequences containing only A, T, C and G (as appropriate)'''
	assert type(string) is str, 'Error, the input has to be a string.'
	string = string.upper()	

	pos_list = []
	for letter in string:
		assert letter in 'NMRWSYKVHDBGATC', 'Error, "%s" is not a valid ambigous nucleotide.' 
		if 'A' == letter:
			pos_list.append(['A'])

		elif 'C' == letter:
			pos_list.append(['C'])

		elif 'T' == letter:
			pos_list.append(['T'])

		elif 'G' == letter:
			pos_list.append(['G'])

		elif 'M' == letter:
			pos_list.append(['A','C'])

		elif 'Y' == letter:
			pos_list.append(['C','T'])

		elif 'K' == letter:
			pos_list.append(['G','T'])

		elif 'S' == letter:
			pos_list.append(['C','G'])

		elif 'W' == letter:
			pos_list.append(['A','T'])

		elif 'R' == letter:
			pos_list.append(['A','G'])

		elif 'H' == letter:
			pos_list.append(['C','T','A'])

		elif 'V' == letter:
			pos_list.append(['C','A','G'])

		elif 'D' == letter:
			pos_list.append(['T','A','G'])

		elif 'B' == letter:
			pos_list.append(['C','T','G'])

		elif 'N' == letter: 
			pos_list.append(['C','T','A','G'])

	return combine(pos_list)


def randomizeSeq(seq):
	'''
	Randomize a given DNA or protein sequence.
	'''
	assert type(seq) == str or type(seq) == unicode, 'Error, input sequence must be a string or unicode'
	seq = list(seq)
	output_list = []
	while seq:
		char = random.choice(seq)
		seq.remove(char)
		output_list.append(char)
	output = ''.join(output_list)
	return output


#add function for fetching DNA from uniprot, ncbi...


############### Identity functions ##########################

def PairIdent(Seq1, Seq2, single_gaps=True):
	'''
	Takes two aligned sequences and returns their percent identity.
	Assumes that Seq1 and Seq2 are sequence strings.
	single_gaps determines whether single gaps should be included or not.
	'''

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

	if single_gaps is True: #include single gaps 
		percent = round(100*(i/l),1) #calculate identity
	elif single_gaps is False: #exclude single gaps
		if n == l: #if gaps same as total length identity is 0
			percent = 0.0
		else:	
			percent = round(100*(i/(l-n)),1) #calculate identity
	return percent




def SingleIdent(align_file, single_gaps=True):
	'''
	Get identities for all proteins compared to one reference sequence in aligned FASTA file.
	Input is the absolute path to a file.
	Assumes that align_file is a fasta-formatted alignment.
	Assumes that reference sequence is the first one.
	Returns a tuple of the results.
	'''

	results = ()
	input_handle = open(align_file, "rU")
	alignment = AlignIO.read(input_handle, "fasta")
	input_handle.close()
	for record2 in range(1, len(alignment)):
		percent = PairIdent(alignment[0].seq, alignment[record2].seq, single_gaps)
		results += (alignment[0].id + ' vs ' + alignment[record2].id+' '+str(percent),)
	return results

	
def AllIdent(align_file, single_gaps=True):
	'''
	Get identities for all combinations of proteins in aligned FASTA file.
	Input is the absolute path to a file.
	Assumes that align_file is a fasta-formatted alignment.
	Returns a tuple of the results.
	'''

	results = ()
	input_handle = open(align_file, "rU")
	alignment = AlignIO.read(input_handle, "fasta")
	input_handle.close()
	for record in range(len(alignment)):
		for record2 in range(record+1, len(alignment)):
			percent = PairIdent(alignment[record].seq, alignment[record2].seq, single_gaps)
			results += (alignment[record].id + ' vs ' + alignment[record2].id+' '+str(percent),)
	return results

####################### End identity functions #################################


####################### Alignment functions ####################################

def alignString(string):


	records_handle = StringIO(string) #turn string into a handle
	tempdata = records_handle.getvalue()
	
	#for separate fasta entries
	muscle_cline = MuscleCommandline()
	stdout, stderr = muscle_cline(stdin=tempdata)
	stdout = fasta.parseString(stdout)

#	#for aligned fasta entries
#	muscle_cline = MuscleCommandline(clw=True)
#	stdout, stderr = muscle_cline(stdin=tempdata)
#	print(stdout)
	
#	#the clustalw-type output can be further formated
#	align = AlignIO.read(StringIO(stdout), "clustal")
#	print(align)

	return stdout
	
def alignList(tlist):
	'''
	Muscle alignment, takes a list of tuples as input and aligns all at once.
	Example input: [('id1', 'CACC'), ('id2', 'CATC'), ('id3', 'TACC')]
	The alignment data is returned as id, dna tuples in a generator.
	'''
	#prepare input as a 'virtual' FASTA file
		
	string = ''
	for entry in tlist:
		if entry[0][0] != '>':
			string += '>%s\n%s\n' % entry
		elif entry[0][0] == '>':
			string += '%s\n%s\n' % entry
		else:
			print('Muscle name error')	
	return alignString(string)

	
def alignFile(filepath):
	'''
	Function for aligning all the records in a file containing fasta sequences.
	The input is the absolute filepath including the filename (as a string).
	The alignment data is returned as id, dna tuples in a generator.
	'''
	string = open(filepath)
	return alignString(string)




########################################################

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
	


def smithWaterman(DNA1, DNA2):
	"""Module aligns two DNA sequences and returns aligned seqA, alignment "bars (|)" and aligned seqB"""
	##got script from Kevin Kwok and modified...
	##http://kevinakwok.tumblr.com/post/18372160357/smith-waterman-algorithm
	seqA = DNA1
	seqB = DNA2
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

