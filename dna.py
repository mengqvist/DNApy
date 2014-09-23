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

import fasta
import string
import random
import genbank
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


def Translate(DNA, table=1):
	"""
	Returns protein sequence from DNA string input.
	The table variable specifies which codon table should be used.
	table defaults to the standard codon table 1
	"""
	assert type(DNA) == str or type(DNA) == unicode, 'Error, input sequence must be a string or unicode'
	codons = CodonTable(table).getCodons()

	protein = []
	DNA = DNA.upper()
	DNA = DNA.replace('\n', '') #get rid of row breaks
	DNA = DNA.replace('U', 'T')

	for i in range(len(DNA)):
		if i%3==0:
			if i+3>len(DNA):
				pass
			elif any(DNA[i:(i+3)] in s for s in codons['F']):
				protein.append('F')
			elif any(DNA[i:(i+3)] in s for s in codons['L']):
				protein.append('L')
			elif any(DNA[i:(i+3)] in s for s in codons['S']):
				protein.append('S')
			elif any(DNA[i:(i+3)] in s for s in codons['Y']):
				protein.append('Y')
			elif any(DNA[i:(i+3)] in s for s in codons['stop']):
				protein.append('*')
			elif any(DNA[i:(i+3)] in s for s in codons['C']):
				protein.append('C')
			elif any(DNA[i:(i+3)] in s for s in codons['W']):
				protein.append('W')
			elif any(DNA[i:(i+3)] in s for s in codons['P']):
				protein.append('P')
			elif any(DNA[i:(i+3)] in s for s in codons['H']):
				protein.append('H')
			elif any(DNA[i:(i+3)] in s for s in codons['E']):
				protein.append('E')
			elif any(DNA[i:(i+3)] in s for s in codons['R']):
				protein.append('R')
			elif any(DNA[i:(i+3)] in s for s in codons['I']):
				protein.append('I')
			elif any(DNA[i:(i+3)] in s for s in codons['M']):
				protein.append('M')
			elif any(DNA[i:(i+3)] in s for s in codons['T']):
				protein.append('T')
			elif any(DNA[i:(i+3)] in s for s in codons['N']):
				protein.append('N')
			elif any(DNA[i:(i+3)] in s for s in codons['K']):
				protein.append('K')
			elif any(DNA[i:(i+3)] in s for s in codons['V']):
				protein.append('V')
			elif any(DNA[i:(i+3)] in s for s in codons['A']):
				protein.append('A')
			elif any(DNA[i:(i+3)] in s for s in codons['D']):
				protein.append('D')
			elif any(DNA[i:(i+3)] in s for s in codons['Q']):
				protein.append('Q')
			elif any(DNA[i:(i+3)] in s for s in codons['G']):
				protein.append('G')
			else:
				raise Error, '"%s" is not a valid codon' % DNA[i:(i+3)]
	return ''.join(protein)	


def TranslateRC(DNA, table=1):
	'''Translate the reverse complement of DNA'''
	assert type(DNA) == str or type(DNA) == unicode, 'Error, input sequence must be a string or unicode'
	DNA = RC(DNA)
	return Translate(DNA, table)


	
def GetCodons(AA, table=1, separate=False):
	'''
	Get the codons for a specified AA. Returns a list of strings.
	The variable table specifies which codon table should be used.
	table defaults to the standard codon table 1
	The separate variable determines whether codons for one amino acid with dissimilar first two nucleotides should be seperated out.
	For example if separate=False the codons for L are 	['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'].
	If separate=True they are split up as L = ['TTA', 'TTG'] and L2 = ['CTT', 'CTC', 'CTA', 'CTG']
	'''
	AA = AA.upper()
	if separate is False:
		assert len(AA) == 1, 'Error, function takes a single amino acid as input'
		assert AA in 'FLSYCWPHERIMTNKVADQG', 'Error, %s is not a valid amino acid' % str(AA)
	
	codons = CodonTable(table).getCodons(separate)
	aacodons = codons[AA]
	return aacodons	
	
def ReverseTranslate(protein, table=1):
	'''
	Translate protein to DNA.
	The input is a protein sequence as a string.
	The output is a list of codons (with ambigous bases) that describe that protein.
	For some amino acids there will be two possible ambigous codons.
	Run the combine() function to convert the list to all possible dna sequences.
	'''
	assert type(protein) == str or type(protein) == unicode, 'Error, input sequence must be a string or unicode'
	dnalist = []
	for aa in protein:
		dnalist.append(GetCodons(aa, table))
	return dnalist


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
	return output #return the list

def Amb(list):
	'''
	This function finds the degenerate nucleotide for a list containing CATG nucleotides.
	Output is a single ambiguous DNA nucleotide as a string.
	Example input is: ['A','T','C','G']. The output for that input is 'N'
	'''
	list = [s.upper() for s in list]
	if all([x in 'A' for x in list]): #test whether each item in a string is present in the list
		output = 'A'

	elif all([x in 'G' for x in list]):
		output = 'G'
		
	elif all([x in 'C' for x in list]): 
		output = 'C'

	elif all([x in 'T' for x in list]):
		output = 'T'

	elif all([x in 'CT' for x in list]): 
		output = 'Y'

	elif all([x in 'GT' for x in list]): 
		output = 'K'

	elif all([x in 'AC' for x in list]):
		output = 'M'

	elif all([x in 'CG' for x in list]): 
		output = 'S'

	elif all([x in 'AT' for x in list]): 
		output = 'W'

	elif all([x in 'AG' for x in list]): 
		output = 'R'

	elif all([x in 'CTA' for x in list]): 
		output = 'H'

	elif all([x in 'CAG' for x in list]): 
		output = 'V'

	elif all([x in 'TAG' for x in list]): 
		output = 'D'		
		
	elif all([x in 'CTG' for x in list]): 
		output = 'B'

	elif all([x in 'CTAG' for x in list]): 
		output = 'N'
	else:
		raise ValueError, 'Error, input must be a list of standard GATC nucleotides.'
	return output
	
	
def MultipleAmb(largelist):
	'''
	This function finds the degenerate nucleotide for a list of lists containing CATG nucleotides.
	Output is a DNA string with the ambigous nucleotides.
	Example input is: [['A','T','C','G'],['A','T','C','G'],['G','T']]. The output for that is 'NNK'
	'''
	output = []

	largelist = [x.upper() for x in largelist] #make the uppercase
	for list in largelist:
		output.append(Amb(list))
	
	return ''.join(output)

	
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

	return combine(pos_list) #call combine function and return the result as a list of strings


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

def NCBIfetch(id, db = 'nucleotide', rettype = 'gb', retmode = 'text'):
	'''
	Retrieve records from NCBI.
	Format of the string is http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=<database>&id=<uid_list>&rettype=<retrieval_type>&retmode=<retrieval_mode>

	For any series of more than 100 requests, do this at weekends or outside USA peak times. This is up to you to obey.
	Use the http://eutils.ncbi.nlm.nih.gov address, not the standard NCBI Web address. Biopython uses this web address.
	Make no more than three requests every seconds (relaxed from at most one request every three seconds in early 2009).
	Use the optional email parameter so the NCBI can contact you if there is a problem. You can either explicitly set this as a parameter with each call to Entrez (e.g. include email="A.N.Other@example.com" in the argument list).
	'''
	import urllib2
	url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=%s&id=%s&rettype=%s&retmode=%s&tool=DNApy&email=martin.engqvist@chalmers.se' % (db, id, rettype, retmode)
	file = urllib2.urlopen(url).read() #this returns the result as a string. I'll need to parse it to get the info out.
	return file

def NCBIfetchGB(id):
	'''
	Same as NCBIfetch, but instead of returning a text string it returns a genbank object.
	'''
	file = NCBIfetch(id).split('\n') #get file and split it
	file = [x+'\n' for x in file] #add back the linebreaks
	x = genbank.gbobject()
	x.readgb(file)
	return x

	
def NCBIblast(blast_type, database, seq): #function for blasting 
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
	

				
				
class CodonTable:
	'''
	A DNA codon object.
	Used to retrieve codon tables and codons for specified codon tables.	
	Pass a valid integer value when instantiating to choose which codon table to use.
	'''
	def __init__(self, number):
		self.code = False
		self.table = False
		self.codons = False
		self.setTable(number) #get the specified codon table (returned as list of strings)
		self.setCodons() #convert the codon table information to codons
		
	
	def setTable(self, number):
		'''
		Find information for specified genetic code and use for downstream methods.
		Method is not intended for direct use.
		'''
		number = int(number)
		
		if number == 1:
			#Genetic Code [1]
			code = "Standard Code (transl_table=1)"
			AAs  = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
			Starts = "---M---------------M---------------M----------------------------"
			Base1  = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
			Base2  = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
			Base3  = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"

		elif number == 2:
			#Genetic Code [2]
			code = "Vertebrate Mitochondrial Code (transl_table=2)"
			AAs  = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG"
			Starts = "--------------------------------MMMM---------------M------------"
			Base1  = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
			Base2  = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
			Base3  = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"

		elif number == 3:
			#Genetic Code [3]
			code = "Yeast Mitochondrial Code (transl_table=3)"
			AAs  = "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
			Starts = "----------------------------------MM----------------------------"
			Base1  = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
			Base2  = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
			Base3  = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"

		elif number == 4:
			#Genetic Code [4]
			code = "Mold, Protozoan, Coelenterate Mitochondrial Code & Mycoplasma/Spiroplasma  Code (transl_table=4)"
			AAs  = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
			Starts = "--MM---------------M------------MMMM---------------M------------"
			Base1  = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
			Base2  = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
			Base3  = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"

		elif number == 5:
			#Genetic Code [5] 
			code = "Invertebrate Mitochondrial Code (transl_table=5)"
			AAs  = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG"
			Starts = "---M----------------------------MMMM---------------M------------"
			Base1  = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
			Base2  = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
			Base3  = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"

		elif number == 6:
			#Genetic Code [6]
			code = "Ciliate, Dasycladacean and Hexamita Nuclear Code (transl_table=6)"
			AAs  = "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
			Starts = "-----------------------------------M----------------------------"
			Base1  = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
			Base2  = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
			Base3  = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"

		elif number == 9:
			#Genetic Code [9]  
			code = "Echinoderm and Flatworm Mitochondrial Code (transl_table=9)"
			AAs  = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG"
			Starts = "-----------------------------------M---------------M------------"
			Base1  = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
			Base2  = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
			Base3  = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"

		elif number == 10:
			#Genetic Code [10]   
			code = "Euplotid Nuclear Code (transl_table=10)"
			AAs  = "FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
			Starts = "-----------------------------------M----------------------------"
			Base1  = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
			Base2  = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
			Base3  = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"

		elif number == 11:
			#Genetic Code [11]
			code = "Bacterial, Archaeal and Plant Plastid Code (transl_table=11)"
			AAs  = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
			Starts = "---M---------------M------------MMMM---------------M------------"
			Base1  = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
			Base2  = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
			Base3  = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"

		elif number == 12:
			#Genetic Code [12]
			code = "Alternative Yeast Nuclear Code (transl_table=12)"
			AAs  = "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
			Starts = "-------------------M---------------M----------------------------"
			Base1  = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
			Base2  = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
			Base3  = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"
																		   
		elif number == 13:
			#Genetic Code [13]
			code = "Ascidian Mitochondrial Code (transl_table=13)"
			AAs  = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG"
			Starts = "---M------------------------------MM---------------M------------"
			Base1  = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
			Base2  = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
			Base3  = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"

		elif number == 14:
			#Genetic Code [14]
			code = "Alternative Flatworm Mitochondrial Code (transl_table=14)"
			AAs  = "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG"
			Starts = "-----------------------------------M----------------------------"
			Base1  = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
			Base2  = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
			Base3  = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"

		elif number == 15:
			#Genetic Code [15]
			code = "Blepharisma Nuclear Code (transl_table=15)"
			AAs  = "FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
			Starts = "-----------------------------------M----------------------------"
			Base1  = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
			Base2  = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
			Base3  = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"

		elif number == 16:
			#Genetic Code [16]
			code = "Chlorophycean Mitochondrial Code (transl_table=16)"
			AAs  = "FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
			Starts = "-----------------------------------M----------------------------"
			Base1  = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
			Base2  = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
			Base3  = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"

		elif number == 21:
			#Genetic Code [21]
			code = "Trematode Mitochondrial Code (transl_table=21)"
			AAs  = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG"
			Starts = "-----------------------------------M---------------M------------"
			Base1  = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
			Base2  = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
			Base3  = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"

		elif number == 22:
			#Genetic Code [22]
			code = "Scenedesmus obliquus mitochondrial"
			AAs  = "FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
			Starts = "-----------------------------------M----------------------------"
			Base1  = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
			Base2  = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
			Base3  = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"

		elif number == 23:
			#Genetic Code [23]
			code = "Thraustochytrium Mitochondrial Code (transl_table=23)"
			AAs  = "FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
			Starts = "--------------------------------M--M---------------M------------"
			Base1  = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
			Base2  = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
			Base3  = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"

		elif number == 24:			
			#Genetic Code [24]
			code = "Pterobranchia mitochondrial code (transl_table=24)"
			AAs  = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG"
			Starts = "---M---------------M---------------M---------------M------------"
			Base1  = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
			Base2  = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
			Base3  = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"

		elif number == 25:			
			#Genetic Code [25]
			code = "Candidate Division SR1 and Gracilibacteria Code (transl_table=25)"
			AAs  = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG"
			Starts = "---M---------------M---------------M----------------------------"
			Base1  = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
			Base2  = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
			Base3  = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"
		else:
			raise ValueError, '%s is not a valid genetic code number' % number
		self.code = code
		self.table = [code, AAs, Starts, Base1, Base2, Base3]
		
	def setCodons(self):
		'''
		Use a predetermined codon table to generate a dictionary of amino acids with their codons.
		Method is not intended for direct use.
		'''
		code, AAs, Starts, Base1, Base2, Base3 = self.getTable()
		codons = {'start':[], 'F':[], 'L':[], 'S':[], 'Y':[], 'C':[], 'W':[], 'P':[], 'H':[], 'E':[], 'R':[], 'I':[], 'M':[], 'T':[], 'N':[], 'K':[], 'V':[], 'A':[], 'D':[], 'Q':[], 'G':[], 'stop':[]}
		for aa, s, b1, b2, b3 in zip(AAs, Starts, Base1, Base2, Base3):
			codon = b1+b2+b3
			
			if aa == '*': #if aa is stop
				codons['stop'].append(codon)
			elif aa in 'FLSYCWPHERIMTNKVADQG':
				codons[aa].append(codon)
			else:
				raise Error, '"%s" is not a valid amino acid' % aa
				
			if s != '-': #if the codon is start
				codons['start'].append(codon)
				
		self.codons = codons

	######## API intended for use #########

	def getCode(self):
		'''
		Return which genetic code is represented.
		The output is a string which specifies the code'
		Method is not intended for direct use.
		'''
		return self.code

		
	def getTable(self):
		'''
		Return the codon table data for the specified genetic code.
		The output is a list of strings.
		'''
		return self.table

	def getCodons(self, separate=False):
		'''
		Returns a dictionary of amino acids with their codons for the specified codon table.
		The separate variable determines whether codons for one amino acid with dissimilar first two nucleotides should be separated out.
		For example if separate=False the codons for L are 	['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'].
		If separate=True they are split up as L = ['TTA', 'TTG'] and L1 = ['CTT', 'CTC', 'CTA', 'CTG']
		'''
		assert type(separate) is bool, 'Error, "separate" must be True or False'
		if separate is False: #returns a dictionary containing the codon table
			return self.codons 
		elif separate is True:
			newdict = {}
			for aa in 'FLSYCWPHERIMTNKVADQG':
				f = lambda x: [codon[0:2] for codon in x] #function to get all first two nucleotides for an aa
				firsttwolist = list(set(f(self.codons[aa]))) #list of all unique first two nucleotides for an aa. For example ['TT', 'CT'] for leucine
#				print('aa', aa)
#				print('ftl', firsttwolist)
				if len(firsttwolist) > 1: #if there is more than one set of the first two nucleotides for this amino acid
					for i in range(len(firsttwolist)):
						if i == 0:
							newaa = aa
						else:
							newaa = aa+str(i) #add a number after the amino acid
						newdict[newaa] = [x for x in self.codons[aa] if x[0:2] in firsttwolist[i]] #add all the codons that match the first two
				else:
					newdict[aa] = self.codons[aa] #
			return newdict				

					
				
	def printTable(self):
		'''
		Print specified codon table.
		'''
		code, AAs, Starts, Base1, Base2, Base3 = self.table
		print('Code   = %s' % code)
		print('AAs    = %s' % AAs)
		print('Starts = %s' % Starts)
		print('Base1  = %s' % Base1)
		print('Base2  = %s' % Base2)
		print('Base3  = %s' % Base3)

	def printCodons(self):	
		'''
		Print codons specified by codon table.
		'''
		codons = self.getCodons()
		print('start = %s' %codons['start'])
		print('stop  = %s' %codons['stop'])
		for aa in 'FLSYCWPHERIMTNKVADQG':
			print('%s     = %s' % (aa, codons[aa]))
		
		
		
