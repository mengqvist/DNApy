#!/usr/bin/env python3


#DNApy is a DNA editor written purely in python.
#The program is intended to be an intuitive and fully featured
#editor for molecular and synthetic biology.
#Enjoy!
#
#copyright (C) 2014-2017  Martin Engqvist |
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

from resources import fasta
import string
import random
from functools import reduce


class BioSeq(object):
	"""
	"""
	def __init__(self, sequence):
		assert type(sequence) == str or type(sequence) == str, 'Error, input sequence must be a string or unicode'
		self.__check_seq(sequence)

	def __repr__(self):
		return self.sequence

	def __str__(self):
		return self.sequence

	def __check_seq(self, sequence):
		"""
		Verify that there are not weird symbols
		"""

		#if it's all ok
		self.sequence = sequence

	def randomize(self):
		'''
		Randomize a given DNA or protein sequence.
		'''
		seq = list(self.sequence)
		output_list = []
		while seq:
			char = random.choice(seq)
			seq.remove(char)
			output_list.append(char)
		output = ''.join(output_list)
		return output

	def length(self):
		"""
		"""
		return length(self.sequence)



class DNA(BioSeq):

	def __init__(self, sequence):
		BioSeq.__init__(self, sequence)

	def type(self):
		"""
		What type sequence is it
		"""
		return 'DNA'

	def clean_dna(self, ambiguous=False, silent=True):
		'''
		Function for cleaning DNA of non-DNA characters.
		The variable ambiguous (bool) determines whether the full set of ambiguous codons are used or not.
		the variable silent determines whether an output is printed every time a non-DNA character is omitted.
		'''
		if ambiguous == False:
			bases = 'GATC'
		elif ambiguous == True:
			bases = 'GATCRYWSMKHBVDN'
		cleaned_seq = ''
		for char in DNA:
			if (char in string.digits) or (char in string.whitespace) or char==('/'or'\\'or ' ' or '  ') or (char.upper() not in bases):
				if silent == False:
					print(('Character "%s" is not a valid DNA character and was deleted' % char))
			else:
				cleaned_seq += char
		return DNA(cleaned_seq)


	def reverse(self):
		"""
		Returns the reverse of a DNA string.
		"""
		return DNA(self.sequence[::-1])  #makes the reverse of the input string


	def complement(self):
		"""
		Returns the complement of a DNA string.
		"""
		complement = {'a':'t', 't':'a', 'c':'g', 'g':'c', 'y':'r', 'r':'y', 'w':'w', 's':'s', 'k':'m', 'm':'k', 'd':'h', 'v':'b', 'h':'d', 'b':'v', 'n':'n',
						'A':'T', 'T':'A', 'C':'G', 'G':'C', 'Y':'R', 'R':'Y', 'W':'W', 'S':'S', 'K':'M', 'M':'K', 'D':'H', 'V':'B', 'H':'D', 'B':'V', 'N':'N'}
		bases = list(self.sequence)
		bases = [complement[base] for base in bases]
		return DNA(''.join(bases))


	def reverse_complement(self):
		"""
		Returns the reverse complement of a DNA string.
		"""
		return self.reverse().complement()

	def transcribe(self):
		"""
		Return RNA version of DNA
		"""
		transl = {'t':'u', 'T':'U'}

		RNA = []
		for base in self.sequence:
			if base in ['T', 't']:
				RNA.append(transl[base])
			else:
				RNA.append(base)

		return ''.join(RNA)

	def translate(self, table=1):
		"""
		Returns protein sequence from DNA string input.
		The table variable specifies which codon table should be used.
		table defaults to the standard codon table 1
		"""
		codons = CodonTable(table).getCodons()

		protein = []
		RNA = self.sequence.upper().transcribe()

		for i in range(0, len(RNA), 3):
			if i+3>len(RNA):
				pass
			elif any(RNA[i:(i+3)] in s for s in codons['F']):
				protein.append('F')
			elif any(RNA[i:(i+3)] in s for s in codons['L']):
				protein.append('L')
			elif any(RNA[i:(i+3)] in s for s in codons['S']):
				protein.append('S')
			elif any(RNA[i:(i+3)] in s for s in codons['Y']):
				protein.append('Y')
			elif any(RNA[i:(i+3)] in s for s in codons['*']):
				protein.append('*')
			elif any(RNA[i:(i+3)] in s for s in codons['C']):
				protein.append('C')
			elif any(RNA[i:(i+3)] in s for s in codons['W']):
				protein.append('W')
			elif any(RNA[i:(i+3)] in s for s in codons['P']):
				protein.append('P')
			elif any(RNA[i:(i+3)] in s for s in codons['H']):
				protein.append('H')
			elif any(RNA[i:(i+3)] in s for s in codons['E']):
				protein.append('E')
			elif any(RNA[i:(i+3)] in s for s in codons['R']):
				protein.append('R')
			elif any(RNA[i:(i+3)] in s for s in codons['I']):
				protein.append('I')
			elif any(RNA[i:(i+3)] in s for s in codons['M']):
				protein.append('M')
			elif any(RNA[i:(i+3)] in s for s in codons['T']):
				protein.append('T')
			elif any(RNA[i:(i+3)] in s for s in codons['N']):
				protein.append('N')
			elif any(RNA[i:(i+3)] in s for s in codons['K']):
				protein.append('K')
			elif any(RNA[i:(i+3)] in s for s in codons['V']):
				protein.append('V')
			elif any(RNA[i:(i+3)] in s for s in codons['A']):
				protein.append('A')
			elif any(RNA[i:(i+3)] in s for s in codons['D']):
				protein.append('D')
			elif any(RNA[i:(i+3)] in s for s in codons['Q']):
				protein.append('Q')
			elif any(RNA[i:(i+3)] in s for s in codons['G']):
				protein.append('G')
			else:
				protein.append('X')
				#raise Exception('"%s" is not a valid codon' % DNA[i:(i+3)])
				# better just ignore it, because else we destroy usability
		return Protein(''.join(protein))

	def translate_reverse_complement(self, table=1):
		'''Translate the reverse complement of DNA'''
		return self.reverse_complement().translate(table)

	def count_bases(self):
		'''
		Counts the number of different bases in a DNA sequence.
		Input should be a string comprising dna bases. These can be a, t, c, g or any of the ambiguous bases.
		Output is a dictionary with DNA base keys and integer values.
		'''
		seq = self.sequence.upper()
		assert all([s in 'GATCRYWSMKHBVDN' for s in seq]) is True, 'Error, unknown dna base(s) %s in sequence: %s' % (str([s for s in seq if s not in 'GATCRYWSMKHBVDN']), seq)

		bases = {'G':seq.count('G'),
		'A':seq.count('A'),
		'T':seq.count('T'),
		'C':seq.count('C'),
		'R':seq.count('R'),
		'Y':seq.count('Y'),
		'W':seq.count('W'),
		'S':seq.count('S'),
		'M':seq.count('M'),
		'K':seq.count('K'),
		'H':seq.count('H'),
		'B':seq.count('B'),
		'P':seq.count('P'),
		'H':seq.count('H'),
		'V':seq.count('V'),
		'D':seq.count('D'),
		'N':seq.count('N')}

		assert sum([s[1] for s in bases.items()]) == len(seq), 'Error, total number of counted bases is not the same as the sequence length'

		return bases


	def count_codons(self):
		'''
		Counts codons in a DNA sequence.
		Input should be a string comprising whole codons.
		No ambiguous codons allowed.
		Output is a dictionary with codon keys and integer values.
		'''
		assert len(self.sequence) % 3 == 0, 'Error, the sequence must be a string comprising whole codons.'

		codons = {'UUU': 0, 'UUC': 0, 'UUA': 0, 'UUG': 0, 'CUU': 0,
					'CUC': 0, 'CUA': 0, 'CUG': 0, 'AUU': 0, 'AUC': 0,
					'AUA': 0, 'AUG': 0, 'GUU': 0, 'GUC': 0, 'GUA': 0,
					'GUG': 0, 'UAU': 0, 'UAC': 0, 'UAA': 0, 'UAG': 0,
					'CAU': 0, 'CAC': 0, 'CAA': 0, 'CAG': 0, 'AAU': 0,
					'AAC': 0, 'AAA': 0, 'AAG': 0, 'GAU': 0, 'GAC': 0,
					'GAA': 0, 'GAG': 0, 'UCU': 0, 'UCC': 0, 'UCA': 0,
					'UCG': 0, 'CCU': 0, 'CCC': 0, 'CCA': 0, 'CCG': 0,
					'ACU': 0, 'ACC': 0, 'ACA': 0, 'ACG': 0, 'GCU': 0,
					'GCC': 0, 'GCA': 0, 'GCG': 0, 'UGU': 0, 'UGC': 0,
					'UGA': 0, 'UGG': 0, 'CGU': 0, 'CGC': 0, 'CGA': 0,
					'CGG': 0, 'AGU': 0, 'AGC': 0, 'AGA': 0, 'AGG': 0,
					'GGU': 0, 'GGC': 0, 'GGA': 0, 'GGG': 0}

		seq = self.sequence.upper().translate()
		for i in range(0, len(seq), 3):
			codon = seq[i:i+3]
			if codon in list(codons.keys()):
				codons[codon] += 1
			else:
				raise ValueError('Codon %s is not valid.' % codon)
		return codons


	def make_codon_freq_table(self, file):
		'''
		Input is a file path.
		Counts the usage of each codon in a FASTA file of DNA sequences.
		Then converts that as codon usage per 1000 codons.
		Good for generating codon tables.
		Output is a dictionary of codon frequencies per 1000 codons and the total number in brackets.
		'''

		num_table = {'UUU': 0, 'UUC': 0, 'UUA': 0, 'UUG': 0, 'CUU': 0,
					'CUC': 0, 'CUA': 0, 'CUG': 0, 'AUU': 0, 'AUC': 0,
					'AUA': 0, 'AUG': 0, 'GUU': 0, 'GUC': 0, 'GUA': 0,
					'GUG': 0, 'UAU': 0, 'UAC': 0, 'UAA': 0, 'UAG': 0,
					'CAU': 0, 'CAC': 0, 'CAA': 0, 'CAG': 0, 'AAU': 0,
					'AAC': 0, 'AAA': 0, 'AAG': 0, 'GAU': 0, 'GAC': 0,
					'GAA': 0, 'GAG': 0, 'UCU': 0, 'UCC': 0, 'UCA': 0,
					'UCG': 0, 'CCU': 0, 'CCC': 0, 'CCA': 0, 'CCG': 0,
					'ACU': 0, 'ACC': 0, 'ACA': 0, 'ACG': 0, 'GCU': 0,
					'GCC': 0, 'GCA': 0, 'GCG': 0, 'UGU': 0, 'UGC': 0,
					'UGA': 0, 'UGG': 0, 'CGU': 0, 'CGC': 0, 'CGA': 0,
					'CGG': 0, 'AGU': 0, 'AGC': 0, 'AGA': 0, 'AGG': 0,
					'GGU': 0, 'GGC': 0, 'GGA': 0, 'GGG': 0}
		records = fasta.parseFile(file)
		for record in records:
			cds = record[1]
			codons = count_codons(cds)
			for key in list(codons.keys()):
				num_table[key] += codons[key]

		#sum codons
		sum = 0.0
		for key in list(num_table.keys()):
			sum += num_table[key]

		#divide each by the sum and multiply by 1000
		freq_table = {}
		for key in list(num_table.keys()):
			freq_table[key] = '%s(%s)' % (1000*(num_table[key]/sum), num_table[key]) #ouput is following format: freq/thousand(number)
		return freq_table








class RNA(BioSeq):

	def __init__(self, sequence):
		BioSeq.__init__(self, sequence)

	def type(self):
		"""
		What type sequence is it
		"""
		return 'RNA'





class Protein(BioSeq):
	def __init__(self, sequence):
		BioSeq.__init__(self, sequence)


	def type(self):
		"""
		What type sequence is it
		"""
		return 'Protein'

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
	#	if separate is False: #what the heck is this if clause for?? Can I remove it?
		assert len(AA) == 1, 'Error, function takes a single amino acid as input'
		assert AA in 'FLSYCWPHERIMTNKVADQG*', 'Error, %s is not a valid amino acid' % str(AA)

		codons = CodonTable(table).getCodons(separate)

		return codons[AA]

	def ReverseTranslate(protein, table=1):
		'''
		Translate protein to DNA.
		The input is a protein sequence as a string.
		The output is a list of codons (with ambigous bases) that describe that protein.
		For some amino acids there will be two possible ambigous codons.
		Run the combine() function to convert the list to all possible dna sequences.
		'''
		assert type(protein) == str or type(protein) == str, 'Error, input sequence must be a string or unicode'
		dnalist = []
		for aa in protein:
			dnalist.append(GetCodons(aa, table))
		return dnalist

	def one_to_three(one_letter):
		'''
		Convert a one letter code amino acid to a three letter code.
		'''
		assert one_letter.upper() in 'FLSYCWPHERIMTNKVADQG*X', 'Error, %s is not a valid amino acid' % one_letter

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
		'*':'***',
		'X':'Unknown'}
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
		assert one_letter.upper() in 'FLSYCWPHERIMTNKVADQG*X', 'Error, %s is not a valid amino acid' % one_letter
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
		'G':'Glycine',
		'X': 'Unknown'}
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
								'glycine',
								'unknown'], ('Error, %s is not a valid amino acid' % full)

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
				'glycine':'G',
				'unknown': 'X'} # to handle N as a nuclec acid
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



class CodonTable(object):
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
			code = "Scenedesmus obliquus mitochondrial (transl_table=22)"
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
			raise ValueError('%s is not a valid genetic code number' % number)
		self.code = code
		self.table = [code, AAs, Starts, Base1, Base2, Base3]


	def setCodons(self):
		'''
		Use a predetermined codon table to generate a dictionary of amino acids with their codons.
		Method is not intended for direct use.
		'''
		code, AAs, Starts, Base1, Base2, Base3 = self.getTable()
		codons = {'start':[], 'F':[], 'L':[], 'S':[], 'Y':[], 'C':[], 'W':[], 'P':[], 'H':[], 'E':[], 'R':[], 'I':[], 'M':[], 'T':[], 'N':[], 'K':[], 'V':[], 'A':[], 'D':[], 'Q':[], 'G':[], '*':[]}
		for aa, s, b1, b2, b3 in zip(AAs, Starts, Base1, Base2, Base3):
			codon = b1+b2+b3

			if aa in 'FLSYCWPHERIMTNKVADQG*':
				codons[aa].append(codon)
			else:
				raise Error('"%s" is not a valid amino acid' % aa)

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
		If separate=True they are split up as L = [['TTA', 'TTG'], ['CTT', 'CTC', 'CTA', 'CTG']]
		'''
		assert type(separate) is bool, 'Error, "separate" must be True or False'
		if separate is False: #returns a dictionary containing the codon table
			return self.codons
		elif separate is True:
			newdict = {}
			for aa in 'FLSYCWPHERIMTNKVADQG*':
				f = lambda x: [codon[0:2] for codon in x] #function to get all first two nucleotides for an aa
				firsttwolist = list(set(f(self.codons[aa]))) #list of all unique first two nucleotides for an aa. For example ['TT', 'CT'] for leucine
#				print('aa', aa)
#				print('ftl', firsttwolist)
				if len(firsttwolist) > 1: #if there is more than one set of the first two nucleotides for this amino acid
					newdict[aa] = []
					for i in range(len(firsttwolist)):
						newdict[aa].append([x for x in self.codons[aa] if x[0:2] in firsttwolist[i]]) #add all the codons that match the first two
				else:
					newdict[aa] = self.codons[aa] #
			return newdict



	def printTable(self):
		'''
		Print specified codon table.
		'''
		code, AAs, Starts, Base1, Base2, Base3 = self.table
		print(('\nCode   = %s' % code))
		print(('AAs    = %s' % AAs))
		print(('Starts = %s' % Starts))
		print(('Base1  = %s' % Base1))
		print(('Base2  = %s' % Base2))
		print(('Base3  = %s' % Base3))

	def printCodons(self):
		'''
		Print codons specified by codon table.
		'''
		codons = self.getCodons()
		print(('start = %s' %codons['start']))
		print(('stop  = %s' %codons['stop']))
		for aa in 'FLSYCWPHERIMTNKVADQG*':
			print(('%s     = %s' % (aa, codons[aa])))




######################
######################




############################################################

def combine(self, input_list, max_total=50000):
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



def listupper(t):
	'''
	Capitalizes all strings in nested lists.
	'''
	if isinstance(t, list):
		return [listupper(s) for s in t]
	elif isinstance(t, str):
		return t.upper()
	else:
		return t


def Amb(list):
	'''
	This function finds the degenerate nucleotide for a list containing CATG nucleotides.
	Output is a single ambiguous DNA nucleotide as a string.
	Example input is: ['A','T','C','G']. The output for that input is 'N'
	'''
	list = listupper(list)

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
		raise ValueError('Error, input must be a list of standard GATC nucleotides.')
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
	'''
	Converts an ambiguous nucleotide sequence to a list of sequences containing only A, T, C and G (as appropriate).
	'''
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




def commonNuc(nuc_list, greedy=False):
	"""
	This function takes a list of lists and finds all degenerate symbols that represent
	at least one nucleotide from each of the lists.
	The variable "greedy" determines whether the algorithm is greedy or not.

	With greedy=False
	An example input is: [['T', 'C', 'A', 'G'], ['T', 'C'], ['T', 'C']].
	T and C are both present in all lists, therefore, both 'T' and 'C' are acceptable returned as ['T', 'C'].

	With greedy=False
	Another example input is: [['G'], ['T'], ['T']].
	In this case either G or T is present in all lists, therefore the only acceptable output is ['K'] (ambiguous nucleotide for G and T).


	With greedy=True
	For the input: [['T', 'C', 'A', 'G'], ['T', 'C'], ['T', 'C']]
	The greedy output includes all degenerate nucleotides that contain the desired regular nucleotides:
	['C', 'T', 'Y', 'K', 'M', 'S', 'W', 'H', 'V', 'D', 'B', 'N']

	With greedy=True
	For the input: [['G'], ['T'], ['T']]
	The greedy output is: ['K', 'D', 'B', 'N']
	"""
	nuc_list = listupper(nuc_list)
	output = []

	if all(['A' in s for s in nuc_list]):
		output.append('A')

	if all(['G' in s for s in nuc_list]):
		output.append('G')

	if all(['C' in s for s in nuc_list]):
		output.append('C')

	if all(['T' in s for s in nuc_list]):
		output.append('T')

	if greedy is False and len(output)>0:
		return output


	if all(['C' in s or 'T' in s for s in nuc_list]):
		output.append('Y')

	if all(['G' in s or 'T' in s for s in nuc_list]):
		output.append('K')

	if all(['A' in s or 'C' in s for s in nuc_list]):
		output.append('M')

	if all(['C' in s or 'G' in s for s in nuc_list]):
		output.append('S')

	if all(['A' in s or 'T' in s for s in nuc_list]):
		output.append('W')

	if all(['A' in s or 'G' in s for s in nuc_list]):
		output.append('R')

	if greedy is False and len(output)>0:
		return output


	if all(['C' in s or 'T' in s or 'A' in s for s in nuc_list]):
		output.append('H')

	if all(['C' in s or 'A' in s or 'G' in s for s in nuc_list]):
		output.append('V')

	if all(['T' in s or 'A' in s or 'G' in s for s in nuc_list]):
		output.append('D')

	if all(['C' in s or 'T' in s or 'G' in s for s in nuc_list]):
		output.append('B')

	if greedy is False and len(output)>0:
		return output


	if all(['C' in s or 'T' in s or 'A' in s or 'G' in s for s in nuc_list]):
		output.append('N')

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
