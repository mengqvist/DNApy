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


#TODO
#finish protein class
#make freq table class and use it to compute/read codon frequency tables. Use that to assign codons in a probabalistic way when reverse-translating protein.
#finish ambDNA class
#finish ambRNA class
#finish codon class

## Sub-classing of the str class and modification of some of the methods ##

class _BioSeq(str):
	"""
	"""
	def __init__(self, sequence, alphabet=''):
		assert type(sequence) == str or type(sequence) == str, 'Error, input sequence must be a string or unicode'
		self.alphabet = alphabet
		self.__check_seq(sequence)
		self._mytype = None

	def __repr__(self):
		return self.sequence

	def __str__(self):
		return self.sequence

	def __iter__(self):
		for item in iter(self.sequence):
			yield type(self)(item)

	def __add__(self, b):
		return type(self)(self.sequence+str(b))

	def __mul__(self, b):
		return type(self)(self.sequence*b)

	def __getitem__(self, index):
		return type(self)(self.sequence[index])

	## Unmodified builtins ##
	#__class__
	#__contains__
	#__delattr__
	#__dir__
	#__doc__
	#__eq__
	#__format__
	#__ge__
	#__getattribute__
	#__getnewargs__
	#__gt__
	#__hash__
	#__init__
	#__le__
	#__len__
	#__lt__
	#__mod__
	#__ne__
	#__new__
	#__reduce__
	#__reduce_ex__
	#__rmod__
	#__rmul__
	#__setattr__
	#__sizeof__
	#__subclasshook__


	def __check_seq(self, sequence):
		"""
		Verify that there are no weird symbols
		"""
		valid_chars = [s for s in sequence if s in self.alphabet]
		diff = len(sequence) - len(valid_chars)
		if diff != 0:
			print('%s characters were non-compliant with %s and were removed from sequence' % (diff, self.alphabet))
		self.sequence = ''.join(valid_chars)

	### Re-define builtins for string ###

	def capitalize(self):
		"""
		"""
		return type(self)(self.sequence.capitalize())

	def casefold(self):
		"""
		"""
		return type(self)(self.sequence.casefold())

	def center(*arg):
		"""
		"""
		raise NotImplementedError

	def encode(*arg):
		"""
		"""
		raise NotImplementedError

	def expandtabs(*arg):
		"""
		"""
		raise NotImplementedError

	def format(*arg):
		"""
		"""
		raise NotImplementedError

	def format_map(*arg):
		"""
		"""
		raise NotImplementedError

	def join(self, string):
		"""
		"""
		return type(self)(self.sequence.join(string))

	def ljust(*arg):
		"""
		"""
		raise NotImplementedError

	def lower(self):
		"""
		"""
		return type(self)(self.sequence.lower())

	def lstrip(*arg):
		"""
		"""
		raise NotImplementedError

	def maketrans(*arg):
		"""
		"""
		raise NotImplementedError

	def partition(self, sep):
		"""
		"""
		return (type(self)(s) for s in self.sequence.partition(sep))

	def replace(self, old, new, count=None):
		"""
		"""
		if count:
			return type(self)(self.sequence.replace(old, new, count))
		else:
			return type(self)(self.sequence.replace(old, new))

	def rjust(*arg):
		"""
		"""
		raise NotImplementedError

	def rpartition(self, sep):
		"""
		"""
		return (type(self)(s) for s in self.sequence.rpartition(sep))

	def rsplit(self, sep=None, maxsplit=-1):
		"""
		"""
		return [type(self)(s) for s in self.sequence.rsplit(sep, maxsplit)]

	def rstrip(*arg):
		"""
		"""
		raise NotImplementedError

	def split(self, sep=None, maxsplit=-1):
		"""
		"""
		return [type(self)(s) for s in self.sequence.split(sep, maxsplit)]

	def splitlines(*arg):
		"""
		"""
		raise NotImplementedError

	def strip(*arg):
		"""
		"""
		raise NotImplementedError

	def swapcase(self):
		"""
		"""
		return type(self)(self.sequence.swapcase())

	def title(*arg):
		"""
		"""
		raise NotImplementedError

	def translate(*arg):
		"""
		"""
		raise NotImplementedError

	def upper(self):
		"""
		"""
		return type(self)(self.sequence.lower())

	def zfill(*arg):
		"""
		"""
		raise NotImplementedError


	## Unmodified methods ##
	#endswith
	#find
	#index
	#isalnum
	#isalpha
	#isdecimal
	#isdigit
	#isidentifier
	#islower
	#isnumeric
	#isprintable
	#isspace
	#istitle
	#isupper
	#rfind
	#rindex
	#startswith


	######################################

	def length(self):
		"""
		"""
		return length(self.sequence)


	def randomize(self):
		'''
		Randomize a given DNA, RNA or protein sequence.
		'''
		seq = list(self.sequence)
		output_list = []
		while seq:
			char = random.choice(seq)
			seq.remove(char)
			output_list.append(char)
		return type(self)(''.join(output_list))


	def clean(self, silent=True):
		'''
		Function for cleaning DNA of non-DNA characters.
		The variable ambiguous (bool) determines whether the full set of ambiguous codons are used or not.
		the variable silent determines whether an output is printed every time a non-DNA character is omitted.
		'''

		cleaned_seq = []
		for char in self.sequence:
			if char not in self.alphabet:
				if silent == False:
					print(('Character "%s" was deleted, it is not a valid character of type "%s"' % (char, self._mytype) ))
			else:
				cleaned_seq += char
		return type(self)(''.join(cleaned_seq))


	def type(self):
		"""
		What type sequence is it
		"""
		return self._mytype


## Common class to be used for DNA and RNA through sub-classing ##

class _NucleotideBaseClass(_BioSeq):

	def __init__(self, sequence, alphabet):
		self.alphabet = 'atcgATCG'
		_BioSeq.__init__(self, sequence, alphabet)


	def reverse(self):
		"""
		Returns the reverse of a DNA or RNA string.
		"""
		return type(self)(self.sequence[::-1])  #makes the reverse of the input string


	def complement(self):
		"""
		Returns the complement of a DNA or RNA string.
		"""
		if 'u' in self.alphabet:
			transl = {'a':'u', 'u':'a', 'c':'g', 'g':'c', 'y':'r', 'r':'y', 'w':'w', 's':'s', 'k':'m', 'm':'k', 'd':'h', 'v':'b', 'h':'d', 'b':'v', 'n':'n',
						'A':'U', 'U':'A', 'C':'G', 'G':'C', 'Y':'R', 'R':'Y', 'W':'W', 'S':'S', 'K':'M', 'M':'K', 'D':'H', 'V':'B', 'H':'D', 'B':'V', 'N':'N'}
		else:
			transl = {'a':'t', 't':'a', 'c':'g', 'g':'c', 'y':'r', 'r':'y', 'w':'w', 's':'s', 'k':'m', 'm':'k', 'd':'h', 'v':'b', 'h':'d', 'b':'v', 'n':'n',
						'A':'T', 'T':'A', 'C':'G', 'G':'C', 'Y':'R', 'R':'Y', 'W':'W', 'S':'S', 'K':'M', 'M':'K', 'D':'H', 'V':'B', 'H':'D', 'B':'V', 'N':'N'}
		bases = [transl[base] for base in self.sequence]
		return type(self)(''.join(bases))


	def reverse_complement(self):
		"""
		Returns the reverse complement of a DNA or RNA string.
		"""
		return self.reverse().complement()


	def transcribe(self):
		"""
		Return RNA version of DNA
		"""
		transl = {'t':'u', 'T':'U'}

		transcript = []
		for base in self.sequence:
			if base in ['T', 't']:
				transcript.append(transl[base])
			else:
				transcript.append(base)

		return RNA(''.join(transcript))


	def reverse_transcribe(self):
		"""
		Return DNA version of RNA
		"""
		transl = {'u':'t', 'U':'T'}

		transcript = []
		for base in self.sequence:
			if base in ['U', 'u']:
				transcript.append(transl[base])
			else:
				transcript.append(base)

		return DNA(''.join(transcript))


	def translate(self, table=1):
		"""
		Returns protein sequence from DNA or RNA string.
		The table variable specifies which codon table should be used.
		table defaults to the standard codon table 1
		"""
		codons = CodonTable(table).get_codons()
		protein = []
		transcript = self.sequence.upper().transcribe()

		for i in range(0, len(transcript), 3):
			codon = transcript[i:(i+3)]
			if i+3>len(transcript): pass
			elif any(codon in s for s in codons['F']): protein.append('F')
			elif any(codon in s for s in codons['L']): protein.append('L')
			elif any(codon in s for s in codons['S']): protein.append('S')
			elif any(codon in s for s in codons['Y']): protein.append('Y')
			elif any(codon in s for s in codons['*']): protein.append('*')
			elif any(codon in s for s in codons['C']): protein.append('C')
			elif any(codon in s for s in codons['W']): protein.append('W')
			elif any(codon in s for s in codons['P']): protein.append('P')
			elif any(codon in s for s in codons['H']): protein.append('H')
			elif any(codon in s for s in codons['E']): protein.append('E')
			elif any(codon in s for s in codons['R']): protein.append('R')
			elif any(codon in s for s in codons['I']): protein.append('I')
			elif any(codon in s for s in codons['M']): protein.append('M')
			elif any(codon in s for s in codons['T']): protein.append('T')
			elif any(codon in s for s in codons['N']): protein.append('N')
			elif any(codon in s for s in codons['K']): protein.append('K')
			elif any(codon in s for s in codons['V']): protein.append('V')
			elif any(codon in s for s in codons['A']): protein.append('A')
			elif any(codon in s for s in codons['D']): protein.append('D')
			elif any(codon in s for s in codons['Q']): protein.append('Q')
			elif any(codon in s for s in codons['G']): protein.append('G')
			elif any(codon in s for s in codons['X']): protein.append('X')
			else: raise ValueError('Codon %s is not valid.' % codon)

		return Protein(''.join(protein))


	def translate_reverse_complement(self, table=1):
		'''Translate the reverse complement of DNA'''
		return self.reverse_complement().translate(table)


	def mass(self):
		"""
		Determine mass of single-stranded sequence in g/mol (da)
		"""
		if 'u' in self.alphabet:
			mass_vals ={'A':347.2, 'C':323.2, 'G':363.2, 'U':324.2}
		 else:
			mass_vals = {'A':331.2, 'C':307.2, 'G':347.2, 'T':322.2}

		 return sum([mass_vals[base.upper()] for base in self.sequence])


	def count_bases(self):
		'''
		Counts the number of different bases in a DNA sequence.
		Input should be a string comprising dna bases. These can be a, t, c, g or any of the ambiguous bases.
		Output is a dictionary with DNA base keys and integer values.
		'''
		base_count = {'G':0, 'A':0, 'T':0, 'C':0,
			'R':0, 'Y':0, 'W':0, 'S':0,
			'M':0, 'K':0, 'H':0, 'B':0,
			'P':0, 'H':0, 'V':0, 'D':0, 'N':0}
		return [base_count[s] += 1 for s in self.sequence.upper()]


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



## Common class to be used for ambDNA and ambRNA through sub-classing ##

class _AmbiguousNucleotideBaseClass(_BioSeq):

	def __init__(self, sequence, alphabet):
		self.alphabet = 'atcgATCG'
		_BioSeq.__init__(self, sequence, alphabet)


	def __combine(self, input_list, max_total=50000):
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


	def listupper(self, t):
		'''
		Capitalizes all strings in nested lists.
		'''
		if isinstance(t, list):
			return [listupper(s) for s in t]
		elif isinstance(t, str):
			return t.upper()
		else:
			return t


	def amb(self, in_list):
		'''
		This function finds the degenerate nucleotide for a list containing CATG nucleotides.
		Output is a single ambiguous DNA nucleotide as a string.
		Example input is: ['A','T','C','G']. The output for that input is 'N'
		'''
		in_list = self.listupper(in_list)

		if all([x in 'A' for x in in_list]):  output = 'A'
		elif all([x in 'G' for x in in_list]): output = 'G'
		elif all([x in 'C' for x in in_list]): output = 'C'
		elif all([x in 'T' for x in in_list]): output = 'T'
		elif all([x in 'CT' for x in in_list]): output = 'Y'
		elif all([x in 'GT' for x in in_list]): output = 'K'
		elif all([x in 'AC' for x in in_list]): output = 'M'
		elif all([x in 'CG' for x in in_list]): output = 'S'
		elif all([x in 'AT' for x in in_list]): output = 'W'
		elif all([x in 'AG' for x in in_list]): output = 'R'
		elif all([x in 'CTA' for x in in_list]): output = 'H'
		elif all([x in 'CAG' for x in in_list]): output = 'V'
		elif all([x in 'TAG' for x in in_list]): output = 'D'
		elif all([x in 'CTG' for x in in_list]): output = 'B'
		elif all([x in 'CTAG' for x in in_list]): output = 'N'
		else: raise ValueError('Error, input must be a list of standard GATC nucleotides.')
		return output


	def multi_amb(self, largelist):
		'''
		This function finds the degenerate nucleotide for a list of lists containing CATG nucleotides.
		Output is a DNA string with the ambigous nucleotides.
		Example input is: [['A','T','C','G'],['A','T','C','G'],['G','T']]. The output for that is 'NNK'
		'''
		output = []

		largelist = [x.upper() for x in largelist] #make the uppercase
		for lst in largelist:
			output.append(amb(lst))

		return ''.join(output)


	def un_amb(self, string):
		'''
		Converts an ambiguous nucleotide sequence to a list of sequences containing only A, T, C and G (as appropriate).
		'''
		assert type(string) is str, 'Error, the input has to be a string.'
		string = string.upper()

		pos_list = []
		for letter in string:
			assert letter in 'NMRWSYKVHDBGATC', 'Error, "%s" is not a valid ambigous nucleotide.'
			if 'A' == letter: pos_list.append(['A'])
			elif 'C' == letter: pos_list.append(['C'])
			elif 'T' == letter: pos_list.append(['T'])
			elif 'G' == letter: pos_list.append(['G'])
			elif 'M' == letter: pos_list.append(['A','C'])
			elif 'Y' == letter: pos_list.append(['C','T'])
			elif 'K' == letter: pos_list.append(['G','T'])
			elif 'S' == letter: pos_list.append(['C','G'])
			elif 'W' == letter: pos_list.append(['A','T'])
			elif 'R' == letter: pos_list.append(['A','G'])
			elif 'H' == letter: pos_list.append(['C','T','A'])
			elif 'V' == letter: pos_list.append(['C','A','G'])
			elif 'D' == letter: pos_list.append(['T','A','G'])
			elif 'B' == letter: pos_list.append(['C','T','G'])
			elif 'N' == letter: pos_list.append(['C','T','A','G'])

		return self.__combine(pos_list) #call combine function and return the result as a list of strings




	def common_nuc(nuc_list, greedy=False):
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

		if all(['A' in s for s in nuc_list]): output.append('A')
		if all(['G' in s for s in nuc_list]): output.append('G')
		if all(['C' in s for s in nuc_list]): output.append('C')
		if all(['T' in s for s in nuc_list]): output.append('T')
		if greedy is False and len(output)>0: return output

		if all(['C' in s or 'T' in s for s in nuc_list]): output.append('Y')
		if all(['G' in s or 'T' in s for s in nuc_list]): output.append('K')
		if all(['A' in s or 'C' in s for s in nuc_list]): output.append('M')
		if all(['C' in s or 'G' in s for s in nuc_list]): output.append('S')
		if all(['A' in s or 'T' in s for s in nuc_list]): output.append('W')
		if all(['A' in s or 'G' in s for s in nuc_list]): output.append('R')
		if greedy is False and len(output)>0: return output

		if all(['C' in s or 'T' in s or 'A' in s for s in nuc_list]): output.append('H')
		if all(['C' in s or 'A' in s or 'G' in s for s in nuc_list]): output.append('V')
		if all(['T' in s or 'A' in s or 'G' in s for s in nuc_list]): output.append('D')
		if all(['C' in s or 'T' in s or 'G' in s for s in nuc_list]): output.append('B')
		if greedy is False and len(output)>0: return output

		if all(['C' in s or 'T' in s or 'A' in s or 'G' in s for s in nuc_list]): output.append('N')
		return output


class DNA(_NucleotideBaseClass):

	def __init__(self, sequence):
		self.alphabet = 'atcgATCG'
		_NucleotideBaseClass.__init__(self, sequence, alphabet)
		self._mytype = 'DNA'


class ambDNA(_AmbiguousNucleotideBaseClass):

	def __init__(self, sequence):
		self.alphabet = 'gatcrywsmkhbvdnGATCRYWSMKHBVDN'
		_AmbiguousNucleotideBaseClass.__init__(self, sequence, alphabet)
		self._mytype = 'ambDNA'


class RNA(_NucleotideBaseClass):

	def __init__(self, sequence):
		self.alphabet = 'aucgAUCG'
		_NucleotideBaseClass.__init__(self, sequence)
		self._mytype = 'RNA'


class ambRNA(_AmbiguousNucleotideBaseClass):

	def __init__(self, sequence):
		self.alphabet = 'gaucrywsmkhbvdnGAUCRYWSMKHBVDN'
		_AmbiguousNucleotideBaseClass.__init__(self, sequence, alphabet)
		self._mytype = 'ambRNA'


class Protein(_BioSeq):
	def __init__(self, sequence):
		self.alphabet = 'flsycwpherimtnkvadqgxFLSYCWPHERIMTNKVADQG*X'
		_BioSeq.__init__(self, sequence)
		self._mytype = 'Protein'


	def get_codons(self, AA, table=1, separate=False):
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
		assert AA in self.alphabet, 'Error, %s is not a valid amino acid' % str(AA)

		codons = CodonTable(table).get_codons(separate)

		return codons[AA]


	def one_to_three(self, one_letter):
		'''
		Convert a one letter code amino acid to a three letter code.
		'''
		AA = {'I':'Ile', 'V':'Val',
		'L':'Leu', 'F':'Phe',
		'C':'Cys', 'M':'Met',
		'A':'Ala', 'G':'Gly',
		'T':'Thr', 'W':'Trp',
		'S':'Ser', 'Y':'Tyr',
		'P':'Pro', 'H':'His',
		'E':'Glu', 'Q':'Gln',
		'D':'Asp', 'N':'Asn',
		'K':'Lys', 'R':'Arg',
		'*':'***', 'X':'Unknown'}
		assert one_letter.upper() in AA.keys(), 'Error, %s is not a valid amino acid' % one_letter
		return AA[one_letter.upper()]


	def three_to_one(self, three_letter):
		'''
		Convert a three letter code amino acid to a one letter code.
		'''
		AA = {'ILE':'I', 'VAL':'V',
		'LEU':'L', 'PHE':'F',
		'CYS':'C', 'MET':'M',
		'ALA':'A', 'GLY':'G',
		'THR':'T', 'TRP':'W',
		'SER':'S', 'TYR':'Y',
		'PRO':'P', 'HIS':'H',
		'GLU':'E', 'GLN':'Q',
		'ASP':'D', 'ASN':'N',
		'LYS':'K', 'ARG':'R',
		'***':'*', 'UNKNOWN':'X'}
		assert three_letter.upper() in AA.keys(), 'Error, %s is not a valid amino acid' % three_letter
		return AA[three_letter.upper()]


	def one_to_full(self, one_letter):
		'''
		Convert one-letter amino acid code to full amino acid name.
		'''
		AA = {'F':'Phenylalanine', 'L':'Leucine',
		'S':'Serine', 'Y':'Tyrosine',
		'*':'Stop', 'C':'Cysteine',
		'W':'Tryptophan', 'P':'Proline',
		'H':'Histidine', 'Q':'Glutamine',
		'R':'Arginine', 'I':'Isoleucine',
		'M':'Methionine', 'T':'Threonine',
		'N':'Asparagine', 'K':'Lysine',
		'V':'Valine', 'A':'Alanine',
		'D':'Aspartic acid', 'E':'Glutamic acid',
		'G':'Glycine', 'X': 'Unknown'}
		assert one_letter.upper() in AA.keys(), 'Error, %s is not a valid amino acid' % one_letter
		return AA[one_letter.upper()]


	def full_to_one(self, full):
		'''
		Convert full amino acid name to one-letter amino acid code.
		'''
		AA = {'phenylalanine':'F', 'leucine':'L',
				'serine':'S', 'tyrosine':'Y',
				'stop':'*', 'cysteine':'C',
				'tryptophan':'W', 'proline':'P',
				'histidine':'H', 'glutamine':'Q',
				'arginine':'R', 'isoleucine':'I',
				'methionine':'M', 'threonine':'T',
				'asparagine':'N', 'lysine':'K',
				'valine':'V', 'alanine':'A',
				'aspartic acid':'D', 'glutamic acid':'E',
				'glycine':'G', 'unknown': 'X'} # to handle N as a nuclec acid
		assert full.lower() in AA.keys(), 'Error, %s is not a valid amino acid' % full
		return AA[full.lower()]


	def three_to_full(self, three_letter):
		'''
		Convert amino acid three letter code to full amino acid names.
		'''
		return one_to_full(three_to_one(three_letter))


	def full_to_three(self, full):
		'''
		Convert full amino acid names to three letter code.
		'''
		return one_to_three(full_to_one(full))


	def count_aa(self, seq):
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


	def reverse_translate(self, prot_seq, table=1):
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
			possible = self.get_codons(AA, table=table, separate=False)
			dna_seq.append(random.choice(possible))

		return RNA(''.join(dna_seq))


	def reverse_translate_ambiguous(self, table=1):
		'''
		Translate protein to RNA.
		The input is a protein sequence as a string.
		The output is a list of codons (with ambigous bases) that describe that protein.
		For some amino acids there will be two possible ambigous codons.
		Run the __combine() function to convert the list to all possible dna sequences.
		'''

		#### This one needs attention! #########

		dnalist = []
		for aa in self.sequence:
			dnalist.append(self.get_codons(aa, table=table, separate=False))
		return dnalist


	def mass(self):
		"""
		"""
		amino_acids = {'A':89, 'R':174,
			 'N':132, 'D':133,
			 'C':121, 'Q':146,
			 'E':147, 'G':75,
			 'H':155, 'I':131,
			 'L':131, 'K':146,
			 'M':149, 'F':165,
			 'P':115, 'S':105,
			 'T:'119, 'W':204,
			 'Y':181, 'V':117}

		 ## fix this ##




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
		self.__set_table(number) #get the specified codon table (returned as list of strings)
		self.__set_codons() #convert the codon table information to codons


	def __set_table(self, number):
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


	def __set_codons(self):
		'''
		Use a predetermined codon table to generate a dictionary of amino acids with their codons.
		Method is not intended for direct use.
		'''
		code, AAs, Starts, Base1, Base2, Base3 = self.get_table()
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

	def get_code(self):
		'''
		Return which genetic code is represented.
		The output is a string which specifies the code'
		Method is not intended for direct use.
		'''
		return self.code


	def get_table(self):
		'''
		Return the codon table data for the specified genetic code.
		The output is a list of strings.
		'''
		return self.table


	def get_codons(self, separate=False):
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


	def print_table(self):
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


	def print_codons(self):
		'''
		Print codons specified by codon table.
		'''
		codons = self.get_codons()
		print(('start = %s' %codons['start']))
		print(('stop  = %s' %codons['stop']))
		for aa in 'FLSYCWPHERIMTNKVADQG*':
			print(('%s     = %s' % (aa, codons[aa])))




######################
######################



class FreqTable(object):

	def __init__(self):
		pass

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
