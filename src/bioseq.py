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

import string
import random
from functools import reduce
import oligo_localizer
import peptide_localizer

#TODO
#finish protein class, What is left?
#use codon freq table to assign codons in a probabalistic way when reverse-translating protein.
#finish ambDNA class
#finish ambRNA class
#finish codon class

## Sub-classing of the str class and modification of some of the methods ##





	def mutate(self, mutationtype, mutationframe, mutation, silent=False):
		'''Mutates a given amino acid or DNA base.
			Mutationtype decides whether AA or DNA.
			Mutationframe decides whether in a certain feture or in DNA.
			Mutation is the actual mutation. This can be a list of mutations or a single mutation.
			Amino acid mutations are in the format: D121E(GAG) were the leading letter connotates the amino acid already present at position.
			The number is the amino acid number.
			The trailing letter is the amino acid that one whishes to introduce at the position.
			The three letters in the bracket is the codon which you wish to use to make the chosen mutation.
			The leading letter and the codon (within brackets) are optional.
			DNA mutataions are in the format: A234C
			The first letter is the base present at the chosen position.
			The numbers desgnate the chosen position.
			The trailing letter dessignates the base you want at that position.
			The leading letter is optional.
			If the 'silent' variable is False, a feature marking the mutation will be added. If True, no feature will be made.
			'''

		#mutation input can be a single mutation or a list of mutations
		if type(mutation) is list: #if it's a list, run the method on each of them
			for mut in mutation:
				self.mutate(mutationtype, mutationframe, mut)
		else:
			assert type(mutation) is str or type(mutation) is str, 'Error, input must be a string or unicode.'
			mutation = mutation.upper()
			if mutationframe == -1: #mutation frame of -1 means entire molecule
				complement = False
			else:
				complement = self.get_feature_complement(mutationframe) #find whether feature is reverse-complement or not

			if mutationtype == 'A': #if amino acid
				leadingAA = ''
				position = ''
				trailingAA = ''
				codon = ''

				#make sure input has the right pattern
				#matches the patterns '121E' or 'D121E' or '121E(GAA)' or 'D121E(GAA)'  (of course not explicitly, posotion, aa and codon can change)
				regular_expression = re.compile(r'''^							#match beginning of string
													[FLSYCWPHERIMTNKVADQG]? 	#zero or one occurances of amino acid
													[1234567890]+				#one or more digits
													[FLSYCWPHERIMTNKVADQG]{1}	#exactly one amino acid
													([(][ATCG]{3}[)])?			#zero or one occurances of brackets with codon inside
													$							#match end of string
													''', re.VERBOSE)
				assert regular_expression.match(mutation) != None, 'Error, the mutation %s is not a valid input.' % mutation


				#assumes the pattern D121E(GAG) (with leading letter and trailing bracket with codon being optional)
				for i in range(0,len(mutation)):
					if i == 0 and mutation[i] in 'FLSYCWPHERIMTNKVADQG': #leading AA if any
						leadingAA = mutation[i]
					elif mutation[i].isdigit(): #for position
						position += mutation[i]
					elif mutation[i] in 'ATCG' and type(position) is int: #getting the codon in brackets
						codon += mutation[i]
					elif mutation[i] in 'FLSYCWPHERIMTNKVADQG': #trailing AA
						trailingAA = mutation[i]
						position = int(position) #important to convert only after the second AA has been found

				#check that position does not exceed feature length
				if mutationframe == -1: #-1 for entire molecule
					length = len(self.GetDNA())/3
				else:
					length = len(self.GetFeatureDNA(mutationframe))/3
				assert position <= length, 'Error, the actual length of feature is %s AA long and is shorter than the specified position %s.' % (str(length), str(position))

				#check that the codon matches the specified amino acid
				if codon == '':
					#assign one...
					codon = dna.GetCodons(trailingAA)[0]
				assert trailingAA == dna.Translate(codon), 'Error, the specified codon %s does not encode the amino acid %s.' % (codon, trailingAA)

				#make sure the position has the AA that is specified in the leading letter
				global_position = self.FindAminoAcid(position, mutationframe) #get the global position (on enire dna) of the mutation
				if complement is True:
					positionAA = dna.TranslateRC(self.GetDNA(global_position[0][0], global_position[0][1])) #find the AA at that position
				if complement is False or mutationframe == -1:
					positionAA = dna.Translate(self.GetDNA(global_position[0][0], global_position[0][1])) #find the AA at that position
				if leadingAA != '':
					assert positionAA == leadingAA, 'Error, position %s has amino acid %s, and not %s as specified.' % (str(position),  positionAA, leadingAA)

				#now make the mutation and add corresponding feature if the 'silent' variable is False
				if complement is False or mutationframe == -1:
					self.changegbsequence(global_position[0][0], global_position[0][1], 'r', codon)
					if silent is False: self.add_feature(key='modified_base', qualifiers=['/note="%s%s%s(%s)"' % (positionAA, position, trailingAA, codon)], location=['%s..%s' % (global_position[0][0], global_position[0][1])], complement=False, join=False, order=False)
					print(('Mutation %s%s%s(%s) performed.' % (positionAA, position, trailingAA, codon)))
				elif complement is True:
					self.changegbsequence(global_position[0][0], global_position[0][1], 'r', dna.reverse_complement(codon))
					if silent is False: self.add_feature(key='modified_base', qualifiers=['/note="%s%s%s(%s)"' % (positionAA, position, trailingAA, codon)], location=['%s..%s' % (global_position[0][0], global_position[0][1])], complement=True, join=False, order=False)
					print(('Mutation %s%s%s(%s) performed.' % (positionAA, position, trailingAA, codon)))
				else:
					raise ValueError


			elif mutationtype == 'D': #if DNA
				leadingnucleotide = ''
				position = ''
				trailingnucleotide = ''

				#make sure input has the right pattern
				#matches the patterns '325T' or 'A325T' (of course not explicitly, position and nucleotide can change)
				regular_expression = re.compile(r'''^							#match beginning of string
													[ATCG]? 					#zero or one occurances of nucleotide
													[1234567890]+				#one or more digits
													[ATCG]{1}					#exactly one nucleotide
													$							#match end of string
													''', re.VERBOSE)
				assert regular_expression.match(mutation) != None, 'Error, the mutation %s is not a valid input.' % mutation


				#assumes the pattern 'A325T' (with leading nucleotide being optional)
				for i in range(0,len(mutation)):
					if i == 0 and mutation[i] in 'ATCG': #leading base if any
						leadingnucleotide = mutation[i]
					elif mutation[i].isdigit(): #for position
						position += mutation[i]
					elif mutation[i] in 'ATCG': #trailing nucleotide
						trailingnucleotide = mutation[i]
						position = int(position) #important to convert only after the second nucleotide has been found

				#check that position does not exceed feature length
				if mutationframe == -1: #-1 for entire molecule
					length = len(self.GetDNA())
				else:
					length = len(self.GetFeatureDNA(mutationframe))
				assert position <= length, 'Error, the actual length of feature is %s nucleotides long and is shorter than the specified position %s.' % (str(length), str(position))

				#make sure the position has the nucleotide that is specified in the leading letter
				global_position = self.FindNucleotide(position, mutationframe) #get the global position (on enire dna) of the mutation
				if complement is True:
					positionnucleotide = dna.reverse_complement(self.GetDNA(global_position[0][0], global_position[0][1])).upper()
				elif complement is False:
					positionnucleotide = self.GetDNA(global_position[0][0], global_position[0][1]).upper()
				if leadingnucleotide != '':
					assert positionnucleotide == leadingnucleotide, 'Error, position %s has the nucleotide %s, and not %s as specified.' % (str(position),  positionnucleotide, leadingnucleotide)

				#now make the mutation and add corresponding feature if the 'silent' variable is False
				if complement is False or mutationframe == -1:
					self.changegbsequence(global_position[0][0], global_position[0][1], 'r', trailingnucleotide)
					if silent is False: self.add_feature(key='modified_base', qualifiers=['/note=%s%s%s' % (positionnucleotide, position, trailingnucleotide)], location=['%s..%s' % (global_position[0][0], global_position[0][1])], complement=False, join=False, order=False)
					print(('Mutation %s%s%s performed.' % (positionnucleotide, position, trailingnucleotide)))
				elif complement is True:
					self.changegbsequence(global_position[0][0], global_position[0][1], 'r', dna.reverse_complement(trailingnucleotide))
					if silent is False: self.add_feature(key='modified_base', qualifiers=['/note=%s%s%s' % (positionnucleotide, position, trailingnucleotide)], location=['%s..%s' % (global_position[0][0], global_position[0][1])], complement=True, join=False, order=False)
					print(('Mutation %s%s%s performed.' % (positionnucleotide, position, trailingnucleotide)))
				else:
					raise ValueError














	def find_nuc(self, searchstring, searchframe=-1, searchRC=False):
		'''Method for finding a certain DNA sequence. Degenerate codons are supported.
		Searchstring should be a integer number or string of DNA characers.
		searchRC is True or False and determines whether the reverse complement should also be searched.
		Searchframe should be an index for a feature. -1 indicates search in entire genbank file
		indeces >-1 are feature indeces'''


		### Need to implement the searchRC variable

		#fix the set_dna_selection functions here
		assert type(searchstring) == str or type(searchstring) == str or type(searchstring) == int, 'Error, search takes a string of DNA or a string of numbers as input.'
		assert type(searchframe) == int and -1<=searchframe<len(self.get_all_features()), 'Error, %s is not a valid argument for searchframe.' % str(searchframe)
		assert type(searchRC) == bool, 'Error, searchRC must be True or False'
		#empty search string
		if searchstring=='':
			print('The searchstring is missing, please check your input')
			return []

		#searching for a position (by number)
		elif type(searchstring) is int: #if search is numbers only
			complement = self.get_feature_complement(searchframe) # is feature complement or not
			if searchframe == -1: #if index is -1, that means search in molecule
				search_hits = [(int(searchstring), int(searchstring))]
				return search_hits

			else: #otherwise search in feature indicated by the index
				feature = self.get_feature(searchframe)
				locations = copy.deepcopy(feature['location']) #get list of all locations
				cleaned_locations = []
				for i in range(len(locations)): # remove all < > and such..
					cleaned_locations.append(self.get_location(locations[i]))

				gaps = [] #for storing gap sizes (between locations)
				for i in range(0,len(cleaned_locations)-1):
					gaps.append(cleaned_locations[i+1][0] - cleaned_locations[i][1]) #subtract end of last location from the beginnning of current

				#now find where the searchstring is located
				if complement == False:
					for i in range(0,len(cleaned_locations)):
						if i == 0:
							start, finish = cleaned_locations[i]
							searchstring += start-1
						else:
							start, finish = cleaned_locations[i]
							searchstring += gaps[i-1]-1
						if start<=searchstring<=finish: #if the chosen number is in the range of the current location
							return [(searchstring,searchstring)]
							break
					print('No matches were found')
					return []

				elif complement == True:
					cleaned_locations = cleaned_locations[::-1] #reverse location list
					for i in range(len(cleaned_locations)):
						cleaned_locations[i] = cleaned_locations[i][::-1]
					gaps = gaps[::-1]
					for i in range(0,len(cleaned_locations)):
						if i == 0:
							start, finish = cleaned_locations[i]
							searchstring = -searchstring
							searchstring += start+1
						else:
							start, finish = cleaned_locations[i]
							searchstring -= gaps[i-1]-1

						#is the chosen number is in the range of the current location
						if start>=searchstring>=finish:
							return [(searchstring,searchstring)]
							break
					return []

		#searching for string, not numbers
		else:
			complement = self.get_feature_complement(searchframe) # is feature complement or not
			#need to test that the characters are valid
			if searchframe == -1: #if index is -1, that means search in molecule
				dna_seq = self.GetDNA()
				search_hits = []
				for match in oligo_localizer.match_oligo(dna_seq, searchstring):
					search_hits.append((match[0], match[1]))
				return search_hits

			else:
				#strategy is to find location of the sting inside the feature DNA and then to map those onto the whole molecule
				DNA = self.GetFeatureDNA(searchframe)
				search_hits = []
				for match in oligo_localizer.match_oligo(DNA, searchstring):
					search_hits.append((match[0], match[1]))

				#now map them one by one to the whole molecule
				re_mapped_search_hits = []
				for i in range(len(search_hits)):
					start = self.FindNucleotide(search_hits[i][0], searchframe)[0][0]
					finish = self.FindNucleotide(search_hits[i][1], searchframe)[0][0]
					re_mapped_search_hits.append((start,finish))

				#if feature is on complement then I need to reverse the list of hits and the tuples inside
				if complement == True:
					re_mapped_search_hits = re_mapped_search_hits[::-1]
					for i in range(len(re_mapped_search_hits)):
						re_mapped_search_hits[i] = re_mapped_search_hits[i][::-1]
				return re_mapped_search_hits  ## need to fix this so that hits spanning a gap get correctly colored #######


	def find_aa(self, searchstring, searchframe, searchRC=False):
		'''Method for finding a certain protein sequence, or position, in the file. Degenerate codons are supported'''
		assert type(searchstring) == str or type(searchstring) == str or type(searchstring) == int, 'Error, search takes a string of DNA or a string of numbers as input.'
		assert type(searchframe) == int and -1<=searchframe<len(self.get_all_features()), 'Error, %s is not a valid argument for searchframe.' % str(searchframe)
		assert type(searchRC) == bool, 'Error, searchRC must be True or False'


		#empty search string
		if searchstring=='':
			print('The searchstring is missing, please check your input')
			return []

		#searching for a position (by number)
		elif type(searchstring) is int: #if search is numbers only
			#get the dna triplet positions for the amino acid
			start = searchstring*3 -2
			finish = searchstring*3
			if searchframe == -1: #if index is -1, that means search in molecule
				search_hits = [(start, finish)]
				return search_hits

			else: #otherwise search in feature indicated by the index
				start = self.FindNucleotide(start, searchframe)[0][0]
				finish = self.FindNucleotide(finish, searchframe)[0][0]

				search_hits = [(start, finish)]
				search_hits[0] = tuple(sorted(search_hits[0]))
				return search_hits

		#searching for an amino acid sequence
		else:
			#strategy is to find location of the sting inside the feature DNA and then to map those onto the whole molecule
			complement = self.get_feature_complement(searchframe) # is feature complement or not
			if searchframe == -1: #if index is -1, that means search in molecule
				DNA = self.GetDNA()
				protein = dna.Translate(DNA)

				#find hits on protein level
				search_hits = []
				for match in peptide_localizer.match_peptide(protein, searchstring):
					search_hits.append((match[0],match[1]))
#				print('hits', search_hits)

				#now map them one by one to the whole molecule
				re_mapped_search_hits = []
				for i in range(len(search_hits)):
					start = self.FindAminoAcid(search_hits[i][0], searchframe)[0][0]
					finish = self.FindAminoAcid(search_hits[i][1], searchframe)[0][0]+2
					re_mapped_search_hits.append((start,finish))
#				print('re-mapped', re_mapped_search_hits)
				return re_mapped_search_hits

			else:
				DNA = self.GetFeatureDNA(searchframe)
				protein = dna.Translate(DNA)

				search_hits = []
				for match in peptide_localizer.match_peptide(protein, searchstring):
					search_hits.append((match[0],match[1]))
#				print('protein hits', search_hits)

				#now map them one by one to the whole molecule
				re_mapped_search_hits = []
				for i in range(len(search_hits)):
					if complement == False:
						start = self.FindAminoAcid(search_hits[i][0], searchframe)[0][0]
						finish = self.FindAminoAcid(search_hits[i][1], searchframe)[0][0]+2
					elif complement == True:
						start = self.FindAminoAcid(search_hits[i][0], searchframe)[0][0]
						finish = self.FindAminoAcid(search_hits[i][1], searchframe)[0][0]+2
					re_mapped_search_hits.append((start,finish))
#				print('re-mapped', re_mapped_search_hits)
				return re_mapped_search_hits  ## need to fix this so that hits spanning a gap get correctly colored #######




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

	def lower(self, start=None, end=None):
		"""
		"""
		self.subsequence(start, end)
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

	def upper(self, start=None, end=None):
		"""
		"""
		self.subsequence(start, end)
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

	def length(self, start=None, end=None):
		"""
		"""
		self.subsequence(start, end)
		return length(self.sequence)


	def randomize(self, start=None, end=None):
		'''
		Randomize a given DNA, RNA or protein sequence.
		'''
		self.subsequence(start, end)
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


	def assert_section(self, start, end):
		'''
		Make sure sequence selection is ok
		'''
		assert (type(start) == int and type(end) == int), 'Function requires two integers.'
		assert start <= end, 'Startingpoint must be before end'


	def subsequence(self, start, end):
		'''
		Get sub-sequence of the current one.
		'''
		if start is None or end is None:
			pass
		else:
			self.assert_selection(start, end)
			self.sequence = self.sequence[start-1, end-1]


	def delete(self, start, end):
		'''Deletes portion of sequence'''
		if start is None or end is None:
			pass
		else:
			self.assert_selection(start, end)
			self.sequence = self.sequence[1, start-1] + self.sequence[end-1, self.length()]


	def insert(self, pos, seq):
		'''
		Insert sequence at specified position
		'''
		seq = type(self)(seq)
		self.sequence = self.sequence[1, pos-1] + seq + self.sequence[pos-1, self.length()]


	def type(self):
		"""
		What type sequence is it
		"""
		return self._mytype


## Common class to be used for DNA and RNA through sub-classing ##

class _NucleotideBaseClass(_BioSeq):

	def __init__(self, sequence, alphabet):
		_BioSeq.__init__(self, sequence, alphabet=alphabet)


	def reverse(self, start=None, end=None):
		"""
		Returns the reverse of a DNA or RNA string.
		"""
		self.subsequence(start, end)
		return type(self)(self.sequence[::-1])  #makes the reverse of the input string


	def complement(self, start=None, end=None):
		"""
		Returns the complement of a DNA or RNA string.
		"""
		self.subsequence(start, end)
		if 'u' in self.alphabet:
			transl = {'a':'u', 'u':'a', 'c':'g', 'g':'c', 'y':'r', 'r':'y', 'w':'w', 's':'s', 'k':'m', 'm':'k', 'd':'h', 'v':'b', 'h':'d', 'b':'v', 'n':'n',
						'A':'U', 'U':'A', 'C':'G', 'G':'C', 'Y':'R', 'R':'Y', 'W':'W', 'S':'S', 'K':'M', 'M':'K', 'D':'H', 'V':'B', 'H':'D', 'B':'V', 'N':'N'}
		else:
			transl = {'a':'t', 't':'a', 'c':'g', 'g':'c', 'y':'r', 'r':'y', 'w':'w', 's':'s', 'k':'m', 'm':'k', 'd':'h', 'v':'b', 'h':'d', 'b':'v', 'n':'n',
						'A':'T', 'T':'A', 'C':'G', 'G':'C', 'Y':'R', 'R':'Y', 'W':'W', 'S':'S', 'K':'M', 'M':'K', 'D':'H', 'V':'B', 'H':'D', 'B':'V', 'N':'N'}
		bases = [transl[base] for base in self.sequence]
		return type(self)(''.join(bases))


	def reverse_complement(self, start=None, end=None):
		"""
		Returns the reverse complement of a DNA or RNA string.
		"""
		self.subsequence(start, end)
		return self.reverse().complement()


	def transcribe(self, start=None, end=None):
		"""
		Return RNA version of DNA
		"""
		self.subsequence(start, end)
		transl = {'t':'u', 'T':'U'}

		transcript = []
		for base in self.sequence:
			if base in ['T', 't']:
				transcript.append(transl[base])
			else:
				transcript.append(base)

		return RNA(''.join(transcript))


	def reverse_transcribe(self, start=None, end=None):
		"""
		Return DNA version of RNA
		"""
		self.subsequence(start, end)
		transl = {'u':'t', 'U':'T'}

		transcript = []
		for base in self.sequence:
			if base in ['U', 'u']:
				transcript.append(transl[base])
			else:
				transcript.append(base)

		return DNA(''.join(transcript))


	def translate(self, table=1, start=None, end=None):
		"""
		Returns protein sequence from DNA or RNA string.
		The table variable specifies which codon table should be used.
		table defaults to the standard codon table 1
		"""
		self.subsequence(start, end)
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


	def translate_reverse_complement(self, table=1, start=None, end=None):
		'''Translate the reverse complement of DNA'''
		self.subsequence(start, end)
		return self.reverse_complement().translate(table)


	def mass(self, start=None, end=None):
		"""
		Determine mass of single-stranded sequence in g/mol (da)
		"""
		self.subsequence(start, end)
		if 'u' in self.alphabet:
			mass_vals ={'A':347.2, 'C':323.2, 'G':363.2, 'U':324.2}
		else:
			mass_vals = {'A':331.2, 'C':307.2, 'G':347.2, 'T':322.2}
		return sum([mass_vals[base.upper()] for base in self.sequence])


	def count_bases(self, start=None, end=None):
		'''
		Counts the number of different bases in a DNA sequence.
		Input should be a string comprising dna bases. These can be a, t, c, g or any of the ambiguous bases.
		Output is a dictionary with DNA base keys and integer values.
		'''
		self.subsequence(start, end)
		base_count = {'G':0, 'A':0, 'T':0, 'C':0,
			'R':0, 'Y':0, 'W':0, 'S':0,
			'M':0, 'K':0, 'H':0, 'B':0,
			'P':0, 'H':0, 'V':0, 'D':0, 'N':0}
		for base in self.sequence.upper():
			base_count[base] += 1
		return base_count


	def count_codons(self, start=None, end=None):
		'''
		Counts codons in a DNA sequence.
		Input should be a string comprising whole codons.
		No ambiguous codons allowed.
		Output is a dictionary with codon keys and integer values.
		'''
		self.subsequence(start, end)

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
		_BioSeq.__init__(self, sequence, alphabet=alphabet)


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
		_NucleotideBaseClass.__init__(self, sequence, alphabet='atcgATCG')
		self._mytype = 'DNA'


class ambDNA(_AmbiguousNucleotideBaseClass):

	def __init__(self, sequence):
		_AmbiguousNucleotideBaseClass.__init__(self, sequence, alphabet='gatcrywsmkhbvdnGATCRYWSMKHBVDN')
		self._mytype = 'ambDNA'


class RNA(_NucleotideBaseClass):

	def __init__(self, sequence):
		_NucleotideBaseClass.__init__(self, sequence, alphabet='aucgAUCG')
		self._mytype = 'RNA'


class ambRNA(_AmbiguousNucleotideBaseClass):

	def __init__(self, sequence):
		_AmbiguousNucleotideBaseClass.__init__(self, sequence, alphabet='gaucrywsmkhbvdnGAUCRYWSMKHBVDN')
		self._mytype = 'ambRNA'


class Protein(_BioSeq):
	def __init__(self, sequence):
		_BioSeq.__init__(self, sequence, alphabet='flsycwpherimtnkvadqgxFLSYCWPHERIMTNKVADQG*X')
		self._mytype = 'Protein'


	def one_to_three(self, one_letter, start=None, end=None):
		'''
		Convert a one letter code amino acid to a three letter code.
		'''
		self.subsequence(start, end)
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


	def three_to_one(self, three_letter, start=None, end=None):
		'''
		Convert a three letter code amino acid to a one letter code.
		'''
		self.subsequence(start, end)
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


	def one_to_full(self, one_letter, start=None, end=None):
		'''
		Convert one-letter amino acid code to full amino acid name.
		'''
		self.subsequence(start, end)
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


	def full_to_one(self, full, start=None, end=None):
		'''
		Convert full amino acid name to one-letter amino acid code.
		'''
		self.subsequence(start, end)
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


	def three_to_full(self, three_letter, start=None, end=None):
		'''
		Convert amino acid three letter code to full amino acid names.
		'''
		self.subsequence(start, end)
		return one_to_full(three_to_one(three_letter))


	def full_to_three(self, full, start=None, end=None):
		'''
		Convert full amino acid names to three letter code.
		'''
		self.subsequence(start, end)
		return one_to_three(full_to_one(full))


	def count_aa(self, start=None, end=None):
		'''
		Count occurrences of all amino acids in sequence.
		The X character for unknown amino acid is allowed.
		Return as dictionary.
		'''
		self.subsequence(start, end)
		aa_count = {'I':0, 'V':0, 'L':0, 'F':0,
		'C':0, 'M':0, 'A':0, 'G':0,
		'T':0, 'W':0, 'S':0, 'Y':0,
		'P':0, 'H':0, 'E':0, 'Q':0,
		'D':0, 'N':0, 'K':0, 'R':0,
		'*':0, 'X':0}
		for aa in self.sequence.upper():
			aa_count[aa] += 1
		return aa_count


	def reverse_translate(self, table=1, start=None, end=None, freq_table=None):
		'''
		Reverse translates protein sequence to DNA sequence using the specified codon table.
		For each amino acid the DNA codon is chosen randomly from those that encode the specified amino acid.
		table defaults to the standard codon table 1

		Input is a protein sequence in one-letter code as a string.
		Output is a DNA sequence as a string.
		table defaults to the standard codon table 1
		'''
		self.subsequence(start, end)
		rna_seq = []
		for AA in prot_seq:
			possible = CodonTable(table).get_codons(separate=False)[AA]
			if freq_table is None:
				rna_seq.append(random.choice(possible))
			else:
				raise NotImplementedError
				#use frequency table to weight probability of codon selection

		return RNA(''.join(rna_seq))


	def reverse_translate_ambiguous(self, start=None, end=None, table=1):
		'''
		Translate protein to RNA.
		The input is a protein sequence as a string.
		The output is a list of codons (with ambigous bases) that describe that protein.
		For some amino acids there will be two possible ambigous codons.
		Run the __combine() function to convert the list to all possible dna sequences.
		'''
		self.subsequence(start, end)
		#### This one needs attention! #########

		rnalist = []
		for aa in self.sequence:
			rnalist.append(CodonTable(table).get_codons(separate=False)[AA])
		return rnalist


	def mass(self, start=None, end=None):
		"""
		"""
		self.subsequence(start, end)
		mass_vals = {'A':89, 'R':174,
			 'N':132, 'D':133,
			 'C':121, 'Q':146,
			 'E':147, 'G':75,
			 'H':155, 'I':131,
			 'L':131, 'K':146,
			 'M':149, 'F':165,
			 'P':115, 'S':105,
			 'T':119, 'W':204,
			 'Y':181, 'V':117,
			 '*':0, 'X':110} #110 is the average amino acid weight
		return sum([mass_vals[aa.upper()] for aa in self.sequence])




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
				codons['start'].append(RNA(codon))

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
		for aa in self.alphabet:
			print(('%s     = %s' % (aa, codons[aa])))






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
