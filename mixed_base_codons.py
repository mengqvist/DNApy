#!/usr/bin/env python


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

from copy import deepcopy
import dna
import re
	
	
#TODO 
#save all the codons that give the same number of off-targets, then compare them in regard to how many times they code each AA

def getinput():
	'''
	This method gets user input if the script is run directly from the command line.
	First a codon table is chosen by entering a valid number. 
	Amino acids are then entered one by one as one-letter code. 
	'''
	###This is where I decide which AA I want
	DesiredAA = [];
	AllNaturalAA = ['F','L','S','Y','C','W','P','H','E','R','I','M','T','N','K','V','A','D','Q','G','*'];
	AA = 'x'
	table = raw_input('Input a number to indicate which codon table to use (standard is 1): ')
	if table is not int:
		table = 1
	while AA != '':
		AA = raw_input('Input single AA in single letter code. If done, press enter: ')
		AA = AA.upper()	
		if AA == '':
			pass	
		elif AA in AllNaturalAA:	
			DesiredAA.append(AA)	
		else:
			print('This is not a valid AA')
	DesiredAA = sorted(list(set(DesiredAA))) #condense list to what is unique
	#print(DesiredAA)
	return DesiredAA, table
	


class AmbigousCodon:
	'''
	Class that holds methods and values for computing the ambiguous codon for a list of amino acids. 
	Alternatively, the class can evaluate an ambigous codon which is provided by the user. 
	Required input is an integer that determines the codon table to use and either a list of desired amino acids in single letter code 
	OR a three-letter codon using the IUPAC Nucleotide ambiguity code (G, A, T, C, R, Y, W, S, M, K, H, B, V, D, N).
	
	The algorithm works as follows:
	If a list of amino acids (as single letter code) is passed to the algorithm;
	All possible regular codons for those amino acids are looked up and returned as a list of lists. -> get_triplets()
	All nucleotides for first, second and third position extracted into separate lists while retaining list structure. -> sumupcodons()
	Separately, for pos 1, 2 and 3, a degenerate nucleotide is found that matches at least one nucleotide in each lists. -> commonNuc()
	These are concatenated and represents the first degenerate triplet.
	
	The degenerate triplet is converted back to real (GATC) codons xxxxx()
	translate to amino acids and check against what you wanted in the first place xxxxxxxx()

	Add more !!!!!!!!
	
	If an ambiguous codon is passed to the algorithm;
	...
	...
	...
	
	
	Regardless of which of the two ways are used to generate the codon object the methods to retrieve information are the same:

	To get the ambiguous codon (as a three-character string):
	codon_object.getTriplet()
	
	To get the target amino acids (as a list of upper case amino acids in single letter code):
	codon_object.getTarget()
	
	To get the off-target amino acids (as a list of upper case amino acids in single letter code):
	codon_object.getOfftarget()

	To get which amino acids can still be chosen without further off-target amino acids (as a list of upper case amino acids in single letter code):	
	codon_object.getPossible()
	
	To get which codon table was used for the computation (as an integer):
	codon_object.getTable()
	'''
	
	def __init__(self, input, table):
		self.setTable(table)
		
		#input can be either a three-nucleotide string or a list of amino acids
		if len(input) == 3 and type(input) == str: #if string i.e. an ambiguous codon
			input = input.upper()
			self.evaluateTriplet(input)
		elif type(input) == list: #if list, i.e. a list of amino acids to evaluate
			input = [s.upper() for s in input]
			self.setTarget(input)
		else:
			raise ValueError, 'The input is not valid'
		return

		
	###############################################################
	#######  Public methods intended for user interaction  ########
	###############################################################
	
	def getTarget(self):
		'''
		Retrieves a list of the target amino acids.
		Output is a list of upper case amino acids in single letter code.
		'''
		return self.target
		
	def getOfftarget(self):
		'''
		Retrieves a list of the off-target amino acids.
		Output is a list of upper case amino acids in single letter code.
		'''
		return self.offtarget
		
	def getPossible(self):
		'''
		Retrieves a list of amino acids still possible without further off-targets.
		Output is a list of upper case amino acids in single letter code.
		'''
		return self.possible
		
	def getTriplet(self):
		'''
		Retrieves the ambiguous codon.
		Output is a three-letter string of upper case characters.
		'''
		return self.triplet	
		
	def getTable(self):
		'''
		Retrieves which genetic code was used.
		Output is an integer.
		'''
		return self.table
		


	################################################################
	######## Methods NOT intended for direct user interaction ######
	################################################################
	
	def sumupcodons(self, DesiredCodons):
		'''
		Takes a list of regular codon lists and does two things; first it splits them up based on the first, second and third position of the triplet,
		then it checks which nucleotide (ambiguous allowed) that matches at least one of the nucleotides from each amino acid.	
		
		Example input, where alanine, cysteine and tyrosine are desired is: [['GCT', 'GCC', 'GCA', 'GCG'], ['TGT', 'TGC'], ['TAT', 'TAC']]
		The objective is to keep the list structure but to make three separate list where each 
		holds all the unique nucleotides for either the first, second or third base of the codons.
		The correct output for the first position would be: [['G'], ['T'], ['T']],
		for the second [['C'], ['G'], ['A']],
		and for the third [['T', 'C', 'A', 'G'], ['T', 'C'], ['T', 'C']].
		These lists can then be passed on to another function for finding a nucleotide that matches at least one of the nucleotides from each amino acid.	
		'''
		allcodon1 = []
		allcodon2 = []
		allcodon3 = []
		for entry in DesiredCodons:
			codon1 = []
			codon2 = []
			codon3 = []
			for codon in entry: #splits up the codons of each AA into triplets
				if codon[0] not in codon1:
					codon1.extend(codon[0]) #Takes first base in the triplet
				if codon[1] not in codon2:
					codon2.extend(codon[1]) #Takes the second base in triplet
				if codon[2] not in codon3:
					codon3.extend(codon[2]) #Takes third base in triplet
			allcodon1.append(codon1) #Adds up all the first bases
			allcodon2.append(codon2) #Adds up all the second bases
			allcodon3.append(codon3) #Adds up all the third bases
		return (allcodon1, allcodon2, allcodon3)


	def extra_list_elements(self, list_A, list_B): 
		'''
		Method for comparing two lists to find which elements are not present in both.
		The elements which are not common to both lists are returned as a list.
		Duplicates are ignored such that if list A has the integer 15 twice and list B only once, this will NOT result in 15 being returned.
		
		This method is used to check which off-target amino acids there are.
		'''
		not_in_B = [s for s in list_A if s not in list_B]
		not_in_A = [s for s in list_B if s not in list_A]
		not_in_both = list(set(not_in_A + not_in_B))
		return not_in_both
	
	
	def flatten_codon_list(self, codon_list):
		'''
		Some of the codons are list of lists (happens when the amino acid has codons at different parts of the codon circle).
		This method flattens the structure one level such that a list of codons remains.
		'''
		output = []
		for pos in codon_list:
			output_len = len(output)
			if isinstance(pos[0], str):
				if output_len == 0:
					output.append([pos])
				else:
					for o in range(output_len):
						output[o].append(pos)
					
			elif isinstance(pos[0], list):
				if output_len == 0:
					for p in pos:
						output.append([p])
				else:
					output.extend(deepcopy(output * (len(pos)-1)))
					for i in range(len(pos)):
						for j in range(output_len):
							output[j+i*output_len].append(pos[i])
		return output

		
	def find_degenerate(self, AA_list):
		'''
		Method for finding an ambiguous codon encoding a list of desired amino acids.
		The method finds the codon with fewest off-target amino acids.
		The input is a list of upper case amino acids in the single-letter code.
		The valid values are: FLSYCWPHERIMTNKVADQG*		
		
		The output is a tuple of the best ambiguous codon and the off-target amino acids.
		The ambiguous codon is a string of three of the following characters: GATCRYWSMKHBVDN
		The off-target amino acids is a list of upper case amino acids in single letter code.
		'''
		#make sure input is OK
		assert all([s in 'FLSYCWPHERIMTNKVADQG*' for s in AA_list]), 'Error, one or more of the amino acids %s are not valid.' % AA_list
		
		#get all codons for chosen amino acids
		regular_triplets = [dna.GetCodons(aa, table=self.getTable(), separate=True) for aa in AA_list]
		
		#some of the codons are list of lists (happens when the amino acid has codons at different parts of the codon circle)
		#I need to flatten this into separate lists with which go on further
		regular_triplets = self.flatten_codon_list(regular_triplets)
		best_score = None
		for codon_list in regular_triplets:
			#get all nucleotides for first, second and third position while retaining list structure		
			first, second, third = self.sumupcodons(codon_list) 
			
			#check which degenerate nucleotide can be used to find at least one match in each of the lists
			possible_triplets = dna.combine([dna.commonNuc(first), dna.commonNuc(second), dna.commonNuc(third)])
			
			#now go through them and see which is best
			for triplet in possible_triplets:
				#convert the triplet back to a list of real codons 
				Realcodons = dna.combine([dna.UnAmb(triplet[0]), dna.UnAmb(triplet[1]), dna.UnAmb(triplet[2])]) #condense the different codons for position 1, 2, and 3 to a list of triplets
			
				#Check which AA these codons code for
				ResultingAA = [dna.Translate(codon, table=self.getTable()) for codon in Realcodons]

				#compare which amino acids were desired with the ones resulting from the degenerate codon
				offtarget = sorted(self.extra_list_elements(AA_list, ResultingAA))
				
				#if there are fewer off-target amino acids with the new codon, keep it
				if len(offtarget) < best_score or best_score == None:
					best_score = len(offtarget)
					good_triplets = []
					good_triplets.append(triplet)
				elif len(offtarget) == best_score:
					good_triplets.append(triplet)
		
		#the saved triplets all have the same number of off-target amino acids, now keep the one with the lowest number of codons for each AA
		best_triplet = None #for storing best degenerate triplet
		best_offtarget = None #for storing the off-target AA of the best triplet
		best_score = None #for storing the length of the off-target list
		for triplet in good_triplets:
			#convert the triplet back to a list of real codons 
			Realcodons = dna.combine([dna.UnAmb(triplet[0]), dna.UnAmb(triplet[1]), dna.UnAmb(triplet[2])]) #condense the different codons for position 1, 2, and 3 to a list of triplets
		
			#save the stats in case there are fewer codons
			if len(Realcodons) < best_score or best_score == None:
				#Check which AA these codons code for
				ResultingAA = [dna.Translate(codon, table=self.getTable()) for codon in Realcodons]

				#compare which amino acids were desired with the ones resulting from the degenerate codon
				offtarget = sorted(self.extra_list_elements(AA_list, ResultingAA))
				
				#save stats
				best_score = len(Realcodons)
				best_triplet = triplet
				best_offtarget = offtarget
				
		return best_triplet, best_offtarget

	
	def next_steps(self):
		'''
		Method for finding which other amino acids can be selected without introducing 
		further (more in number) off-target ones.
		The output is a list of upper case amino acids in single letter code.
		'''

		possibleAA = []
		targetAA = self.getTarget()
		unusedAA = self.extra_list_elements(targetAA, list('FLSYCWPHERIMTNKVADQG*'))
		
		if len(unusedAA)>0:
			for AA in unusedAA:
				temptargetAA = targetAA[:]
				temptargetAA.append(AA)
				triplet, offtarget = self.find_degenerate(temptargetAA)
				if len(offtarget) <= len(self.getOfftarget()):
					possibleAA.append(AA)

		return sorted(possibleAA)


	def evaluateTriplet(self, amb_codon):
		'''
		Evaluate the ambiguous codon by computing which amino acids it codes for.
		The input is a string, three letters long and comprising only IUPAC Nucleotide ambiguity code.
		The valid values is any combination of three of the following: GATCRYWSMKHBVDN		
		'''
		#make sure input is OK
		assert type(amb_codon) is str and len(amb_codon) == 3, 'Error, the ambiguous codon must be a string three characters long.'
		m = re.match('^[GATCRYWSMKHBVDN]{3}$', amb_codon)
		assert m != None, 'Error, the codon %s is not valid. It may only use the chracters GATCRYWSMKHBVDN.' % amb_codon
		
		#compute target amino acids and set variables
		self.target = list(set([dna.Translate(s, self.getTable()) for s in dna.UnAmb(amb_codon)]))
		self.setTriplet(amb_codon)
		self.setOfftarget([])

		#see which other amino acids are possible without further off-target
		self.setPossible(self.next_steps())

		
	def setTarget(self, AA_list):
		'''
		Set which target amino acids are desired.
		The input is a list of upper case amino acids in the single-letter code.	
		The valid values are: FLSYCWPHERIMTNKVADQG*	
		'''
		#make sure input is OK
		assert all([s in 'FLSYCWPHERIMTNKVADQG*' for s in AA_list]), 'Error, one or more of the amino acids %s are not valid.' % AA_list
		self.target = AA_list

		#compute the triplet and the off-target AA
		triplet, offtarget = self.find_degenerate(self.getTarget())
		self.setTriplet(triplet)
		self.setOfftarget(offtarget)

		#see which other amino acids are possible without further off-target
		self.setPossible(self.next_steps())

	
	def setOfftarget(self, AA_list):
		'''
		Set which amino acids can still be chosen in the next step without further off-target amino acids.
		The input is a list of upper case amino acids in the single-letter code.	
		The valid values are: FLSYCWPHERIMTNKVADQG*
		'''
		assert all([s in 'FLSYCWPHERIMTNKVADQG*' for s in AA_list]), 'Error, one or more of the amino acids %s are not valid.' % AA_list
		self.offtarget = AA_list

		
	def setPossible(self, AA_list):
		'''
		Set which amino acids can still be chosen in the next step without further off-target amino acids.
		The input is a list of upper case amino acids in the single-letter code.
		The valid values are: FLSYCWPHERIMTNKVADQG*
		'''
		assert all([s in 'FLSYCWPHERIMTNKVADQG*' for s in AA_list]), 'Error, one or more of the amino acids %s are not valid.' % AA_list		
		self.possible = AA_list
	
	
	def setTriplet(self, amb_codon):
		'''
		Set the ambiguous codon.
		The input is a string, three letters long and comprising only IUPAC Nucleotide ambiguity code.
		The valid values is any combination of three of the following: GATCRYWSMKHBVDN
		'''
		assert type(amb_codon) is str and len(amb_codon) == 3, 'Error, the ambiguous codon must be a string three characters long.'
		m = re.match('^[GATCRYWSMKHBVDN]{3}$', amb_codon)
		assert m != None, 'Error, the codon %s is not valid. It may only use the chracters GATCRYWSMKHBVDN.' % amb_codon
		self.triplet = amb_codon
		
		
	def setTable(self, table):
		'''
		Set which codon table to use.
		The input is an integer.
		The valid values are: 1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23, 24, 25
		'''
		table = int(table)
		assert table in [1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23, 24, 25], 'Error, %s is an invalid codon table.' % table
		self.table = table
	
	################################################################		

	
		
if __name__ == '__main__':
	AA, table = getinput()
	codon_object = AmbigousCodon(AA, table)
	print("mixed base codon: %s" % codon_object.getTriplet())
	print("target AA: %s" % codon_object.getTarget())
	print("off-target AA: %s" % codon_object.getOfftarget())
	print("AA that can be added w/o off-targets: %s" % codon_object.getPossible())
