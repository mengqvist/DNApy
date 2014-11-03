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

	

def getinput():
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
	DesiredAA = list(set(DesiredAA)) #condense list to what is unique
	#print(DesiredAA)
	return DesiredAA, table
	

def find_multiple_codons(target_AA):	
	'''For difficult desired AA, find a combination of codons that will give you no off-targets'''
	target_AA.sort()
	symbols = 'ATGCRYMKSWHBVDN'
	two_codons = []
	three_codons = []
	four_codons = []
	five_codons = []
	for i in target_AA:
		for j in target_AA:
			if i == j:
				pass
			for k in target_AA:
				for l in target_AA:	
					
					
					two_codons1.append([i])
					
					three_codons = []
					four_codons = []
					five_codons = []

				codons = []
				mixed_base_codon = i+j+k
				codons.append(mixed_base_codon)

				#change!
				AA_list = translate(dna.UnAmb(mixed_base_codon)) #wrong usage, have to take a base at the time!!!
				codon, target, offtarget = get_codon_for_chosen_AA(AA_list)
				triplet, result = get_triplets(AA_list)
				target = result[0].sort()
				offtarget = result[1]
				if len(offtarget) == 0 and target == target_AA: #if no offtarget aa are present and all target aa are present
					return codons
				elif len(offtarget) == 0 and len(target) < len(target_AA): #codon matches some of the aa, with no offtargets
					remaining_AA = []
					for aa in target:
						if aa in target_AA:
							pass
						else:	
							remaining_AA.append(aa)
					triplet, result = get_triplets(remaining_AA)
					target = result[0].sort()
					offtarget = result[1]
					if len(offtarget) == 0 and target == remaining_AA: #if no offtarget aa are present and all target aa are present
						codons.append(triplet)
						return codons


					

#	dna.UnAmb()

#	codons = []
#	triplet, target, offtarget = get_codon_for_chosen_AA(target_AA)
#	if len(offtarget) == 0:
#		return triplet
#	else:
#		print('New round')
#		for i in range(len(target_AA)):
#			test_AA = target_AA[:]
#			del test_AA[i]
#			codons = find_multiple_codons(test_AA[:]) + find_multiple_codons(target_AA[i])
#			if codon1 != None and codon2 != None:
#				print(codon1)
#				print(codon2)
#				print('')
#				break
#			elif i == len(target_AA)-1:
#				print('Not found')
				


class AmbigousCodon:
	'''
	Class that holds methods and values for computing the ambigous codon for a list of amino acids.
	Required input is a list of desired amino acids in single letter code and 
	an integer that determines the codon table to use.
	
	The algorithm works as follows:
	A list of amino acids (as single letter code) is passed to the algorithm.
	All possible regular codons for those amino acids are looked up and returned as a list of lists -> get_triplets()
	All nucleotides for first, second and third position extracted into seperate lists while retaining list structure -> sumupcodons()
	Separately, for pos 1, 2 and 3, a degenerate nucleotide is found that matches at least one nucleotide in each lists -> commonNuc()
	These are concatenated and represents the first degenerate triplet
	The degenerate triplet is converted back to real codons xxxxx()
	translate to amino acids and check against what you wanted in the first place xxxxxxxx()

	Add more !!!!!!!!
	'''
	def __init__(self, AA_list, table):
		self.setTable(table)
		AA_list = [s.upper() for s in AA_list]
		self.setTarget(AA_list)
		print('Desired AA: %s' % AA_list)

		#compute the triplet and the off-target AA
		triplet, offtarget = self.find_degenerate(self.getTarget())
		self.setTriplet(triplet)
		self.setOfftarget(offtarget)

		#see which other can be added without further off-target
		possible = self.next_steps()
		self.setPossible(possible)
		return

		
	######## Public methods intended for user interaction #########
	def getTarget(self):
		'''
		
		'''
		return self.target
		
	def getOfftarget(self):
		'''
		
		'''
		return self.offtarget
		
	def getPossible(self):
		'''
		
		'''
		return self.possible
		
	def getTriplet(self):
		'''
		
		'''
		return self.triplet	
		
	def getTable(self):
		'''
		
		'''
		return self.table
		
	################################################################


	
	######## Methods NOT intended for direct user interaction ######
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



	def chosenvsresulting(self, DesiredAA, AAlist): 
		"""
		Function for checking a list vs another to find which are present in both and which are not.
		"""
		TargetAA = []
		OffTargetAA = []
		for entry in AAlist:
			if any(entry in s for s in DesiredAA):	
				TargetAA.append(entry)	
			else:
				OffTargetAA.append(entry)
		return OffTargetAA
	
	
	def flatten_codon_list(self, codon_list):
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
	
		#get all codons for chosen amino acids
		regular_triplets = [dna.GetCodons(aa, table=self.table, separate=True) for aa in AA_list]
#		print('Codons of desired AA: %s' % regular_triplets)
		
		#some of the codons are list of lists (happens when the amino acid has codons at different parts of the codon circle)
		#I need to flatten this into separate lists with which go on further
		regular_triplets = self.flatten_codon_list(regular_triplets)

		best_triplet = None #for storing best degenerate triplet
		best_offtarget = None #for storing the off-target AA of the best triplet
		best_score = None #for storing the length of the off-target list
		
		for codon_list in regular_triplets:
			#get all nucleotides for first, second and third position while retaining list structure		
			first, second, third = self.sumupcodons(codon_list) 
	#		print('First: %s' % first)
	#		print('Second: %s' % second)
	#		print('Third: %s' % third)
			
			#check which degenerate nucleotide can be used to find at least one match in each of the lists
			possible_triplets = dna.combine([dna.commonNuc(first), dna.commonNuc(second), dna.commonNuc(third)])
			
			#now go through them and see which is best
			for triplet in possible_triplets:
				#convert the triplet back to a list of real codons 
				Realcodons = dna.combine([dna.UnAmb(triplet[0]), dna.UnAmb(triplet[1]), dna.UnAmb(triplet[2])]) #condense the different codons for position 1, 2, and 3 to a list of triplets
	#			print('Codons when converting back degenerate triplet: %s' % Realcodons)
			
				#Check which AA these codons code for
				ResultingAA = [dna.Translate(codon, table=self.table) for codon in Realcodons]
	#			print('AA encoded by degenerate codon: %s' % ResultingAA)

				#compare which amino acids were desired with the ones resulting from the degenerate codon
				offtarget = self.chosenvsresulting(AA_list, ResultingAA)

				if len(offtarget) < best_score or best_score == None:
					best_triplet = triplet
					best_offtarget = offtarget
					best_score = len(best_offtarget)
#				print('Best triplet: %s' % best_triplet)
#				print('Offtarget: %s' % best_score)
		
		return best_triplet, best_offtarget

	
	def next_steps(self):
		"""Function for finding which other amino acids can be selected without introducing further (more in number) off-target ones"""

		possibleAA = []
		targetAA = self.getTarget()
		AllNaturalAA = dna.CodonTable(self.getTable()).getCodons().keys()
		if 'start' in AllNaturalAA:
			AllNaturalAA.remove('start')
		if 'stop' in AllNaturalAA:
			AllNaturalAA.remove('stop')
			AllNaturalAA.append('*')
		
		unusedAA = self.chosenvsresulting(targetAA, AllNaturalAA)

		
		if len(unusedAA)>0:
			for AA in unusedAA:
				temptargetAA = targetAA[:]
				temptargetAA.append(AA)
				triplet, offtarget = self.find_degenerate(temptargetAA)
				if len(offtarget) <= len(self.getOfftarget()):
					possibleAA.append(AA)
		return possibleAA

		
	def setTarget(self, AA_list):
		for s in AA_list:
			assert s in 'FLSYCWPHERIMTNKVADQG*'
		self.target = AA_list
	
	def setOfftarget(self, AA_list):
		self.offtarget = AA_list
	
	def setPossible(self, AA_list):
		self.possible = AA_list
	
	def setTriplet(self, triplet_string):
		self.triplet = triplet_string
		
	def setTable(self, table):
		self.table = table
	################################################################		

	
		
if __name__ == '__main__': #if script is run by itself and not loaded	
	AA, table = getinput()
	codon_object = AmbigousCodon(AA, table)
	print("mixed base codon: %s" % codon_object.getTriplet())
	print("target AA: %s" % codon_object.getTarget())
	print("off-target AA: %s" % codon_object.getOfftarget())
	print("AA that can be add w/o off-targets: %s" % codon_object.getPossible())
