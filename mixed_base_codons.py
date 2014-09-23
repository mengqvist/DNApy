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


#TODO 
#this code is awful, I'll need to re-write it so that others, and myself can understand it..

from copy import deepcopy
import dna



## Explanation of algorithm ##
#a list of amino acids (as single letter code) is passed to the algorithm



DesiredAA = [] #list that holds the desired amino acids

def sumupcodons(DesiredCodons):
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
	#print('sumcodons', DesiredCodons)
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
		allcodon1.append(codon1) #Adds upp all the first bases
		allcodon2.append(codon2) #Adds upp all the second bases
		allcodon3.append(codon3) #Adds upp all the third bases
	return (allcodon1, allcodon2, allcodon3)



		
def commonNuc(nuc_list): 
	"""
	This function takes a list of lists and finds the degenerate symbol that represents at least one nucleotide from each of the lists.

	An example input is: [['T', 'C', 'A', 'G'], ['T', 'C'], ['T', 'C']].
	T and C are both present in all lists, therefore, both 'T' and 'C' are acceptable return values.
	In this case the nucleotide (T or C) first tested by the algorithm is returned as the output.
	Typically these are tested in the order 'A' 'G' 'C' 'T'.

	Another example input is: [['G'], ['T'], ['T']].
	In this case G or T is present in all lists, therefore the only acceptable output is 'K' (ambiguous nucleotide for G and T). 
	"""
	nuc_list = [s.upper() for s in nuc_list]
	
	print('list to degerate', list)
	
	for nuc in 'AGCTYKMSWRHVDBN'
	[]
	
	for x in list:	#x is the list of codons for a certain AA, it cycles through all chosen AA
		if 'A' in x:
			var = 'A'
		else:
			var = ''
			break
	if var != '': return var

	for x in list:
		if 'G' in x:
			var = 'G'
		else:
			var = ''
			break
	if var != '': return var

	for x in list:		
		if 'C' in x:
			var = 'C'
		else:
			var = ''
			break
	if var != '': return var

	for x in list:
		if 'T' in x:
			var = 'T'
		else:
			var = ''
			break
	if var != '': return var

	for x in list:
		if 'C' in x or 'T' in x:
			var = 'Y'
		else:
			var = ''
			break
	if var != '': return var

	for x in list:
		if 'G' in x or 'T' in x:
			var = 'K'
		else:
			var = ''
			break
	if var != '': return var

	for x in list:	
		if 'A' in x or 'C' in x:
			var = 'M'
		else:
			var = ''
			break
	if var != '': return var

	for x in list:	
		if 'C' in x or 'G' in x:
			var = 'S'	
		else:
			var = ''
			break
	if var != '': return var

	for x in list:	
		if 'A' in x or 'T' in x:
			var = 'W'
		else:
			var = ''
			break
	if var != '': return var

	for x in list:
		if 'A' in x or 'G' in x:
			var = 'R'
		else:
			var = ''
			break
	if var != '': return var

	for x in list:
		if 'C' in x or 'T' in x or 'A' in x:
			var = 'H'
		else:
			var = ''
			break
	if var != '': return var

	for x in list:
		if 'C' in x or 'A' in x or 'G' in x:
			var = 'V'
		else:
			var = ''
			break
	if var != '': return var

	for x in list:
		if 'T' in x or 'A' in x or 'G' in x:
			var = 'D'
		else:
			var = ''
			break
	if var != '': return var

	for x in list:
		if 'C' in x or 'T' in x or 'G' in x:
			var = 'B'
		else:
			var = ''
			break
	if var != '': return var

	for x in list:
		if 'C' in x or 'T' in x or 'A' in x or 'G' in x: 
			var = 'N'
		else: 
			var = 'X'
			break
	if var != '': return var

	return var;


def count_codon_list(codonlist):
	"""
	Takes a list of codons and counts how many of each AA are coded for. Returns a dictionary.
	"""
	#needs fixing
	F = ['TTT', 'TTC'];
	L = ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'];
	S = [ 'AGC', 'AGT', 'TCT', 'TCC', 'TCA', 'TCG'];
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

	ResultingAA = {'F':0,'L':0,'S':0,'Y':0,'C':0,'W':0,'P':0,'H':0,'E':0,'R':0,'I':0,'M':0,'T':0,'N':0,'K':0,'V':0,'A':0,'D':0,'Q':0,'G':0, 'stop':0}
	for i in range(0, len(codonlist)):
		if any(codonlist[i] in s for s in F):
			ResultingAA['F'] += 1
		elif any(codonlist[i] in s for s in L):
			ResultingAA['L'] += 1
		elif any(codonlist[i] in s for s in S):
			ResultingAA['S'] += 1
		elif any(codonlist[i] in s for s in Y):
			ResultingAA['Y'] += 1
		elif any(codonlist[i] in s for s in stop):
			ResultingAA['stop'] += 1
		elif any(codonlist[i] in s for s in C):
			ResultingAA['C'] += 1
		elif any(codonlist[i] in s for s in W):
			ResultingAA['W'] += 1
		elif any(codonlist[i] in s for s in P):
			ResultingAA['P'] += 1
		elif any(codonlist[i] in s for s in H):
			ResultingAA['H'] += 1
		elif any(codonlist[i] in s for s in E):
			ResultingAA['E'] += 1
		elif any(codonlist[i] in s for s in R):
			ResultingAA['R'] += 1
		elif any(codonlist[i] in s for s in I):
			ResultingAA['I'] += 1
		elif any(codonlist[i] in s for s in M):
			ResultingAA['M'] += 1
		elif any(codonlist[i] in s for s in T):
			ResultingAA['T'] += 1
		elif any(codonlist[i] in s for s in N):
			ResultingAA['N'] += 1
		elif any(codonlist[i] in s for s in K):
			ResultingAA['K'] += 1
		elif any(codonlist[i] in s for s in V):
			ResultingAA['V'] += 1
		elif any(codonlist[i] in s for s in A):
			ResultingAA['A'] += 1
		elif any(codonlist[i] in s for s in D):
			ResultingAA['D'] += 1
		elif any(codonlist[i] in s for s in Q):
			ResultingAA['Q'] += 1
		elif any(codonlist[i] in s for s in G):
			ResultingAA['G'] += 1
		else:
			print('Error, this is not an amino acid codon')
	return ResultingAA	



def chosenvsresulting(DesiredAA, AAlist): 
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
	return (TargetAA, OffTargetAA)

	

def getinput():
	###This is where I decide which AA I want
	DesiredAA = [];
	AllNaturalAA = ['F','L','S','Y','C','W','P','H','E','R','I','M','T','N','K','V','A','D','Q','G'];
	AA = 'x'
	table = raw_input('Input a number to indicate which codon table to use (standard is 1): ')
	while AA != '':
		AA = raw_input('Input single AA in single letter code. If done, press enter: ')
		AA = AA.upper()	
		if AA == '':
			pass	
		elif any(AA in s for s in AllNaturalAA):	
			DesiredAA.append(AA)	
		else:
			print('This is not a valid AA')
	#print(DesiredAA)
	return DesiredAA, table
	
	

def get_triplets(chosenAA):	
	'''
	Takes a list of chosen amino acids (in single letter code) and returns a list of all possible triplet codons.
	'''
	
	regular_triplets = []

	#get non-ambigous codons for the amino acids
	for aa in chosenAA:
		regular_triplets.append(dna.GetCodons(aa, separate=True))

	######## this is where I need to use the codon table #############
	first, second, third = sumupcodons(regular_triplets) #takes the triplets and splits them up into their first, second and third positions

	degenerate1 = dna.Amb(first) #Gets degenerate codon that represents all bases at position 1
	print('degenerate1', degenerate1)
	degenerate2 = dna.Amb(second) #Gets degenerate codon that represents all bases at position 2
	print('degenerate2', degenerate2)
	degenerate3 = dna.Amb(third) #Gets degenerate codon that represents all bases at position 3
	print('degenerate3', degenerate3)

	if degenerate1 == 'N' and degenerate2 == 'N' and degenerate3 == 'N': 
		triplet = degenerate1 + degenerate2 + 'K' #for this special case I want NNK, not NNN
	else: 
		triplet = degenerate1 + degenerate2 + degenerate3 # summing up
	
	#self.setTriplet(triplet)

	
	
	################ split here ################
	
	
	###now I just need to convert this to a list of real codons and then check to which aa they match
	Realcodons = dna.combine([dna.UnAmb(triplet[0]), dna.UnAmb(triplet[1]), dna.UnAmb(triplet[2])]) #condense the different codons for position 1, 2, and 3 to a list of triplets
	
	ResultingAA = []
	for codon in Realcodons:
		ResultingAA.append(dna.Translate(codon)) #Check which AA these codons code for

	
	TargetAA, OffTargetAA = chosenvsresulting(chosenAA, ResultingAA) #Check which of these AA were desired and which were not

	return (triplet, TargetAA, OffTargetAA)

def test_alternate(DesiredAA, AA, triplet, result):
	'''Sees whether alternate codons causes less off-target amino acids'''
	alternateAA = []
	alternateAA.append(DesiredAA)
	for n in range(len(AA)): #make all combinations of codons
		tempAA = [x[:] for x in alternateAA] #creates a deep copy which is needed to copy nested lists
		for i in range(len(tempAA)):
			try:
				index = tempAA[i].index(AA[n])
				tempAA[i][index] = '%s2' % AA[n]	
			except:
				pass
		alternateAA.extend(tempAA)
	
	for entry in alternateAA: #test which combination gives least off-target hits
		testset = entry
		temptriplet, tempresult = get_triplets(testset)
		if len(tempresult[1]) < len(result[1]):
			DesiredAA = testset
			result = tempresult
			triplet = temptriplet
	return (triplet, result)

def next_steps(targetAA, offtargetAA):
	"""Function for finding which other amino acids can be selected without introducion further off-target ones"""
	possibleAA = [] #for storing which ones are possible
	AllNaturalAA = dna.CodonTable(1).getCodons(separate=True).keys()
	print('natural', AllNaturalAA)
	unusedAA = [] #get the unused amino acids

	for AA in AllNaturalAA:
		if AA in targetAA:
			pass
		else:
			unusedAA.append(AA)

	if len(unusedAA)>0:
		for AA in unusedAA : #check which ones should be added to the new testset
			TrueFalseList = []
			testAA = deepcopy(targetAA)
			testAA.append(AA)
			new_codon, new_target, new_offtarget = get_codon_for_chosen_AA(testAA)   
			for i in range(len(new_offtarget)):
				if (new_offtarget[i] in offtargetAA) == True:
					TrueFalseList.append(True)
				elif (new_offtarget[i] in offtargetAA) == False:
					TrueFalseList.append(False)
			if False in TrueFalseList:
				pass
			else:
				possibleAA.append(AA)
	#print('possibleAA', possibleAA)

	#change S2 to S and L2 to L
#	new_possibleAA = []
#	for AA in possibleAA:
#		if AA == 'S2' and ('S' in new_possibleAA) == False:
#			new_possibleAA.append('S')
#		elif AA == 'L2' and ('L' in new_possibleAA) == False:
#			new_possibleAA.append('L')
#		else:
#			new_possibleAA.append(AA)
	return possibleAA




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
				

def get_codon_for_chosen_AA(DesiredAA):	
	'''
	
	'''
	#DesiredAA containst list of amino acids.
	triplet, target, offtarget = get_triplets(DesiredAA)
	AA = ''
	if ('S' in DesiredAA):
		AA += 'S'
	if ('R' in DesiredAA):
		AA += 'R'
	if ('L' in DesiredAA):
		AA += 'L'
	#triplet, target, offtarget = test_alternate(DesiredAA, AA, triplet, result) #fix this!!!!!!!!!!!!!!!!!!!!!


#	print('Degenerate codon: ', triplet)
#	print('For the chosen AA: ', result[0])
#	print('And the off-target AA: ', result[1])
	
	return (triplet, target, offtarget) #codon, target, offtarget
	


class AmbigousCodon:
	'''
	Class that holds methods and values for computing the ambigous codon for a list of amino acids.
	Required input is a list of desired amino acids in single letter code and 
	an integer that determines the codon table to use.
	'''
	def __init__(self, AA_list, table):
		self.setTable(table)
		AA_list = [s.upper() for s in AA_list]
		self.setTarget(AA_list)

		#compute triplet and offtarget based on input
		triplet, target, offtarget = get_codon_for_chosen_AA(AA)		
		
		#now find which aa are possible without adding  off-target hits
		possibleAA = next_steps(target, offtarget)

	#	if len(offtarget) != 0: #if there are offtargets, find a combination of codons that will give you no off-targets
	#		multiple_codons = find_multiple_codons(target)
		
		return
	
	def setTarget(self, AA_list):
		for s in AA_list:
			assert s in 'FLSYCWPHERIMTNKVADQG'
		self.target = AA_list
	def getTarget(self):
		return self.target
	
	def setOfftarget(self, AA_list):
		self.offtarget = AA_list
	def getOfftarget():
		return self.offtarget
	
	def setPossible(self, AA_list):
		self.possible = AA_list
	def getPossible(self):
		return self.possible
	
	def setTriplet(self, triplet_string):
		self.triplet = triplet_string
	def getTriplet(self):
		return self.triplet
		
	def setTable(self, table):
		self.tableobject = dna.CodonTable(table)
	def getTable(self):
		return self.tableobject

		
if __name__ == '__main__': #if script is run by itself and not loaded	
	AA, table = getinput()
	codon_object = AmbigousCodon(AA, table)
	print("mixed base codon: %s" % codon_object.getTriplet())
	print("target AA: %s" % codon_object.getTarget())
	print("off-target AA: %s" % codon_object.getOfftarget())
	print("AA that can be add w/o off-targets: %s" % codon_object.getPossible())
