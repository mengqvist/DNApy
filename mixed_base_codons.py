#!/usr/bin/env python


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


#TODO 
#this code is awful, I'll need to re-write it so that others, and myself can understand it..

from copy import deepcopy

DesiredAA = [] #list that holds the desired amino acids


def AAtocodons(AAlist, option):
	"""Fetches codons for a given AA. Optionally pass 'y' if a restricted codon table without rare codons should be used"""
	codonlist = []
	if option == 'y' or option == 'Y':
		F = ['TTT', 'TTC'];
		L = ['TTA', 'TTG'];
		L2 = ['CTT', 'CTC', 'CTG'];
		S = ['TCT', 'TCC', 'TCA', 'TCG'];
		S2 = ['AGC', 'AGT']
		Y = ['TAT', 'TAC'];
		stop = ['TAA', 'TAG', 'TGA'];
		C = ['TGT', 'TGC'];
		W = ['TGG'];
		P = ['CCT', 'CCA', 'CCG'];
		H = ['CAT', 'CAC'];
		E = ['GAA', 'GAG'];
		R = ['CGT', 'CGC'];
		I = ['ATT', 'ATC'];
		M = ['ATG'];
		T = ['ACG', 'ACC', 'ACA', 'ACT'];
		N = ['AAT', 'AAC'];
		K = ['AAA', 'AAG'];
		V = ['GTT', 'GTA', 'GTG', 'GTC'];
		A = ['GCG', 'GCC', 'GCA', 'GCT'];
		D = ['GAT', 'GAC'];
		Q = ['CAG', 'CAA'];
		G = ['GGT', 'GGC', 'GGG']	
	else:
		F = ['TTT', 'TTC'];
		L = ['TTA', 'TTG'];
		L2 = ['CTT', 'CTC', 'CTA', 'CTG'];
		S = ['TCT', 'TCC', 'TCA', 'TCG'];
		S2 = ['AGC', 'AGT'];
		Y = ['TAT', 'TAC'];
		stop = ['TAA', 'TAG', 'TGA'];
		C = ['TGT', 'TGC'];
		W = ['TGG'];
		P = ['CCT', 'CCA', 'CCG', 'CCC'];
		H = ['CAT', 'CAC'];
		E = ['GAA', 'GAG'];
		R = ['CGT', 'CGA', 'CGG', 'CGC'];
		R2 = ['AGG', 'AGA'];		
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
	
	for i in range(0, len(AAlist)):
		codonlist.append(eval(AAlist[i]))  #eval is key for converting the string AA's to codons
	return codonlist


def sumupcodons(DesiredCodons):
	allcodon1 = [];
	allcodon2 = [];
	allcodon3 = [];
	for entry in DesiredCodons:     #splits up the list in segments containing codons for each AA
		codon1 = [];
		codon2 = [];
		codon3 = [];

		for codon in entry:      #splits up the codons of each AA into triplets
			if codon[0:1] in codon1: 
				pass
			else:
				codon1.append(codon[0:1]) #Takes first base in the triplet
	
			if codon[1:2] in codon2: 
				nothing = 1 #do nothing
			else:
				codon2.append(codon[1:2]) #Takes the second base in triplet

			if codon[2:3] in codon3: 
				nothing = 1 #do nothing
			else:
				codon3.append(codon[2:3]) #Takes third base in triplet
		
		allcodon1.append(codon1) #Adds upp all the first bases
		allcodon2.append(codon2) #Adds upp all the second bases
		allcodon3.append(codon3) #Adds upp all the third bases
	#print('1', allcodon1)
	#print('2', allcodon2)
	#print('3', allcodon3)
	return (allcodon1, allcodon2, allcodon3)


def degenerate(list): 
	"""This function finds the degenerate symbol for a certain base"""
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


def checkdegenerate(codon):  
	"""Converts a degenerate codon to a list of regular (ATCG) codons"""
	var = []	
	if 'A' in codon:
		var.append('A')

	elif 'C' in codon:
		var.append('C')

	elif 'T' in codon:
		var.append('T')

	elif 'G' in codon:
		var.append('G')

	elif 'M' in codon:
		var.append('A')
		var.append('C')

	elif 'Y' in codon:
		var.append('C')
		var.append('T')

	elif 'K' in codon:
		var.append('G')
		var.append('T')

	elif 'S' in codon:
		var.append('C')
		var.append('G')

	elif 'W' in codon:
		var.append('A')
		var.append('T')

	elif 'R' in codon:
		var.append('A')
		var.append('G')

	elif 'H' in codon:
		var.append('C')
		var.append('T')
		var.append('A')

	elif 'V' in codon:
		var.append('C')
		var.append('A')
		var.append('G')

	elif 'D' in codon:
		var.append('T')
		var.append('A')
		var.append('G')

	elif 'B' in codon:
		var.append('C')
		var.append('T')
		var.append('G')

	elif 'N' in codon: 
		var.append('C')
		var.append('T')
		var.append('A')
		var.append('G')
	else: 
		var.append('X')

	return var;


def combinelists(list1, list2, list3):
	"""Function makes every combination of list1, list2, and list3 to make actual codons out of them, but retains their internal order."""
	some_var = []
	for i in range(0, len(list1)):
		some_var.append(list1[i])

	some_var2 = []
	for i in range(0, len(list2)):
		for n in range(0, len(some_var)):
			some_var2.append(some_var[n] + list2[i])

	some_var3 = []
	for i in range(0, len(list3)):
		for n in range(0, len(some_var2)):
			some_var3.append(some_var2[n] + list3[i])
	return some_var3

def count_codon_list(codonlist):
	"""Takes a list of codons and counts how many of each AA are coded for. Returns a dictionary."""
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

def translatecodonlist(codonlist):
	"""Takes a list of codons and returns the amino acids encoded by them"""
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

	ResultingAA = []
	for i in range(0, len(codonlist)):
		if any(codonlist[i] in s for s in F) and any('F' in s for s in ResultingAA) == False:
			if any('F' in s for s in ResultingAA) == False:
				ResultingAA.append('F')
		elif any(codonlist[i] in s for s in L):
			if any('L' in s for s in ResultingAA) == False:
				ResultingAA.append('L')
		elif any(codonlist[i] in s for s in S):
			if any('S' in s for s in ResultingAA) == False:
				ResultingAA.append('S')
		elif any(codonlist[i] in s for s in Y):
			if any('Y' in s for s in ResultingAA) == False:
				ResultingAA.append('Y')
		elif any(codonlist[i] in s for s in stop):
			if any('stop' in s for s in ResultingAA) == False:	
				ResultingAA.append('stop')
		elif any(codonlist[i] in s for s in C):
			if any('C' in s for s in ResultingAA) == False:
				ResultingAA.append('C')
		elif any(codonlist[i] in s for s in W):
			if any('W' in s for s in ResultingAA) == False:
				ResultingAA.append('W')
		elif any(codonlist[i] in s for s in P):
			if any('P' in s for s in ResultingAA) == False:
				ResultingAA.append('P')
		elif any(codonlist[i] in s for s in H):
			if any('H' in s for s in ResultingAA) == False:
				ResultingAA.append('H')
		elif any(codonlist[i] in s for s in E):
			if any('E' in s for s in ResultingAA) == False:
				ResultingAA.append('E')
		elif any(codonlist[i] in s for s in R):
			if any('R' in s for s in ResultingAA) == False:
				ResultingAA.append('R')
		elif any(codonlist[i] in s for s in I):
			if any('I' in s for s in ResultingAA) == False:	
				ResultingAA.append('I')
		elif any(codonlist[i] in s for s in M):
			if any('M' in s for s in ResultingAA) == False:
				ResultingAA.append('M')
		elif any(codonlist[i] in s for s in T):
			if any('T' in s for s in ResultingAA) == False:
				ResultingAA.append('T')
		elif any(codonlist[i] in s for s in N):
			if any('N' in s for s in ResultingAA) == False:
				ResultingAA.append('N')
		elif any(codonlist[i] in s for s in K):
			if any('K' in s for s in ResultingAA) == False:
				ResultingAA.append('K')
		elif any(codonlist[i] in s for s in V):
			if any('V' in s for s in ResultingAA) == False:
				ResultingAA.append('V')
		elif any(codonlist[i] in s for s in A):
			if any('A' in s for s in ResultingAA) == False:
				ResultingAA.append('A')
		elif any(codonlist[i] in s for s in D):
			if any('D' in s for s in ResultingAA) == False:	
				ResultingAA.append('D')
		elif any(codonlist[i] in s for s in Q):
			if any('Q' in s for s in ResultingAA) == False:
				ResultingAA.append('Q')
		elif any(codonlist[i] in s for s in G):
			if any('G' in s for s in ResultingAA) == False:
				ResultingAA.append('G')
		else:
			ResultingAA.append('?')
	return ResultingAA


def chosenvsresulting(DesiredAA, AAlist): 
	"""Function for checking a list vs another to find which are present in both and which are not"""
	TargetAA = []
	OffTargetAA = []
	for entry in AAlist:
		if any(entry in s for s in DesiredAA):	
			TargetAA.append(entry)	
		else:
			OffTargetAA.append(entry)
	return (TargetAA, OffTargetAA)


def getinput():
	input_var = raw_input('Would you like to exclude the rare codons CGA, CGG, AGA, AGG, CUA, AUA, GGA, and CCC (y/n)? ')
	###This is where I decide which AA I want
	DesiredAA = [];
	AllNaturalAA = ['F','L','S','Y','C','W','P','H','E','R','I','M','T','N','K','V','A','D','Q','G'];
	AA = 'x'
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
	return DesiredAA
	
	

def evaluate(chosenAA, codon_restriction):	
	DesiredCodons = AAtocodons(chosenAA, codon_restriction) #get triplet codons for the desired AA
	allcodons = sumupcodons(DesiredCodons) #takes the triplets and splits them up into their first, second and third positions

	degenerate1 = degenerate(allcodons[0]) #Gets degenerate codon that represents all bases at position 1
	degenerate2 = degenerate(allcodons[1]) #Gets degenerate codon that represents all bases at position 2
	degenerate3 = degenerate(allcodons[2]) #Gets degenerate codon that represents all bases at position 3


	if degenerate1[0] == 'N' and degenerate2[0] == 'N' and degenerate3[0] == 'N': 
		triplet = degenerate1[0] + degenerate2[0] + 'K' #for this special case I want NNK, not NNN
	else: 
		triplet = degenerate1[0] + degenerate2[0] + degenerate3[0] # summing up
	##triplet holds the degenerate triplet codon


	###now I just need to convert this to a list of real codons and then check to which aa they match
	Realcodons = combinelists(checkdegenerate(triplet[0]), checkdegenerate(triplet[1]), checkdegenerate(triplet[2])) #condense the different codons for position 1, 2, and 3 to a list of triplets
	ResultingAA = translatecodonlist(Realcodons) #Check which AA these codons code for
	result = chosenvsresulting(chosenAA, ResultingAA) #Check which of these AA were desired and which were not
	return (triplet, result)

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
		temptriplet, tempresult = evaluate(testset, 'n')
		if len(tempresult[1]) < len(result[1]):
			DesiredAA = testset
			result = tempresult
			triplet = temptriplet
	return (triplet, result)

def next_steps(targetAA, offtargetAA):
	"""Function for finding which other amino acids can be selected without introducion further off-target ones"""
	possibleAA = [] #for storing which ones are possible
	AllNaturalAA = ['F','L', 'L2','S', 'S2', 'Y','C','W','P','H','E','R', 'R2','I','M','T','N','K','V','A','D','Q','G'];
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

def get_codon_for_chosen_AA(AA):	
	#DesiredAA containst list of amino acids and also the option wheter codons should be excluded or not.
	DesiredAA = AA
	codon_restriction = 'n'	

	triplet, result = evaluate(DesiredAA, codon_restriction)
	AA = ''
	if ('S' in DesiredAA):
		AA += 'S'
	if ('R' in DesiredAA):
		AA += 'R'
	if ('L' in DesiredAA):
		AA += 'L'
	triplet, result = test_alternate(DesiredAA, AA, triplet, result)


#	print('Degenerate codon: ', triplet)
#	print('For the chosen AA: ', result[0])
#	print('And the off-target AA: ', result[1])
	
	return (triplet, result[0], result[1]) #codon, target, offtarget


def find_multiple_codons(target_AA):	
	'''For difficult desired AA, find a combination of codons that will give you no off-targets'''
	target_AA.sort()
	codon_restriction = 'n'
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
				AA_list = translate(checkdegenerate(mixed_base_codon)) #get AA for given codon
				codon, target, offtarget = get_codon_for_chosen_AA(AA_list)
				triplet, result = evaluate(AA_list, codon_restriction)
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
					triplet, result = evaluate(remaining_AA, codon_restriction)
					target = result[0].sort()
					offtarget = result[1]
					if len(offtarget) == 0 and target == remaining_AA: #if no offtarget aa are present and all target aa are present
						codons.append(triplet)
						return codons


					

#	checkdegenerate()

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
				

def run(AA):
	triplet, target, offtarget = get_codon_for_chosen_AA(AA)


	#now find which aa are possible without adding  off-target hits
	possibleAA = next_steps(target, offtarget)

#	if len(offtarget) != 0: #if there are offtargets, find a combination of codons that will give you no off-targets
#		multiple_codons = find_multiple_codons(target)

	return (triplet, target, offtarget, possibleAA) #codon, target, offtarget, possibleAA


if __name__ == '__main__': #if script is run by itself and not loaded	
	AA = getinput()
	print("mixed base codon: " + str(run(AA)[0]))
	print("target AA: " + str(run(AA)[1]))
	print("off-target AA: " + str(run(AA)[2]))
	print("AA that can be add w/o off-targets: " + str(run(AA)[3]))
