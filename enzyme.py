#!/usr/bin/python


#DNApy is a DNA editor written purely in python. 
#The program is intended to be an intuitive, fully featured, 
#extendable, editor for molecular and synthetic biology.  
#Enjoy!
#
#copyright (C) 2014  Martin Engqvist 
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
#Get source code at: https://github.com/0b0bby0/DNApy
#
#

#
# This file should contain classes to enable enzyme functions
# as of now it contains the following classes:
# 
# restrictionEnzyme()		store info about an restriction enzyme type II
# initRestriction()			class to store and update all avialible restriction enzymes
#




from base_class import DNApyBaseClass
import wx
import re
import string
import genbank
import math
import collections

# class to make restriktion enzymes more usefull
# name = name of enzyme
# pattern = recognition site
# length = length of pattern
# ncuts = number of cuts made by enzyme
#         Zero represents unknown
# blunt = true if blunt end cut, false if sticky
# c1 = First 5' cut
# c2 = First 3' cut
# c3 = Second 5' cut
# c4 = Second 3' cut



##########################################################################
# Class restrictionEnzyme()
class restrictionEnzyme():
	def __init__(self, name, pattern, length, ncuts,blunt,c51,c31,c52,c32,regex):
		self.name 		= str(name)
		self.pattern 	= str(pattern)
		self.length 	= int(length)
		self.ncuts 		= int(ncuts)
		self.blunt		= int(blunt)
		self.c51		= int(c51)
		self.c31		= int(c31)
		self.c52		= int(c52)
		self.c32		= int(c32)
		
		self.regex		= regex
		
		self.restrictionSites	= []	# empty list
		
	def addRestrictionSite(self,start,end,cut51,cut52,dnaMatch):
		self.restrictionSites.append([self.name,start,end,cut51,cut52,dnaMatch])
		
	def resetRestrictionSite(self):
		self.restrictionSites = []
# End of Class restrictionEnzyme()
##########################################################################



##########################################################################
# Class initRestriction()
# 
# store and update info about restriction enzymes
#
class initRestriction():
	def __init__(self):
		'''Class to be loaded at every file change. It then evaluates the restrictionsites
		and saves every Info. This enables fast access from every part of the software to 
		restriction sites''' 
		
		# on init get the current DNA
		self.oldDna 	= ''
		self.currentDNA = genbank.gb.gbfile["dna"]
		
		# load enzymes
		self.loadEnzymes()
		
		# add the other stuff
		self.reloadEnzymes()
			
	def reloadEnzymes(self):
		'''should be called if you want to reload the enzymes'''

		self.currentDNA = genbank.gb.gbfile["dna"]
		
		# just if the dna had changed, reload the restriction sites
		if self.currentDNA != self.oldDna:
			
			# updates cuts:
			self.findRestrictionSites()
			
			self.oldDna 			= genbank.gb.gbfile["dna"]
		
		
		
		

	def loadEnzymes(self):
		''' load the restrictionenzymes from the emboss file, once. '''	
		self.enzymeObj			= collections.OrderedDict()
		
		self.allEnzymes			= []


		# for this we have enzymes in folder /resources
		with open('resources/emboss_e.txt') as f:
			for line in f:	
				# each line is one enzyme, except the header
				if ( line[:1] != "#" ):	
		
					# split the line
					lineparts = re.split(r'\t+', line)	
			
			
			

					# form the regexp:
					#Code	Meaning			Etymology	Complement	Opposite
					#A	A			Adenosine	T	B
					#T/U	T			Thymidine/Uridine	A	V
					#G	G			Guanine	C	H
					#C	C			Cytidine	G	D
					#K	G or T			Keto	M	M
					#M	A or C			Amino	K	K
					#R	A or G			Purine	Y	Y
					#Y	C or T			Pyrimidine	R	R
					#S	C or G			Strong	S	W
					#W	A or T			Weak	W	S
					#B	C or G or T		not A (B comes after A)	V	A
					#V	A or C or G		not T/U (V comes after U)	B	T/U
					#H	A or C or T		not G (H comes after G)	D	G
					#D	A or G or T		not C (D comes after C)	H	C
					#X/N	G or A or T or C	any	N	.
					#.	not G or A or T or C	.	N
					#-	gap of indeterminate length	
					regpattern = lineparts[1]
					regpattern = string.replace(regpattern,"U","T")
					regpattern = string.replace(regpattern,"K","(G|T)")
					regpattern = string.replace(regpattern,"M","(A|C)")
					regpattern = string.replace(regpattern,"R","(A|G)")
					regpattern = string.replace(regpattern,"Y","(C|T)")
					regpattern = string.replace(regpattern,"S","(C|G)")
					regpattern = string.replace(regpattern,"W","(A|T)")
					regpattern = string.replace(regpattern,"B","[CGT]")
					regpattern = string.replace(regpattern,"V","[ACG]")
					regpattern = string.replace(regpattern,"H","[ACT]")
					regpattern = string.replace(regpattern,"D","[AGT]")
					regpattern = string.replace(regpattern,"N","[AGTC]")
					regpattern = string.replace(regpattern,"X","[AGTC]")

			
				
					# use new object instead:
					name 		= lineparts[0]
					pattern		= lineparts[1]
					length		= lineparts[2]
					ncuts		= int(lineparts[3])
					blunt		= lineparts[4]
					c51			= int(lineparts[5])
					c31			= int(lineparts[6])
					c52			= int(lineparts[7])
					c32			= int(lineparts[8])
					
					
					regex		= re.compile(regpattern, re.IGNORECASE)
				
					self.enzymeObj[name] = restrictionEnzyme(name, pattern, length, ncuts,blunt,c51,c31,c52,c32,regex)
					self.allEnzymes.append(name)

		return self.enzymeObj
	
	
	
	# load and find all restriction sites	
	def findRestrictionSites(self):


		dnaseq      = self.currentDNA
		if dnaseq != None:
		
			# circular dna
			# we just inspect the region +-100 
			# to find enzyme, which cut near 0 in a ciclic plasmid
			if genbank.gb.gbfile['locus']['topology'] == 'circular':
				circularDnaHelper = dnaseq[:100]  # helper of 200bp from 0 to 100
			else:
				circularDnaHelper = ''
					
			wholeDNA2Inspect = '%s%s' % (dnaseq, circularDnaHelper)		     # dna width circular helper added

			# loop all the selected enzymes
			for enzyme in self.enzymeObj:
				
				# reset restrictionsites first:
				self.enzymeObj[enzyme].resetRestrictionSite()
				
				restrictionsitesList = [] # variable to return the restriction sites
				r         	= self.enzymeObj[enzyme].regex
				# get the cut position
				offset1     = self.enzymeObj[enzyme].c51	# first 5' cut
				offset2     = self.enzymeObj[enzyme].c52 	# second 5' cut if any
				# handle the cuts
				iterator    = r.finditer(wholeDNA2Inspect)      # find in dnaseq and circularDnaHelper
				for match in iterator:
					newEnz = []
					modulo = len(dnaseq)		# all positions are modulo len(dnaseq)
					# add offset to match position
					name 		= str(enzyme)
					start		= (match.start()+1) % modulo
					end 		= match.end() % modulo
					cut51		= (match.start() + offset1) % modulo	# cuts after this base

					if offset2 == 0: # if it just cuts once, we'll say None to the second cut
						cut52	= None
					else:
						ncut52	= (match.start() + offset2) % modulo

					dnaMatch 	= wholeDNA2Inspect[match.start():match.end()]
				
					# save the new restriction site if its new
					if start not in restrictionsitesList:
						restrictionsitesList.append(start)
						
					# now add the sites to the enzyme
					self.enzymeObj[enzyme].addRestrictionSite(start,end,cut51,cut52,dnaMatch)
		
			
			
			return restrictionsitesList
# End of Class initRestriction()
##########################################################################


