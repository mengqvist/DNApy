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
#fix header parsing
#fix undo/redo
#fix qualifier buttons
#integrate the feature selection variable (from the genbank file)

import dna
from copy import deepcopy
import pyperclip
import oligo_localizer



class gbobject():
	"""Class that reads a genbank file (.gb) and has functions to edit its features and DNA sequence"""
	def __init__(self):
		self.gbfile = {} #this variable stores the whole genbank file
		self.filepath = '' #for keeping track of the file being edited
		self.clipboard = {}
		self.clipboard['dna'] = ''
		self.clipboard['features'] = []
		self.search_hits = []	# variable for storing a list of search hits
		self.dna_selection = (0, 0)	 #variable for storing current DNA selection
		self.feature_selection = False

	def treat_input_line(self, tempstr):
		'''Function for parsing a string containing feature information into the correct data format'''
		templist = tempstr.split('  ')

		templist[:] = [x for x in templist if x != ''] #remove empty entries
		
		#to remove any \r and \n newline characters at the end
		for i in range(len(templist)): 
			templist[i] = templist[i].rstrip('\r\n')

		#remove single whitespace in front
		for i in range(len(templist)): 
			if templist[i][0] == ' ':
				templist[i] = templist[i][1:]
	
		#to deal with feature descriptions or amino acid sequences that break over several lines 
		done = False
		i = 2
		while done != True:
			listlen = len(templist)
			if templist[i][0] != '/': 
				templist[i-1] += templist[i]
				del templist[i]
				i = 1
			if i >= listlen-1:
				done = True
			i += 1
		return templist

	def readgb(self, filepath):
		"""Function takes self.filepath to .gb file and extracts the header, features and DNA sequence"""
		a = open(filepath, 'r') #open it for reading
		gbfile = a.read()
		a.close()
		
		#split the file
		headerandfeatures, dna = gbfile.split('ORIGIN') #origin is the word in the file that seperates features from dna
		header, features = headerandfeatures.split('FEATURES')
		header = header[0:-1] #removes the \n at the end		
		
		#get the header.....
		self.gbfile['header'] = header
		
		
		#get the DNA and clean it from numbering and spaces
		DNA = '' 
		for i in range(len(dna)):
			if dna[i].lower() == 'a' or dna[i].lower() == 't' or dna[i].lower() == 'c' or dna[i].lower() == 'g':
				DNA = DNA + dna[i]
		self.gbfile['dna'] = DNA
		
		#get the features
		## need to add single base support!!!! ##
		
		featurelist2 = []
		featurelist = features.split('\n')
		templist = []
		tempstr = ''
		
		if featurelist[-1] == '': del featurelist[-1] #last entry tends to be empty, if so, remove
		
		for line in range(1, len(featurelist)):
			if ('..' in featurelist[line] and tempstr != '') == True:
				if tempstr[-1] == ',': #to protect against numberings that extend over several rows
					tempstr += featurelist[line]

				elif line+1 == len(featurelist):
					tempstr += featurelist[line]
			
				featurelist2.append(self.treat_input_line(tempstr))
				tempstr = featurelist[line]
				
			elif '..' in featurelist[line] and tempstr == '': #first feature
				tempstr = featurelist[line]

			else:						#in-between features
				tempstr += featurelist[line]

		#catch final entry
		featurelist2.append(self.treat_input_line(tempstr))
		

		#now arrange features into the correct data structure
		features = []
		for i in range(len(featurelist2)):
			Feature = {}
			
			#get key
			Feature['key'] = featurelist2[i][0] #append type of feature
			
			#get whether complement or not, join or not, order or not
			if 'complement' in featurelist2[i][1]:
				Feature['complement'] = True
			else:
				Feature['complement'] = False

			if 'join' in featurelist2[i][1]:
				Feature['join'] = True
			else:
				Feature['join'] = False
				
			if 'order' in featurelist2[i][1]:
				Feature['order'] = True
			else:
				Feature['order'] = False
								
			#get location
			tempsites = []
			commasplit = featurelist2[i][1].split(',')
			for entry in commasplit:
				tempstr = ''
				for n in range(len(entry)):
					if (entry[n] == '1' or
						entry[n] == '2' or
						entry[n] == '3' or
						entry[n] == '4' or
						entry[n] == '5' or
						entry[n] == '6' or
						entry[n] == '7' or
						entry[n] == '8' or
						entry[n] == '9' or
						entry[n] == '0' or
						entry[n] == '<' or
						entry[n] == '>' or
						entry[n] == '.'):			   
						tempstr += entry[n]
				tempsites.append(tempstr)
			Feature['location'] = tempsites
					
		
			#add qualifiers		
			templist = []
			for n in range(2, len(featurelist2[i])): #get all other tags
				templist.append(featurelist2[i][n])
			Feature['qualifiers'] = templist
			
			#add dictionary to list
			features.append(Feature)
			
		self.gbfile['features'] = features
		self.gbfile['filepath'] = filepath
		self.clutter = self.ApEandVNTI_clutter() #check for Vector NTI and ApE clutter and store result

	
	def makegb(self):
		'''Method that creates a new, empty genbank file'''
		self.gbfile['features'] = []
		self.gbfile['dna'] = ''
		self.gbfile['header'] = ''
		self.gbfile['filepath'] = ''


##### Get and Set methods #####
	
	#DNA#
	def get_dna_selection(self):
		'''Method for getting which DNA range is currently selected'''
		return self.dna_selection[0], self.dna_selection[1]

	def set_dna_selection(self, selection):
		'''Method for selecting a certain DNA range'''
		#input needs to be a touple of two values
		self.dna_selection = selection
		start = selection[0]
		finish = selection[1]
		print('Selection from %s to %s') % (start, finish)


	#Feature#
	def get_feature(self, index):
		"""Returns the entire feature from a certain index"""
		try:
			return self.gbfile['features'][index]
		except:
			print('This is not a valid index')
			return False

	def get_feature_index(self, feature):
		'''Get the index of a feature in the self.gbfile data structure'''
		if len(self.gbfile['features']) == 0:
			print('Error, no index found')
			return False
		else:
			for i in range(len(self.gbfile['features'])):
				if self.gbfile['features'][i]['key'] != feature['key']: continue
				if self.gbfile['features'][i]['location'] != feature['location']: continue
				if self.gbfile['features'][i]['qualifiers'][0] != feature['qualifiers'][0]: continue
				if self.gbfile['features'][i]['complement'] == feature['complement']: #it is == here
					return i
			print('Error, no index found')
			return False		

	def get_feature_selection(self):
		'''Get index of currently selected feature, if any'''
		return self.feature_selection 

	def set_feature_selection(self, index):
		'''Set currently selected feature'''
		#this is currently independent from DNA selection
		self.feature_selection = index
		#add logic to find first and last position for feature and make DNA selection match.
		print('Feature "%s" selected') % (self.get_feature_label(self.feature_selection))


	def get_feature_label(self, index):
		"""This method extracts the first qualifier and returns that as a label"""
		try:
			return self.gbfile['features'][index]['qualifiers'][0].split('=')[1]
		except:
			print('This is not a valid index')
			return False
				


	def get_feature_type(self, index):
		'''Get feature type (key) for feature with given index'''
		try:
			return self.gbfile['features'][index]['key']
		except:
			print('This is not a valid index')
			return False

	def set_feature_type(self, feature, newkey):
		'''Changes feature type of the feature passed to method'''
		index = self.get_feature_index(feature)
		if index is False:
			print('Error, no index found')
		else:
			self.gbfile['features'][index]['key'] = newkey

	def get_feature_complement(self, index):
		'''Get whether a feature is on leading or complement DNA strand'''
		try:
			return self.gbfile['features'][index]['complement']
		except:
			print('This is not a valid index')
			return False

	def set_feature_complement(self, feature, complement):
		'''Changes whether a feature is on leading or complement DNA strand'''
		index = self.get_feature_index(feature)
		if index is False:
			print('Error, no index found')
		else:
			self.gbfile['features'][index]['complement'] = complement 


	def get_feature_join(self, index):
		'''Get whether a feature with multiple locations should be joined or not'''
		try:
			return self.gbfile['features'][index]['join']
		except:
			print('This is not a valid index')
			return False

	def set_feature_join(self, feature, join):
		'''Change whether a feature with multiple locations should be joined or not'''
		index = self.get_feature_index(feature)
		if index is False:
			print('Error, no index found')
		else:
			self.gbfile['features'][index]['join'] = join 


	def get_feature_order(self, index):
		'''Get whether a feature with multiple locations should be indicated as having the specified order or not'''
		try:
			return self.gbfile['features'][index]['order']
		except:
			print('This is not a valid index')
			return False

	def set_feature_order(self, feature, order):
		'''Change whether a feature with multiple locations should be indicated as having the specified order or not'''
		index = self.get_feature_index(feature)
		if index is False:
			print('Error, no index found')
		else:
			self.gbfile['features'][index]['order'] = order 


	def get_feature_location(self, index):
		'''Gets all locations for a feature'''
		try:
			return self.gbfile['features'][index]['location']
		except:
			print('This is not a valid index')
			return False

	def set_feature_location(self, feature, newlocation):
		'''Sets all location for a feature'''
		index = self.get_feature_index(feature)
		if index is False:
			print('Error, no index found')
		else:
			self.gbfile['features'][index]['location'] = newlocation

	def get_location(self, entry):
		'''Returns start and end location for an entry of a location list'''
		tempentry = ''
		for n in range(len(entry)):
			if entry[n] != '<' and entry[n] != '>': tempentry += entry[n]
		start, finish = tempentry.split('..')
		return int(start), int(finish)


	def remove_location(self, index, number):
		'''Removes locaiton in self.gbfile['features'][index]['location'][number]'''
		del self.gbfile['features'][index]['location'][number]
		if len(self.gbfile['features'][index]['location']) == 0: # if no locations are left for that feature, delete feature
			del self.gbfile['features'][index] 



##### DNA modification methods #####

	def uppercase(self):
		'''Change DNA selection to uppercase characters'''
		start, finish = self.get_dna_selection()
		string = self.get_dna()[start:finish]
		self.changegbsequence(start, finish, 'r', string.upper())

	def lowercase(self):
		'''Change DNA selection to lowercase characters'''
		start, finish = self.get_dna_selection()
		string = self.get_dna()[start:finish]
		self.changegbsequence(start, finish, 'r', string.lower())


	def reverse_complement_clipboard(self):	
		'''Reverse-complements the DNA and all features in clipboard'''
		self.clipboard['dna'] = dna.reversecomplement(self.clipboard['dna']) #change dna sequence
		pyperclip.copy(self.clipboard['dna'])	
		for i in range(len(self.clipboard['features'])): #checks self.allgbfeatures to match dna change	
			if self.clipboard['features'][i]['complement'] == True: self.clipboard['features'][i]['complement'] = False
			elif self.clipboard['features'][i]['complement'] == False: self.clipboard['features'][i]['complement'] = True

			for n in range(len(self.clipboard['features'][i]['location'])):
				start, finish = self.get_location(self.clipboard['features'][i]['location'][n])
				featurelength = finish - start    #how much sequence in feature?
				trail = len(self.clipboard['dna']) - finish     # how much sequence after feature?

				self.clipboard['features'][i]['location'][n] = self.add_or_subtract_to_locations(self.clipboard['features'][i]['location'][n], -finish+(trail+featurelength+1), 'f')	
				self.clipboard['features'][i]['location'][n] = self.add_or_subtract_to_locations(self.clipboard['features'][i]['location'][n], -start+trail+1, 's')
			self.clipboard['features'][i]['location'].reverse() #reverse order of list elements				


	def reverse_complement_selection(self):
		'''Reverse-complements current DNA selection'''
		start, finish = self.get_dna_selection()
		if start != finish: #must be a selection
			self.copy()
			self.delete()
			self.reverse_complement_clipboard()		
			self.paste()




	def delete(self):
		'''Deletes current DNA selection'''
		start, finish = self.get_dna_selection()
		if start != finish: #must be a selection
			deletedsequence = self.get_dna()[start:finish]
			self.changegbsequence(start+1, finish+1, 'd', deletedsequence)
		self.set_dna_selection((start, start))

	def cut(self):
		'''Cut current DNA selection and place it in clipboard together with any features present on that DNA'''
		start, finish = self.get_dna_selection()
		if start != finish: #must be a selection
			self.copy()
			self.delete()

	def cut_reverse_complement(self):
		start, finish = self.get_dna_selection()
		if start != finish: #must be a selection
			self.copy()
			self.delete()
			self.reverse_complement_clipboard()

	def paste(self):
		'''Paste DNA present in clipboard and any features present on that DNA'''
		start, finish = self.get_dna_selection()
		if start != finish: #If a selection, remove sequence
			self.delete()

		temp_clipboard = deepcopy(self.clipboard) #creates a deep copy which is needed to copy nested lists
		if temp_clipboard['dna'] != pyperclip.paste(): #if internal and system clipboard is not same then system clipboard takes presidence
			print('internal clipboard override')
			temp_clipboard['dna'] = pyperclip.paste()

		DNA = str(temp_clipboard['dna'])
		self.changegbsequence(start+1, start+1, 'i', DNA) #change dna sequence	
		for i in range(len(temp_clipboard['features'])): #add features from clipboard
			self.paste_feature(temp_clipboard['features'][i], start)

	def paste_reverse_complement(self):
		'''Paste reverse complement of DNA in clipboard'''
		self.reverse_complement_clipboard()
		self.paste()
		self.reverse_complement_clipboard() #change it back

	def copy(self):
		'''Copy DNA and all the features for a certain selection'''
		start, finish = self.get_dna_selection()
		if start != finish: #must be a selection
			pyperclip.copy(self.get_dna()[start:finish]) #copy dna to system clipboard (in case I want to paste it somwhere else)
			self.clipboard = {}
			self.clipboard['dna'] = self.get_dna()[start:finish] #copy to internal clipboard
			self.clipboard['features'] = []
			self.allgbfeatures_templist = deepcopy(self.gbfile['features'])
			for i in range(len(self.gbfile['features'])): #checks to match dna change
				if len(self.gbfile['features'][i]['location']) == 1:
					featurestart, featurefinish = self.get_location(self.gbfile['features'][i]['location'])
				else:
					n = 0
					featurestart = self.get_location(self.gbfile['features'][i]['location'][n])[0]
					n = len(self.gbfile['features'][i]['location'])-1
					featurefinish = self.get_location(self.gbfile['features'][i]['location'][n])[1]
				
				if start<featurestart<=featurefinish<=finish: #if change encompasses whole feature
					self.clipboard['features'].append(self.allgbfeatures_templist[i])
					for n in range(len(self.gbfile['features'][i]['location'])):
						newlocation = self.add_or_subtract_to_locations(self.gbfile['features'][i]['location'][n], -start, 'b')
						self.clipboard['features'][-1]['location'][n] = newlocation

	def copy_reverse_complement(self):
		''' '''
		self.copy()
		self.reverse_complement_clipboard()




#### Feature modification methods ####


	def paste_feature(self, feature, insertlocation):
		'''Paste feature from clipboard into an existing genbankfile'''
		#adjust the clipboard location to the insert location and append
		feature = deepcopy(feature)
		for n in range(len(feature['location'])):
			feature['location'][n] = self.add_or_subtract_to_locations(feature['location'][n], insertlocation, 'b')
				
		self.gbfile['features'].append(feature)


	def sort_features():
		'''sort features based on their starting point'''
		#not at all done, needs major changes

		#add the feature at the correct location
#		featurestart, featurefinish = self.get_location(feature['location'][0])
#		for i in range(len(self.gbfile['features'])):
#			start, finish = self.get_location(self.gbfile['features'][i]['location'][0])
#			if featurestart > start:
#				self.gbfile['features'].insert(i, feature)
#				break	



	def add_feature(self, key, qualifiers, location, complement, join, order):
		"""Method adds a new feature to the genbank file"""
		featuretypes = ["modified_base", "variation", "enhancer", "promoter", "-35_signal", "-10_signal", "CAAT_signal", "TATA_signal", "RBS", "5'UTR", "CDS", "gene", "exon", "intron", "3'UTR", "terminator", "polyA_site", "rep_origin", "primer_bind", "protein_bind", "misc_binding", "mRNA", "prim_transcript", "precursor_RNA", "5'clip", "3'clip", "polyA_signal", "GC_signal", "attenuator", "misc_signal", "sig_peptide", "transit_peptide", "mat_peptide", "STS", "unsure", "conflict", "misc_difference", "old_sequence", "LTR", "repeat_region", "repeat_unit", "satellite", "mRNA", "rRNA", "tRNA", "scRNA", "snRNA", "snoRNA", "misc_RNA", "source", "misc_feature", "misc_binding", "misc_recomb", "misc_structure", "iDNA", "stem_loop", "D-loop", "C_region", "D_segment", "J_segment", "N_region", "S_region", "V_region", "V_segment"]
		feature = {}

		if key in featuretypes: 
			feature['key'] = key
		else: 
			print('Key error')
			return False

		if type(qualifiers) == 'list': 
			feature['qualifiers'] = qualifiers
		else:
			print('Qualifiers error')
			return  False
		
		if type(location) == 'list':	#need more checks here to make sure the numbers are ok
			feature['location'] = location
		else:
			print('Location error')
			return False

		if complement == True or complement == False:
			feature['complement'] = complement
		else:
			print('Complement error')
			return False

		if join == True or join == False:
			feature['join'] = join
		else: 
			print('Join error')
			return False
		
		if order == True or order == False:
			feature['order'] = order
		else:
			print('Order error')
			return False

		self.gbfile['features'].append(feature) #change append to sth that works for dicts

		
	def remove_feature(self, feature):
		"""Function removes the feature that is passed to it from the genbank file"""
		position = self.get_feature_index(feature)
		
		if position is False:
			print('feature identify error')
		else:
			del self.gbfile['features'][position]


	def move_feature(self, feature, upordown):
		'''Moves a feature one step up or down the list (up defined as closer to the beginning)'''
		index = self.get_feature_index(feature)
		if upordown == 'u' and index != 0:
			self.gbfile['features'][index-1], self.gbfile['features'][index] = self.gbfile['features'][index], self.gbfile['features'][index-1]

		elif upordown == 'd' and index != len(self.gbfile['features'])-1:
			self.gbfile['features'][index+1], self.gbfile['features'][index] = self.gbfile['features'][index], self.gbfile['features'][index+1]


	def add_qualifier(self, feature, newqualifier):
		'''Adds qualifier tag to existing feature'''
		index = self.get_feature_index(feature)
		if index is False:
			print('Error, no index found')
		else:
			self.gbfile['features'][index]['qualifiers'].append(newqualifier)  #change append to sth that works for dicts
		
		
	def remove_qualifier(self, feature, number):
		'''Removes a qualifier tag from an existing feature'''
		index = self.get_feature_index(feature)
		if index is False:
			print('Error, no index found')
		else:
			del self.gbfile['features'][index]['qualifiers'][number]


	def move_qualifier(self, feature, number, upordown):
		'''Moves a qualifier one step up or down the list (up defined as closer to the beginning)'''
		index = self.get_feature_index(feature)
		if index is False:
			print('Error, no index found')
		else:
			if upordown == 'u' and number != 0:
				self.gbfile['features'][index]['qualifiers'][number-1], self.gbfile['features'][index]['qualifiers'][number] = self.gbfile['features'][index]['qualifiers'][number], self.gbfile['features'][index]['qualifiers'][number-1]

			elif upordown == 'd' and number != len(self.gbfile['features'][index]['qualifiers'])-1:
				self.gbfile['features'][index]['qualifiers'][number+1], self.gbfile['features'][index]['qualifiers'][number] = self.gbfile['features'][index]['qualifiers'][number], self.gbfile['features'][index]['qualifiers'][number+1]


	def ApEandVNTI_clutter(self):
		'''Find out whether there is clutter from Vector NTI or ApE in the genbank file'''
		for i in range(len(self.gbfile['features'])):
			for n in range(len(self.gbfile['features'][i]['qualifiers'])):
				if 'ApEinfo' in self.gbfile['features'][i]['qualifiers'][n] or 'vntifkey' in self.gbfile['features'][i]['qualifiers'][n] :
					return True
		return False


	def clean_clutter(self):
		'''Method for removing ApE- and Vector NTI-specific codes in the qualifiers.'''
		deletionlist = []
		for i in range(len(self.gbfile['features'])):
			for n in range(len(self.gbfile['features'][i]['qualifiers'])):
				if 'ApEinfo' in self.gbfile['features'][i]['qualifiers'][n] or 'vntifkey' in self.gbfile['features'][i]['qualifiers'][n] :
					deletionlist.append((i, n))
#				elif ''\'' in self.gbfile['features'][i]['qualifiers'][n]:
#					self.gbfile['features'][i]['qualifiers'][n] = self.gbfile['features'][i]['qualifiers'][n].replace(''\'', ' ')

		#remove qualifiers based on deletionlist
		while len(deletionlist)>0:
			index, number = deletionlist[-1]
			del self.gbfile['features'][index]['qualifiers'][number]
			del deletionlist[-1]


	def changegbfeatureid(self, oldfeatureid, newfeatureid):
		"""Function changes the ID for a certain feature in self.allgbfeatures"""
		for i in range(len(self.gbfile['features'])):
			if self.gbfile['features'][i]['qualifiers'][0] == '/label='+oldfeatureid:
				self.gbfile['features'][i]['qualifiers'][0] = '/label='+newfeatureid

	def add_or_subtract_to_locations(self, location, changenumber, option):
		'''Given a single location and a number, it adds or subtracts that number to the location'''						
		#get location numbers
		templocations = ''
		brackets = ''					
		for x in range(len(location)):
			if location[x] != '<' and location[x] != '>': templocations += location[x]
			else: brackets += location[x]
		start, finish = templocations.split('..')
		start = int(start)
		finish = int(finish)
		#change them
		if option == 's': #s for start
			start = start + changenumber

		elif option == 'f': #f for finish
			finish = finish + changenumber

		elif option == 'b': #b for both
			start = start + changenumber
			finish = finish + changenumber

		else:
			print('Error, wrong option in add_or_subtract_to_locations')	

		#now add brackets back
		if len(brackets) == 0: pass					
		elif brackets[0] == '<': start = '<' + start
		elif brackets[0] == '>' or brackets[1] == '>': finish = '>' + finish
		location = str(start) + '..' + str(finish)		
		return location		

	
	def get_feature_dna(self, index):
		'''Retrieves the dna sequence for a feature specified by its index'''
		DNA = ''
		for entry in self.gbfile['features'][index]['location']:
			start, finish = self.get_location(entry)
			
			if self.gbfile['features'][index]['join'] == False and self.gbfile['features'][index]['order'] == False:
				DNA = self.gbfile['dna'][start+1:finish+1]
			
			elif self.gbfile['features'][index]['join'] == True:
				DNA += self.gbfile['dna'][start+1:finish+1]
				
			elif self.gbfile['features'][index]['order'] == True: #I probably need to change the reversecomplement function to not remove /n
				DNA += self.gbfile['dna'][start+1:finish+1] + '\n'
		
		if self.gbfile['features'][index]['complement'] == True:
			DNA = dna.reversecomplement(DNA)
		return DNA

	def get_dna(self):
		'''Get the entire DNA sequence from the self.gbfile'''
		return self.gbfile['dna']
	
	def get_filepath(self):
		'''Get the self.filepath for the opened file'''
		return self.gbfile['filepath']
	
	def update_filepath(self, new_path):
		'''Update the self.filepath where a file is saved'''
		self.gbfile['filepath'] = new_path

	def get_features_for_pos(self, position):
		'''For a given DNA position, return any feature entry that spanns that position'''
		features = []
		for feature in self.gbfile['features']:
			for entry in feature['location']: #there may be several starts and finishes in a feature...
				start, finish = self.get_location(entry)
				start -= 1
				finish -= 1

				if start <= position <= finish:
					features.append(feature)
		return features		


	def get_featurename_for_pos(self, position):
		'''For a given dna position, which features are there?'''		
		featurename = ''
		features = self.get_features_for_pos(position)
		for feature in features:
			if featurename != '': featurename += ', '
			featurename += feature['qualifiers'][0].split('=')[1] #removes the '/label=' part from qualifier
		if featurename == '': featurename = 'None'
		return featurename


	def list_features(self):
		'''List all features as a string output'''
		featurelist = ''		
		for i in range(len(self.gbfile['features'])):
			feature = self.gbfile['features'][i]
			if feature['complement'] == True:
				complement = 'complement'
			else:
				complement = 'leading'
			currentfeature = ('>[%s] %s at %s on %s strand\n' % (str(i), feature['qualifiers'][0], feature['location'], complement))
			print(currentfeature)
			print('\n')
			featurelist += currentfeature
		return featurelist


#### Find methods ####

	def find_dna(self, searchstring):
		'''Method for finding a certain DNA sequence in the file. Degenerate codons are supported'''
		dna_seq = self.get_dna()
		dna_seq = dna_seq.upper()
		oligo = searchstring
		oligo = oligo.upper()
		dna_seq = list(dna_seq)
		dna_seq = oligo_localizer.cleaner(dna_seq)

		if oligo=='':
			print 'The searchstring is missing, please check your input'
			self.search_hits = []
		else:
			if oligo_localizer.match_oligo(dna_seq,oligo)!=[]:
				print('The following matches were found: ')
				self.search_hits = []
				for match in oligo_localizer.match_oligo(dna_seq,oligo):
					lm=len(match[2])
					print('from %s to %s %s' % (match[0],match[1]+lm,match[2]))
					self.search_hits.append((match[0],match[1]+lm))
				self.set_dna_selection(self.search_hits[0])
			else:
				print('Sorry, no matches were found')			

	def find_protein():
		'''Method for finding a certain protein sequence in the file. Degenerate codons are supported'''
		pass

	def find_previous(self):
		'''Switch to the previous search hit'''
		start, finish = self.get_dna_selection()
		for i in range(len(self.search_hits)):
			if start < self.search_hits[0][0]:
				self.set_dna_selection(self.search_hits[-1])
				break
			elif start <= self.search_hits[i][0]:
				self.set_dna_selection(self.search_hits[i-1])
				break

	def find_next(self):
		'''Switch to the next search hit'''
		start, finish = self.get_dna_selection()
		for i in range(len(self.search_hits)):
			if start < self.search_hits[i][0]:
				self.set_dna_selection(self.search_hits[i])
				break
			elif i == len(self.search_hits)-1:
				self.set_dna_selection(self.search_hits[0])
				break



#################################

	def mutate():
		pass


	def get_location(self, location):
		'''Takes a location entry and extracts the start and end numbers'''	
		templocation = ''
		for n in range(len(location)): #for each character
			if location[n] != '<' and location[n] != '>': templocation += location[n]
		start, finish = templocation.split('..')
		return int(start), int(finish)


	def get_all_feature_positions(self):
		'''Get type, complement, start and finish for each feature'''		
		positionlist = []
		
		for feature in self.gbfile['features']:
			Key = feature['key']
			Complement = feature['complement']
			for entry in feature['location']:
				start, finish = self.get_location(entry)
				start = int(start)-1
				finish = int(finish)
				positionlist.append([Key, Complement, start, finish])
		return positionlist


# does this work? is it needed?	
	def changegbfeaturepos(self, featureid, newstart, newend, complement, changetype):
		"""Function changes the position of a certain feature in self.allgbfeatures"""
		for i in range(len(self.gbfile['features'])):
			if (('/label=' in featureid) == False and self.gbfile['features'][i]['qualifiers'][0] == '/label=' + featureid) or self.gbfile['features'][i]['qualifiers'][0] == featureid:
				n = 0
				while n <= len(self.gbfile['features'][i]['location']):
					
					#info for current location
					start, finish = get_locations(self.gbfile['features'][i]['qualifiers'][n])
					startbracket, finishbracket = self.gbfile['features'][i]['qualifiers'][n].split('..')
										
					if firstbracket[0] == '<': 
						firstbracket = '<'
					else:
						firstbracket = ''
					
					if secondbracket[0] == '>': 
						secondbracket = '>'
					else:
						secondbracket = ''					
					
					#info for next location
					if n == len(self.gbfile['features'][i]['location']):
						pass
					else:
						nextstart, nextfinish = get_locations(self.gbfile['features'][i]['qualifiers'][n+1])
						nextstartbracket, nextfinishbracket = self.gbfile['features'][i]['qualifiers'][n].split('..')
						if nextstartbracket[0] == '<': 
							nextstartbracket = '<'
						else:
							nextstartbracket = ''
					
						if nextfinishbracket[0] == '>': 
							nextfinishbracket = '>'
						else:
							nextfinishbracket = ''

					
					#if two existing locations overlap, join them
					if start<=nextstart<=finish<=nextfinish:
						self.gbfile['features'][i]['qualifiers'][n] = '%s%s..%s%s' %(startbracket, start, nextfinishbracket, nextfinish)
						del self.gbfile['features'][i]['qualifiers'][n+1]
						n = -1
						
					#add to feature locations
					if changetype == 'a':
						if start<=newstart<finish<=newfinish: #if it overlaps with finish of current location
							self.gbfile['features'][i]['qualifiers'][n] = '%s%s..%s%s' %(firstbracket, start, secondbracket, newfinish)
							n = -1
							
						elif newstart<=start<newfinish<=finish: # if it overlaps with the start of current location
							self.gbfile['features'][i]['qualifiers'][n] = '%s%s..%s%s' %(firstbracket, newstart, secondbracket, finish)
							n = -1
						
						elif start<finish<newstart<newfinish and newstart<newfinish<nextstart<nextfinish: #if it is btween two locations
							self.gbfile['features'][i]['qualifiers'].insert(n+1, '%s..%s' %(newstart, newfinish))
							n = -1
							
		
					#remove coverage from feature locations
					if changetype == 'r':
						pass
						#specify sequence to remove, shorten features if they overlap, if overlap is complete, remove
		
					#move feature locations
					if changetype == 'm':
						pass
						#move all of the locations in a sequence up or down
		

					n += 1
				

#				self.gbfile['features'][i][complement] = complement
				
		
	def changegbsequence(self, changestart, changeend, changetype, change):
		"""Function changes the dna sequence of a .gb file and modifies the feature positions accordingly"""
		#this does not yet handle split features... need to fix that!!!

		if changetype == 'r': #replacement
			self.gbfile['dna'] = self.gbfile['dna'][:changestart] + change + self.gbfile['dna'][changestart+len(change):] #is this correct???
			#need to add feature modifications here!!!!!!
					
		elif changetype == 'i': #insertion
			olddnalength = len(self.gbfile['dna']) #for changing header
			self.gbfile['dna'] = self.gbfile['dna'][:changestart-1] + change + self.gbfile['dna'][changestart-1:]
			self.gbfile['header'] = self.gbfile['header'].replace('%s bp' % olddnalength, '%s bp' % len(self.gbfile['dna'])) #changing header
			for i in range(len(self.gbfile['features'])): #change features already present
				for n in range(len(self.gbfile['features'][i]['location'])):
					start, finish = self.get_location(self.gbfile['features'][i]['location'][n])
					if start<changestart<=finish: #if change is within the feature
						self.gbfile['features'][i]['location'][n] = self.add_or_subtract_to_locations(self.gbfile['features'][i]['location'][n], len(change), 'f')
					elif changestart<=start<=finish: #if change is before feature
						self.gbfile['features'][i]['location'][n] = self.add_or_subtract_to_locations(self.gbfile['features'][i]['location'][n], len(change), 'b')
		
		
		elif changetype == 'd': #deletion
			deletionlist = []
			olddnalength = len(self.gbfile['dna']) #for changing header
			self.gbfile['dna'] = self.gbfile['dna'][:changestart-1] + self.gbfile['dna'][changeend-1:]
			self.gbfile['header'] = self.gbfile['header'].replace('%s bp' % olddnalength, '%s bp' % len(self.gbfile['dna'])) #changing header
			for i in range(len(self.gbfile['features'])): #modifies self.allgbfeatures to match dna change
				for n in range(len(self.gbfile['features'][i]['location'])):
					start, finish = self.get_location(self.gbfile['features'][i]['location'][n])
					if i >= len(self.gbfile['features']):
						break
					if start<=changestart<=int(changeend)<=finish: #if change is within the feature, change finish
						self.gbfile['features'][i]['location'][n] = self.add_or_subtract_to_locations(self.gbfile['features'][i]['location'][n], -len(change), 'f')
					elif changestart<=start<=changeend<=finish: #if change encompasses start, change start and finish
						self.gbfile['features'][i]['location'][n] = self.add_or_subtract_to_locations(self.gbfile['features'][i]['location'][n], changeend-start, 's')						
						self.gbfile['features'][i]['location'][n] = self.add_or_subtract_to_locations(self.gbfile['features'][i]['location'][n], -len(change), 'b')
					elif start<changestart<=finish<=changeend: #if change encompasses finish, change finish
						self.gbfile['features'][i]['location'][n] = self.add_or_subtract_to_locations(self.gbfile['features'][i]['location'][n], -(finish-changestart-1), 'f')	
					elif changestart<=start<=finish<=changeend or changestart==start<=finish==changefinish: #if change encompasses whole feature, add to deletion list
						deletionlist.append(deepcopy((i, n)))
					elif changestart<=changeend<=start<=finish: #if change is before feature, change start and finish
						self.gbfile['features'][i]['location'][n] = self.add_or_subtract_to_locations(self.gbfile['features'][i]['location'][n], -len(change), 'b')				
			#execute deletions (if any)
			while len(deletionlist)>0:
				index, number = deletionlist[-1]
				self.remove_location(index, number)
				del deletionlist[-1]
		else:
			print('%s is not a valid argument for changetype' % changtype)


	def make_gbstring(self):
		'''Prepare data stored in gbfile into one string. Used for saving and for displayig gbfile'''
		string = ''
		string += (self.gbfile['header']+'\n')
		string += ('FEATURES             Location/Qualifiers\n')
		for entry in self.gbfile['features']:
			
			locations = entry['location'][0]
			if len(entry['location']) > 1:
				for i in range(1, len(entry['location'])):
					locations += ',' + entry['location'][i]

				
			if entry['complement'] == True and entry['join'] == True and entry['order'] == False:
				locations = 'complement(join(%s))' % locations
			
			elif entry['complement'] == True and entry['join'] == False and entry['order'] == True:
				locations = 'complement(order(%s))' % locations
			
			elif entry['complement'] == True and entry['join'] == False and entry['order'] == False:
				locations = 'complement(%s)' % locations
			
			
			elif entry['complement'] == False and entry['join'] == True and entry['order'] == False:
				locations = 'join(%s)' % locations
			
			elif entry['complement'] == False and entry['join'] == False and entry['order'] == True:
				locations = 'order(%s)' % locations
			
			elif entry['complement'] == False and entry['join'] == False and entry['order'] == False:
				locations = '%s' % locations			

			string += (('     %s%s%s\n') % (entry['key'], ' '*(16-len(str(entry['key']))), locations))

			for i in range(len(entry['qualifiers'])):
				string += (('                     %s\n') % entry['qualifiers'][i])

		string += ('ORIGIN\n')
		for i in range(0, len(self.gbfile['dna']), 60):
			string += ('%s%d %s %s %s %s %s %s\n' % (' '*(9-len(str(i))), i+1, self.gbfile['dna'][i:i+10], self.gbfile['dna'][i+10:i+20],self.gbfile['dna'][i+20:i+30],self.gbfile['dna'][i+30:i+40],self.gbfile['dna'][i+40:i+50],self.gbfile['dna'][i+50:i+60]))
		string += ('//')		

		return string

	def write_file(self, filepath):
		"""Function writes data stored in header, featruelist and dna to a .gb file"""

		#need to add conditions in case header and features are empty
		string = self.make_gbstring()
		a = open(filepath, 'w')
		a.write(string)
		a.close()




	def protein_mw(self, evt, protein):
		if ('*' in protein[:-1]) == True:
			pass
		else:
			pass
			#add function for calculating MW

	def align_proteins(self, evt):
		pass
		
	def align_proteins_fasta(self, evt):
		pass
	
	def protein_identsim(self, evt):
		pass

	def mutate_positions(self, evt):
		pass
	
	def codon_optimize():
		pass
	
	def protein_to_dna():
		pass
	
	def minilib_design():
		pass
	
	def find_barcode():
		pass

	def insert_barcode():
		pass
	
	def read_barcode():
		pass	
	
	def protein_pattern_find():
		pass


	
	def align_dna_fasta():
		pass	
	
	def restriction_enzyme():
		pass
		#don't forget methylation toggle
#		'''--G  T-C-G-A-C--
#		--|          |--
#		--C-A-G-C-T  G--'''
	
	def seq_analysis():
		#analyze sequencing results
		pass
	
	def gbfile_from_mut():
		pass
		
	def blast_dna():
		pass
	
	def find_hairpins():
		pass
	
	def codon_usage():
		pass
	
	def find_bad_codons():
		pass
	
	def primer_design():
		pass
	
	def primer_database():
		pass	


def open_file(filepath):
	'''Function that makes a gbobject with a specified genbank file'''
	global gb
	gb = gbobject()
	gb.readgb(filepath)


def new_file():
	'''Function that makes a gbobject with an empty genbank file'''
	global gb
	gb = gbobject()
	gb.makegb()


#if __name__ == '__main__': #if script is run by itself and not loaded	
	#import sysargw
	#get the imput...
	#open_file(input)
	#return gb

