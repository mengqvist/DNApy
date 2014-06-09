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
#fix search and mutate
#match features with qualifiers (mandatory and optional)
#add a function that checks that everything is ok
#add a way of actually adding a new feature...
#make changes to how genbank handles /qualifier=xyz, the '=' is not always there...

import dna
from copy import deepcopy
import pyperclip
import oligo_localizer

#for asserts and tests
import types
import unittest

#class feature(object):
#	"""A featue object class that defines the key, location and qualifiers of a feature and methods to interact with these.
#	Data structure is as follows:
#	{key:string #feature key
#		location:list #list of locations on DNA for feature
#		qualifiers:list #list of qualifiers attached to the fature
#		complement:bool #is feature on complement strand or not
#		join:bool #should feature locations be joined or not
#		order:bool #are feature locations in a certain order or not
#		}""" 

#	def __init__(self, inittype, initlocation, initqualifiers, initcomplement, initjoin, initorder):
#		self.SetType(inittype) #the type or "key" of the feature
#		self.SetLocation(initlocation)
#		self.SetQualifiers(initqualifiers)
#		self.SetComplement(initcomplement)
#		self.SetJoin(initjoin)
#		self.SetOrder(initorder)

#	def GetType(self):
#		'''Returns the feature type (its "key")'''
#		return self.type

#	def SetType(self, newtype):
#		'''Sets the feature type (its "key"). Input is a string and must match a feature type as specified by the genbank format.'''
#		assert type(newtype) == str, 'Error, %s is not a string' % str(newtype)
#		assert newtype in ["modified_base", "variation", "enhancer", "promoter", "-35_signal", "-10_signal", "CAAT_signal", "TATA_signal", "RBS", "5'UTR", "CDS", "gene", "exon", "intron", "3'UTR", "terminator", "polyA_site", "rep_origin", "primer_bind", "protein_bind", "misc_binding", "mRNA", "prim_transcript", "precursor_RNA", "5'clip", "3'clip", "polyA_signal", "GC_signal", "attenuator", "misc_signal", "sig_peptide", "transit_peptide", "mat_peptide", "STS", "unsure", "conflict", "misc_difference", "old_sequence", "LTR", "repeat_region", "repeat_unit", "satellite", "mRNA", "rRNA", "tRNA", "scRNA", "snRNA", "snoRNA", "misc_RNA", "source", "misc_feature", "misc_binding", "misc_recomb", "misc_structure", "iDNA", "stem_loop", "D-loop", "C_region", "D_segment", "J_segment", "N_region", "S_region", "V_region", "V_segment"], 'Error, %s is not a valid feature type' % newtype

#		#assert that qualifiers are ok for this one

#		self.type = newtype

#	def GetQualifiers(self):
#		'''Returns a list of qualifiers belonging to the feature.'''
#		return self.qualifiers

#	def SetQualifiers(self, newqualifiers):
#		'''Takes a list of strings and sets which qualifiers belong to the feature'''
#		assert type(newqualifiers) == list, 'Error, %s is not a list.' % str(newqualifiers)
#		for entry in newqualifiers:
#			assert type(entry) == str, 'Error, entry %s in qualifiers is not a string.' % entry

#		#assert that qualifiers are valid for feature type

#	def GetLocations(self):
#		'''Returns a list of locations belonging to the feature.'''
#		return self.location

#	def SetLocations(self, newlocations):
#		'''Takes a list of strings and sets the locations for the feature.'''
#		assert type(newlocations) == list, 'Error, %s is not a list.' % str(newlocations)
#		for entry in newlocations:
#			assert type(entry) == str, 'Error, entry %s in location is not a string.' % entry
#		self.location = newlocations

#	def GetOrder(self):
#		'''Returns a boolean of whether locations are in a certain order or not.'''
#		return self.order

#	def SetOrder(self, neworder):
#		'''Takes a boolean and sets whether locations are in a certain order or not.'''
#		assert type(neworder) == bool, 'Error, %s is not a boolean.' % str(neworder)
#		self.order = neworder	
#		
#	def GetJoin(self):
#		'''Returns a boolian of whether locations should be joined or not.'''
#		return self.join

#	def SetJoin(self, newjoin):
#		'''Takes a boolean and sets whether locations should be joined or not.'''
#		assert type(newjoin) == bool, 'Error, %s is not a boolean.' % str(newjoin)
#		self.join = newjoin

#	def GetComplement(self):
#		'''Returns a boolian of whether feature is on complement strand or not.'''
#		return self.complement

#	def SetComplement(self, newcomplement):
#		'''Takes a boolean and sets whtehr feature is on complement strand or not.'''
#		assert type(newcomplement) == bool, 'Error, %s is not a boolean.' % str(newcomplement)
#		self.complement = newcomplement



class gbobject(object):
	"""Class that reads a genbank file (.gb) and has functions to edit its features and DNA sequence"""
	def __init__(self, filepath = None):
		self.gbfile = {} #this variable stores the whole genbank file
		self.filepath = '' #for keeping track of the file being edited
		self.clipboard = {}
		self.clipboard['dna'] = ''
		self.clipboard['features'] = []

		self.search_hits = []	# variable for storing a list of search hits


		self.file_versions = (deepcopy(self.gbfile),) #stores version of the file
		self.file_version_index = 0 #stores at which index the current file is
		if filepath == None:
			self.gbfile['features'] = []
			self.gbfile['dna'] = ''
			self.gbfile['header'] = ''
			self.gbfile['filepath'] = ''
		else:
			self.readgb(filepath)

###############################


	def treat_input_line(self, tempstr):
		'''Function for parsing a string containing genbank feature information into the correct data format'''
		assert type(tempstr) is types.StringType, "Error parsing genbank line. Input is not a string: %s" % str(tempstr)

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
		assert (type(filepath) == types.StringType or type(filepath) is types.UnicodeType) , "Error opening genbank file. Filepath is not a string: %s" % str(filepath)
#		unittest.assertTrue(type(filepath) is types.StringType or type(filepath) is types.UnicodeType)
		
		try:
			a = open(filepath, 'r') #open it for reading
			gbfile = a.read()
			a.close()

		except IOError:
			print('Error opening file: %s' % str(filepath))

		else:
		
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

	


############# undo and redo functions #################

	def get_file_version(self):
		'''Get the current file version'''
		return self.file_versions[self.get_file_version_index()]

	def add_file_version(self):
		'''Add another file version to the version list.
			This should be added after the current version and should delete all later versions if there are any.'''
		index = self.get_file_version_index()
		if index == len(self.file_versions)-1: #if the current version is the last one
			print('last version.......')
			self.file_versions += (deepcopy(self.gbfile),)
		else:
			self.file_versions = self.file_versions[:index]+(deepcopy(self.gbfile),)
		self.set_file_version_index(index+1)
		print('index', self.get_file_version_index())

	def get_file_version_index(self):
		'''Get the index of the current file version'''
		return self.file_version_index

	def set_file_version_index(self, index):
		'''Set the index of the current file version'''
		assert type(index) == int, 'Error, the index %s is not an integer.' % str(index)
		self.file_version_index = index

	def Undo(self):
		'''Function for undoing user action (such as deleting or adding dna or features).
			Function intended for use by the user.'''
		if self.get_file_version_index <= 0:
			print("Can't undo there are no previous versions")
		else:
			old_index = self.get_file_version_index()
			new_index = old_index-1
			self.set_file_version_index(new_index)
			self.gbfile = self.get_file_version()
			print('index after undo', self.get_file_version_index())

	def Redo(self):
		'''Function for redoing user action (such as deleting or adding dna or features).
			Function intended for use by the user.'''
		pass


####################################################################


##### Get and Set methods #####
	



	#Feature#
	def get_all_features(self):
		return self.gbfile['features']

	def get_feature(self, index):
		"""Returns the entire feature from a certain index"""
		num_features = len(self.get_all_features()) #number of features already present
		if index == -1:
			index = num_features-1
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


	def get_feature_label(self, index):
		"""This method extracts the first qualifier and returns that as a label.
			Index should be an integer."""
		assert type(index) == int, "Error, index must be an integer."
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
			if entry[n] != '<' and entry[n] != '>': 
				tempentry += entry[n]
		start, finish = tempentry.split('..')
		return int(start), int(finish)

	def GetFirstLastLocation(self, feature):
		'''Returns ultimate first and ultimate last position of a feature, regardless of how many pieces it is broken into'''
		index = self.get_feature_index(feature)
		if index is False:
			print('Error, no index found')
		else:
			locations = self.gbfile['features'][index]['location']
			start = locations[0].split('..')[0]
			finish = locations[-1].split('..')[1]
			print(start, finish)
			return int(start), int(finish)

	def remove_location(self, index, number):
		'''Removes locaiton in self.gbfile['features'][index]['location'][number]'''
		del self.gbfile['features'][index]['location'][number]
		if len(self.gbfile['features'][index]['location']) == 0: # if no locations are left for that feature, delete feature
			del self.gbfile['features'][index] 

	def IsValidLocation(self, locationlist):
		'''Takes a location list and tests whether it is valid'''
		result = True		
		try:
			assert type(locationlist) == types.ListType
			for location in locationlist:
				start, finish = location.split('..')
				start = int(start)
				finish = int(finish)
				assert (start == 0 and finish == 0) == False
				assert start <= finish
				assert finish <= len(self.GetDNA())
		except:
			result = False
		return result

	def get_qualifiers(self, index):
		'''Returns all qualifiers for a feature'''
		try:
			return self.gbfile['features'][index]['qualifiers']
		except:
			raise ValueError('This is not a valid index')

	def get_qualifier(self, index, number):
		'''Returns specified qualifier for specified feature'''
		try:
			return self.gbfile['features'][index]['qualifiers'][number][1:].split('=')
		except:
			raise ValueError('Index or number is not valid')

	def set_qualifier(self, index, number, qualifier, tag):
		'''Sets the qualifier and tag for a given qualifier'''
		assert type(index) is types.IntType, "Index is not an integer: %s" % str(index)
		assert type(number) is types.IntType, "Number is not an integer: %s" % str(number)
		assert type(qualifier) is types.StringType, "Qualifier is not a string: %s" % str(qualifier)
		assert type(tag) is types.StringType, "Tag is not a string: %s" % str(tag)
		try:
			self.gbfile['features'][index]['qualifiers'][number] = '/%s=%s' % (qualifier, tag)
		except:
			raise IOError('Error setting qualifier')


##### DNA modification methods #####

	def Upper(self, start, finish):
		'''Change DNA selection to uppercase characters'''
		assert (type(start) == int and type(finish) == int), 'Function requires two integers.'
		assert start <= finish, 'Startingpoint must be before finish'
		string = self.GetDNA(start, finish)
		self.changegbsequence(start, finish, 'r', string.upper())
		self.add_file_version()

	def Lower(self, start, finish):
		'''Change DNA selection to lowercase characters'''
		assert (type(start) == int and type(finish) == int), 'Function requires two integers.'
		assert start <= finish, 'Startingpoint must be before finish'
		string = self.GetDNA(start, finish)
		self.changegbsequence(start, finish, 'r', string.lower())
		self.add_file_version()

	def reverse_complement_clipboard(self):	
		'''Reverse-complements the DNA and all features in clipboard'''
		self.clipboard['dna'] = dna.RC(self.clipboard['dna']) #change dna sequence
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


	def RCselection(self, start, finish):
		'''Reverse-complements current DNA selection'''
		assert (type(start) == int and type(finish) == int), 'Function requires two integers.'
		assert start <= finish, 'Startingpoint must be before finish'
		self.Copy(start, finish)
		self.reverse_complement_clipboard()	
		deletedsequence = self.GetDNA(start, finish)
		self.changegbsequence(start, finish, 'd', deletedsequence)
		self.Paste(start)

	def Delete(self, start, finish, visible=True):
		'''Deletes current DNA selection.
			Start and finish should be integers.
			The optional variable 'hidden' can be set to True or False. 
			If set to True, it is a hidden deletion that does not trigger other events.
			If set to False, it does trigger other events.'''
		assert (type(start) == int and type(finish) == int), 'Function requires two integers.'
		assert start <= finish, 'Startingpoint must be before finish'
		deletedsequence = self.GetDNA(start, finish)
		self.changegbsequence(start, finish, 'd', deletedsequence)
		if visible == True:
			self.add_file_version()

	def Cut(self, start, finish):
		'''Cuts selection and place it in clipboard together with any features present on that DNA'''
		assert (type(start) == int and type(finish) == int), 'Function requires two integers.'
		assert start <= finish, 'Startingpoint must be before finish'
		self.Copy(start, finish)
		deletedsequence = self.GetDNA(start, finish)
		self.changegbsequence(start, finish, 'd', deletedsequence)
		self.add_file_version()

	def CutRC(self, start, finish):
		'''Cuts the reverese-complement of a selection and place it in clipboard together with any features present on that DNA'''
		assert (type(start) == int and type(finish) == int), 'Function requires two integers.'
		assert start <= finish, 'Startingpoint must be before finish'
		self.Cut(start, finish)
		self.reverse_complement_clipboard()

	def Paste(self, ip):
		'''Paste DNA present in clipboard and any features present on that DNA'''
		assert type(ip) == int, 'The insertion point must be an integer.'

		temp_clipboard = deepcopy(self.clipboard) #creates a deep copy which is needed to copy nested lists
		if temp_clipboard['dna'] != pyperclip.paste(): #if internal and system clipboard is not same then system clipboard takes presidence
			print('internal clipboard override')
			temp_clipboard['dna'] = pyperclip.paste()

		DNA = str(temp_clipboard['dna'])
		self.changegbsequence(ip, ip, 'i', DNA) #change dna sequence	
		for i in range(len(temp_clipboard['features'])): #add features from clipboard
			self.paste_feature(temp_clipboard['features'][i], ip)
		self.add_file_version()

	def PasteRC(self, ip):
		'''Paste reverse complement of DNA in clipboard'''
		assert type(ip) == int, 'The insertion point must be an integer.'
		self.reverse_complement_clipboard()
		self.Paste(ip)
		self.reverse_complement_clipboard() #change it back


	def Copy(self, start, finish):
		'''Copy DNA and all the features for a certain selection'''
		assert (type(start) == int and type(finish) == int), 'Function requires two integers.'
		assert start <= finish, 'Startingpoint must be before finish'
		pyperclip.copy(self.GetDNA(start, finish)) #copy dna to system clipboard (in case I want to paste it somwhere else)
		self.clipboard = {}
		self.clipboard['dna'] = self.GetDNA(start, finish) #copy to internal clipboard
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
			
			if (start<=featurestart and featurefinish<=finish) == True: #if change encompasses whole feature
				self.clipboard['features'].append(self.allgbfeatures_templist[i])
				for n in range(len(self.gbfile['features'][i]['location'])):
					newlocation = self.add_or_subtract_to_locations(self.gbfile['features'][i]['location'][n], -start, 'b')
					self.clipboard['features'][-1]['location'][n] = newlocation

	def CopyRC(self, start, finish):
		'''Copy the reverse complement of DNA and all the features for a certain selection'''
		assert (type(start) == int and type(finish) == int), 'Function requires two integers.'
		assert start <= finish, 'Startingpoint must be before finish'
		self.Copy(start, finish)
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

		if type(qualifiers) == list: 
			feature['qualifiers'] = qualifiers
		else:
			print('Qualifiers error')
			return  False
		
		if type(location) == list:	#need more checks here to make sure the numbers are ok
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
			DNA = dna.RC(DNA)
		return DNA

	def GetDNA(self, start=1, finish=-1):
		'''Get the entire DNA sequence from the self.gbfile'''
		if (start == 1 and finish == -1) == True:
			return self.gbfile['dna']
		else:	
			assert (type(start) == int and type(finish) == int), 'Function requires two integers.'
			assert start <= finish, 'Starting point must be before finish.'
			assert start > 0 and start <= len(self.gbfile['dna']), 'Starting point must be between 1 and the the DNA length.'
			assert finish > 0 and finish <= len(self.gbfile['dna']), 'Finish must be between 1 and the the DNA length.'
			return self.gbfile['dna'][start-1:finish]
	
	def get_filepath(self):
		'''Get the self.filepath for the opened file'''
		return self.gbfile['filepath']
	
	def set_filepath(self, new_path):
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


	def ListFeatures(self):
		'''List all features as a string output'''
		featurelist = ''		
		for i in range(len(self.gbfile['features'])):
			feature = self.gbfile['features'][i]
			if feature['complement'] == True:
				complement = 'complement'
			else:
				complement = 'leading'
			currentfeature = '>[%s] %s at %s on %s strand' % (str(i), feature['qualifiers'][0], feature['location'], complement)
			print(currentfeature)
#			print('\n')
#			featurelist += currentfeature
#		return featurelist


#### Find methods ####

	def FindNucleotide(self, searchstring, searchframe):
		'''Method for finding a certain DNA sequence. Degenerate codons are supported.
Searchstring should be a string of numbers or DNA characers and searchframe should be an index for a feature. 
-1 indicates search in entire genbank file
indeces >-1 are feature indeces'''

		#fix the set_dna_selection functions here
	
		assert type(searchstring) == str or type(searchstring) == unicode
		assert type(searchframe) == int

		if searchstring=='':
			print 'The searchstring is missing, please check your input'
			return None

		elif searchstring.isdigit(): #if search is numbers only
			if searchframe == -1: #if index is -1, that means search in molecule
				search_hits = [(int(searchstring)-1, int(searchstring))]
				return search_hits
#				self.set_dna_selection(self.search_hits[0])
			else: #otherwise search in feature indicated by the index
				complement = self.get_feature_complement(searchframe) # is feature complement or not
				feature = self.get_feature(searchframe)
				first, last = self.GetFirstLastLocation(feature)

				#add check ensuring that you cannot search outsite of feature bounds	
				if complement == True:
					search_hits = [(last+1 - int(searchstring)-1, last+1 - int(searchstring))]
					return search_hits
#					self.set_dna_selection(self.search_hits[0])					
				elif complement == False:
					search_hits = [(first-1 + int(searchstring)-1, first-1 + int(searchstring))]
					return search_hits
#					self.set_dna_selection(self.search_hits[0])					
				
		else: #if searchstring is not numbers
			#need to test that the characters are valid
			if searchframe == -1: #if index is -1, that means search in molecule
				dna_seq = self.GetDNA()
				dna_seq = dna_seq.upper()
				searchstring = searchstring.upper()
				dna_seq = list(dna_seq)
				dna_seq = oligo_localizer.cleaner(dna_seq)
				if oligo_localizer.match_oligo(dna_seq, searchstring)!=[]:
					search_hits = []
					for match in oligo_localizer.match_oligo(dna_seq, searchstring):
						lm=len(match[2])
						search_hits.append((match[0],match[1]+lm))
					return search_hits
#					self.set_dna_selection(self.search_hits[0])
				else:
					print('No matches were found')			
					return None

	def FindAminoAcid(self, searchsting, searchframe):
		'''Method for finding a certain protein sequence, or position, in the file. Degenerate codons are supported'''
		pass
		return None


	def FindFeature(self, searchstring):
		'''Method for finding a certain feature in a genbank file'''
		if searchstring=='':
			print 'The searchstring is missing, please check your input'
			return None
		else:
			search_hits = []
			hits = []
			features = self.get_all_features()
			for i in range(len(features)):
				qualifiers = self.get_qualifiers(i)
				for qualifier in qualifiers:
					if searchstring in qualifier.split('=')[1] and (features[i] in hits) == False:
						hits.append(features[i])

		if len(hits) == 0:
			print('Sorry, no matches were found')
			return None
		else:
			for feature in hits:
				start, finish = self.GetFirstLastLocation(feature)
				search_hits.append((start-1, finish))
		search_hits = sorted(search_hits)	
		return search_hits	
#		self.set_dna_selection(self.search_hits[0])


	def find_previous(self):
		'''Switch to the previous search hit'''
		start, finish = self.get_dna_selection()
		for i in range(len(self.search_hits)):
			print('start', start)
			print('hits', self.search_hits[0][0])
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
		"""Function changes the dna sequence of a .gb file and modifies the feature positions accordingly."""
		#assumes that changestart and changeend are list positions
#		changestart -= 1 #convert back from DNA numbering to string numbering
#		changeend += 1
		print(change)
		print(changestart, changeend)
		if changetype == 'r': #replacement
			self.gbfile['dna'] = self.gbfile['dna'][:changestart-1] + change + self.gbfile['dna'][changestart-1+len(change):] #is this correct???
			
					
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
			self.gbfile['dna'] = self.gbfile['dna'][:changestart-1] + self.gbfile['dna'][changeend:]
			self.gbfile['header'] = self.gbfile['header'].replace('%s bp' % olddnalength, '%s bp' % len(self.gbfile['dna'])) #changing header
			for i in range(len(self.gbfile['features'])): #modifies self.allgbfeatures to match dna change
				for n in range(len(self.gbfile['features'][i]['location'])):
					start, finish = self.get_location(self.gbfile['features'][i]['location'][n])
					if i >= len(self.gbfile['features']):
						break
					if (start<=changestart and changeend<finish) or (start<changestart and changeend<=finish): #if change is within the feature, change finish
						self.gbfile['features'][i]['location'][n] = self.add_or_subtract_to_locations(self.gbfile['features'][i]['location'][n], -len(change), 'f') #finish
						print('within feature')
					elif changestart<start and start<=changeend<finish: #if change encompasses start, change start and finish
						self.gbfile['features'][i]['location'][n] = self.add_or_subtract_to_locations(self.gbfile['features'][i]['location'][n], changeend+1-start, 's')	#start					
						self.gbfile['features'][i]['location'][n] = self.add_or_subtract_to_locations(self.gbfile['features'][i]['location'][n], -len(change), 'b') #both
						print('encompass start')
					elif start<changestart<=finish and finish<changeend: #if change encompasses finish, change finish
						self.gbfile['features'][i]['location'][n] = self.add_or_subtract_to_locations(self.gbfile['features'][i]['location'][n], -(finish-changestart), 'f')	
						print('encompass finish')
					elif changestart<=start and finish<changeend: #if change encompasses whole feature, add to deletion list
						deletionlist.append(deepcopy((i, n)))
						print('encompass all')
					elif changestart<start and changeend<start: #if change is before feature, change start and finish
						self.gbfile['features'][i]['location'][n] = self.add_or_subtract_to_locations(self.gbfile['features'][i]['location'][n], -len(change), 'b')				
						print('before start')
			#execute deletions (if any)
			while len(deletionlist)>0:
				index, number = deletionlist[-1]
				self.remove_location(index, number)
				del deletionlist[-1]
		else:
			print('%s is not a valid argument for changetype' % changetype)


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





#if __name__ == '__main__': #if script is run by itself and not loaded	
#	import sys
#	assert (len(sys.argv) > 0 and len(sys.argv) <= 2), 'Error, this script takes zero or one arguments.'

#	if len(sys.argv) == 1:
#		print('No file specified. Creting new genbank file.')
#		#then create file

#	elif len(sys.argv) == 2:
#		print('Opening %s' % str(sys.argv[1]))
#		#then open file




