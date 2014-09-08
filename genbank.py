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
#improve spead at which large genbank files are loaded
#fix header parsing
#match features with qualifiers (mandatory and optional)
#add a function that checks that everything is ok
#make changes to how genbank handles /qualifier=xyz, the '=' is not always there...

import dna
import copy
import pyperclip
import oligo_localizer
import peptide_localizer
import re
import sys

class feature(object):
	"""A feature object class that defines the key, location and qualifiers of a feature and methods to interact with these.
	Data structure is as follows:
	{key:string #feature key
		location:list #list of locations on DNA for feature
		qualifiers:list #list of qualifiers attached to the fature
		complement:bool #is feature on complement strand or not
		join:bool #should feature locations be joined or not
		order:bool #are feature locations in a certain order or not
		}""" 

	def __init__(self, inittype, initlocation, initqualifiers, initcomplement, initjoin, initorder):
		self.SetType(inittype) #the type or "key" of the feature
		self.SetLocation(initlocation)
		self.SetQualifiers(initqualifiers)
		self.SetComplement(initcomplement)
		self.SetJoin(initjoin)
		self.SetOrder(initorder)

	def GetType(self):
		'''Returns the feature type (its "key")'''
		return self.type

	def SetType(self, newtype):
		'''Sets the feature type (its "key"). Input is a string and must match a feature type as specified by the genbank format.'''
		assert type(newtype) == str, 'Error, %s is not a string' % str(newtype)
		assert newtype in ["modified_base", "variation", "enhancer", "promoter", "-35_signal", "-10_signal", "CAAT_signal", "TATA_signal", "RBS", "5'UTR", "CDS", "gene", "exon", "intron", "3'UTR", "terminator", "polyA_site", "rep_origin", "primer_bind", "protein_bind", "misc_binding", "mRNA","ncRNA", "prim_transcript", "precursor_RNA", "5'clip", "3'clip", "polyA_signal", "GC_signal", "attenuator", "misc_signal", "sig_peptide", "transit_peptide", "mat_peptide", "STS", "unsure", "conflict", "misc_difference", "old_sequence", "LTR", "repeat_region", "repeat_unit", "satellite", "mRNA", "rRNA", "tRNA", "scRNA", "snRNA", "snoRNA", "misc_RNA", "source", "misc_feature", "misc_binding", "misc_recomb", "misc_structure", "iDNA", "stem_loop", "D-loop", "C_region", "D_segment", "J_segment", "N_region", "S_region", "V_region", "V_segment"], 'Error, %s is not a valid feature type' % newtype

		#assert that qualifiers are ok for this one

		self.type = newtype

	def GetQualifiers(self):
		'''Returns a list of qualifiers belonging to the feature.'''
		return self.qualifiers

	def SetQualifiers(self, newqualifiers):
		'''Takes a list of strings and sets which qualifiers belong to the feature'''
		assert type(newqualifiers) == list, 'Error, %s is not a list.' % str(newqualifiers)
		for entry in newqualifiers:
			assert type(entry) == str, 'Error, entry %s in qualifiers is not a string.' % entry

		#assert that qualifiers are valid for feature type

	def GetLocations(self):
		'''Returns a list of locations belonging to the feature.'''
		return self.location

	def SetLocations(self, newlocations):
		'''Takes a list of strings and sets the locations for the feature.'''
		assert type(newlocations) == list, 'Error, %s is not a list.' % str(newlocations)
		for entry in newlocations:
			assert type(entry) == str, 'Error, entry %s in location is not a string.' % entry
		self.location = newlocations

	def GetOrder(self):
		'''Returns a boolean of whether locations are in a certain order or not.'''
		return self.order

	def SetOrder(self, neworder):
		'''Takes a boolean and sets whether locations are in a certain order or not.'''
		assert type(neworder) == bool, 'Error, %s is not a boolean.' % str(neworder)
		self.order = neworder	
		
	def GetJoin(self):
		'''Returns a boolian of whether locations should be joined or not.'''
		return self.join

	def SetJoin(self, newjoin):
		'''Takes a boolean and sets whether locations should be joined or not.'''
		assert type(newjoin) == bool, 'Error, %s is not a boolean.' % str(newjoin)
		self.join = newjoin

	def GetComplement(self):
		'''Returns a boolian of whether feature is on complement strand or not.'''
		return self.complement

	def SetComplement(self, newcomplement):
		'''Takes a boolean and sets whtehr feature is on complement strand or not.'''
		assert type(newcomplement) == bool, 'Error, %s is not a boolean.' % str(newcomplement)
		self.complement = newcomplement



class gbobject(object):
	"""Class that reads a genbank file (.gb) and has functions to edit its features and DNA sequence"""
	def __init__(self, filepath = None):
		self.gbfile = {} #this variable stores the whole genbank file
		self.clipboard = {}
		self.clipboard['dna'] = ''
		self.clipboard['features'] = []
		self.featuretypes = ["modified_base", "variation", "enhancer", "promoter", "-35_signal", "-10_signal", "CAAT_signal", "TATA_signal", "RBS", "5'UTR", "CDS", "gene", "exon", "intron", "3'UTR", "terminator", "polyA_site", "rep_origin", "primer_bind", "protein_bind", "misc_binding", "mRNA", "prim_transcript", "precursor_RNA", "5'clip", "3'clip", "polyA_signal", "GC_signal", "attenuator", "misc_signal", "sig_peptide", "transit_peptide", "mat_peptide", "STS", "unsure", "conflict", "misc_difference", "old_sequence", "LTR", "repeat_region", "repeat_unit", "satellite", "mRNA", "rRNA", "tRNA", "scRNA", "snRNA", "snoRNA", "misc_RNA", "source", "misc_feature", "misc_binding", "misc_recomb", "misc_structure", "iDNA", "stem_loop", "D-loop", "C_region", "D_segment", "J_segment", "N_region", "S_region", "V_region", "V_segment"]


		self.search_hits = []	# variable for storing a list of search hits

		self.file_versions = (copy.deepcopy(self.gbfile),) #stores version of the file
		self.file_version_index = 0 #stores at which index the current file is

		#space constants for where info starts
		self.space_header = 12*' '
		self.space_feature = 21*' '
		self.space_DNA = 10*' '
		
		#set up header info
		self.gbfile['header'] = {}
		self.gbfile['header']['locus'] = []
		self.gbfile['header']['definition'] = ''
		self.gbfile['header']['accession'] = ''
		self.gbfile['header']['version'] = []
		self.gbfile['header']['dblink'] = []
		self.gbfile['header']['keywords'] = []
		self.gbfile['header']['source'] = []
		self.gbfile['header']['references'] = []
		self.gbfile['header']['comments'] = []

		#set up feature info
		self.gbfile['features'] = []
		
		#set up dna info
		self.gbfile['dna'] = ''
		
		#set up file info
		self.gbfile['filepath'] = filepath
		self.fileName = 'New DNA' #name of the file/plasmid name

		if filepath == None:
			pass
		else:
			platform = sys.platform
			if platform == 'win32':
				filepath = filepath.replace('\\', '/')
			self.readgb(filepath)

###############################


	def readgb(self, filepath):
		"""Function takes self.filepath to .gb file and extracts the header, features and DNA sequence"""
		assert type(filepath) == str or type(filepath) == unicode , "Error opening genbank file. Filepath is not a string: %s" % str(filepath)
		
		self.fileName = filepath.split('/')[-1] #update the fileName variable
		file_read = False #indicate whether the file was read
		
		#open it and read line by line (very memory efficient)
		with open(filepath) as infile:
			#to keep track where in the document I currently am. 
			#in the header section both sect_feature and sect_origin are false
			#in the feature section sect_feature is True and the other False
			#in the DNA section both are True
			sect_feature = False
			sect_origin = False
			whole_line = '' #many things are broken over several lines, use this variable to collect them.
			for line in infile:
				if line[0:8] == 'FEATURES': #indicates start of feature section
					self.parse_header_line(whole_line)
					sect_feature = True
					whole_line = ''

				elif line[0:6] == 'ORIGIN': #indicates start of DNA section
					self.parse_feature_line(whole_line)
					sect_origin = True
					self.gbfile['dna'] = [] #temporarily set this up as a list. This is only while parsing the DNA and it will be converted to a string after.
					whole_line = ''

				elif line[0:2] == '//' and sect_feature is True and sect_origin is True: #this indicates the end of the file
					self.parse_dna_line(whole_line)
					self.gbfile['dna'] = ''.join(self.gbfile['dna']) #while parsing the DNA was kept as a list, now make it a string
					file_read = True
					
				elif sect_feature is False and sect_origin is False:
					if line[0:1] == ' ': #this line is a continuation of a previous header section entry
						whole_line += line
					elif whole_line == '': #start of a new entry
						whole_line = line
					else: #the entry is over and needs to be parsed
						self.parse_header_line(whole_line)
						whole_line = line

				elif sect_feature is True and sect_origin is False:
					if line[0:21] == self.space_feature: #this line is a continuation of a previous feature entry
						whole_line += line
					elif whole_line == '': #start of a new entry
						whole_line = line
					else: #the entry is over and needs to be parsed
						self.parse_feature_line(whole_line)
						whole_line = line
						
				elif sect_feature is True and sect_origin is True:
					self.parse_dna_line(line)

		assert file_read is True, 'There was an unknown error reading the file.'
		self.add_file_version()

		
	def parse_header_line(self, line):
		'''
		Parses a header line into the correct data format.
		The term 'line' is used loosly and can actually be comprised of several lines that together describe a header entry subpart.
		'''
		
		print('header', line)
		if 'LOCUS' in line[0:12]:
			pass
		elif 'DEFINITION' in line[0:12]:
			pass
		elif 'ACCESSION' in line[0:12]:
			pass
		elif 'VERSION' in line[0:12]:
			pass
		elif 'DBLINK' in  line[0:12]:
			pass
		elif 'KEYWORDS' in line[0:12]:
			pass
		elif 'SOURCE' in line[0:12]:
			pass
		elif 'REFERENCE' in line[0:12]:
			pass
		elif 'COMMENTS' in line[0:12]:
			pass

			
	def parse_feature_line(self, inputline):
		'''
		Takes a feature line consisting info for an entire feature and parses them into the DNApy data format.
		This method only works if only one feature is to be parsed at a time.
		'''
		assert type(inputline) == str, 'Error, the input has to be a string.'
		
		#first I need to go through the input and see if there are any info that wraps over several lines
		linelist = inputline.split('\n')
		wholelinelist = [] #for storing the complete un-wrapped lines
		whole_line = ''
		for line in linelist:
			if whole_line == '': #to catch the first line
				whole_line = line		
			elif line[0:21] == 21*' ' and line[21] != '/': #continuation of line which is broken on several lines
				whole_line += line
			elif line == '': #catch empty lines (that should not be there in any case) and move to next line
				pass
			elif line[0:21] == 21*' ' and line[21] == '/': #new qualifier line
				wholelinelist.append(whole_line)
				whole_line = line
		wholelinelist.append(whole_line) #catch the last one

		
		#now parse the whole lines into the correct data structure
		feature = {}
		feature['qualifiers'] = []
		for line in wholelinelist:
			assert type(line) is str, "Error parsing genbank line. Input line is not a string, it is %s:" % str(type(line))
		
			#the line either starts with 5 blanks and then a feature key (feature key line), or it start with 21 blanks and then / (qualifier line)
			if line[0:5] == 5*' ' and (line[6] in 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ-') == True:
				key = line[5:20].rstrip('\t\n\x0b\x0c\r ')
				location = line[21:].rstrip('\t\n\x0b\x0c\r ')
				feature['key'] = key
				feature['location'], feature['complement'], feature['join'], feature['order'] = self.parse_location(location)

			elif line[0:21] == '                     ' and line[21] == '/':
				qualifier = line[21:].rstrip('\t\n\x0b\x0c\r ')
				feature['qualifiers'].append(qualifier)

			else:
				print('error', line)		

		self.gbfile['features'].append(copy.deepcopy(feature)) #add feature to the global data structure
		self.clutter = self.ApEandVNTI_clutter() #check the stored features for Vector NTI and ApE clutter and store result

		

	def parse_location(self, locationstring):
		'''get whether complement or not, join or not, order or not'''
		assert type(locationstring) == str, "Error, locationstring must be a string"
		if 'complement' in locationstring:
			complement = True
		else:
			complement = False

		if 'join' in locationstring:
			join = True
		else:
			join = False
	
		if 'order' in locationstring:
			order = True
		else:
			order = False
						
		#get location numbering
		## need to add single base support!!!! ##
		tempsites = []
		commasplit = locationstring.split(',')
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
		location = tempsites
			
		return location, complement, join, order

	def parse_dna_line(self, line):
		#clean DNA from numbering and spaces and append it to a list. This list will later be converted to a string.
		self.gbfile['dna'].append(dna.CleanDNA(line, silent=True) )
	
	
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



############# undo and redo functions #################

	def get_file_version(self):
		'''Get the current file version'''
		return self.file_versions[self.get_file_version_index()]

	def add_file_version(self):
		'''Add another file version to the version list.
			This should be added after the current version and should delete all later versions if there are any.'''
		index = self.get_file_version_index()
		if index == len(self.file_versions)-1: #if the current version is the last one
			self.file_versions += (copy.deepcopy(self.gbfile),)
		else:
			self.file_versions = self.file_versions[0:index]+(copy.deepcopy(self.gbfile),)
		self.set_file_version_index(index+1)

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
		if self.get_file_version_index() <= 0:
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
		if self.get_file_version_index() >= len(self.file_versions)-1:
			print("Can't redo, there are no later versions")
		else:
			old_index = self.get_file_version_index()
			new_index = old_index+1
			self.set_file_version_index(new_index)
			self.gbfile = self.get_file_version()
			print('index after redo', self.get_file_version_index())


####################################################################


##### Get and Set methods #####
	



	#Feature#
	def get_all_features(self):
		return self.gbfile['features']

	def get_feature(self, index):
		"""Returns the entire feature from a certain index"""
		assert -1<=index<len(self.get_all_features()), 'This is not a valid index'
		return self.gbfile['features'][index]

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
			self.add_file_version()

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
			self.add_file_version() 


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
			self.add_file_version()


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
			self.add_file_version()


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
			self.add_file_version()

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
			assert type(locationlist) == list
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
		assert type(index) is int, "Index is not an integer: %s" % str(index)
		assert type(number) is int, "Number is not an integer: %s" % str(number)
		assert type(qualifier) is str, "Qualifier is not a string: %s" % str(qualifier)
		assert type(tag) is str, "Tag is not a string: %s" % str(tag)
		try:
			self.gbfile['features'][index]['qualifiers'][number] = '/%s=%s' % (qualifier, tag)
			self.add_file_version()
		except:
			raise IOError('Error setting qualifier')


##### DNA modification methods #####

	def Upper(self, start=1, finish=-1):
		'''Change DNA selection to uppercase characters.
			Start and finish are optional arguments (integers). If left out change applies to the entire DNA molecule.
			If specified, start and finish determines the range for the change.'''
		if finish == -1:
			finish = len(self.GetDNA())
		assert (type(start) == int and type(finish) == int), 'Function requires two integers.'
		assert start <= finish, 'Startingpoint must be before finish'
		string = self.GetDNA(start, finish)
		self.changegbsequence(start, finish, 'r', string.upper())
		self.add_file_version()

	def Lower(self, start=1, finish=-1):
		'''Change DNA selection to lowercase characters.
			Start and finish are optional arguments (integers). If left out change applies to the entire DNA molecule.
			If specified, start and finish determines the range for the change.'''
		if finish == -1:
			finish = len(self.GetDNA())
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
			if self.clipboard['features'][i]['complement'] == True: 
				self.clipboard['features'][i]['complement'] = False
			elif self.clipboard['features'][i]['complement'] == False: 
				self.clipboard['features'][i]['complement'] = True

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
		self.Delete(start, finish, visible=False)
		self.Paste(start)

	def Delete(self, start, finish, visible=True):
		'''Deletes current DNA selection.
			Start and finish should be integers.
			The optional variable 'hidden' can be set to True or False. 
			If set to True, no file versions are added to the undo/redo record.
			If set to False, it does add file versions to the undo/redo record.'''
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


	def Paste(self, ip, DNA=None):
		'''Paste DNA present in clipboard and any features present on that DNA
			If a string is passed to DNA then this will over-ride anything that is present in the clipboard'''
		assert type(ip) == int, 'The insertion point must be an integer.'

		if DNA == None:
			system_clip = re.sub(r'\s+', '', pyperclip.paste()) #remove whitespace from item in system clipboard
			if self.clipboard['dna'] != system_clip: #if internal and system clipboard is not same, then system clipboard takes presidence
				self.clipboard['dna'] = system_clip
				self.clipboard['features'] = []
			DNA = copy.copy(self.clipboard['dna'])
			self.changegbsequence(ip, ip, 'i', DNA) #change dna sequence	
			for i in range(len(self.clipboard['features'])): #add features from clipboard
				self.paste_feature(self.clipboard['features'][i], ip-1)
		else:
			self.changegbsequence(ip, ip, 'i', DNA) #change dna sequence
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
		self.allgbfeatures_templist = copy.deepcopy(self.gbfile['features'])
		for i in range(len(self.gbfile['features'])): #checks to match dna change
			if len(self.gbfile['features'][i]['location']) == 1:
				featurestart, featurefinish = self.get_location(self.gbfile['features'][i]['location'][0])
			else:
				n = 0
				featurestart = self.get_location(self.gbfile['features'][i]['location'][n])[0]
				n = len(self.gbfile['features'][i]['location'])-1
				featurefinish = self.get_location(self.gbfile['features'][i]['location'][n])[1]
			
			if (start<=featurestart and featurefinish<=finish) == True: #if change encompasses whole feature
				self.clipboard['features'].append(self.allgbfeatures_templist[i])
				for n in range(len(self.gbfile['features'][i]['location'])):
					newlocation = self.add_or_subtract_to_locations(self.gbfile['features'][i]['location'][n], -start+1, 'b')
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
		feature = copy.deepcopy(feature)
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
		feature = {}

		if key in self.featuretypes: 
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
		self.add_file_version()
		
	def remove_feature(self, feature):
		"""Function removes the feature that is passed to it from the genbank file"""
		position = self.get_feature_index(feature)
		
		if position is False:
			print('feature identify error')
		else:
			del self.gbfile['features'][position]
		self.add_file_version()

	def move_feature(self, feature, upordown):
		'''Moves a feature one step up or down the list (up defined as closer to the beginning)'''
		index = self.get_feature_index(feature)
		if upordown == 'u' and index != 0:
			self.gbfile['features'][index-1], self.gbfile['features'][index] = self.gbfile['features'][index], self.gbfile['features'][index-1]

		elif upordown == 'd' and index != len(self.gbfile['features'])-1:
			self.gbfile['features'][index+1], self.gbfile['features'][index] = self.gbfile['features'][index], self.gbfile['features'][index+1]
		self.add_file_version()

	def add_qualifier(self, feature, newqualifier):
		'''Adds qualifier tag to existing feature'''
		index = self.get_feature_index(feature)
		if index is False:
			print('Error, no index found')
		else:
			self.gbfile['features'][index]['qualifiers'].append(newqualifier)  #change append to sth that works for dicts
		self.add_file_version()		
		
	def remove_qualifier(self, feature, number):
		'''Removes a qualifier tag from an existing feature'''
		index = self.get_feature_index(feature)
		if index is False:
			print('Error, no index found')
		else:
			del self.gbfile['features'][index]['qualifiers'][number]
		self.add_file_version()

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
		self.add_file_version()

	def ApEandVNTI_clutter(self):
		'''Find out whether there is clutter from Vector NTI or ApE in the genbank file'''
		for i in range(len(self.gbfile['features'])):
			for n in range(len(self.gbfile['features'][i]['qualifiers'])):
				if 'ApEinfo' in self.gbfile['features'][i]['qualifiers'][n] or 'vntifkey' in self.gbfile['features'][i]['qualifiers'][n] :
					return True
		return False


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
		if len(brackets) == 0: 
			pass					
		elif brackets[0] == '<': 
			start = '<' + str(start)
		elif brackets[0] == '>' or brackets[1] == '>': 
			finish = '>' + str(finish)
		location = str(start) + '..' + str(finish)		
		return location		

	
	def GetFeatureDNA(self, index):
		'''Retrieves the dna sequence for a feature specified by its index'''
		assert type(index) == int, "Index must be an integer."
		assert (index>=0 and index<len(self.gbfile['features'])), "The index must be between 0 and %s" % str(len(self.gbfile['features'])-1)
		DNA = ''
		for entry in self.gbfile['features'][index]['location']:
			start, finish = self.get_location(entry)
			DNA += self.GetDNA(start, finish)
		
		if self.gbfile['features'][index]['complement'] == True:
			DNA = dna.RC(DNA)
		return DNA


	def GetDNA(self, start=1, finish=-1):
		'''Get the entire DNA sequence from the self.gbfile'''
		if (start == 1 and finish == -1) == True:
			return self.gbfile['dna']
		else:
			if finish == -1:
				finish = len(self.GetDNA())
			assert (type(start) == int and type(finish) == int), 'Function requires two integers.'
			assert start <= finish, 'Starting point must be before finish.'
			assert start > 0 and start <= len(self.gbfile['dna']), 'Starting point must be between 1 and the the DNA length.'
			assert finish > 0 and finish <= len(self.gbfile['dna']), 'Finish must be between 1 and the the DNA length.'
			return self.gbfile['dna'][start-1:finish]
	
	def GetFilepath(self):
		'''Get the self.filepath for the opened file'''
		return self.gbfile['filepath']
	
	def SetFilepath(self, new_path):
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
		featurelist = []		
		for i in range(len(self.gbfile['features'])):
			feature = self.gbfile['features'][i]
			if feature['complement'] == True:
				complement = 'complement'
			else:
				complement = 'leading'
			currentfeature = '>[%s] %s at %s on %s strand' % (str(i), feature['qualifiers'][0], feature['location'], complement)
			print(currentfeature)
#			print('\n')
			featurelist.append(currentfeature)
		return featurelist


################## Find methods ###################################

	def FindNucleotide(self, searchstring, searchframe=-1, searchRC=False):
		'''Method for finding a certain DNA sequence. Degenerate codons are supported.
Searchstring should be a integer number or string of DNA characers.
searchRC is True or False and determines whether the reverse complement should also be searched. 
Searchframe should be an index for a feature. -1 indicates search in entire genbank file
indeces >-1 are feature indeces'''


### Need to implement the searchRC variable

		#fix the set_dna_selection functions here
		assert type(searchstring) == str or type(searchstring) == unicode or type(searchstring) == int, 'Error, search takes a string of DNA or a string of numbers as input.'
		assert type(searchframe) == int and -1<=searchframe<len(self.get_all_features()), 'Error, %s is not a valid argument for searchframe.' % str(searchframe)
		assert type(searchRC) == bool, 'Error, searchRC must be True or False'
		#empty search string
		if searchstring=='':
			print 'The searchstring is missing, please check your input'
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


	def FindAminoAcid(self, searchstring, searchframe, searchRC=False):
		'''Method for finding a certain protein sequence, or position, in the file. Degenerate codons are supported'''
		assert type(searchstring) == str or type(searchstring) == unicode or type(searchstring) == int, 'Error, search takes a string of DNA or a string of numbers as input.'
		assert type(searchframe) == int and -1<=searchframe<len(self.get_all_features()), 'Error, %s is not a valid argument for searchframe.' % str(searchframe)
		assert type(searchRC) == bool, 'Error, searchRC must be True or False'


		#empty search string
		if searchstring=='':
			print 'The searchstring is missing, please check your input'
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


	def FindFeature(self, searchstring):
		'''Method for finding a certain feature in a genbank file'''
		assert type(searchstring) == str or type(searchstring) == unicode or type(searchstring) == int, 'Error, search takes a string of DNA or a string of numbers as input.'

		if searchstring=='':
			print 'The searchstring is missing, please check your input'
			return []
		elif type(searchstring) is int:
			#assert that the integer is within the number of features
			#return the position for that feature
			raise NotImplementedError 
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
			return []
		else:
			for feature in hits:
				start, finish = self.GetFirstLastLocation(feature)
				search_hits.append((start-1, finish))
		search_hits = sorted(search_hits)
		return search_hits	



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
			assert type(mutation) is str or type(mutation) is unicode, 'Error, input must be a string or unicode.'
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
					print('Mutation %s%s%s(%s) performed.' % (positionAA, position, trailingAA, codon))
				elif complement is True:
					self.changegbsequence(global_position[0][0], global_position[0][1], 'r', dna.RC(codon))
					if silent is False: self.add_feature(key='modified_base', qualifiers=['/note="%s%s%s(%s)"' % (positionAA, position, trailingAA, codon)], location=['%s..%s' % (global_position[0][0], global_position[0][1])], complement=True, join=False, order=False)
					print('Mutation %s%s%s(%s) performed.' % (positionAA, position, trailingAA, codon))
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
					positionnucleotide = dna.RC(self.GetDNA(global_position[0][0], global_position[0][1])).upper()
				elif complement is False:
					positionnucleotide = self.GetDNA(global_position[0][0], global_position[0][1]).upper()
				if leadingnucleotide != '':
					assert positionnucleotide == leadingnucleotide, 'Error, position %s has the nucleotide %s, and not %s as specified.' % (str(position),  positionnucleotide, leadingnucleotide)

				#now make the mutation and add corresponding feature if the 'silent' variable is False
				if complement is False or mutationframe == -1:
					self.changegbsequence(global_position[0][0], global_position[0][1], 'r', trailingnucleotide)
					if silent is False: self.add_feature(key='modified_base', qualifiers=['/note=%s%s%s' % (positionnucleotide, position, trailingnucleotide)], location=['%s..%s' % (global_position[0][0], global_position[0][1])], complement=False, join=False, order=False)
					print('Mutation %s%s%s performed.' % (positionnucleotide, position, trailingnucleotide))
				elif complement is True:
					self.changegbsequence(global_position[0][0], global_position[0][1], 'r', dna.RC(trailingnucleotide))
					if silent is False: self.add_feature(key='modified_base', qualifiers=['/note=%s%s%s' % (positionnucleotide, position, trailingnucleotide)], location=['%s..%s' % (global_position[0][0], global_position[0][1])], complement=True, join=False, order=False)
					print('Mutation %s%s%s performed.' % (positionnucleotide, position, trailingnucleotide))
				else:
					raise ValueError


#################################
#################################



	def get_location(self, location):
		'''Takes a location entry and extracts the start and end numbers'''	
		#This needs a serious update to deal with more exotic arrangements
		templocation = ''
		for n in range(0,len(location)): #for each character
			if location[n] != '<' and location[n] != '>': 
				templocation += location[n]
		start, finish = templocation.split('..')
		return int(start), int(finish)


	def get_all_feature_positions(self):
		'''Get type, complement, start and finish for each feature'''		
		positionlist = []
		
		for i in range(0,len(self.gbfile['features'])):
			Key = self.gbfile['features'][i]['key']
			Complement = self.gbfile['features'][i]['complement']
			name = self.gbfile['features'][i]['qualifiers'][0].split('=')[1]
			for entry in self.gbfile['features'][i]['location']:
				start, finish = self.get_location(entry)
				start = int(start)
				finish = int(finish)
				positionlist.append([Key, Complement, start, finish, name, i])
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
		if changetype == 'r': #replacement. This method does NOT modify feture positions. Use with caution.
			self.gbfile['dna'] = self.gbfile['dna'][:changestart-1] + change + self.gbfile['dna'][changeend:]
					
		elif changetype == 'i': #insertion
			olddnalength = len(self.gbfile['dna']) #for changing header
			self.gbfile['dna'] = self.gbfile['dna'][:changestart-1] + change + self.gbfile['dna'][changestart-1:]
			self.gbfile['header'] = self.gbfile['header'].replace('%s bp' % olddnalength, '%s bp' % len(self.gbfile['dna'])) #changing header
			for i in range(len(self.gbfile['features'])): #change features already present
				for n in range(len(self.gbfile['features'][i]['location'])):
					start, finish = self.get_location(self.gbfile['features'][i]['location'][n])
					if start<changestart<=finish: #if change is within the feature
						self.gbfile['features'][i]['location'][n] = self.add_or_subtract_to_locations(self.gbfile['features'][i]['location'][n], len(change), 'f')
					elif changestart<=start: #if change is before feature
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
						#print('within feature')
					elif changestart<start and start<=changeend<finish: #if change encompasses start, change start and finish
						self.gbfile['features'][i]['location'][n] = self.add_or_subtract_to_locations(self.gbfile['features'][i]['location'][n], changeend+1-start, 's')	#start					
						self.gbfile['features'][i]['location'][n] = self.add_or_subtract_to_locations(self.gbfile['features'][i]['location'][n], -len(change), 'b') #both
						#print('encompass start')
					elif start<changestart<=finish and finish<changeend: #if change encompasses finish, change finish
						self.gbfile['features'][i]['location'][n] = self.add_or_subtract_to_locations(self.gbfile['features'][i]['location'][n], -(finish-changestart)-1, 'f')	 #finish
						#print('encompass finish')
					elif changestart<=start and finish<=changeend: #if change encompasses whole feature, add to deletion list
						deletionlist.append(copy.deepcopy((i, n)))
						#print('encompass all')
					elif changestart<start and changeend<start: #if change is before feature, change start and finish
						self.gbfile['features'][i]['location'][n] = self.add_or_subtract_to_locations(self.gbfile['features'][i]['location'][n], -len(change), 'b')	 #both			
						#print('before start')
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

	def Save(self, filepath=None):
		"""Function writes data stored in header, featruelist and dna to a .gb file"""
		if filepath==None:
			filepath = self.GetFilepath()

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




