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




import dna
from copy import deepcopy
import pyperclip



class gbobject():
	"""Class that reads a genbank file (.gb) and has functions to edit its features and DNA sequence"""
	def __init__(self):
		self.gbfile = {} #this variable stores the whole genbank file
		self.allgbfeatures = [] #will contain all features in list format
		self.filepath = '' #for keeping track of the file being edited
		
		
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
#		print(featurelist)
		templist = []
		tempstr = ''
		
		if featurelist[-1] == '': del featurelist[-1] #last entry tends to be empty, if so, remove
		
		for line in range(1, len(featurelist)):
			if ('..' in featurelist[line] and tempstr != '') == True:
				if tempstr[-1] == ',': #to protect against numberings that extend over several rows
					tempstr += featurelist[line]

				elif line+1 == len(featurelist):
					tempstr += featurelist[line]


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
			
			
				featurelist2.append(templist)
				tempstr = featurelist[line]
				
			elif '..' in featurelist[line] and tempstr == '': #first feature
				tempstr = featurelist[line]

			else:						#in-between features
				tempstr += featurelist[line]
		


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
		
	def get_location(self, entry):
		'''Returns start and end location for an entry of a location list'''
		tempentry = ''
		for n in range(len(entry)):
			if entry[n] != '<' and entry[n] != '>': tempentry += entry[n]
		start, finish = tempentry.split('..')
		return int(start), int(finish)


			


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


	def add_empty_feature(self):
		"""Function adds a feature to self.allgbfeatures"""
		feature = {}
		feature['key'] = ""
		feature['qualifiers'] = ['/label=empty']
		feature['location'] = []
		feature['complement'] = False
		self.gbfile['features'].append(feature) #change append to sth that works for dicts

	
	def identify_feature(self, feature):
		'''Used to find the position of a feature in the self.gbfile data structure'''
		if len(self.gbfile['features']) == 0:
			return False
		else:
			for i in range(len(self.gbfile['features'])):
				if self.gbfile['features'][i]['key'] != feature['key']: continue
				if self.gbfile['features'][i]['location'] != feature['location']: continue
				if self.gbfile['features'][i]['qualifiers'][0] != feature['qualifiers'][0]: continue
				if self.gbfile['features'][i]['complement'] == feature['complement']: #it is == here
					return i
			return False

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
			
		
	def remove_feature(self, feature):
		"""Function removes a feature from self.allgbfeatures based on its ID (which is the first qualifier)"""
		position = self.identify_feature(feature)
		
		if position is False:
			print('feature identify error')
		else:
			del self.gbfile['features'][position]


	def add_qualifier(self, feature, newqualifier):
		'''Adds qualifier tag to existing feature'''
		index = self.identify_feature(feature)
		if index is False:
			print('Error, no index found')
		else:
			self.gbfile['features'][index]['qualifiers'].append(newqualifier)  #change append to sth that works for dicts
		
		
	def remove_qualifier(self, feature, number):
		'''Removes a qualifier tag from an existing feature'''
		index = self.identify_feature(feature)
		if index is False:
			print('Error, no index found')
		else:
			del self.gbfile['features'][index]['qualifiers'][number]


	def change_feature_type(self, feature, newkey):
		'''Changes feature type of the feature passed to method'''
		index = self.identify_feature(feature)
		if index is False:
			print('Error, no index found')
		else:
			self.gbfile['features'][index]['key'] = newkey

	def change_feature_complement(self, feature, complement):
		'''Changes whether a feature is on leading or complement DNA strand'''
		index = self.identify_feature(feature)
		if index is False:
			print('Error, no index found')
		else:
			self.gbfile['features'][index]['complement'] = complement 

	def changegbfeatureid(self, oldfeatureid, newfeatureid):
		"""Function changes the ID for a certain feature in self.allgbfeatures"""
		for i in range(len(self.gbfile['features'])):
			if self.gbfile['features'][i]['qualifiers'][0] == '/label='+oldfeatureid:
				self.gbfile['features'][i]['qualifiers'][0] = '/label='+newfeatureid


	def copy(self, copystart, copyend):
		'''Copy DNA and all the features for a certain selection'''
		self.clipboard = {}
		self.clipboard['dna'] = self.get_dna()[copystart:copyend] #copy to internal clipboard
		self.clipboard['features'] = []
		self.allgbfeatures_templist = deepcopy(self.gbfile['features'])

		for i in range(len(self.gbfile['features'])): #checks to match dna change
	
			if len(self.gbfile['features'][i]['location']) == 1:
				start, finish = self.get_location(self.gbfile['features'][i]['location'])
			else:
				n = 0
				start = self.get_location(self.gbfile['features'][i]['location'][n])[0]
				n = len(self.gbfile['features'][i]['location'])-1
				finish = self.get_location(self.gbfile['features'][i]['location'][n])[1]
				
			if copystart<start<=finish<=copyend: #if change encompasses whole feature
				self.clipboard['features'].append(self.allgbfeatures_templist[i])
				for n in range(len(self.gbfile['features'][i]['location'])):
					newlocation = self.add_or_subtract_to_locations(self.gbfile['features'][i]['location'][n], -copystart, 'b')
					self.clipboard['features'][-1]['location'][n] = newlocation
	

	def paste(self, pastestart):
		'''Makes an insertion and updates dna and features'''
		temp_clipboard = deepcopy(self.clipboard) #creates a deep copy which is needed to copy nested lists
		
		if temp_clipboard['dna'] != pyperclip.paste(): #if internal and system clipboard is not same then system clipboard takes presidence
			print('internal clipboard override')
			temp_clipboard['dna'] = pyperclip.paste()

		DNA = str(temp_clipboard['dna'])
		self.changegbsequence(pastestart+1, pastestart+1, 'i', DNA) #change dna sequence	

		for i in range(len(temp_clipboard['features'])): #add features from clipboard
			self.paste_feature(temp_clipboard['features'][i], pastestart)

#		if len(self.gb.self.allgbfeatures) > 1: #make sure that all features have unique names
#			self.check_for_unique_feature_names(self.gb.self.allgbfeatures)	


	def remove_location(self, index, number):
		'''Removes locaiton in self.gbfile['features'][index]['location'][number]'''
		del self.gbfile['features'][index]['location'][number]
		if len(self.gbfile['features'][index]['location']) == 0: # if no locations are left for that feature, delete feature
			del self.gbfile['features'][index] 


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

	
	def get_feature_dna(self, featureid):
		'''Retrieves the dna sequence for a specified feature'''
		DNA = ''
		for i in range(len(self.gbfile['features'])):
			if (('/label=' in featureid) == False and self.gbfile['features'][i]['qualifiers'][0] == '/label=' + featureid) or self.gbfile['features'][i]['qualifiers'][0] == featureid:
				for entry in self.gbfile['features'][i]['location']:
					start, finish = self.get_location(entry)
					
					if self.gbfile['features'][i]['join'] == False and self.gbfile['features'][i]['order'] == False:
						DNA = self.gbfile['dna'][start+1:finish+1]
					
					elif self.gbfile['features'][i]['join'] == True:
						DNA += self.gbfile['dna'][start+1:finish+1]
						
					elif self.gbfile['features'][i]['order'] == True: #I probably need to change the reversecomplement function to not remove /n
						DNA += self.gbfile['dna'][start+1:finish+1] + '\n'
				
				if self.gbfile['features'][i]['complement'] == True:
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

	def get_feature_for_pos(self, position):
		'''For a given dna position, which features are there?'''		
		Feature = ''
		for feature in self.gbfile['features']:
			for entry in feature['location']: #there may be several starts and finishes in a feature...
				start, finish = self.get_location(entry)
				start -= 1
				finish -= 1

				if start <= position <= finish:
					if Feature != '': Feature += ', '
					Feature += feature['qualifiers'][0].split('=')[1] #removes the '/label=' part from qualifier

		if Feature == '': Feature = 'None'
		return Feature

	def list_features(self):
		'''List all features as a string output'''
		featurelist = ''		
		for i in range(len(self.gbfile['features'])):
			feature = self.gbfile['features'][i]
			if feature['complement'] == True:
				complement = 'complement'
			else:
				complement = 'leading'
			featurelist += ('>%s at %s on %s strand\n' % (feature['qualifiers'][0], feature['location'], complement))
		return featurelist

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
			#print(entry)
			
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




gb = gbobject()





