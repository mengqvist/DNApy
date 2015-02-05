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
#fix the line shortening in the header section
#fix file saving
#add single base support to location parsing
#add methods to modify header sections

#match features with mandatory and optional qualifiers
#add a function that checks that everything is ok. I.e. that it conforms to the genbank format.
#make changes to how qualifiers are parsed. For example /qualifier=xyz, the '=' is not always there...

import dna
import copy
import pyperclip
import oligo_localizer
import peptide_localizer
import re
import sys
import wx	# for richCopy
import json	# for richCopy


#the feature class is not currently used. 
#Something I started working on and may or may not continue.
class feature(object):
	"""
	A feature object class that defines the key, location and qualifiers of a feature and methods to interact with these.
	Data structure is as follows:
	{key:string #feature key
		location:list #list of locations on DNA for feature
		qualifiers:list #list of qualifiers attached to the fature
		complement:bool #is feature on complement strand or not
		join:bool #should feature locations be joined or not
		order:bool #are feature locations in a certain order or not
		}
	""" 

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
		self.clipboard = {}
		self.clipboard['dna'] = ''
		self.clipboard['features'] = []
		self.featuretypes = ["modified_base", "variation", "enhancer", "promoter", "-35_signal", "-10_signal", "CAAT_signal", "TATA_signal", "RBS", "5'UTR", "CDS", "gene", "exon", "intron", "3'UTR", "terminator", "polyA_site", "rep_origin", "primer_bind", "protein_bind", "misc_binding", "mRNA", "prim_transcript", "precursor_RNA", "5'clip", "3'clip", "polyA_signal", "GC_signal", "attenuator", "misc_signal", "sig_peptide", "transit_peptide", "mat_peptide", "STS", "unsure", "conflict", "misc_difference", "old_sequence", "LTR", "repeat_region", "repeat_unit", "satellite", "mRNA", "rRNA", "tRNA", "scRNA", "snRNA", "snoRNA", "misc_RNA", "source", "misc_feature", "misc_binding", "misc_recomb", "misc_structure", "iDNA", "stem_loop", "D-loop", "C_region", "D_segment", "J_segment", "N_region", "S_region", "V_region", "V_segment"]


		self.search_hits = []	# variable for storing a list of search hits
		
		#set up data structure
		self.gbfile = {} #this variable stores the whole genbank file
		self.gbfile['locus'] = {}
		self.gbfile['locus']['name'] = None
		self.gbfile['locus']['length'] = None
		self.gbfile['locus']['type'] = None
		self.gbfile['locus']['topology'] = None
		self.gbfile['locus']['division'] = None
		self.gbfile['locus']['date'] = None		
		self.gbfile['definition'] = None
		self.gbfile['accession'] = None
		self.gbfile['version'] = None
		self.gbfile['gi'] = None
		self.gbfile['dblink'] = None
		self.gbfile['keywords'] = None
		self.gbfile['segment'] = None		
		self.gbfile['source'] = None
		self.gbfile['organism'] = None
		self.gbfile['references'] = None
		self.gbfile['comments'] = None
		self.gbfile['dbsource'] = None	#not sure this is a valid keyword
		self.gbfile['primary'] = None	#not sure this is a valid keyword	
		self.gbfile['features'] = None
		self.gbfile['origin'] = None
		self.gbfile['contig'] = None		
		self.gbfile['dna'] = None

		self.gbfile['filepath'] = filepath

		self.fileName = 'New DNA' #name of the file/plasmid name

		
		## compile regular expressions used for parsing ##
		self._re_locus = re.compile(r'''
				LOCUS \s+? ([a-zA-Z0-9:=_-]+?) \s+?												#match name
				([0-9]+?)[ ](?:bp|aa) \s+														#match length
				((?:ss-|ds-|ms-)*(?:NA|DNA|RNA|tRNA|rRNA|mRNA|uRNA))* \s+						#match type, zero or one time
				(linear|circular)* \s+															#match topology, zero or one time
				(PRI|ROD|MAM|VRT|INV|PLN|BCT|VRL|PHG|SYN|UNA|EST|PAT|STS|GSS|HTG|HTC|ENV)* \s+  #match division, zero or one time
				([0-9]{2}-(?:JAN|FEB|MAR|APR|MAY|JUN|JUL|AUG|SEP|OCT|NOV|DEC){1}-[0-9]{4})*		#match date, zero or one time
				''', re.VERBOSE)

		self._re_version = re.compile(r'''
					VERSION \s+ ([-_.a-zA-Z0-9]+)* \s+ 		#match version, zero or one time
					(?:GI:([0-9]+))*							#match gi, zero or one time
					''', re.VERBOSE)

		#for undo/redo
		self.file_versions = [] #stores version of the file
		self.file_version_index = 0 #stores at which index the current file is					
		
		if filepath == None:
			pass
		else:
			platform = sys.platform
			if platform == 'win32':
				filepath = filepath.replace('\\', '/')
			self.opengb(filepath)

		
					
					
###############################

	def opengb(self, filepath):
		'''
		Function opens a genbank file.
		'''
		assert type(filepath) == str or type(filepath) == unicode , "Error opening genbank file. Filepath is not a string: %s" % str(filepath)
		self.fileName = filepath.split('/')[-1] #update the fileName variable
		
		#open it and read line by line (very memory efficient)
		with open(filepath) as infile:
			self.readgb(infile)

			
	def readgb(self, infile):
		"""
		Method takes a genbank file and parses it into a manageable data structure.
		"""

		line_list = [] #for collecting the dna
		for line in infile:
			line = line.replace('\r', '') #important for removing linux \r newline character

			if re.match('^[ \t\n]+$', line): #get rid of blank lines
				pass

			elif line[0] == ' ' and line[0:10] != '        1 ': #a line that starts with a space is a continuation of the previous line (and belongs to the previous keyword). I have to have a special case '        1 ' for where the DNA sequence starts.
				line_list.append(line)

			elif (line[0] != ' ' or line[0:10] == '        1 ') and len(line_list) == 0: #a line that does not start with a space marks a new keyword
				line_list.append(line)
				
			elif (line[0] != ' ' or line[0:10] == '        1 ') and len(line_list) > 0: #a line that does not start with a space marks a new keyword
				self.parse(''.join(line_list)) 
				line_list = []
				line_list.append(line)
			
		#add the loaded file state to the undo/redo list
#		self.add_file_version()

		
	def parse(self, line):
		'''
		Method for parsing a genbank line. The term line is used very flexibly and indicates all the lines belonging to a certain keyword.
		It matches certain keywords to the beginning of the line and parses accordingly.
		'''
		if 'LOCUS' in line[0:12]:
			#LOCUS A short mnemonic name for the entry, chosen to suggest the sequence's definition. Mandatory keyword/exactly one record.
			
			#match line with re
			m = re.match(self._re_locus, line)
			
			#assign values
			self.gbfile['locus']['name'] = m.group(1) #locus id
			self.gbfile['locus']['length'] = m.group(2) #sequence length
			self.gbfile['locus']['type'] = m.group(3) #type of sequence
			self.gbfile['locus']['topology'] = m.group(4) #linear or circular
			self.gbfile['locus']['division'] = m.group(5) #which GenBank division the sequence belongs to
			self.gbfile['locus']['date'] = m.group(6) #modification date
			
		elif 'DEFINITION' in line[0:12]:
			#A concise description of the sequence. Mandatory keyword/one or more records.
			line_list = line[12:].split('\n') #split by line
			self.gbfile['definition'] = ''.join([re.sub(' +', ' ', x) for x in line_list]) #remove spaces longer than one and empty entries
			
		elif 'ACCESSION' in line[0:12]:
			#The primary accession number is a unique, unchanging identifier assigned to each GenBank sequence record. (Please use this identifier when citing information from GenBank.) Mandatory keyword/one or more records.
			self.gbfile['accession'] = line[12:].strip('\t\n\x0b\x0c\r ')
			
		elif 'VERSION' in line[0:12]:
			#A compound identifier consisting of the primary accession number and a numeric version number associated with the current version of the sequence data in the record. This is optionally followed by an integer identifier (a "GI") assigned to the sequence by NCBI. Mandatory keyword/exactly one record.
			m = re.match(self._re_version, line)
			
			self.gbfile['version'] = m.group(1)
			self.gbfile['gi'] = m.group(2)
		
		elif 'NID' in line[0:3]:
			#An alternative method of presenting the NCBI GI identifier (described above).
			#NOTE: The NID linetype is obsolete and was removed from the GenBank flatfile format in December 1999.
			if self.gbfile['gi'] == None: #if the gi has not yet been assigned
				self.gbfile['gi'] = ''.join(line[0:12].split('\n')).strip()
			else:
				pass
			
		elif 'PROJECT' in line[0:12]:
			#The identifier of a project (such as a Genome Sequencing Project) to which a GenBank sequence record belongs. Optional keyword/one or more records.
			#NOTE: The PROJECT linetype is obsolete and was removed from the GenBank flatfile format after Release 171.0 in April 2009.
			self.gbfile['project'] = ''.join(line[12:].split('\n')).strip()
			
		elif 'DBLINK' in  line[0:12]:
			#Provides cross-references to resources that support the existence a sequence record, such as the Project Database and the NCBI Trace Assembly Archive. Optional keyword/one or more records.
			line_list = line[12:].split('\n') #split by line only
			self.gbfile['dblink'] = ' '.join([s.strip() for s in line[12:].split('\n')])

		elif 'KEYWORDS' in line[0:12]:
			#Short phrases describing gene products and other information about an entry. Mandatory keyword in all annotated entries/one or more records.
			line_list = re.split('[\n ]', line[12:]) #split by line and space
			self.gbfile['keywords'] = ' '.join([re.sub(' +', ' ', x) for x in line_list if x != '']) #remove spaces longer than one, and also empty entries

		elif 'SEGMENT' in line[0:12]:
			#Information on the order in which this entry appears in a series of discontinuous sequences from the same molecule. Optional keyword (only in segmented entries)/exactly one record.
			self.gbfile['segment'] = line[12:].strip()
			
		elif 'SOURCE' in line[0:6]:
			#Common name of the organism or the name most frequently used in the literature. Mandatory keyword in all annotated entries/one or more records/includes one subkeyword.
			#contains the sub-class ORGANISM, Formal scientific name of the organism (first line) and taxonomic classification levels (second and subsequent lines). Mandatory subkeyword in all annotated entries/two or more records.
			
			line_list = line[12:].split('\n') #split by line break
			index = next(line_list.index(x) for x in line_list if 'ORGANISM' in x)  #find the index of 'ORGANISM'
			
			#the free-format organism information (which may stretch over several lines)
			self.gbfile['source'] = ' '.join([re.sub('[ .\n]+', ' ', x) for x in line_list[0:index]]) #the first entry (the source line)
	
			#The formal scientific name for the source organism
			self.gbfile['organism'] = []
			self.gbfile['organism'].append(line_list[index][12:]) #this is the 'organism' line
			
			#here comes the taxonomy
			taxon = ''.join(line_list[index+1:]).split(';') #join the rest of the entries and split on the ; delimiter
			self.gbfile['organism'].extend([re.sub('[ .]+', '', x.replace('[','').replace(']','')) for x in taxon]) #now remove spaces and punctuations and add the taxonomy to the global data structure

		elif 'REFERENCE' in line[0:12]:
			#Citations for all articles containing data reported in this entry. Includes seven subkeywords and may repeat. Mandatory keyword/one or more records.
			#AUTHORS - Lists the authors of the citation. Optional subkeyword/one or more records.
			#CONSRTM - Lists the collective names of consortiums associated with the citation (eg, International Human Genome Sequencing Consortium), rather than individual author names. Optional subkeyword/one or more records.
			#TITLE - Full title of citation. Optional subkeyword (present in all but unpublished citations)/one or more records.
			#JOURNAL - Lists the journal name, volume, year, and page numbers of the citation. Mandatory subkeyword/one or more records.
			#MEDLINE - Provides the Medline unique identifier for a citation. Optional subkeyword/one record.
			#NOTE: The MEDLINE linetype is obsolete and was removed
			#PUBMED - Provides the PubMed unique identifier for a citation. Optional subkeyword/one record.
			#REMARK	- Specifies the relevance of a citation to an entry. Optional subkeyword/one or more records. from the GenBank flatfile format in April 2005.
			if self.gbfile['references'] == None:
				self.gbfile['references'] = []
			current_ref = {}
			current_ref['reference'] = None
			current_ref['authors'] = None
			current_ref['consrtm'] = None
			current_ref['title'] = None
			current_ref['journal'] = None
			current_ref['medline'] = None
			current_ref['pubmed'] = None
			current_ref['remark'] = None
			
			ref_lines = re.split('\n(?!            )', line) #match newline only if it is not matched by 12 whitespace characters
			
			for l in ref_lines:
				if 'REFERENCE' in l[0:12]:
					current_ref['reference'] = ' '.join([s.strip() for s in l[12:].split('\n')])			
				elif 'AUTHORS' in l[0:12]:
					current_ref['authors'] = ' '.join([s.strip() for s in l[12:].split('\n')])
				elif 'CONSRTM' in l[0:12]:
					current_ref['consrtm'] = ' '.join([s.strip() for s in l[12:].split('\n')])
				elif 'TITLE' in l[0:12]:
					current_ref['title'] = ' '.join([s.strip() for s in l[12:].split('\n')])
				elif 'JOURNAL' in l[0:12]:
					current_ref['journal'] = ' '.join([s.strip() for s in l[12:].split('\n')])
				elif 'MEDLINE' in l[0:12]:
					current_ref['medline'] = ' '.join([s.strip() for s in l[12:].split('\n')])
				elif 'NOTE:' in l[0:12]:
					pass
				elif 'PUBMED' in l[0:12]:
					current_ref['pubmed'] = ' '.join([s.strip() for s in l[12:].split('\n')])
				elif 'REMARK' in l[0:12]:
					current_ref['remark'] = ' '.join([s.strip() for s in l[12:].split('\n')])
			
			self.gbfile['references'].append(current_ref)
			

		elif 'COMMENT' in line[0:12]:
			#Cross-references to other sequence entries, comparisons to other collections, notes of changes in LOCUS names, and other remarks. Optional keyword/one or more records/may include blank records.
			if self.gbfile['comments'] == None:
				self.gbfile['comments'] = []
			line_list = line[12:].split('\n') #split by line
			self.gbfile['comments'].append(''.join([re.sub(' +', ' ', x) for x in line_list if x != '']))  #remove spaces longer than one and empty entries
			
		elif 'DBSOURCE' in  line[0:12]: #should I really include this one?
			self.gbfile['dbsource'] = ''.join(line[12:].split('\n')).strip()
			
		elif 'PRIMARY' in line[0:12]: #should I really include this one?
			self.gbfile['primary'] = ''.join(line[12:].split('\n')).strip()
			#see gi_625194262
		
		elif 'FEATURES' in line[0:12]:
			#Table containing information on portions of the sequence that code for proteins and RNA molecules and information on experimentally determined sites of biological significance. Optional keyword/one or more records.
			feature_lines = re.split('\n(?!                     )', line) #match newline only if it is not matched by 21 whitespace characters
			if self.gbfile['features'] == None:
				self.gbfile['features'] = []
			
			for l in feature_lines:
				if l != '':
					self.parse_feature_line(l)
			
		elif 'BASE COUNT' in line[0:12]:
			#Summary of the number of occurrences of each basepair code (a, c, t, g, and other) in the sequence. Optional keyword/exactly one record. 
			#NOTE: The BASE COUNT linetype is obsolete and was removed from the GenBank flatfile format in October 2003.
			pass
			
		elif 'ORIGIN' in line[0:12]:
			#Specification of how the first base of the reported sequence is operationally located within the genome. Where possible, this includes its location within a larger genetic map. Mandatory keyword/exactly one record.
			self.gbfile['origin'] = line[12:].strip()
			
		elif '        1 ' in line[0:10]:
			#the DNA/RNA/protein sequence
			self.gbfile['dna'] = ''.join([re.match('[a-zA-Z]', b).group(0) for b in line if re.match('[a-zA-Z]', b) is not None]) #make the DNA string while skipping anything that is not a-z

		elif 'CONTIG' in line[0:12]:
			#This linetype provides information about how individual sequence
			#records can be combined to form larger-scale biological objects, such as
			#chromosomes or complete genomes. Rather than presenting actual sequence
			#data, a special join() statement on the CONTIG line provides the accession
			#numbers and basepair ranges of the underlying records which comprise the object.
			#It is an alternative to providing a sequence.
			self.gbfile['contig'] = ''.join(line[12:].split('\n')).strip()
		
		else:
			raise ValueError, 'Unparsed line: "%s" in record %s' % (line, self.gbfile['locus']['name'])

		#check the stored features for Vector NTI and ApE clutter and store result
		self.clutter = self.ApEandVNTI_clutter() 
			
			
			
	def parse_feature_line(self, inputline):
		'''
		Takes a feature line consisting info for an entire feature and parses them into the DNApy data format.
		'''
		assert type(inputline) == str, 'Error, the input has to be a string.'
		if 'FEATURES' in inputline[0:12]:
			pass
		else:
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

				if line[0:5] == 5*' ' and (line[5] in 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ35') == True:
					key = line[5:20].rstrip('\t\n\x0b\x0c\r ')
					location = line[21:].rstrip('\t\n\x0b\x0c\r ')
					feature["key"] = key
					feature['location'], feature['complement'], feature['join'], feature['order'] = self.parse_location(location)

				elif line[0:21] == '                     ' and line[21] == '/':
					if '/translation' in line:
						qualifier = re.sub('[ \n]+', '', line[21:]) #remove newline characters and all spaces
					else:
						qualifier = re.sub('[ \n]+', ' ', line[21:]) #remove newline characters and spaces longer than 1
					feature['qualifiers'].append(qualifier)

				else:
					print('error parsing feature line', line)	


			self.gbfile['features'].append(copy.deepcopy(feature)) #add feature to the global data structure
		

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
				if entry[n] in '0123456789<>.':			   
					tempstr += entry[n]
			tempsites.append(tempstr)
		location = tempsites
			
		return location, complement, join, order

		
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
		'''
		Get the current file version
		'''
		return self.file_versions[self.get_file_version_index()]

	def add_file_version(self):
		'''
		Add another file version to the version list.
		This should be added after the current version and should delete all later versions if there are any.
		'''
		index = self.get_file_version_index()
		if len(self.file_versions) == 0: #the first version of the file
			self.file_versions += (copy.deepcopy(self.gbfile),)
		elif index == len(self.file_versions)-1: #if the current version is the last one
			self.file_versions += (copy.deepcopy(self.gbfile),)
			self.set_file_version_index(index+1)
		else:
			self.file_versions = self.file_versions[0:index]+(copy.deepcopy(self.gbfile),)
			self.set_file_version_index(index+1)

	def get_file_version_index(self):
		'''
		Get the index of the current file version
		'''
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
		'''
		Returns ultimate first and ultimate last position of a feature, 
		regardless of how many pieces it is broken into.
		'''
		index = self.get_feature_index(feature)
		if index is False:
			print('Error, no index found')
		else:
			locations = self.gbfile['features'][index]['location']
			start = locations[0].split('..')[0]
			finish = locations[-1].split('..')[1]
#			print(start, finish)
			return int(start), int(finish)

	def remove_location(self, index, number):
		'''Removes locaiton in self.gbfile['features'][index]['location'][number]'''
		del self.gbfile['features'][index]['location'][number]
		if len(self.gbfile['features'][index]['location']) == 0: # if no locations are left for that feature, delete feature
			self.remove_feature(self.gbfile['features'][index])

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


	def RichCopy(self, start, finish):
		'''a method to copy not only dna to clipboard but exchange features and dna 
		   beetween two windows
		   this depreciates the Copy method from this file
		'''
		assert (type(start) == int and type(finish) == int), 'Function requires two integers.'

		dna = self.GetDNA(start, finish)

		dnapyClipboard             = {}
		dnapyClipboard['dna']      = self.GetDNA(start, finish) #copy to internal clipboard
		dnapyClipboard['features'] = []
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
				dnapyClipboard['features'].append(self.allgbfeatures_templist[i])
				for n in range(len(self.gbfile['features'][i]['location'])):
					newlocation = self.add_or_subtract_to_locations(self.gbfile['features'][i]['location'][n], -start+1, 'b')
					dnapyClipboard['features'][-1]['location'][n] = newlocation

		# save clipboard to json string
		dnapyClipboard     = json.dumps(dnapyClipboard)
		
		self.richClipboard = wx.DataObjectComposite()
		
		text               = wx.TextDataObject(dna)

		dnapy              = wx.CustomDataObject("application/DNApy")
		dnapy.SetData(dnapyClipboard)
		
		self.richClipboard.Add(text, True)	# add clear dna text as preffered object
		self.richClipboard.Add(dnapy)		# add rich dna info as json
		

	
		# save to the clipboard
		if wx.TheClipboard.Open():
			wx.TheClipboard.SetData(self.richClipboard)
			wx.TheClipboard.Close()
		
		return True
		
		

	def RichPaste(self, ip):
		'''Paste DNA present in clipboard and any features present on that DNA
			If a string is passed to DNA then this will over-ride anything that is present in the clipboard
			this depreciates the Paste method from this file
			'''
		assert type(ip) == int, 'The insertion point must be an integer.'



		self.textClipboard      = wx.TextDataObject()
		self.dnaClipboard       = wx.CustomDataObject("application/DNApy")

		# retrive the clipboard content
		if wx.TheClipboard.Open():
			containsDNApy  = wx.TheClipboard.GetData(self.dnaClipboard)
			containsText   = wx.TheClipboard.GetData(self.textClipboard)
			wx.TheClipboard.Close()
			
			if containsDNApy:
				pasteRich = True
				# get the json object
				pasteDNA  = json.loads(self.dnaClipboard.GetData())
				
			elif containsText:
				pasteRich = False
				# paste the pure text
				paste     = self.textClipboard.GetDataHere()
			else:
				# else there is nothing to paste
				pasteRich = False
				paste     = ''
		
		# here we paste, mostly the same as Paste() 
		if pasteRich:

			DNA = pasteDNA['dna']
			self.changegbsequence(ip, ip, 'i', DNA) #change dna sequence	
			for i in range(len(pasteDNA['features'])): #add features from clipboard
				self.paste_feature(pasteDNA['features'][i], ip-1)
		else:
			DNA = paste
			self.changegbsequence(ip, ip, 'i', DNA) #change dna sequence
		self.add_file_version()
			

		return True


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
		
		if self.gbfile['features'] == None:
			self.gbfile['features'] = [feature]
		else:
			self.gbfile['features'].append(feature) #change append to sth that works for dicts
		self.add_file_version()
		
		
	def remove_feature(self, feature):
		"""
		Function removes the feature that is passed to it from the genbank file.
		"""
		position = self.get_feature_index(feature)
		
		if position is False:
			print('feature identify error')
		else:
			del self.gbfile['features'][position]
			
		if len(self.gbfile['features']) == 0:
			self.gbfile['features'] = None
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
		if self.gbfile['features'] == None:
			return False
		else:
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
		if self.gbfile['dna'] == None:
			return None
		elif start == 1 and finish == -1:
			return self.gbfile['dna']
		else:
			if finish == -1:
				finish = len(self.GetDNA())
			
			assert (type(start) == int and type(finish) == int), 'Function requires two integers.'
			assert start <= finish, 'Starting point must be before finish.'
			assert start > 0 and start <= len(self.gbfile['dna']), 'Starting point must be between 1 and the the DNA length. %s is not' % start
			assert finish > 0 and finish <= len(self.gbfile['dna']), 'Finish must be between 1 and the the DNA length. %s is not' % finish
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
		try:
			start, finish = templocation.split('..')
		except:
			start = templocation
			finish = templocation
		return int(start), int(finish)


	def get_all_feature_positions(self):
		'''Get type, complement, start and finish for each feature'''		
		if self.gbfile['features'] == None:
			return None
		else:
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
		if changetype == 'r': #replacement. This method does NOT modify feature positions. Use with caution.
			self.gbfile['dna'] = self.gbfile['dna'][:changestart-1] + change + self.gbfile['dna'][changeend:]
					
		elif changetype == 'i': #insertion
			if self.gbfile['dna'] == '' or self.gbfile['dna'] == None: #empty dna
				self.gbfile['dna'] = change
			else:
				olddnalength = len(self.gbfile['dna']) #for changing header
				self.gbfile['dna'] = self.gbfile['dna'][:changestart-1] + change + self.gbfile['dna'][changestart-1:]
				self.gbfile['locus']['length'] = len(self.gbfile['dna']) #changing header
				if self.gbfile['features'] != None: #change features already present
					for i in range(len(self.gbfile['features'])): 
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
			self.gbfile['locus']['length'] = len(self.gbfile['dna']) #changing header
			if self.gbfile['features'] != None: #change features already present
				for i in range(len(self.gbfile['features'])):
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

			
	def check_line(self, line, line_type):
		'''
		Evaluates a line (for the genbank file) and makes sure it is maximum 79 characters.
		If the line is longer it finds a good place to break it and returns a modified line.
		Line type is either 'header', 'locations' or 'feature' and indicates from which part of the genbank file the line came.
		'''
		mod_line = ''
		if line_type == 'header':
			if len(line) <= 79: #add short lines directly
				mod_line = line
			else: #split long lines
				words = line[12:].split()
				part_line = ''
				for word in words:
					if part_line == '':
						part_line += line[:12] + word
					elif len(part_line + ' ' + word) <= 67:
						part_line += ' ' + word
					elif len(part_line + ' ' + word) > 67:
						mod_line += part_line + '\n'
						part_line = ' '*12 + word
					else:
						raise ValueError
				mod_line += part_line #catch the last part

				
		elif line_type == 'locations':
			if len(line) <= 58: #add short lines directly
					mod_line = line
			else: #split long lines
				locations = line.split(',')
				part_line = ''
				for location in locations:
					if part_line == '':
						part_line += location
					elif len(part_line + ',' + location) <= 58:
						part_line += ',' + location
					elif len(part_line + ' ' + location) > 58:
						if mod_line == '':
							mod_line += part_line + ',' + '\n'
						else:
							mod_line += ' '*21 + part_line + ',' + '\n'
						part_line = location
					else:
						raise ValueError
				mod_line += ' '*21 + part_line #catch the last part
			
		elif line_type == 'feature':
			if '/translation=' in line:		#it's a sequence, so just break at 58
				for i in range(0, len(line), 58):
					mod_line += ' '*21 + line[0+i:58+i] + '\n'
				mod_line = mod_line[:-1] #remove the last \n
			else: 	#it's regular text, so match line breaks at spaces
				if len(line) <= 58: #add short lines directly
					mod_line = ' '*21 + line
				else: #split long lines
					words = line.split()
					part_line = ''
					for word in words:
						if part_line == '':
							part_line += word
						elif len(part_line + ' ' + word) <= 58:
							part_line += ' ' + word
						elif len(part_line + ' ' + word) > 58:
							mod_line += ' '*21 + part_line + '\n'
							part_line = word
						else:
							raise ValueError
					mod_line += ' '*21 + part_line #catch the last part
		else:
			raise ValueError
		
		return mod_line
	
	
	def make_gbstring(self):
		'''
		Prepare data stored in data structure into one "genbank format" string. 
		Used for saving to flatfile in .gb format and for displaying gbfile.
		'''
		output = [] #collect the output line by line

		## add the LOCUS line ##
		#Positions  Contents
		#---------  --------
		#01-05      'LOCUS'
		#06-12      spaces
		#13-28      Locus name
		#29-29      space
		#30-40      Length of sequence, right-justified
		#41-41      space
		#42-43      bp
		#44-44      space
		#45-47      spaces, ss- (single-stranded), ds- (double-stranded), or
		#           ms- (mixed-stranded)
		#48-53      NA, DNA, RNA, tRNA (transfer RNA), rRNA (ribosomal RNA), 
		#           mRNA (messenger RNA), uRNA (small nuclear RNA).
		#           Left justified.
		#54-55      space
		#56-63      'linear' followed by two spaces, or 'circular'
		#64-64      space
		#65-67      The division code
		#68-68      space
		#69-79      Date, in the form dd-MMM-yyyy (e.g., 15-MAR-1991)
		
		line = ''
		line += 'LOCUS' + ' '*(12-len('LOCUS'))
		
		name = ''
		if self.gbfile['locus']['name'] != None:
			name = name = self.gbfile['locus']['name']
		line += name + ' '*(16-len(name)) + ' '
		
		length = ''
		if self.gbfile['locus']['length'] != None:
			length = str(self.gbfile['locus']['length'])
		line += ' '*(11-len(length)) + length + ' ' + 'bp' + ' '
		
		type = ''
		start = ''
		if self.gbfile['locus']['type'] != None:
			type = self.gbfile['locus']['type']
			if type[2] == '-':
				start = type[:2]
				type = type[2:]
			else:
				start = '   '
		line += start + type + ' '*(6-len(type)) + '  '
			
		topology = ''
		if self.gbfile['locus']['topology'] != None:
			topology = self.gbfile['locus']['topology']
		line += topology + ' '*(8-len(topology)) + ' '
		
		division = ''
		if self.gbfile['locus']['division'] != None:
			division = self.gbfile['locus']['division']
		line += division + ' '
		
		date = ''
		if self.gbfile['locus']['date'] != None:
			date = self.gbfile['locus']['date']
		line += date + ' '*(11-len(date))
		
		output.append(line)
		
			
		## add definition line ##
		if self.gbfile['definition'] != None:
			line = 'DEFINITION' + ' '*(12-len('DEFINITION')) + self.gbfile['definition']
			line = self.check_line(line, 'header')
			output.append(line)
			
		## add accession line ##
		if self.gbfile['accession'] != None:
			line = 'ACCESSION' + ' '*(12-len('ACCESSION')) + self.gbfile['accession']
			output.append(line)
		
		## add version line ##
		if self.gbfile['version'] != None:
			line = 'VERSION' + ' '*(12-len('VERSION')) + self.gbfile['version']
			if self.gbfile['gi'] != None:
				line += '  GI:' + self.gbfile['gi']
			output.append(line)
		
		## add dblink ##
		if self.gbfile['dblink'] != None:
			line = 'DBLINK' + ' '*(12-len('DBLINK')) + self.gbfile['dblink']
			line = self.check_line(line, 'header')
			output.append(line)
		
		## add keywords line ##
		if self.gbfile['keywords'] != None:
			line = 'KEYWORDS' + ' '*(12-len('KEYWORDS')) + self.gbfile['keywords']
			line = self.check_line(line, 'header')
			output.append(line)
		
		## add segment line ##
		if self.gbfile['segment'] != None:
			line = 'SEGMENT' + ' '*(12-len('SEGMENT')) + self.gbfile['segment']
			line = self.check_line(line, 'header')
			output.append(line)
		
		## add source line ##		
		if self.gbfile['source'] != None:
			line = 'SOURCE' + ' '*(12-len('SOURCE')) + self.gbfile['source']
			output.append(line)
		
		## add organism ##
		if self.gbfile['organism'] != None:
			line = '  ORGANISM' + ' '*(12-len('  ORGANISM')) + self.gbfile['organism'][0]
			output.append(line)
			line = ' '*12
			for entry in self.gbfile['organism'][1:]:
				if entry != self.gbfile['organism'][-1]:
					line += entry + '; '
				else:
					line += entry + '.'
			line = self.check_line(line, 'header')
			output.append(line)
		
		## add references ##
		#The REFERENCE field consists of five parts: the keyword REFERENCE, and
		#the subkeywords AUTHORS, TITLE (optional), JOURNAL, MEDLINE (optional),
		#PUBMED (optional), and REMARK (optional).
		if self.gbfile['references'] != None:
			for reference in self.gbfile['references']:
				line = 'REFERENCE' + ' '*(12-len('REFERENCE'))
				if reference['reference'] != None:
					line += reference['reference']
				line = self.check_line(line, 'header')
				output.append(line)
						
				if reference['authors'] != None:
					line = '  AUTHORS' + ' '*(12-len('  AUTHORS')) + reference['authors']
				line = self.check_line(line, 'header')
				output.append(line)
				
				if reference['consrtm'] != None: #optional
					line = '  CONSRTM' + ' '*(12-len('  CONSRTM')) + reference['consrtm']
					line = self.check_line(line, 'header')
					output.append(line)
				
				if reference['title'] != None: #optional
					line = '  TITLE' + ' '*(12-len('  TITLE')) + reference['title']
					line = self.check_line(line, 'header')
					output.append(line)
					
				line = '  JOURNAL' + ' '*(12-len('  JOURNAL'))
				if reference['journal'] != None:
					line += reference['journal']
					line = self.check_line(line, 'header')
					output.append(line)
					
				if reference['medline'] != None: #optional
					line = '   MEDLINE' + ' '*(12-len('   MEDLINE')) + reference['medline']
					line = self.check_line(line, 'header')
					output.append(line)
					
				if reference['pubmed'] != None: #optional
					line = '   PUBMED' + ' '*(12-len('   PUBMED')) + reference['pubmed']
					line = self.check_line(line, 'header')
					output.append(line)
					
				if reference['remark'] != None: #optional
					line = '  REMARK' + ' '*(12-len('  REMARK')) + reference['remark']
					line = self.check_line(line, 'header')
					output.append(line)
		
		## add comments ##
		if self.gbfile['comments'] != None:
			for comment in self.gbfile['comments']:
				line = 'COMMENT' + ' '*(12-len('COMMENT')) + comment
				line = self.check_line(line, 'header')
				output.append(line)
		
		## add dbsource ##
		if self.gbfile['dbsource'] != None:
			line = 'DBSOURCE' + ' '*(12-len('DBSOURCE')) + self.gbfile['dbsource']
			line = self.check_line(line, 'header')
			output.append(line)
		
		## add primary ##
		if self.gbfile['primary'] != None:
			line = 'PRIMARY' + ' '*(12-len('PRIMARY')) + self.gbfile['primary']
			line = self.check_line(line, 'header')
			output.append(line)	
		
		## add features ##
		output.append('FEATURES             Location/Qualifiers')
		if self.gbfile['features'] != None:
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

				line = self.check_line(locations, 'locations')
				output.append(('     %s%s%s') % (entry['key'], ' '*(16-len(str(entry['key']))), line))

				
				for i in range(len(entry['qualifiers'])):
					line = self.check_line(entry['qualifiers'][i], 'feature')
					output.append(line)
		
		## add origin line ##
		origin = ''
		if self.gbfile['origin'] != None:
			origin = self.gbfile['origin']
		line = 'ORIGIN' + ' '*(12-len('ORIGIN')) + origin
		output.append(line)
		
		## add contig ##
		if self.gbfile['contig'] != None:
			line = 'CONTIG' + ' '*(12-len('CONTIG')) + self.gbfile['contig'] #need to fix for when it exceeds 80 characters
			output.append(line)
		
		## add sequence ##		
		if self.gbfile['dna'] != None:
			for i in range(0, len(self.gbfile['dna']), 60):
				output.append('%s%d %s %s %s %s %s %s' % (' '*(9-len(str(i))), i+1, self.gbfile['dna'][i:i+10], self.gbfile['dna'][i+10:i+20],self.gbfile['dna'][i+20:i+30],self.gbfile['dna'][i+30:i+40],self.gbfile['dna'][i+40:i+50],self.gbfile['dna'][i+50:i+60]))
		
		## add the last // ##
		output.append('//')		

		return '\n'.join(output)

		
	def Save(self, filepath=None):
		"""Function writes data stored in header, feature list and dna to a .gb file"""
		if filepath==None: #if no path is provided, attemt to get one stored in the gb data format
			filepath = self.GetFilepath()
		
		if filepath == None: #if it's still none
			return False
		else:
			#need to add conditions in case header and features are empty
			string = self.make_gbstring()
			a = open(filepath, 'w')
			a.write(string)
			a.close()
			return True




	def protein_mw(self, evt, protein):
		if ('*' in protein[:-1]) == True:
			pass
		else:
			pass
			#add function for calculating MW

	
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




