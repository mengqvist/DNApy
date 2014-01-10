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
#fix gb so that new file can be made using just dna
#fix header parsing
#rig catg keys so that when keys are pressed the dna actually goes into the file (and backspace and del keys too)
#when loading file, ask whether it should be cleaned from ApEinfo tags and those from vectorNTI
#clean up functions, organize them and take away everything that is not needed
#Change MyFrame to a panel instead. Move menus and toolbars out into a seperate frame. This is needed for later integration of the other parts.
#Fix "new document" 
#fix undo/redo


from copy import deepcopy
from string import *
from os import access,listdir
import sys, os
import wx
import wx.richtext as rt

from wxPython.lib.buttons import *
from wxPython.lib.colourselect import *

import pyperclip

import output
import dna
import genbank
#import seqfiles
#import featureeditor

files={}   #dictionary with all configuration files

files['default_dir'] = os.path.abspath(os.path.dirname(sys.argv[0]))+"/"
files['default_dir']=replace(files['default_dir'], "\\", "/")
files['default_dir']=replace(files['default_dir'], "library.zip", "")

variables=files['default_dir']+"variables"   ##path to the file of the global variables
settings=files['default_dir']+"settings"   ##path to the file of the global settings

execfile(variables) #gets all the pre-assigned variables
execfile(settings) #gets all the pre-assigned settings

#remove this and make it get the filename proper
tabtext = 'replace!'



########### class for text ######################
class MyFrame(wx.Frame):
	def __init__(self, parent, id, title):
		wx.Frame.__init__(self, parent, id, title, wx.DefaultPosition, wx.Size(1000, 800))
		
		#create toolbars
		self.__generate_toolbar()
		
		
		#create Menu Bar
		self.create_menu()		

		#create splitter
		splitter = wx.SplitterWindow(self, 0, style=wx.SP_3D)	
		
		
		#create statusbar
		self.statusbar = self.CreateStatusBar(2)
		self.statusbar.SetStatusStyles(styles=[wx.SB_FLAT, wx.SB_FLAT])

		#create dna view panel
		self.gbviewer = output.create(splitter, style=wx.VSCROLL|wx.HSCROLL); #create DNA window
		self.gbviewer.SetEditable(False)	
		
		#create output
		self.output = output.create(splitter, style=wx.VSCROLL|wx.HSCROLL); #create output window

			
		splitter.SplitHorizontally(self.gbviewer, self.output, sashPosition=500)

		sizer = wx.BoxSizer(wx.VERTICAL)
		sizer.Add(self.frame_1_toolbar, 0, wx.EXPAND)
		sizer.Add(self.frame_2_toolbar, 0, wx.EXPAND)
		sizer.Add(splitter, -1, wx.EXPAND)
		self.SetSizer(sizer)	
		
		self.Centre()

		

############## Info methods

	def info(self, event):
		string = '''
This file is part of DNApy. DNApy is a DNA editor written purely in python. 
The program is intended to be an intuitive, fully featured, 
extendable, editor for molecular and synthetic biology.  
Enjoy!

Copyright (C) 2014  Martin Engqvist | 

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
LICENSE:
This file is part of DNApy.
DNApy is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.
 
DNApy is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Library General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software Foundation,
Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Get source code at: https://github.com/0b0bby0/DNApy

'''
		self.output.SetInsertionPointEnd()
		self.output.write(string+'\n', 'Text')
		self.output.ShowPosition(self.output.GetLastPosition()) 
		
	def IUPAC_codes(self, event):
		'''Prints info about the IUPAC codes for DNA and amino acids'''
		string = '''
Nucleotide code 	Base
A 				Adenine
C 				Cytosine
G 				Guanine
T (or U) 			Thymine (or Uracil)
R 				A or G
Y 				C or T
S 				G or C
W 				A or T
K 				G or T
M 				A or C
B 				C or G or T
D 				A or G or T
H 				A or C or T
V 				A or C or G
N 				any base
. or - 				gap

Amino acid 	Three letter 	Amino acid
A 			Ala 			Alanine
C 			Cys 			Cysteine
D 			Asp 			Aspartic Acid
E 			Glu 			Glutamic Acid
F 			Phe 			Phenylalanine
G 			Gly 			Glycine
H 			His 			Histidine
I 			Ile 			Isoleucine
K 			Lys 			Lysine
L 			Leu 			Leucine
M 			Met 			Methionine
N 			Asn 			Asparagine
P 			Pro 			Proline
Q 			Gln 			Glutamine
R 			Arg 			Arginine
S 			Ser 			Serine
T 			Thr 			Threonine
V 			Val 			Valine
W 			Trp 			Tryptophan
Y 			Tyr 			Tyrosine
'''
		self.output.write(string+'\n', 'Text')

	
	
	def codon_table(self, event):
		'''Prints the standard codon table for translating DNA to protein'''
		string = '''
Put Table here
'''
		self.output.write(string+'\n', 'Text')
 


#########################################

	def dna_output(self, featurelist):
		'''t'''

		#generate output
#		tabtext = str(self.gbviewer.GetPageText(self.gbviewer.GetSelection()))
		DNA = featurelist[0]
		self.output.write('%s | DNA in clipboard, %d bp' % (tabtext, len(DNA))+'\n', 'File')
		self.output.write(DNA+'\n', 'DNA')
		if len(featurelist) > 1:
			self.output.write('With features: ', 'Text')
			for i in range(1, len(featurelist)):
				feature = featurelist[i]
				self.output.write(('>%s "%s", ' % (feature['key'], feature['qualifiers'][0].split('=')[1])), 'Text')
			self.output.write('\n', 'Text')
		self.output.ShowPosition(self.output.GetLastPosition()) 
		
	def new_document(self, evt):
		'''Create new gb file'''
		self.gb=genbank.gbobject() #make new gb in panel	

	def open_file(self, evt):
		'''Function for opening file'''
		self.dir_to_open = '/home/martin/Python_files/DNApy/'
		dlg = wx.FileDialog( self, style=wx.OPEN|wx.FILE_MUST_EXIST,   defaultDir=self.dir_to_open ,wildcard='GenBank files (*.gb)|*|Any file (*)|*')
		dlg.ShowModal()
		fileName = dlg.GetFilename()
		all_path=dlg.GetPath()
		dire=dlg.GetDirectory()
		dlg.Destroy()
		if(fileName == None or fileName == "" ):
			return1
		
		name, extension = fileName.split('.')
		if extension == 'gb':
			self.gb=genbank.gbobject() #make new gb in panel
			self.gb.readgb(all_path) #read the file
			self.gbviewer.SetValue(self.gb.get_dna()) #put dna in box
			self.paint_features() #update features		
		else:
			print("error, not a gb file")		

		self.Bind(wx.EVT_UPDATE_UI, self.update_statusbar)
		wx.EVT_CLOSE(self, self.OnCloseWindow)
		self.Bind(wx.EVT_KEY_DOWN, self.OnKeyPress)
#		wx.EVT_KEY_DOWN(self, self.OnKeyPress)
	
	def save_all(self, evt):
		pass
	
	def save_file(self, evt):
		'''Function for saving file'''
#		try:


		self.gb.write_file(self.gb.get_filepath())
		
#			self.OnUpdateUI(1) #update statusbar
#		except:
#			self.save_as_file("")
		

	def save_as_file(self, evt):
		'''Function for saving file as'''
		filepath = self.gb.get_filepath()	
		for i in range(len(filepath)): #get directory for file
			if filepath[i] == '/':
				dire = filepath[0:i+1]
		#get save dialog
		dlg = wx.FileDialog( self, style=wx.SAVE | wx.OVERWRITE_PROMPT,defaultDir=dire,wildcard='TXT files (*)|*|Any file (*)|*')
		dlg.ShowModal()
		all_path=dlg.GetPath()
		fileName=dlg.GetFilename()
		dire=dlg.GetDirectory()
		dlg.Destroy()
		if (fileName == None or fileName == ""):
			return
		else:
			#try:
			self.gb.update_filepath(all_path)
			self.save_file("")
			#except:
			#	error_window(7, self)



	def close_single(self, evt):
		pass
	def close_all(self, evt):
		pass

	def search_up(self, evt):
		'''Search upwards'''
		self.search(2)
	
	def search_down(self, evt):
		'''Search downwards'''
		self.search(1)
	
	def search(self, searchtype):
		'''This is the function for the actual search'''
		if(searchtype==2):
			end=0
			start=self.gbviewer.GetSelection()[1]-1 #selection?
			if start == -3:
				start=self.gbviewer.GetInsertionPoint()-1 #or single char position?
			text = self.gbviewer.GetValue() #get document text
			searchword = self.search_word.GetValue() #get searchword
			for i in range(len(text)):
				if searchword == text[start-i-len(searchword):start-i]: #go through text until searchword is found
					self.gbviewer.SetSelection(start-i-len(searchword),start-i)
					break
				
		else:
			start=self.gbviewer.GetSelection()[0]+1
			if start == -3:
				start=self.gbviewer.GetInsertionPoint()+1
			text = self.gbviewer.GetValue()
			searchword = self.search_word.GetValue()
			for i in range(len(text)):
				if searchword == text[start+i:start+i+len(searchword)]:
					self.gbviewer.SetSelection(start+i, start+i+len(searchword))
					break

	def search_tool(self):
		'''Generates the search textbox'''
		self.search_word=wx.TextCtrl(self.frame_1_toolbar, wx.ID_ANY, "")
		self.frame_1_toolbar.AddControl(self.search_word)		
	


	def quit(self, evt):
		'''Function for quiting program'''
		print("close")
#    		self.close_all("")
       		self.Close()
       		#self.Destroy()		


	def OnCloseWindow(self, evt):
		'''not currently used'''
#		self.close_all("")
#		foo=self.GetSize()  ###except for the window size of file 
#		if(self.IsMaximized()==0):
#			file=open(files['size'], "w")
#			file.write(str(foo[0])+"\n"+str(foo[1]))
#			file.close()
		self.Close()
		#self.Destroy()


################ genbank methods ###############
	def OnKeyPress(self, evt):
		'''Checks which key is pressed and if one of them is A, T, C or G inserts the base into the file'''
		#this is not working
		print('ok')
		print(evt)

		keycode = event.GetUnicodeKey()
#		if keycode == :
#			insert...
	
		evt.Skip()

	def select_all(self, evt):
		pass

		
	def Uppercase(self, evt):
		'''Change section to uppercase'''
		start, end = self.gbviewer.GetSelection()
		string = self.gb.get_dna()[start:end]
		self.gb.changegbsequence(start, end, 'r', string.upper())
		
		#update viewer
		self.gbviewer.SetValue(self.gb.get_dna())
		self.paint_features()
		self.gbviewer.SetSelection(start, end)
		
	def Lowercase(self, evt):
		'''Change section to lowercase'''
		start, end = self.gbviewer.GetSelection()
		string = self.gb.get_dna()[start:end]
		self.gb.changegbsequence(start, end, 'r', string.lower())
		
		#update viewer
		self.gbviewer.SetValue(self.gb.get_dna())
		self.paint_features()
		self.gbviewer.SetSelection(start, end)


	def Undo(self, evt):
		pass
		#I should store each genbank file as a deep copy in a list, that way I can get them back by undoing/redoing

	def Redo(self, evt):
		pass

	def RevComp_sel(self, evt):
		'''Funciton to reverse-complement the current selection'''
		cutstart, cutend = self.gbviewer.GetSelection()
		if cutstart != -2 and cutend != -2: #must be a selection

			#make changes
			self.Copy("")
			self.delete_selection()
			self.gb.reverse_complement_clipboard()		
			self.gb.paste(cutstart)
			
			#realize changes
			self.gbviewer.SetValue(self.gb.get_dna()) #put dna in box
			self.paint_features()
			self.gbviewer.SetInsertionPoint(cutstart)




	def delete_selection(self): #done
		'''Deletes a selection and updates dna and features'''
		start, end = self.gbviewer.GetSelection()
		if start != -2 and end != -2: #must be a selection
			deletedsequence = self.gb.get_dna()[start:end]
			self.gb.changegbsequence(start+1, end+1, 'd', deletedsequence)

			self.gbviewer.SetValue(self.gb.get_dna()) #put dna in box
			self.paint_features()
			self.gbviewer.SetInsertionPoint(start)

	def Cut(self, evt):
		'''t'''
		cutstart, cutend = self.gbviewer.GetSelection()
		if cutstart != -2 and cutend != -2: #must be a selection
			self.Copy("")
			self.delete_selection()
			
			#update viewer
			self.gbviewer.SetValue(self.gb.get_dna()) #put dna in box
			self.paint_features()
			self.gbviewer.SetInsertionPoint(cutstart)
#			self.dna_output(self.gb.clipboard)
			
	def Cut_RevComp(self, evt):
		'''Cut reverse complement of DNA'''
		cutstart, cutend = self.gbviewer.GetSelection()
		if cutstart != -2 and cutend != -2: #must be a selection
			self.Copy("")
			self.delete_selection()
			self.gb.reverse_complement_clipboard()
			
			#update viewer
			self.gbviewer.SetValue(self.gb.get_dna()) #put dna in box
			self.paint_features()
			self.gbviewer.SetInsertionPoint(cutstart)	
#			self.dna_output(self.gb.clipboard)

		
		
	def Paste(self, evt):
		'''Paste DNA'''
		control = wx.Window.FindFocus() #which field is selected?
		
		if control == self.search_word: #the searchbox
			self.search_word.SetValue(pyperclip.paste())

		elif control == self.output: #the outputpanel
			pass
			
		elif control == self.gbviewer: #the main dna window

			#add filter to remove non-dna characters
		
			#figure out whether to over-write something
			pastestart, pasteend = self.gbviewer.GetSelection()
			if pastestart == -2 and pasteend == -2: #not a selection
				pastestart = self.gbviewer.GetInsertionPoint()
				pasteend = False

			#if paste is from raw sequence
			#if.....
		
		
			#if paste is from copied/cut genbank dna
			#elif.....
			if pasteend == False:
				self.gb.paste(pastestart)

#			else: #if a selection
#				start, end =self.tab_list[self.current_tab].GetSelection()
#				self.delete_selection()  
#				self.tab_list[self.current_tab].SetInsertionPoint(start)   #this generates probelems with the features. Need to sepereate delete and insert somehow
#				self.gb.paste("")

			#update viwer
			self.gbviewer.SetValue(self.gb.get_dna()) #put dna in box
			self.paint_features()
			self.gbviewer.SetInsertionPoint(pastestart)
		
			#generate output
			self.output.write('Sequence pasted'+'\n', 'Text')

	
	def Paste_RevComp(self, evt):
		'''Paste reverse complement of DNA'''
		self.gb.reverse_complement_clipboard()
		self.Paste("")
		self.gb.reverse_complement_clipboard() #change it back
		#generate output
		self.output.write('Reverse complement of sequence pasted'+'\n', 'Text')


		
	def Copy(self, evt):
		'''Copy DNA and features into clipboard'''
		control = wx.Window.FindFocus() #which field is selected?

		if control == self.search_word: #the searchbox
			start, finish = self.search_word.GetSelection()
			if start != -2 and finish != -2: #must be a selection
				pyperclip.copy(self.search_word.GetValue()[start:finish])

		elif control == self.output: #the outputpanel
			start, finish = self.output.GetSelection()
			if start != -2 and finish != -2: #must be a selection
				pyperclip.copy(self.output.GetValue()[start:finish])
			
		elif control == self.gbviewer: #the main dna window	
			start, finish = self.gbviewer.GetSelection()
			if start != -2 and finish != -2: #must be a selection
				pyperclip.copy(self.gb.get_dna()[start:finish]) #copy dna to system clipboard (in case I want to paste it somwhere else)
				self.gb.copy(start, finish)
#				self.dna_output(self.gb.clipboard)

	
	def Copy_RevComp(self, evt):
		'''Copy reverse complement of DNA'''
		copystart, copyend = self.gbviewer.GetSelection()
		self.Copy("")
		self.gb.reverse_complement_clipboard()
#		self.dna_output(self.gb.clipboard)
		







########## change this #########################################		

	def check_for_unique_feature_names(self, featurelist):
		'''Function for ensuring that all features have unique names'''
		for i in range(1, len(featurelist)):
			for n in range(1, len(featurelist)):
				if i == n:
					pass
				elif featurelist[n][4] == featurelist[i][4]:
					if ('_copy' in featurelist[n][4]) == True:
						feature, number = featurelist[n][4].split('_copy')
					else:
						feature = featurelist[n][4]
						number = 0
					self.gb.allgbfeatures[n][4] = '%s_copy%s' % (feature, str(int(number)+1))








######### Functions for updating text and text highlighting #############
	def remove_styling(self):
		'''Removes background styling for whole document'''
		start = 0
		finish = self.gbviewer.GetLastPosition() #last position in file
		self.attr = rt.RichTextAttr()
		self.attr.SetFlags(wx.TEXT_ATTR_BACKGROUND_COLOUR) #do I need this flag?
		self.attr.SetBackgroundColour('#FFFFFF')
		self.gbviewer.SetStyleEx(rt.RichTextRange(start, finish), self.attr)
	
	def get_feature_color(self, featuretype, complement):
		'''Takes single feature and finds its color, if any'''
		#piece of code to check if there are assigned colors to the features
		
		try: 
			color = eval(featuretype)
							
		except:
			color = grey	
		
		
		if complement == True: self.gbviewer.current_highlight_color = color['rv']
		else: self.gbviewer.current_highlight_color = color['fw']
			


	
	def paint_features(self):
		'''For changing background color of text ranges'''
		self.remove_styling() #first remove old styles
		
		#returns a list of lists [[featuretype1, complement1, start1, end1], [featuretype2, complement2, start2, end2].....] 
		featurelist = self.gb.get_all_feature_positions()
		for entry in featurelist:
			featuretype, complement, start, finish = entry

			self.get_feature_color(featuretype, complement)
	
			color = self.gbviewer.current_highlight_color #get color

			#do the painting
			self.attr = rt.RichTextAttr()
			self.attr.SetFlags(wx.TEXT_ATTR_BACKGROUND_COLOUR) #do I need this flag?
			self.attr.SetBackgroundColour(color)
			self.gbviewer.SetStyleEx(rt.RichTextRange(start, finish), self.attr)
		self.update_text("")

	def update_text(self, ev):
		self.update_featurebuttons()		
		self.update_statusbar("")
		
#		self.tab_list[self.current_tab].modify=1
#		self.current_tab=self.gbviewer.GetSelection()
#		if self.tab_list[self.current_tab].CanUndo(): #if undo is available
#			self.menubar.Enable(9, 1)
#			self.frame_1_toolbar.EnableTool(513, 1) #undo
#		else:
#			self.menubar.Enable(9, 0)
#			self.frame_1_toolbar.EnableTool(513, 0) #undo
#			
#		if self.tab_list[self.current_tab].CanRedo(): #if redo is available
#			self.menubar.Enable(10, 1)
#			self.frame_1_toolbar.EnableTool(514, 1) #redo
#		else:
#			self.menubar.Enable(10, 0)
#			self.frame_1_toolbar.EnableTool(514, 0) #redo
#		self.OnUpdateUI("")	

	def update_featurebuttons(self):
		'''Mehod for generating featurebuttons in the top toolbar'''
			
		try:			
			#make new ones from featurelist
			features = ['DNA']
			featurelist = self.gb.allgbfeatures
			for entry in featurelist:
				feature = entry[4].split('=')[1] #remove the /label= part
				features.append(feature)
						
			#update the 
			#self.featurebox......
			
#			n = 0
#			for entry in featurelist:
#				n += 1
#				label = entry[4][7:]
#				self.button = wx.Button(self.frame_1_toolbar_vertical, id=n, size=(200,30), label=label)
#				self.button.name = label
#				self.button.Bind(wx.EVT_BUTTON, self.feature_dialog) 
				
				#leftklick should select feature
				#rightklick should bring up menu with: Edit, mutate, report, hide, delete
				
				
	#			self.button.Bind(wx.EVT_ENTER_WINDOW, self.OnButtonEnter) #for some reason these don't play nice with wx.EVT_BUTTON
	#			self.button.Bind(wx.EVT_LEAVE_WINDOW, self.OnButtonLeave)
#				self.frame_1_toolbar_vertical.AddControl(self.button)
					
				#setting color
#				self.get_feature_color(entry)
#				color = self.tab_list[self.current_tab].current_highlight_color #get color
#				self.button.SetBackgroundColour('#FCFCFC')		
#				self.frame_1_toolbar_vertical.Refresh()
		except:
			pass


#####################################################





####### Advanced DNA functions #######
#	def align_dna(self, evt):
#		'''Pass a list of dictionaries with name and dna to m_align to align sequences'''
#		#this works sort of ok. Needs improvement though...
#		
#		dna = str(self.gb.get_dna())
#		
#		self.seqlist.append(dict(name='Reference', dna=dna))
#		print(self.seqlist)
#		result = seqfiles.M_align(self.seqlist)
#		self.output.write(result, 'text')
	
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
		
####### Protein functions #######
	def translate_output(self, protein, DNA, info):
		'''Generate output in the output.panel'''
#		tabtext = str(self.gbviewer.GetPageText(self.gbviewer.GetSelection()))
		self.output.SetInsertionPointEnd()
		self.output.write('%s | Translate %s ' % (tabtext, info), 'File')
		self.output.write(('%d AA from %d bases, %d bases left untranslated' % (len(protein), len(DNA), len(DNA)%3))+'\n', 'Text')
		self.output.write(protein, 'Protein')
		self.output.Newline()
		self.output.ShowPosition(self.output.GetLastPosition()) 
	
	def translate_selection(self, evt):
		'''Translate selected DNA'''
		DNA = self.gbviewer.GetStringSelection()
		protein = dna.translate(DNA)
		self.translate_output(protein, DNA, 'leading strand')
		
	def translate_selection_reverse_complement(self, evt):
		'''Translate reverse-complement of selected DNA'''
		DNA = self.gbviewer.GetStringSelection()
		protein = dna.translate(dna.reversecomplement(DNA))
		self.translate_output(protein, DNA, 'complement strand')
				
	def translate_feature(self, evt):
		'''Translate specified feature'''
		feature = self.gb.allgbfeatures[2]
		DNA = self.gb.getdnaforgbfeature(feature[4])
		protein = dna.translate(DNA)
		self.translate_output(protein, DNA, 'feature "%s"' % feature[4][7:])
	
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
	
	def insert_barcode():
		pass
	
	def read_barcode():
		pass	
	
	def protein_pattern_find():
		pass
	
####### Feature functions #######

	
#	def features_output(self, feature, info):
#		#generate output
#		tabtext = str(self.gbviewer.GetPageText(self.gbviewer.GetSelection()))
#		self.output.SetInsertionPointEnd() 
#		self.output.write('%s | List features' % tabtext, 'File')
#		self.output.write(('>%s "%s" from %s to %s on %s strand %s' % (feature[0], feature[4][7:], feature[1], feature[2], feature[3], info)), 'Text')
#		self.output.Newline()
#		self.output.ShowPosition(self.output.GetLastPosition()) 

		
	def remove_feature(self, feature):
		self.gb.removegbfeature(feature)
		self.features_output(allgbfeatures[2], 'removed')
		
	def new_feature():
		pass


	def list_features(self, evt):
		'''List all features in output panel'''
#		tabtext = str(self.gbviewer.GetPageText(self.gbviewer.GetSelection()))
		tabtext = 'Replace!'
		self.output.SetInsertionPointEnd() 
		self.output.write('%s | List features\n' % tabtext, 'File')
		featurelist = self.gb.list_features()
		self.output.write(featurelist, 'Text')
		self.output.ShowPosition(self.output.GetLastPosition()) 

	
		






################ other functions ###############################


	def OnCloseWindow(self, e):
		self.close_all("")
		foo=self.GetSize()  ###except for the window size of file 
		if(self.IsMaximized()==0):
			file=open(files['size'], "w")
			file.write(str(foo[0])+"\n"+str(foo[1]))
			file.close()
		self.Destroy()
		self.paint_features() #refresh everything

	def update_statusbar(self, evt):
		'''Updates statusbar'''
		#this stuff is for the statusbar
#		if len(self.tab_list) == 0:
#			string = 'File unmodified'
#		elif self.tab_list[self.current_tab].modify==0:
#			string = 'File unmodified'
#		elif self.tab_list[self.current_tab].modify==1:
#			string = 'File not yet saved'
		
		mposition, Feature = self.mouse_position("") #get mouse position
		
		
		try:
			Position = str(mposition+1)
		except:
			Position = ""
		
		try:
			Feature = str(Feature)
		except:
			Feature = ""
		
		try:		
			SelectionFrom, SelectionTo = (str(self.gbviewer.GetSelection()[0]+1), str(self.gbviewer.GetSelection()[1]))
			if SelectionFrom == '-1' and SelectionTo == '-2': #no selection if true
				SelectionFrom, SelectionTo = ("0", "0")
		except:
			SelectionFrom, SelectionTo = ("0", "0")
		try:	
			Length = str(self.gbviewer.GetSelection()[1] - self.gbviewer.GetSelection()[0])
		except:
			Length = ""


		self.SetStatusText('Position: %s      Feature: %s' % (Position, Feature), 0) #text in first field
		
		if float(Length)/3 == 1: #if one triplet is selected, show the AA
			AA = ': %s' % dna.translate(self.gbviewer.GetStringSelection())
		else:
			AA = ''
			
		self.SetStatusText('Selection: %s to %s,   %s bp,   %.1f AA%s' % (SelectionFrom, SelectionTo, Length, float(Length)/3, AA), 1) #text in second field

		
	def mouse_position(self, event):
		'''Get which features are at a given position'''		
		xyposition = self.gbviewer.ScreenToClient(wx.GetMousePosition())
		#event.GetPosition() #this can be used if it controlled by an event..
		mposition = self.gbviewer.HitTest(xyposition)[1] #1 for the character num, 0 for linenum

		#which feature corresponds to this pos?
		Feature = self.gb.get_feature_for_pos(mposition)
		return mposition, Feature


	def view_gbfile():
		pass


###############################################################


	def __generate_toolbar(self):
		'''For generating toolbar with buttons'''
				
		self.frame_1_toolbar = wx.ToolBar(self, wx.ID_ANY, style=wx.TB_HORIZONTAL|wx.TB_FLAT|wx.TB_DOCKABLE)

   		#syntax for toolbar
   		#AddLabelTool(self, id, label, bitmap, bmpDisabled, kind, shortHelp, longHelp, clientData) 
   		

   		
   		#New Document
   		self.frame_1_toolbar.AddLabelTool(500, "New Document", wx.Bitmap(files['default_dir']+"/icon/new.png", wx.BITMAP_TYPE_ANY), wx.NullBitmap, wx.ITEM_NORMAL, 'New File', "New File") #last one is the one displayed in status bar
   		wx.EVT_TOOL(self, 500, self.new_document)
		#Open File
   		self.frame_1_toolbar.AddLabelTool(501, "Open File", wx.Bitmap(files['default_dir']+"/icon/open.png", wx.BITMAP_TYPE_ANY), wx.NullBitmap, wx.ITEM_NORMAL, 'Load File', 'Load File')
   		wx.EVT_TOOL(self, 501, self.open_file)
		#Save current file
   		self.frame_1_toolbar.AddLabelTool(502, "Save current file", wx.Bitmap(files['default_dir']+"/icon/save.png", wx.BITMAP_TYPE_ANY), wx.NullBitmap, wx.ITEM_NORMAL, 'Save File', 'Save File')
   		wx.EVT_TOOL(self, 502, self.save_file)
		#Save as
   		self.frame_1_toolbar.AddLabelTool(503, "Save as", wx.Bitmap(files['default_dir']+"/icon/saveas.png", wx.BITMAP_TYPE_ANY), wx.NullBitmap, wx.ITEM_NORMAL, 'Save File As', 'Save File As')
   		wx.EVT_TOOL(self, 503, self.save_as_file)
		#Cut
   		self.frame_1_toolbar.AddLabelTool(504, "Cut", wx.Bitmap(files['default_dir']+"/icon/cut.png", wx.BITMAP_TYPE_ANY), wx.NullBitmap, wx.ITEM_NORMAL, 'Cut', 'Cut')
   		wx.EVT_TOOL(self, 504, self.Cut)
		#Copy
   		self.frame_1_toolbar.AddLabelTool(505, "Copy", wx.Bitmap(files['default_dir']+"/icon/copy.png", wx.BITMAP_TYPE_ANY), wx.NullBitmap, wx.ITEM_NORMAL, 'Copy', 'Copy')
   		wx.EVT_TOOL(self, 505, self.Copy)
		#Paste
   		self.frame_1_toolbar.AddLabelTool(506, "Paste", wx.Bitmap(files['default_dir']+"/icon/paste.png", wx.BITMAP_TYPE_ANY), wx.NullBitmap, wx.ITEM_NORMAL, 'Paste', 'Paste')
   		wx.EVT_TOOL(self, 506, self.Paste)
   		#Undo
   		self.frame_1_toolbar.AddLabelTool(513, "Undo", wx.Bitmap(files['default_dir']+"/icon/undo.png", wx.BITMAP_TYPE_ANY), wx.NullBitmap, wx.ITEM_NORMAL, 'Undo', 'Undo')
   		wx.EVT_TOOL(self, 513, self.Undo)   
   		#Redo
   		self.frame_1_toolbar.AddLabelTool(514, "Redo", wx.Bitmap(files['default_dir']+"/icon/redo.png", wx.BITMAP_TYPE_ANY), wx.NullBitmap, wx.ITEM_NORMAL, 'Redo', 'Redo')
   		wx.EVT_TOOL(self, 514, self.Redo) 
		#Search Upward
   		self.frame_1_toolbar.AddLabelTool(507, "Search Upward", wx.Bitmap(files['default_dir']+"/icon/up.png", wx.BITMAP_TYPE_ANY), wx.NullBitmap, wx.ITEM_NORMAL, 'Search Upward', 'Search Upward')
   		wx.EVT_TOOL(self, 507, self.search_up)
		#Search window
		self.search_tool()
		#Search Downward
   		self.frame_1_toolbar.AddLabelTool(509, "Search Downward", wx.Bitmap(files['default_dir']+"/icon/down.png", wx.BITMAP_TYPE_ANY), wx.NullBitmap, wx.ITEM_NORMAL, 'Search Downward', 'Search Downward')
   		wx.EVT_TOOL(self, 509, self.search_down)
		#Print current window
#   		self.frame_1_toolbar.AddLabelTool(510, "Print current window", wx.Bitmap(files['default_dir']+"/icon/print.png", wx.BITMAP_TYPE_ANY), wx.NullBitmap, wx.ITEM_NORMAL, 'Print Current Window', 'Print Current Window')
 #  		wx.EVT_TOOL(self, 510, self.print_setup)
		
		self.frame_1_toolbar.Realize()
		
		
		##### Toolbar 2 #####
		
		self.frame_2_toolbar = wx.ToolBar(self, wx.ID_ANY, style=wx.TB_HORIZONTAL|wx.TB_FLAT|wx.TB_DOCKABLE)
			
		#Select or Mutate
		self.selormut = wx.ComboBox(self.frame_2_toolbar, id=wx.ID_ANY, size=(85, 28), choices=['Select', 'Mutate'], style=wx.CB_READONLY)
		self.frame_2_toolbar.AddControl(self.selormut)		
		self.selormut.SetSelection(0)
		
		#nucleotide or amino acid
		self.nucleotideoraminoacid = wx.ComboBox(self.frame_2_toolbar, id=wx.ID_ANY, size=(120, 28), choices=['Nucleotide', 'Amino Acid', 'Feature'], style=wx.CB_READONLY)
		self.frame_2_toolbar.AddControl(self.nucleotideoraminoacid)
		self.nucleotideoraminoacid.SetSelection(0)
		
		#'position'
		self.positionbox=wx.TextCtrl(self.frame_2_toolbar, id=wx.ID_ANY, size=(70, 25), value="position")
		self.frame_2_toolbar.AddControl(self.positionbox)
		self.positionbox.SetEditable(False)	
		
		
		#'pos #'
		self.posnum=wx.TextCtrl(self.frame_2_toolbar, id=wx.ID_ANY, size=(90, 25), value="")
		self.frame_2_toolbar.AddControl(self.posnum)
		
		
		#'in'
		self.inbox=wx.TextCtrl(self.frame_2_toolbar, id=wx.ID_ANY, size=(25, 25), value="in")
		self.frame_2_toolbar.AddControl(self.inbox)
		self.inbox.SetEditable(False)
				
		#featurebox
		self.featurebox = wx.ComboBox(self.frame_2_toolbar, id=wx.ID_ANY, size=wx.DefaultSize, choices=['File', 'feature1'], style=wx.CB_READONLY)
		self.frame_2_toolbar.AddControl(self.featurebox)
		
		#'go'

		self.frame_2_toolbar.Realize()






	def create_menu(self):     #method for creating menu
		self.menubar = wx.MenuBar()
		fileitem = wx.Menu()
			
		#new document
		fileitem.Append(1, "New\tCtrl+Q", "New Document")
		wx.EVT_MENU(self, 1, self.new_document)

		#open document
		fileitem.Append(2, "Open\tCtrl+O", "Open Document")
		wx.EVT_MENU(self, 2, self.open_file)
		fileitem.AppendSeparator()

		#save document
		fileitem.Append(3, "Save\tCtrl+S", "Save current document")
		wx.EVT_MENU(self, 3, self.save_file)

		#save document as
		fileitem.Append(4, "Save as", "Save a copy of current document")
		wx.EVT_MENU(self, 4, self.save_as_file)

		#save all
		fileitem.Append(5, "Save all", "Save all open tabs")
		wx.EVT_MENU(self, 5, self.save_all)
		fileitem.AppendSeparator()

		#close single
		fileitem.Append(5, "Close", "Close current document")
		wx.EVT_MENU(self, 5, self.close_single)

		#close all
		fileitem.Append(6, "Close all", "Close all tabs")
		wx.EVT_MENU(self, 6, self.close_all)
		fileitem.AppendSeparator()

		#quit
		fileitem.Append(7, "Exit", "Exit program")
		wx.EVT_MENU(self, 7, self.quit)

		self.menubar.Append(fileitem, "&File")

		######################### For 'Edit DNA' menu item #############################################
		self.edit = wx.Menu()
		#undo
		self.edit.Append(9, "Undo\tCtrl+Z", "Undo")
		wx.EVT_MENU(self, 9, self.Undo)

		#redo
		self.edit.Append(10, "Redo\tCtrl+Y", "Redo")
		wx.EVT_MENU(self, 10, self.Redo)
		self.edit.AppendSeparator() #________________________devider

		#cut
		self.edit.Append(11, "Cut\tCtrl+X", "Cut selected DNA")
		wx.EVT_MENU(self,11, self.Cut)

		#copy
		self.edit.Append(12, "Copy\tCtrl+C", "Copy selected DNA")
		wx.EVT_MENU(self, 12, self.Copy)

		#paste
		self.edit.Append(13, "Paste\tCtrl+V", "Paste DNA")
		wx.EVT_MENU(self, 13, self.Paste)
		
		#cut reverse complement
		self.edit.Append(111, "Cut Rev-Comp\tCtrl+Shift+X", "Cut reverse-complement of selected DNA")
		wx.EVT_MENU(self,111, self.Cut_RevComp)

		#copy reverse complement
		self.edit.Append(121, "Copy Rev-Comp\tCtrl+Shift+C", "Copy reverse-complement of selected DNA")
		wx.EVT_MENU(self, 121, self.Copy_RevComp)		
		
		#paste reverse complement
		self.edit.Append(131, "Paste Rev-Comp\tCtrl+Shift+V", "Paste reverse-complement of DNA")
		wx.EVT_MENU(self, 131, self.Paste_RevComp)

		#reverse-complement selection
		self.edit.Append(141, "Rev-Comp selection\tCtrl+Shift+R", "Reverse-complement seleected DNA")
		wx.EVT_MENU(self,141, self.RevComp_sel)
		self.edit.AppendSeparator() #________________________devider
		
		#select all
		self.edit.Append(14, "Select all", "Select all text")
		wx.EVT_MENU(self, 14, self.select_all)
		self.edit.AppendSeparator() #________________________devider

		#uppercase
		self.edit.Append(34, "Uppercase\tCtrl+W", "Convert selected text to uppercase")
		wx.EVT_MENU(self, 34, self.Uppercase)

		#lowercase
		self.edit.Append(35, "Lowercase\tCtrl+E", "Convert selected text to lowercase")
		wx.EVT_MENU(self, 35, self.Lowercase)
		self.edit.AppendSeparator() #________________________devider

		self.menubar.Append(self.edit, "Edit DNA")
		
	
	
		######## Features menu item ########
		self.features = wx.Menu()
		
		self.features.Append(40, "List Features", "List all features in file")
		wx.EVT_MENU(self,40, self.list_features)
		
#		self.features.Append(41, "Edit Features", "Edit features in file")
#		wx.EVT_MENU(self,41, self.edit_features)		
		
		self.menubar.Append(self.features, "Features")
		


		######## Protein menu item #########
		self.protein = wx.Menu()
		
		#translate
		self.protein.Append(30, "Translate\tCtrl+T", "Translate DNA to protein")
		wx.EVT_MENU(self,30, self.translate_selection)

		#translate reverse complement
		self.protein.Append(31, "Translate Rev-Comp\tCtrl+Shift+T", "Translate DNA to protein")
		wx.EVT_MENU(self,31, self.translate_selection_reverse_complement)

		#translate feature
		self.protein.Append(32, "Translate feature", "Translate DNA feature to protein")
		wx.EVT_MENU(self,32, self.translate_feature)

		self.menubar.Append(self.protein, "Protein")



		########## For 'Help' menu item #############
		self.help = wx.Menu()
		#about program
		self.help.Append(21, "About", "About this program")
		wx.EVT_MENU(self, 21, self.info)

		#print IUPAC codes for dna and AA
		self.help.Append(22, "IUPAC codes", "IUPAC codes of DNA and amino acids")
		wx.EVT_MENU(self, 22, self.IUPAC_codes)


		self.menubar.Append(self.help, "Help")

		self.SetMenuBar(self.menubar)

##### main loop
class MyApp(wx.App):
    def OnInit(self):
        frame = MyFrame(None, -1, "DNApy - GBviewer")
        frame.Show(True)
        self.SetTopWindow(frame)
        return True
        
if __name__ == '__main__': #if script is run by itself and not loaded	
	app = MyApp(0)
	app.MainLoop()


