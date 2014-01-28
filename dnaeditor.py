#!/usr/bin/env python


#This file is part of DNApy. DNApy is a DNA editor written purely in python. 
#The program is intended to be an intuitive, fully featured, 
#extendable, editor for molecular and synthetic biology.  
#Enjoy!
#
#Copyright (C) 2014  Martin K. M. Engqvist | 
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
class MyPanel(wx.Panel):
	def __init__(self, parent):
		wx.Panel.__init__(self, parent)
		
#		self.Bind(wx.EVT_KEY_DOWN, self.OnKeyPress)	

		#create splitter
		splitter = wx.SplitterWindow(self, 0, style=wx.SP_3D)	
		
		
		#create statusbar
#		self.statusbar = self.CreateStatusBar(2)
#		self.statusbar.SetStatusStyles(styles=[wx.SB_FLAT, wx.SB_FLAT])

		#create line count panel
#		self.linecount = output.create(self, style=wx.VSCROLL|wx.HSCROLL|wx.BORDER_NONE); #create DNA window
#		self.linecount.SetEditable(False)



		#create dna view panel
		self.gbviewer = output.create(self, style=wx.VSCROLL|wx.HSCROLL|wx.BORDER_NONE); #create DNA window
		self.gbviewer.SetEditable(False)	



		sizer = wx.BoxSizer(wx.HORIZONTAL)
#		sizer.Add(self.frame_1_toolbar, 0, wx.EXPAND)
#		sizer.Add(self.frame_2_toolbar, 0, wx.EXPAND)
#		sizer.Add(splitter, -1, wx.EXPAND)
#		sizer.Add(self.linecount, proportion=1, flag=wx.EXPAND)
		sizer.Add(self.gbviewer, proportion=10, flag=wx.EXPAND)
		self.SetSizer(sizer)	
		
		self.Centre()


##########################
	def make_outputpopup(self):
		'''Creates a popup window in which output can be printed'''
		self.outputframe = wx.Frame(None, title="Output Panel") # creation of a Frame with a title
		self.output = output.create(self.outputframe, style=wx.VSCROLL|wx.HSCROLL) # creation of a richtextctrl in the frame



	def changeme(self):
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
		self.make_outputpopup()		
		self.output.write(string+'\n', 'Text')
		self.outputframe.Show()
		
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

		self.make_outputpopup()		
		self.output.write(string+'\n', 'Text')
		self.outputframe.Show()


	
	def codon_table(self, event):
		'''Prints the standard codon table for translating DNA to protein'''
		string = '''
Put Table here
'''
		self.make_outputpopup()	
		self.output.write(string+'\n', 'Text')
		self.outputframe.Show()
 


#########################################

	def dna_output(self, featurelist):
		'''Prints output to the output panel'''
		self.make_outputpopup()	
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
		self.outputframe.Show()


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
#		self.search_word=wx.TextCtrl(self.frame_1_toolbar, wx.ID_ANY, "")
#		self.frame_1_toolbar.AddControl(self.search_word)		
		pass	





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
		string = genbank.gb.get_dna()[start:end]
		genbank.gb.changegbsequence(start, end, 'r', string.upper())
		
		#update viewer
		self.gbviewer.SetValue(genbank.gb.get_dna())
		self.updateUI()
		self.gbviewer.SetSelection(start, end)
		
	def Lowercase(self, evt):
		'''Change section to lowercase'''
		start, end = self.gbviewer.GetSelection()
		string = genbank.gb.get_dna()[start:end]
		genbank.gb.changegbsequence(start, end, 'r', string.lower())
		
		#update viewer
		self.gbviewer.SetValue(genbank.gb.get_dna())
		self.updateUI()
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
			genbank.gb.reverse_complement_clipboard()		
			genbank.gb.paste(cutstart)
			
			#realize changes
			self.gbviewer.SetValue(genbank.gb.get_dna()) #put dna in box
			self.updateUI()
			self.gbviewer.SetInsertionPoint(cutstart)


	def delete_selection(self): #done
		'''Deletes a selection and updates dna and features'''
		start, end = self.gbviewer.GetSelection()
		if start != -2 and end != -2: #must be a selection
			deletedsequence = genbank.gb.get_dna()[start:end]
			genbank.gb.changegbsequence(start+1, end+1, 'd', deletedsequence)

			self.gbviewer.SetValue(genbank.gb.get_dna()) #put dna in box
			self.updateUI()
			self.gbviewer.SetInsertionPoint(start)


	def Cut(self, evt):
		'''Cut DNA and store it in clipboard together with any features present on that DNA'''
		cutstart, cutend = self.gbviewer.GetSelection()
		if cutstart != -2 and cutend != -2: #must be a selection
			self.Copy("")
			self.delete_selection()
			
			#update viewer
			self.gbviewer.SetValue(genbank.gb.get_dna()) #put dna in box
			self.updateUI()
			self.gbviewer.SetInsertionPoint(cutstart)
			self.gbviewer.ShowPosition(cutstart) 
#			self.dna_output(genbank.gb.clipboard)
			
			
	def Cut_RevComp(self, evt):
		'''Cut reverse complement of DNA'''
		cutstart, cutend = self.gbviewer.GetSelection()
		if cutstart != -2 and cutend != -2: #must be a selection
			self.Copy("")
			self.delete_selection()
			genbank.gb.reverse_complement_clipboard()
			
			#update viewer
			self.gbviewer.SetValue(genbank.gb.get_dna()) #put dna in box
			self.updateUI()
			self.gbviewer.SetInsertionPoint(cutstart)	
#			self.dna_output(genbank.gb.clipboard)

		
		
	def Paste(self, evt):
		'''Paste DNA and any features'''
		control = wx.Window.FindFocus() #which field is selected?
		
#		if control == self.search_word: #the searchbox
#			self.search_word.SetValue(pyperclip.paste())

#		if control == self.output: #the outputpanel
#			pass
			
		if control == self.gbviewer: #the main dna window

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
				genbank.gb.paste(pastestart)

#			else: #if a selection
#				start, end =self.tab_list[self.current_tab].GetSelection()
#				self.delete_selection()  
#				self.tab_list[self.current_tab].SetInsertionPoint(start)   #this generates probelems with the features. Need to sepereate delete and insert somehow
#				genbank.gb.paste("")

			#update viwer
			self.gbviewer.SetValue(genbank.gb.get_dna()) #put dna in box
			self.updateUI()
			self.gbviewer.SetInsertionPoint(pastestart)
		

	
	def Paste_RevComp(self, evt):
		'''Paste reverse complement of DNA'''
		genbank.gb.reverse_complement_clipboard()
		self.Paste("")
		genbank.gb.reverse_complement_clipboard() #change it back


		
	def Copy(self, evt):
		'''Copy DNA and features into clipboard'''
		control = wx.Window.FindFocus() #which field is selected?

#		if control == self.search_word: #the searchbox
#			start, finish = self.search_word.GetSelection()
#			if start != -2 and finish != -2: #must be a selection
#				pyperclip.copy(self.search_word.GetValue()[start:finish])

#		if control == self.output: #the outputpanel
#			start, finish = self.output.GetSelection()
#			if start != -2 and finish != -2: #must be a selection
#				pyperclip.copy(self.output.GetValue()[start:finish])
			
		if control == self.gbviewer: #the main dna window	
			start, finish = self.gbviewer.GetSelection()
			if start != -2 and finish != -2: #must be a selection
				pyperclip.copy(genbank.gb.get_dna()[start:finish]) #copy dna to system clipboard (in case I want to paste it somwhere else)
				genbank.gb.copy(start, finish)
#				self.dna_output(genbank.gb.clipboard)

	
	def Copy_RevComp(self, evt):
		'''Copy reverse complement of DNA'''
		copystart, copyend = self.gbviewer.GetSelection()
		self.Copy("")
		genbank.gb.reverse_complement_clipboard()
#		self.dna_output(genbank.gb.clipboard)
		







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
					genbank.gb.allgbfeatures[n][4] = '%s_copy%s' % (feature, str(int(number)+1))



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
			


	
	def updateUI(self):
		'''For changing background color of text ranges'''



		self.remove_styling() #first remove old styles
		
		#returns a list of lists [[featuretype1, complement1, start1, end1], [featuretype2, complement2, start2, end2].....] 
		featurelist = genbank.gb.get_all_feature_positions()
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

		#use this to get the first line characters 
		##develop it##
		#set insertion point to beginning... Make it update on resize
#		while self.gbviewer.MoveDown() == True:
#			self.gbviewer.MoveDown()
#			self.linecount.write(str(self.gbviewer.GetInsertionPoint()+1)+'\n', 'Text')

		##############

	def update_text(self, ev):
		self.update_featurebuttons()		
#		self.update_statusbar("")
		
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
			featurelist = genbank.gb.allgbfeatures
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
#		dna = str(genbank.gb.get_dna())
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
		self.make_outputpopup()			
		self.output.write('%s | Translate %s\n' % (tabtext, info), 'File')
		self.output.write(('%d AA from %d bases, %d bases left untranslated' % (len(protein), len(DNA), len(DNA)%3))+'\n', 'Text')
		self.output.write(protein, 'Protein')
		self.outputframe.Show()
	
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

#update this one...
	def translate_feature(self, evt):
		'''Translate specified feature'''
		feature = genbank.gb.allgbfeatures[2]
		DNA = genbank.gb.getdnaforgbfeature(feature[4])
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

		
	def list_features(self, evt):
		'''List all features in output panel'''
		self.make_outputpopup()
#		tabtext = str(self.gbviewer.GetPageText(self.gbviewer.GetSelection()))
		tabtext = 'Replace!'
		self.output.write('%s | List features\n' % tabtext, 'File')
		featurelist = genbank.gb.list_features()
		self.output.write(featurelist, 'Text')
		self.outputframe.Show()
	


################ other functions ###############################


	def OnCloseWindow(self, e):
		self.close_all("")
		foo=self.GetSize()  ###except for the window size of file 
		if(self.IsMaximized()==0):
			file=open(files['size'], "w")
			file.write(str(foo[0])+"\n"+str(foo[1]))
			file.close()
		self.Destroy()
		self.updateUI() #refresh everything

		
	def mouse_position(self, event):
		'''Get which features are at a given position'''		
		xyposition = self.gbviewer.ScreenToClient(wx.GetMousePosition())
		#event.GetPosition() #this can be used if it controlled by an event..
		mposition = self.gbviewer.HitTest(xyposition)[1] #1 for the character num, 0 for linenum

		#which feature corresponds to this pos?
		Feature = genbank.gb.get_feature_for_pos(mposition)
		return mposition, Feature


###############################################################





if __name__ == '__main__': #if script is run by itself and not loaded	
	app = wx.App() # creation of the wx.App object (initialisation of the wxpython toolkit)



	#remove this and make it get the filename proper
	tabtext = 'replace!'

	frame = wx.Frame(None, title="DNA editor") # creation of a Frame with a title
	frame.dnaeditor = MyPanel(frame) # creation of a richtextctrl in the frame
	frame.Show() # frames are invisible by default so we use Show() to make them visible
	frame.dnaeditor.open_file("")
	app.MainLoop() # here the app enters a loop waiting for user input


