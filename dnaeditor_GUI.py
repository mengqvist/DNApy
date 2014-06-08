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
import features_GUI as features

#fix selection so that the statusbar matches the real



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
		

		#create splitter
		splitter = wx.SplitterWindow(self, 0, style=wx.SP_3D)	
		
		#create line count panel
#		self.linecount = output.create(self, style=wx.VSCROLL|wx.HSCROLL|wx.BORDER_NONE); #create DNA window
#		self.linecount.SetEditable(False)

		#create feature list view
		splitter1 = wx.SplitterWindow(self, 0, style=wx.SP_3D)	
		self.feature_list = features.FeatureCreate(splitter1, id=wx.ID_ANY, editor=False)
		

		#create dna view panel
		self.gbviewer = rt.RichTextCtrl(splitter1, style=wx.VSCROLL|wx.HSCROLL|wx.BORDER_NONE)

#		self.gbviewer = output.create(splitter1, style=wx.VSCROLL|wx.HSCROLL|wx.BORDER_NONE); #create DNA window
		self.gbviewer.SetEditable(False)	
		font = wx.Font(pointSize=10, family=wx.FONTFAMILY_DEFAULT, style=wx.FONTSTYLE_NORMAL, weight=wx.FONTWEIGHT_NORMAL, underline=False, faceName='Source Code Pro', encoding=wx.FONTENCODING_DEFAULT) #could also use Inconsolata
		self.gbviewer.SetFont(font)
		self.gbviewer.Bind(wx.EVT_CHAR, self.OnKeyPress) #any key press
		self.gbviewer.Bind(rt.EVT_RICHTEXT_SELECTION_CHANGED, self.OnSelChange)#these don't work!
		self.gbviewer.Bind(rt.EVT_RICHTEXT_SELECTION_CHANGED, self.OnSelChange)



		splitter1.SplitHorizontally(self.feature_list, self.gbviewer,sashPosition=-(windowsize[1]-270))

		sizer = wx.BoxSizer(wx.HORIZONTAL)
#		sizer.Add(splitter, -1, wx.EXPAND)
#		sizer.Add(self.linecount, proportion=1, flag=wx.EXPAND)
		sizer.Add(splitter1, -1, wx.EXPAND)
#		sizer.Add(self.feature_list, proportion=1, flag=wx.EXPAND)
#		sizer.Add(self.gbviewer, proportion=1, flag=wx.EXPAND)
		self.SetSizer(sizer)	
		
		self.Centre()

	def OnKeyPress(self, evt):
		key = evt.GetUniChar()
#		print(key)
#		shift = evt.ShiftDown()
		if key in [97, 65, 116, 84, 99, 67, 103, 71]: # [a, A, t, T, c, C, g, G']
			start, finish = self.gbviewer.GetSelection()
			if (start == -2 and finish == -2) == False: # if a selection, delete it
				self.GetTopLevelParent().gb.changegbsequence(start, finish, 'r', chr(key))
				#self.GetTopLevelParent().gb.Delete(start+1, finish+1)
			else:
				start = self.gbviewer.GetInsertionPoint()
				self.GetTopLevelParent().gb.changegbsequence(start, start, 'i', chr(key))
			self.updateUI()
			self.gbviewer.SetInsertionPoint(start+1)
			self.gbviewer.ShowPosition(start+1) 


		if key == 8: #backspace
			start, finish = self.gbviewer.GetSelection()
			if (start == -2 and finish == -2) == False: # if a selection, delete it
				self.GetTopLevelParent().gb.Delete(start+1, finish)
				self.updateUI()
				self.gbviewer.SetInsertionPoint(start)
				self.gbviewer.ShowPosition(start) 
			else:
				start = self.gbviewer.GetInsertionPoint()
				self.GetTopLevelParent().gb.Delete(start-1, start-1)
				self.updateUI()
				self.GetTopLevelParent().updateUndoRedo() #for updating the undo and redo buttons in the menu
				self.gbviewer.SetInsertionPoint(start-1)
				self.gbviewer.ShowPosition(start-1) 

		if key == 127: #delete
			start, finish = self.gbviewer.GetSelection()
			if (start == -2 and finish == -2) == False: # if a selection, delete it
				self.GetTopLevelParent().gb.Delete(start+1, finish)
			else:
				start = self.gbviewer.GetInsertionPoint()
				self.GetTopLevelParent().gb.Delete(start, start)
			self.updateUI()
			self.GetTopLevelParent().updateUndoRedo() #for updating the undo and redo buttons in the menu
			self.gbviewer.SetInsertionPoint(start)
			self.gbviewer.ShowPosition(start) 

	def OnSelChange(self, evt):
		print('sel changed')
##########################
	def make_outputpopup(self):
		'''Creates a popup window in which output can be printed'''
		self.outputframe = wx.Frame(None, title="Output Panel") # creation of a Frame with a title
		self.output = output.create(self.outputframe, style=wx.VSCROLL|wx.HSCROLL) # creation of a richtextctrl in the frame
		



#########################################

#	def dna_output(self, featurelist):
#		'''Prints output to the output panel'''
#		self.make_outputpopup()	
#		tabtext = str(self.gbviewer.GetPageText(self.gbviewer.GetSelection()))
#		DNA = featurelist[0]
#		self.output.write('%s | DNA in clipboard, %d bp' % (tabtext, len(DNA))+'\n', 'File')
#		self.output.write(DNA+'\n', 'DNA')
#		if len(featurelist) > 1:
#			self.output.write('With features: ', 'Text')
#			for i in range(1, len(featurelist)):
#				feature = featurelist[i]
#				self.output.write(('>%s "%s", ' % (feature['key'], feature['qualifiers'][0].split('=')[1])), 'Text')
#			self.output.write('\n', 'Text')
#		self.outputframe.Show()




################ genbank methods ###############


	def match_selection(self):
		'''Checks whether the dnaview selection matches that stored in the selection variable in genbank and if not, updates it'''
		viewerstart, viewerend = self.gbviewer.GetSelection()
		if viewerstart == -2 and viewerend == -2: # if not a selection
			viewerstart = self.gbviewer.GetInsertionPoint()
			selection = (viewerstart+1, viewerend)
		else:
			selection = (viewerstart+1, viewerend)
		self.GetTopLevelParent().set_dna_selection(selection)
		print(selection)



	def select_all(self):
		pass

	def uppercase(self):
		'''Change selection to uppercase'''
		self.match_selection()
		start, finish = self.GetTopLevelParent().get_dna_selection()
		self.GetTopLevelParent().gb.Upper(start, finish)
		self.updateUI()
		
	def lowercase(self):
		'''Change selection to lowercase'''
		self.match_selection()
		start, finish = self.GetTopLevelParent().get_dna_selection()
		self.GetTopLevelParent().gb.Lower(start, finish)
		self.updateUI() 

	def reverse_complement_selection(self):
		'''Reverse-complement current selection'''
		self.match_selection()
		start, finish = self.GetTopLevelParent().get_dna_selection()
		self.GetTopLevelParent().gb.RCselection(start, finish)
		self.updateUI()

#	def delete(self):
#		'''Deletes a selection and updates dna and features'''
#		self.match_selection()
#		start, finish = self.GetTopLevelParent().get_dna_selection()
#		self.GetTopLevelParent().gb.Delete(start, finish)
#		self.updateUI()
#		self.gbviewer.SetInsertionPoint(start-1)


	def cut(self):
		'''Cut DNA and store it in clipboard together with any features present on that DNA'''
		self.match_selection()
		start, finish = self.GetTopLevelParent().get_dna_selection()
		self.GetTopLevelParent().gb.Cut(start, finish)
		self.updateUI()
		self.gbviewer.SetInsertionPoint(start-1)
			
	def cut_reverse_complement(self):
		'''Cut reverse complement of DNA and store it in clipboard together with any features present on that DNA'''
		self.match_selection()
		start, finish = self.GetTopLevelParent().get_dna_selection()
		self.GetTopLevelParent().gb.CutRC(start, finish)
		self.updateUI()
		self.gbviewer.SetInsertionPoint(start-1)

	def paste(self):
		'''Paste DNA and any features present on that DNA'''
		self.match_selection()
		start, finish = self.GetTopLevelParent().get_dna_selection()
		if start != finish: #If a selection, remove sequence
			self.GetTopLevelParent().gb.Delete(start, finish, visible=False)
		self.GetTopLevelParent().gb.Paste(start)
		self.updateUI() # need to fix so that pasted dna is still highlighted
		self.gbviewer.SetInsertionPoint(start-1)
		
	def paste_reverse_complement(self):
		'''Paste reverse complement of DNA and any features present on that DNA'''
		self.match_selection()
		start, finish = self.GetTopLevelParent().get_dna_selection()
		if start != finish: #If a selection, remove sequence
			self.GetTopLevelParent().gb.Delete(start, finish, visible=False)
		self.GetTopLevelParent().gb.PasteRC(start)
		self.updateUI()
		self.gbviewer.SetInsertionPoint(start-1)

	def copy(self):
		'''Copy DNA and features into clipboard'''
		self.match_selection()
		start, finish = self.GetTopLevelParent().get_dna_selection()
		self.GetTopLevelParent().gb.Copy(start, finish)

	def copy_reverse_complement(self):
		'''Copy reverse complement of DNA'''
		self.match_selection()
		start, finish = self.GetTopLevelParent().get_dna_selection()
		self.GetTopLevelParent().gb.CopyRC(start, finish)


#######################################################


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
		self.gbviewer.SetValue(self.GetTopLevelParent().gb.GetDNA()) #put the DNA in

		#match selection to previous one
#		start, finish = self.GetTopLevelParent().get_dna_selection()
#		if start != finish: self.gbviewer.SetSelection(start, finish)
#		elif start == finish: self.gbviewer.SetInsertionPoint(start-1)
#		self.gbviewer.ShowPosition(start) 

		
		#returns a list of lists [[featuretype1, complement1, start1, end1], [featuretype2, complement2, start2, end2].....] 
		featurelist = self.GetTopLevelParent().gb.get_all_feature_positions()
		for entry in featurelist:
			featuretype, complement, start, finish = entry
			self.get_feature_color(featuretype, complement)
	
			color = self.gbviewer.current_highlight_color #get color

			#do the painting
			self.attr = rt.RichTextAttr()
			self.attr.SetFlags(wx.TEXT_ATTR_BACKGROUND_COLOUR) #do I need this flag?
			self.attr.SetBackgroundColour(color)
			self.gbviewer.SetStyleEx(rt.RichTextRange(start, finish), self.attr)

		
		#color in search hits if any are present
		search_hits = self.GetTopLevelParent().gb.search_hits
		if len(search_hits) != 0:
			color = '#ffff00'
			for i in range(len(search_hits)):
				start = int(search_hits[i][0])
				finish = int(search_hits[i][1])
				
				#do the painting
				self.attr = rt.RichTextAttr()
				self.attr.SetFlags(wx.TEXT_ATTR_BACKGROUND_COLOUR) #do I need this flag?
				self.attr.SetBackgroundColour(color)
				self.gbviewer.SetStyleEx(rt.RichTextRange(start, finish), self.attr)

#I need to think carefully about where to place these...
		#realize the current selection in the DNA editor
		start, finish = self.GetTopLevelParent().get_dna_selection()
		self.gbviewer.SetSelection(start, finish) #update the graphical selection
		self.gbviewer.ShowPosition(start) 

		#update feature list
		self.feature_list.updateUI()

		#use this to get the first line characters 
		##develop it##
		#set insertion point to beginning... Make it update on resize
#		while self.gbviewer.MoveDown() == True:
#			self.gbviewer.MoveDown()
#			self.linecount.write(str(self.gbviewer.GetInsertionPoint()+1)+'\n', 'Text')

		##############



####### Advanced DNA functions #######
#	def align_dna(self, evt):
#		'''Pass a list of dictionaries with name and dna to m_align to align sequences'''
#		#this works sort of ok. Needs improvement though...
#		
#		dna = str(self.GetTopLevelParent().gb.GetDNA())
#		
#		self.seqlist.append(dict(name='Reference', dna=dna))
#		print(self.seqlist)
#		result = seqfiles.M_align(self.seqlist)
#		self.output.write(result, 'text')

		
####### Protein functions #######
	def translate_output(self, protein, DNA, info):
		'''Generate output in the output.panel'''
#		tabtext = str(self.gbviewer.GetPageText(self.gbviewer.GetSelection()))
		self.make_outputpopup()			
		self.output.write('%s | Translate %s\n' % (tabtext, info), 'File')
		self.output.write(('%d AA from %d bases, %d bases left untranslated' % (len(protein), len(DNA), len(DNA)%3))+'\n', 'Text')
		self.output.write(protein, 'Protein')
		self.outputframe.Show()
	
	def translate_selection(self):
		'''Translate selected DNA'''
		DNA = self.gbviewer.GetStringSelection()
		protein = dna.translate(DNA)
		self.translate_output(protein, DNA, 'leading strand')
		
	def translate_selection_reverse_complement(self):
		'''Translate reverse-complement of selected DNA'''
		DNA = self.gbviewer.GetStringSelection()
		protein = dna.translate(dna.revcomp(DNA))
		self.translate_output(protein, DNA, 'complement strand')

#update this one...
	def translate_feature(self):
		'''Translate specified feature'''
		feature = self.GetTopLevelParent().gb.allgbfeatures[2]
		DNA = self.GetTopLevelParent().gb.getdnaforgbfeature(feature[4])
		protein = dna.translate(DNA)
		self.translate_output(protein, DNA, 'feature "%s"' % feature[4][7:])
	

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
#		mposition += 6
#		print(mposition)
		#which feature corresponds to this pos?
		Feature = self.GetTopLevelParent().gb.get_featurename_for_pos(mposition)
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


