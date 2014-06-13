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
import string
from os import access,listdir
import sys, os
import wx
import wx.stc
from wx.lib.pubsub import Publisher as pub

import pyperclip

import output
import dna
import genbank

import features_GUI

#fix selection so that the statusbar matches the real



files={}   #dictionary with all configuration files

files['default_dir'] = os.path.abspath(os.path.dirname(sys.argv[0]))+"/"
files['default_dir']=string.replace(files['default_dir'], "\\", "/")
files['default_dir']=string.replace(files['default_dir'], "library.zip", "")

variables=files['default_dir']+"variables"   ##path to the file of the global variables
settings=files['default_dir']+"settings"   ##path to the file of the global settings

execfile(variables) #gets all the pre-assigned variables
execfile(settings) #gets all the pre-assigned settings




############# Set up custom stc lexer class ##########################

class BaseLexer(object):
	"""Defines simple interface for custom lexer objects"""
	def __init__(self):
		super(BaseLexer, self).__init__()

	def StyleText(self, event):
		raise NotImplementedError

class DNALexer(BaseLexer):
	"""Lexer to highlight features in a genbank file"""
	# Define some style IDs
	STC_STYLE_DEFAULT = 0
	STC_STYLE_RED_FW = 1
	STC_STYLE_RED_RV = 2
	STC_STYLE_ORANGE_FW = 3
	STC_STYLE_ORANGE_RV = 4
	STC_STYLE_YELLOW_FW = 5
	STC_STYLE_YELLOW_RV = 6
	STC_STYLE_GREEN_FW = 7
	STC_STYLE_GREEN_RV = 8
	STC_STYLE_CYAN_FW = 9
	STC_STYLE_CYAN_RV = 10
	STC_STYLE_BLUE_FW = 11
	STC_STYLE_BLUE_RV = 12
	STC_STYLE_PURPLE_FW = 13
	STC_STYLE_PURPLE_RV = 14
	STC_STYLE_MAGENTA_FW = 15
	STC_STYLE_MAGENTA_RV = 16
	STC_STYLE_BURGUNDY_FW = 17
	STC_STYLE_BURGUNDY_RV = 18
	STC_STYLE_GREY_FW = 19
	STC_STYLE_GREY_RV = 20


	def __init__(self):
		super(DNALexer, self).__init__()




	def StyleText(self, event):
		"""Handle the EVT_STC_STYLENEEDED event"""
		stc = event.GetEventObject()

		featurelist = genbank.gb.get_all_feature_positions()
		for entry in featurelist:

			featuretype, complement, start, finish = entry
			featuretype = featuretype.replace('-', 'a') #for -10 and -35 region
			featuretype = featuretype.replace("5'", "a5") #for 5' features
			featuretype = featuretype.replace("3'", "a3") #for 5' features
			color = eval(featuretype)['color'] #get the color of feature (as string)
			assert type(color) == str

			stc.StartStyling(start, 0x1f)
			if color == 'red' and complement == False: 
				style = DNALexer.STC_STYLE_RED_FW
			elif color == 'red' and complement == True: 
				style = DNALexer.STC_STYLE_RED_RV
			elif color == 'orange' and complement == False:
				style = DNALexer.STC_STYLE_ORANGE_FW
			elif color == 'orange' and complement == True:
				style = DNALexer.STC_STYLE_ORANGE_RV
			elif color == 'yellow' and complement == False:
				style = DNALexer.STC_STYLE_YELLOW_FW
			elif color == 'yellow' and complement == True:
				style = DNALexer.STC_STYLE_YELLOW_RV
			elif color == 'green' and complement == False:
				style = DNALexer.STC_STYLE_GREEN_FW
			elif color == 'green' and complement == True:
				style = DNALexer.STC_STYLE_GREEN_RV
			elif color == 'cyan' and complement == False:
				style = DNALexer.STC_STYLE_CYAN_FW
			elif color == 'cyan' and complement == True:
				style = DNALexer.STC_STYLE_CYAN_RV
			elif color == 'blue' and complement == False:
				style = DNALexer.STC_STYLE_BLUE_FW
			elif color == 'blue' and complement == True:
				style = DNALexer.STC_STYLE_BLUE_RV
			elif color == 'purple' and complement == False:
				style = DNALexer.STC_STYLE_PURPLE_FW
			elif color == 'purple' and complement == True:
				style = DNALexer.STC_STYLE_PURPLE_RV
			elif color == 'magenta' and complement == False:
				style = DNALexer.STC_STYLE_MAGENTA_FW
			elif color == 'magenta' and complement == True:
				style = DNALexer.STC_STYLE_MAGENTA_RV
			elif color == 'burgundy' and complement == False:
				style = DNALexer.STC_STYLE_BURGUNDY_FW
			elif color == 'burgundy' and complement == True:
				style = DNALexer.STC_STYLE_BURGUNDY_RV
			elif color == 'grey' and complement == False:
				style = DNALexer.STC_STYLE_GREY_FW
			elif color == 'grey' and complement == True:
				style = DNALexer.STC_STYLE_GREY_RV
			else: 
				print('nonononononon')
				style = DNALexer.STC_STYLE_DEFAULT

			length = finish-start
			stc.SetStyling(length, style)


class CustomSTC(wx.stc.StyledTextCtrl):
	def __init__(self, *args, **kwargs):
		super(CustomSTC, self).__init__(*args, **kwargs)

		# Attributes
		self.custlex = None

		# Event Handlers
		self.Bind(wx.stc.EVT_STC_STYLENEEDED, self.OnStyle)



#	def SetDNA(self, text):
#		self.SetText('')
#		rows = len(text)//100
#		for i in range(0,rows-1):
#			self.AppendText('%s\n' % text[i*100:(i+1)*100])
#		print('SetDNA')

	def OnStyle(self, event):
		# Delegate to custom lexer object if one exists
		if self.custlex:
			self.custlex.StyleText(event)
		else:
			event.Skip()

	def SetLexer(self, lexerid, lexer=None):
		"""Overrides StyledTextCtrl.SetLexer
		Adds optional param to pass in custom container
		lexer object.
		"""
		self.custlex = lexer
		super(CustomSTC, self).SetLexer(lexerid)


################## Done with lexer class ######################




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
		self.feature_list = features_GUI.FeatureList(splitter1, id=wx.ID_ANY)
		

		#create dna view panel
		self.stc = CustomSTC(splitter1)
		self.stc.SetWrapMode(wx.stc.STC_WRAP_WORD) #enable word wrap
		self.stc.SetLayoutCache(wx.stc.STC_CACHE_DOCUMENT) #cache layout calculations and only draw when something changes

		#set lexer styles
		style = DNALexer.STC_STYLE_DEFAULT
		self.stc.StyleSetSpec(style, "fore:#000000,face:Mono,size:10")
		style = DNALexer.STC_STYLE_RED_FW
		self.stc.StyleSetSpec(style, "fore:#000000,back:%s,face:Mono,size:10" % red['fw'])
		style = DNALexer.STC_STYLE_RED_RV
		self.stc.StyleSetSpec(style, "fore:#000000,back:%s,face:Mono,size:10" % red['rv'])
		style = DNALexer.STC_STYLE_ORANGE_FW
		self.stc.StyleSetSpec(style, "fore:#000000,back:%s,face:Mono,size:10" % orange['fw'])
		style = DNALexer.STC_STYLE_ORANGE_RV
		self.stc.StyleSetSpec(style, "fore:#000000,back:%s,face:Mono,size:10" % orange['rv'])
		style = DNALexer.STC_STYLE_YELLOW_FW
		self.stc.StyleSetSpec(style, "fore:#000000,back:%s,face:Mono,size:10" % yellow['fw'])
		style = DNALexer.STC_STYLE_YELLOW_RV
		self.stc.StyleSetSpec(style, "fore:#000000,back:%s,face:Mono,size:10" % yellow['rv'])
		style = DNALexer.STC_STYLE_GREEN_FW
		self.stc.StyleSetSpec(style, "fore:#000000,back:%s,face:Mono,size:10" % green['fw'])
		style = DNALexer.STC_STYLE_GREEN_RV
		self.stc.StyleSetSpec(style, "fore:#000000,back:%s,face:Mono,size:10" % green['rv'])
		style = DNALexer.STC_STYLE_CYAN_FW
		self.stc.StyleSetSpec(style, "fore:#000000,back:%s,face:Mono,size:10" % cyan['fw'])
		style = DNALexer.STC_STYLE_CYAN_RV
		self.stc.StyleSetSpec(style, "fore:#000000,back:%s,face:Mono,size:10" % cyan['rv'])
		style = DNALexer.STC_STYLE_BLUE_FW
		self.stc.StyleSetSpec(style, "fore:#000000,back:%s,face:Mono,size:10" % blue['fw'])
		style = DNALexer.STC_STYLE_BLUE_RV
		self.stc.StyleSetSpec(style, "fore:#000000,back:%s,face:Mono,size:10" % blue['rv'])
		style = DNALexer.STC_STYLE_PURPLE_FW
		self.stc.StyleSetSpec(style, "fore:#000000,back:%s,face:Mono,size:10" % purple['fw'])
		style = DNALexer.STC_STYLE_PURPLE_RV
		self.stc.StyleSetSpec(style, "fore:#000000,back:%s,face:Mono,size:10" % purple['rv'])
		style = DNALexer.STC_STYLE_MAGENTA_FW
		self.stc.StyleSetSpec(style, "fore:#000000,back:%s,face:Mono,size:10" % magenta['fw'])
		style = DNALexer.STC_STYLE_MAGENTA_RV
		self.stc.StyleSetSpec(style, "fore:#000000,back:%s,face:Mono,size:10" % magenta['rv'])
		style = DNALexer.STC_STYLE_BURGUNDY_FW
		self.stc.StyleSetSpec(style, "fore:#000000,back:%s,face:Mono,size:10" % burgundy['fw'])
		style = DNALexer.STC_STYLE_BURGUNDY_RV
		self.stc.StyleSetSpec(style, "fore:#000000,back:%s,face:Mono,size:10" % burgundy['rv'])
		style = DNALexer.STC_STYLE_GREY_FW
		self.stc.StyleSetSpec(style, "fore:#000000,back:%s,face:Mono,size:10" % grey['fw'])
		style = DNALexer.STC_STYLE_GREY_RV
		self.stc.StyleSetSpec(style, "fore:#000000,back:%s,face:Mono,size:10" % grey['rv'])

		self.stc.SetLexer(wx.stc.STC_LEX_CONTAINER, DNALexer())

		self.stc.Bind(wx.EVT_KEY_DOWN, self.OnKeyPress) #This is important for controlling the input into the editor

		splitter1.SplitHorizontally(self.feature_list, self.stc,sashPosition=-(windowsize[1]-270))

		sizer = wx.BoxSizer(wx.HORIZONTAL)
		sizer.Add(splitter1, -1, wx.EXPAND)
		self.SetSizer(sizer)	
		
		self.Centre()

		#subscribe to messages sent from feature edit module
		pub.subscribe(self.listen_to_updateUI, 'dna_edit_updateUI')

	def listen_to_updateUI(self, msg):
		'''Method for listening updateUI requests'''
		self.updateUI()
		

	def OnKeyPress(self, evt):
		key = evt.GetUniChar()
#		print(key)
		shift = evt.ShiftDown() # is shift down?
#		print(shift)
		if key in [97, 65, 116, 84, 99, 67, 103, 71]: # [a, A, t, T, c, C, g, G']
			start, finish = self.stc.GetSelection()
			if start != finish: # if a selection, delete it, then paste
				genbank.gb.Delete(start+1, finish, visible=False)
			if shift == True:
				genbank.gb.Paste(start+1, chr(key).upper())
			elif shift == False:
				genbank.gb.Paste(start+1, chr(key).lower())
			self.updateUI()
			self.GetTopLevelParent().updateUndoRedo()
			self.stc.SetSelection(start+1, start+1)


		elif key == 8: #backspace
			start, finish = self.stc.GetSelection()
			if start != finish: # if a selection, delete it
				genbank.gb.Delete(start+1, finish)
				self.updateUI()
				self.GetTopLevelParent().updateUndoRedo()
				self.stc.SetSelection(start+1, start+1)
			else:
				genbank.gb.Delete(start, start)
				self.updateUI()
				self.GetTopLevelParent().updateUndoRedo() #for updating the undo and redo buttons in the menu
				self.stc.SetSelection(start-1, start-1)

		elif key == 127: #delete
			start, finish = self.stc.GetSelection()
			if start != finish: # if a selection, delete it
				genbank.gb.Delete(start+1, finish)
			else:
				genbank.gb.Delete(start+1, start+1)
			self.updateUI()
			self.GetTopLevelParent().updateUndoRedo() #for updating the undo and redo buttons in the menu
			self.stc.SetSelection(start, start)

		elif key == 314 and shift == False: #left
			self.stc.CharLeft()

		elif key == 314 and shift == True: #left + select
			self.stc.CharLeftExtend()

		elif key == 316 and shift == False: #right
			self.stc.CharRight()

		elif key == 316 and shift == True: #right + select
			self.stc.CharRightExtend()

		elif key == 315 and shift == False: #up
			self.stc.LineUp()

		elif key == 315 and shift == True: #up + select
			self.stc.LineUpExtend()

		elif key == 317 and shift == False: #down
			self.stc.LineDown()

		elif key == 317 and shift == True: #down + select
			self.stc.LineDownExtend()

##########################
	def make_outputpopup(self):
		'''Creates a popup window in which output can be printed'''
		self.outputframe = wx.Frame(None, title="Output Panel") # creation of a Frame with a title
		self.output = output.create(self.outputframe, style=wx.VSCROLL|wx.HSCROLL) # creation of a richtextctrl in the frame
		



#########################################

#	def dna_output(self, featurelist):
#		'''Prints output to the output panel'''
#		self.make_outputpopup()	
#		tabtext = str(self.stc.GetPageText(self.stc.GetSelection()))
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


	def get_selection(self):
		'''Gets the text editor selection and adjusts it to DNA locations'''
		start, finish = self.stc.GetSelection()
		print(start, finish)
		if start == finish: # if not a selection
			selection = (start+1, finish+1)
		else:
			selection = (start+1, finish)
		return selection

	def uppercase(self):
		'''Change selection to uppercase'''
		start, finish = self.get_selection()
		genbank.gb.Upper(start, finish)
		self.updateUI()
		self.stc.SetSelection(start-1, finish)
	
	def lowercase(self):
		'''Change selection to lowercase'''
		start, finish = self.get_selection()
		genbank.gb.Lower(start, finish)
		self.updateUI() 
		self.stc.SetSelection(start-1, finish)

	def reverse_complement_selection(self):
		'''Reverse-complement current selection'''
		start, finish = self.get_selection()
		genbank.gb.RCselection(start, finish)
		self.updateUI()
		self.stc.SetSelection(start-1, finish)

	def delete(self):
		'''Deletes a selection and updates dna and features'''
		start, finish = self.get_selection()
		genbank.gb.Delete(start, finish)
		self.updateUI()
		self.stc.SetSelection(start-1, start-1)

	def cut(self):
		'''Cut DNA and store it in clipboard together with any features present on that DNA'''
		start, finish = self.get_selection()
		genbank.gb.Cut(start, finish)
		self.updateUI()
		self.stc.SetSelection(start-1, start-1)
			
	def cut_reverse_complement(self):
		'''Cut reverse complement of DNA and store it in clipboard together with any features present on that DNA'''
		start, finish = self.get_selection()
		genbank.gb.CutRC(start, finish)
		self.updateUI()
		self.stc.SetSelection(start-1, start-1)

	def paste(self):
		'''Paste DNA and any features present on that DNA'''
		start, finish = self.get_selection()
		if start != finish: #If a selection, remove sequence
			genbank.gb.Delete(start, finish, visible=False)
		genbank.gb.Paste(start)
		self.updateUI() 
		self.stc.SetSelection(start-1, start-1)
		
	def paste_reverse_complement(self):
		'''Paste reverse complement of DNA and any features present on that DNA'''
		start, finish = self.get_selection()
		if start != finish: #If a selection, remove sequence
			genbank.gb.Delete(start, finish, visible=False)
		genbank.gb.PasteRC(start)
		self.updateUI()
		self.stc.SetSelection(start-1, start-1)

	def copy(self):
		'''Copy DNA and features into clipboard'''
		start, finish = self.get_selection()
		genbank.gb.Copy(start, finish)

	def copy_reverse_complement(self):
		'''Copy reverse complement of DNA'''
		start, finish = self.get_selection()
		genbank.gb.CopyRC(start, finish)


#######################################################


######### Functions for updating text and text highlighting #############
#	def remove_styling(self):
#		'''Removes background styling for whole document'''
#		start = 0
#		finish = self.stc.GetLastPosition() #last position in file
#		self.attr = rt.RichTextAttr()
#		self.attr.SetFlags(wx.TEXT_ATTR_BACKGROUND_COLOUR) #do I need this flag?
#		self.attr.SetBackgroundColour('#FFFFFF')
#		self.stc.SetStyleEx(rt.RichTextRange(start, finish), self.attr)
	
#	def get_feature_lexer(self, featuretype, complement):
#		'''Takes single feature and finds its lexer, if any'''
#		#piece of code to check if there are assigned colors to the features




	
	def updateUI(self):
		'''For changing background color of text ranges'''
#		self.remove_styling() #first remove old styles
		self.stc.SetText(genbank.gb.GetDNA()) #put the DNA in

		#match selection to previous one
#		
#		if start != finish: self.stc.SetSelection(start, finish)
#		elif start == finish: self.stc.SetInsertionPoint(start-1)
#		self.stc.ShowPosition(start) 

	
		#returns a list of lists [[featuretype1, complement1, start1, end1], [featuretype2, complement2, start2, end2].....] 
#		featurelist = genbank.gb.get_all_feature_positions()
#		for entry in featurelist:
#			featuretype, complement, start, finish = entry
#			self.get_feature_color(featuretype, complement)

#			color = self.stc.current_highlight_color #get color

#			#do the painting
#			self.attr = rt.RichTextAttr()
#			self.attr.SetFlags(wx.TEXT_ATTR_BACKGROUND_COLOUR) #do I need this flag?
#			self.attr.SetBackgroundColour(color)
#			self.stc.SetStyleEx(rt.RichTextRange(start, finish), self.attr)

	
#		#color in search hits if any are present
#		search_hits = genbank.gb.search_hits
#		if len(search_hits) != 0:
#			color = '#ffff00'
#			for i in range(len(search_hits)):
#				start = int(search_hits[i][0])
#				finish = int(search_hits[i][1])
#				#do the painting
#				self.attr = rt.RichTextAttr()
#				self.attr.SetFlags(wx.TEXT_ATTR_BACKGROUND_COLOUR) #do I need this flag?
#				self.attr.SetBackgroundColour(color)
#				self.stc.SetStyleEx(rt.RichTextRange(start, finish), self.attr)

#I need to think carefully about where to place these...
		#realize the current selection in the DNA editor
#		
#		self.stc.SetSelection(start, finish) #update the graphical selection
#		self.stc.ShowPosition(start) 

		#update feature list
		self.feature_list.updateUI()

		#use this to get the first line characters 
		##develop it##
		#set insertion point to beginning... Make it update on resize
#		while self.stc.MoveDown() == True:
#			self.stc.MoveDown()
#			self.linecount.write(str(self.stc.GetInsertionPoint()+1)+'\n', 'Text')

		##############



####### Advanced DNA functions #######
#	def align_dna(self, evt):
#		'''Pass a list of dictionaries with name and dna to m_align to align sequences'''
#		#this works sort of ok. Needs improvement though...
#		
#		dna = str(genbank.gb.GetDNA())
#		
#		self.seqlist.append(dict(name='Reference', dna=dna))
#		print(self.seqlist)
#		result = seqfiles.M_align(self.seqlist)
#		self.output.write(result, 'text')

		
####### Protein functions #######
	def translate_output(self, protein, DNA, info):
		'''Generate output in the output.panel'''
#		tabtext = str(self.stc.GetPageText(self.stc.GetSelection()))
		self.make_outputpopup()			
		self.output.write('%s | Translate %s\n' % (tabtext, info), 'File')
		self.output.write(('%d AA from %d bases, %d bases left untranslated' % (len(protein), len(DNA), len(DNA)%3))+'\n', 'Text')
		self.output.write(protein, 'Protein')
		self.outputframe.Show()
	
	def translate_selection(self):
		'''Translate selected DNA'''
		DNA = self.stc.GetStringSelection()
		protein = dna.Translate(DNA)
		self.translate_output(protein, DNA, 'leading strand')
		
	def translate_selection_reverse_complement(self):
		'''Translate reverse-complement of selected DNA'''
		DNA = self.stc.GetStringSelection()
		protein = dna.Translate(dna.RC(DNA))
		self.translate_output(protein, DNA, 'complement strand')

#update this one...
	def translate_feature(self):
		'''Translate specified feature'''
		feature = genbank.gb.allgbfeatures[2]
		DNA = genbank.gb.getdnaforgbfeature(feature[4])
		protein = dna.Translate(DNA)
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

		
#	def mouse_position(self, event):
#		'''Get which features are at a given position'''		
#		xposition, yposition = self.stc.ScreenToClient(wx.GetMousePosition())
#		if xposition > 0 and yposition > 0:

#			point = wx.Point(xposition, yposition)
#			#event.GetPosition() #this can be used if it controlled by an event..
#			mposition = self.stc.HitTest((1,1))

#			print(mposition)
#			print(type(mposition))
#			#which feature corresponds to this pos?
#			Feature = genbank.gb.get_featurename_for_pos(mposition)
#			return mposition, Feature
#		else:
#			return None, None

###############################################################
if __name__ == '__main__': #if script is run by itself and not loaded	
	app = wx.App() # creation of the wx.App object (initialisation of the wxpython toolkit)

	frame = wx.Frame(None, title="DNA editor") # creation of a Frame with a title
	frame.dnaeditor_GUI = MyPanel(frame) # creation of a richtextctrl in the frame
	frame.Show() # frames are invisible by default so we use Show() to make them visible
	frame.dnaeditor_GUI.open_file("")
	app.MainLoop() # here the app enters a loop waiting for user input


