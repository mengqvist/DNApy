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



from wx.lib.pubsub import setupkwargs #this line not required in wxPython2.9.
 	                                  #See documentation for more detail
from wx.lib.pubsub import pub

from copy import deepcopy
import string
from os import access,listdir


import sys, os
import string

import genbank
import wx.stc

import pyperclip
import copy

import output
import dna
from base_class import DNApyBaseClass
import featurelist_GUI
import plasmid_GUI

#TODO
#fix statusbar in general
#fix save option
#fix undo/redo on windows




files={}   #list with all configuration files
files['default_dir'] = os.path.abspath(os.path.dirname(sys.argv[0]))+"/"
files['default_dir']=string.replace(files['default_dir'], "\\", "/")
files['default_dir']=string.replace(files['default_dir'], "library.zip", "")
settings=files['default_dir']+"settings"   ##path to the file of the global settings
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
	STC_STYLE_SEARCH_HITS = 21
	STC_STYLE_ENZYMES = 22
	STC_STYLE_SELECTION = 23
	

	def __init__(self):
		super(DNALexer, self).__init__()




	def StyleText(self, event):
		"""Handle the EVT_STC_STYLENEEDED event"""
		stc = event.GetEventObject()

		featurelist = genbank.gb.get_all_feature_positions()
		if featurelist == None: #no features
			pass
		else:
			for entry in featurelist:

				featuretype, complement, start, finish, name, index = entry
				featuretype = featuretype.replace('-', 'a') #for -10 and -35 region
				featuretype = featuretype.replace("5'", "a5") #for 5' features
				featuretype = featuretype.replace("3'", "a3") #for 5' features
				color = eval(featuretype)['color'] #get the color of feature (as string)
				assert type(color) == str

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
					style = DNALexer.STC_STYLE_DEFAULT
				start -= 1
				stc.StartStyling(start, 0x1f)
				length = finish-start
				if length == 0:
					stc.StartStyling(start-1, 0x1f)
					length = 1
				stc.SetStyling(length, style)

		#color restriction sites
		for enzyme in genbank.restriction_sites:

			sites = genbank.restriction_sites[enzyme].restrictionSites

			if len(sites) > 0:
				for site in sites:
					start  = site[1]
					finish = site[2]
					style = DNALexer.STC_STYLE_ENZYMES
					if start < finish:  # enzyme cuts in the strain
						stc.StartStyling(start-1, 0x1f)
						length = finish-(start-1)
						if length == 0:
							stc.StartStyling(start-1, 0x1f)
							length = 1
						stc.SetStyling(length, style)
					else: # enzyme is exactly throughout the 0
						# the beginning:
						stc.StartStyling(0, 0x1f)
						length = finish
						stc.SetStyling(length, style)

						# and the end:
						stc.StartStyling(start-1, 0x1f)
						length = len(genbank.gb.gbfile['dna']) - start + 1
						stc.SetStyling(length, style)


		#color search hits
		if genbank.search_hits != None:
			for hit in genbank.search_hits:
				##add logic here for hits that are broken up
				start, finish = hit
				style = DNALexer.STC_STYLE_SEARCH_HITS
				stc.StartStyling(start-1, 0x1f)
				length = finish-(start-1)
				if length == 0:
					stc.StartStyling(start-1, 0x1f)
					length = 1
				stc.SetStyling(length, style)



#		stc.IndicatorSetStyle(2, wx.stc.STC_INDIC_PLAIN)
#		stc.IndicatorSetForeground(0, wx.RED)
#		stc.StartStyling(1, wx.stc.STC_INDICS_MASK)
#		stc.SetStyling(10, wx.stc.STC_INDIC2_MASK)

		#PAINT SELECTION
#		start, finish = genbank.dna_selection
#		print('start, finish', start, finish)
#		style = DNALexer.STC_STYLE_SELECTION
#		stc.StartStyling(start-1, 0x1f)
#		length = finish-(start-1)
#		if length == 0:
#			stc.StartStyling(start-1, 0x1f)
#			length = 1
#		stc.SetStyling(length, style)

class CustomSTC(wx.stc.StyledTextCtrl):
	def __init__(self, *args, **kwargs):
		super(CustomSTC, self).__init__(*args, **kwargs)

		# Attributes
		self.custlex = None

		# Event Handlers
		self.Bind(wx.stc.EVT_STC_STYLENEEDED, self.OnStyle)


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
class TextEdit(DNApyBaseClass):
	def __init__(self, parent, id):
		wx.Panel.__init__(self, parent)
		

		
#		#create line count panel
##		self.linecount = output.create(self, style=wx.VSCROLL|wx.HSCROLL|wx.BORDER_NONE); #create DNA window
##		self.linecount.SetEditable(False)

		#create dna view panel
		self.stc = CustomSTC(self)
		self.stc.SetWrapMode(wx.stc.STC_WRAP_CHAR) #enable word wrap
		self.stc.SetLayoutCache(wx.stc.STC_CACHE_DOCUMENT) #cache layout calculations and only draw when something changes


		self.stc.Bind(wx.EVT_LEFT_DOWN, self.OnLeftDown)
		self.stc.Bind(wx.EVT_LEFT_UP, self.OnLeftUp)
		self.stc.Bind(wx.EVT_RIGHT_UP, self.OnRightUp)
		self.stc.Bind(wx.EVT_MOTION, self.OnMotion)

		#set lexer styles
		style = DNALexer.STC_STYLE_DEFAULT
		self.stc.StyleSetSpec(style, "fore:#000000,back:#FFFFFF,face:Mono,size:10")
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
		style = DNALexer.STC_STYLE_SEARCH_HITS
		self.stc.StyleSetSpec(style, "fore:#000000,back:#CCFFOO,face:Mono,bold,size:10")
		style = DNALexer.STC_STYLE_SELECTION
		self.stc.StyleSetSpec(style, "fore:#F6F036,back:#CCFF33,face:Mono,bold,size:10")

		style = DNALexer.STC_STYLE_ENZYMES
		self.stc.StyleSetSpec(style, "fore:#000000,back:#FFOO33,face:Mono,size:10")


		self.stc.SetLexer(wx.stc.STC_LEX_CONTAINER, DNALexer())

		self.stc.Bind(wx.EVT_KEY_DOWN, self.OnKeyPress) #This is important for controlling the input into the editor

		sizer = wx.BoxSizer(wx.HORIZONTAL)
		sizer.Add(item=self.stc, proportion=-1, flag=wx.EXPAND)
		self.SetSizer(sizer)	
		
		self.Centre()


####### Modify methods from base calss to fit current needs #########

	def update_globalUI(self):
		'''
		Method should be modified as to update other panels in response to changes in own panel.
		'''
		MSG_CHANGE_TEXT = "change.text"
		pub.sendMessage(MSG_CHANGE_TEXT, text="DNA view says update!")

	
	def update_ownUI(self):
		'''
		For changing background color of text ranges and updating selection.
		'''
		sequence = genbank.gb.GetDNA()
		if sequence == None:
			sequence = ''	
		self.stc.SetText(sequence) #put the DNA in the editor
		start, finish = genbank.dna_selection
		if finish == -1: #a caret insertion (and no selection). Will actually result in a selection where start is one larger than finish
			self.stc.GotoPos(start-1) #set caret insertion
		else:
			self.stc.SetSelection(start-1, finish) #update selection


######################################################

	def OnLeftDown(self, event):
		event.Skip() #very important to make the event propagate and fulfill its original function

	def OnLeftUp(self, event):
		self.set_dna_selection() #update the varable keeping track of DNA selection
		event.Skip() #very important to make the event propagate and fulfill its original function		

	def OnMotion(self, event):
		#if event.Dragging() and event.LeftIsDown():
		#	self.set_dna_selection() #update the varable keeping track of DNA selection
		event.Skip() #very important to make the event propagate and fulfill its original function

	def OnRightUp(self, event):
		print('dna right up')
		event.Skip() #very important to make the event propagate and fulfill its original function

#####################################################################

	def set_dna_selection(self):
		'''
		Updates DNA selection.
		'''
		selection = self.get_selection()
		genbank.dna_selection = selection
		self.update_globalUI()
		#self.update_ownUI()


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
			self.update_ownUI()
			self.update_globalUI()
			self.stc.SetSelection(start+1, start+1)


		elif key == 8: #backspace
			start, finish = self.stc.GetSelection()
			if start != finish: # if a selection, delete it
				genbank.gb.Delete(start+1, finish)
				self.update_ownUI()
				self.update_globalUI()
				self.stc.SetSelection(start+1, start+1)
			else:
				genbank.gb.Delete(start, start)
				self.update_ownUI()
				self.update_globalUI()
				self.stc.SetSelection(start-1, start-1)

		elif key == 127: #delete
			start, finish = self.stc.GetSelection()
			if start != finish: # if a selection, delete it
				genbank.gb.Delete(start+1, finish)
			else:
				genbank.gb.Delete(start+1, start+1)
			self.update_ownUI()
			self.update_globalUI()
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

		self.set_dna_selection() #update the varable keeping track of DNA selection

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
	def select_all(self):
		'''Select the entire dna sequence'''
		start = 0
		finish = len(genbank.gb.GetDNA())
		self.stc.SetSelection(start, finish)
		self.set_dna_selection()	
		self.update_ownUI()
		self.update_globalUI()

	def get_selection(self):
		'''Gets the text editor selection and adjusts it to DNA locations.'''
		start, finish = self.stc.GetSelection()
		if start == finish: #not a selection
			finish = -1
		elif start > finish: #the selection was made backwards
			start, finish = finish, start
		
		selection = (start+1, finish)
#		print('selection', selection)
		return selection

	def uppercase(self):
		'''Change selection to uppercase'''
		start, finish = self.get_selection()
		if finish == -1:
			raise ValueError, 'Cannot modify an empty selection'
		else:
			genbank.gb.Upper(start, finish)
			self.update_ownUI()
			self.stc.SetSelection(start-1, finish)
	
	def lowercase(self):
		'''Change selection to lowercase'''
		start, finish = self.get_selection()
		if finish == -1:
			raise ValueError, 'Cannot modify an empty selection'
		else:
			genbank.gb.Lower(start, finish)
			self.update_ownUI() 
			self.stc.SetSelection(start-1, finish)

	def reverse_complement_selection(self):
		'''Reverse-complement current selection'''
		start, finish = self.get_selection()
		if finish == -1:
			raise ValueError, 'Cannot modify an empty selection'
		else:
			genbank.gb.RCselection(start, finish)
			self.update_ownUI()
			self.update_globalUI()
			self.stc.SetSelection(start-1, finish)

	def delete(self):
		'''Deletes a selection and updates dna and features'''
		start, finish = self.get_selection()
		if finish == -1:
			raise ValueError, 'Cannot delete an empty selection'
		else:
			genbank.gb.Delete(start, finish)
			self.update_ownUI()
			self.update_globalUI()
			self.stc.SetSelection(start-1, start-1)

	def cut(self):
		'''Cut DNA and store it in clipboard together with any features present on that DNA'''
		start, finish = self.get_selection()
		if finish == -1:
			raise ValueError, 'Cannot cut an empty selection'
		else:
			genbank.gb.Cut(start, finish)
			self.update_ownUI()
			self.update_globalUI()
			self.stc.SetSelection(start-1, start-1)
			
	def cut_reverse_complement(self):
		'''Cut reverse complement of DNA and store it in clipboard together with any features present on that DNA'''
		start, finish = self.get_selection()
		if finish == -1:
			raise ValueError, 'Cannot cut an empty selection'
		else:
			genbank.gb.CutRC(start, finish)
			self.update_ownUI()
			self.update_globalUI()
			self.stc.SetSelection(start-1, start-1)

	def paste(self):
		'''Paste DNA and any features present on that DNA'''
		start, finish = self.get_selection()
		if finish == -1:
			pass
		else: #If a selection, remove sequence
			genbank.gb.Delete(start, finish, visible=False)
		#genbank.gb.Paste(start)
		genbank.gb.RichPaste(start)
		self.update_ownUI() 
		self.update_globalUI()
		self.stc.SetSelection(start-1, start-1)
		
	def paste_reverse_complement(self):
		'''Paste reverse complement of DNA and any features present on that DNA'''
		start, finish = self.get_selection()
		if finish == -1:
			pass
		else: #If a selection, remove sequence
			genbank.gb.Delete(start, finish, visible=False)
		genbank.gb.PasteRC(start)
		self.update_ownUI()
		self.update_globalUI()
		self.stc.SetSelection(start-1, start-1)

	def copy(self):
		'''Copy DNA and features into clipboard'''
		start, finish = self.get_selection()
		if finish == -1:
			raise ValueError, 'Cannot copy an empty selection'
		else:
			#genbank.gb.Copy(start, finish)
		
			# try the new copy:
			genbank.gb.RichCopy(start, finish)

	def copy_reverse_complement(self):
		'''Copy reverse complement of DNA'''
		start, finish = self.get_selection()
		if finish == -1:
			raise ValueError, 'Cannot copy an empty selection'
		else:
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
		self.output.write('Translate %s\n' % (info), 'File')
		self.output.write(('%d AA from %d bases, %d bases left untranslated' % (len(protein), len(DNA), len(DNA)%3))+'\n', 'Text')
		self.output.write(protein, 'Protein')
		self.outputframe.Show()
	
	def translate_selection(self):
		'''Translate selected DNA'''
		start, finish = self.get_selection()
		if finish == -1:
			raise ValueError, 'Cannot translate an empty selection'
		else:
			DNA = genbank.gb.GetDNA(start, finish)
			protein = dna.Translate(DNA)
			self.translate_output(protein, DNA, 'leading strand')
		
	def translate_selection_reverse_complement(self):
		'''Translate reverse-complement of selected DNA'''
		start, finish = self.get_selection()
		if finish == -1:
			raise ValueError, 'Cannot translate an empty selection'
		else:
			DNA = genbank.gb.GetDNA(start, finish)
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
		self.update_ownUI() #refresh everything

		
	def mouse_position(self, event):
		'''Get which features are at a given position'''		
		xposition, yposition = self.stc.ScreenToClient(wx.GetMousePosition())
		if xposition > 1 and yposition > 1:
			print('x,y', xposition, yposition)
			mposition = self.stc.CharPositionFromPoint(xposition, yposition)
			print(mposition)

#			#which feature corresponds to this pos?
			Feature = genbank.gb.get_featurename_for_pos(mposition)
			return mposition, Feature
		else:
			return None, None




			
			
######################################
######################################

		
##### main loop
class MyApp(wx.App):
	def OnInit(self):
		frame = wx.Frame(None, -1, title="DNA Edit", size=(700,600))
		panel = DNAedit(frame, -1)
		frame.Centre()
		frame.Show(True)
		self.SetTopWindow(frame)
		return True



if __name__ == '__main__': #if script is run by itself and not loaded	

	files={}   #list with all configuration files
	files['default_dir'] = os.path.abspath(os.path.dirname(sys.argv[0]))+"/"
	files['default_dir']=string.replace(files['default_dir'], "\\", "/")
	files['default_dir']=string.replace(files['default_dir'], "library.zip", "")
	settings=files['default_dir']+"settings"   ##path to the file of the global settings
	execfile(settings) #gets all the pre-assigned settings

	genbank.dna_selection = (1, 1)	 #variable for storing current DNA selection
	genbank.feature_selection = False #variable for storing current feature selection
	genbank.search_hits = []
	
	import sys
	assert len(sys.argv) == 2, 'Error, this script requires a path to a genbank file as an argument.'
	print('Opening %s' % str(sys.argv[1]))

	genbank.gb = genbank.gbobject(str(sys.argv[1])) #make a genbank object and read file


	app = MyApp(0)
	app.MainLoop()
