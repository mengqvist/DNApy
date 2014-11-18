#!/usr/bin/env python


#This file is part of DNApy. DNApy is a DNA editor written purely in python. 
#The program is intended to be an intuitive, fully featured, 
#extendible, editor for molecular and synthetic biology.  
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



import wx
#import wx.lib.graphics
import math
from base_class import DNApyBaseDrawingClass
import mixed_base_codons as mbc
import dna

#class Button(DNApyBaseDrawingClass):
#	def __init__(self, parent, id, xpos, ypos):
#		wx.Panel.__init__(self, parent, id, size=(30, 30), style=wx.SUNKEN_BORDER)
#		self.parent = parent
#		self.xpos = xpos
#		self.ypos = ypos
#		self.font = wx.Font(9, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL,
#			wx.FONTWEIGHT_NORMAL, False, 'Courier 10 Pitch')

#	def Draw(self, dc):
#		dc.SetBackground(wx.Brush("Black"))
#		dc.Clear() # make sure you clear the bitmap!
#		self.gcdc = wx.GCDC(dc)
#		self.gcdc.SetPen(wx.Pen(colour='#666666', width=0))
	



class Button:
	def __init__(self, parent, id=-1, size=(10,10), label='x'):
		'''A custom class that draws on the parent GCDC.
		The idea is to add some of the functionality and abstraction of a regular button.
		This approach seems easer than adding a bona fide button on top of a gcdc.'''

		self.parent = parent
		self.id = id
		self.size=size
		self.label=label

		self.no_hl_pen = wx.Pen(colour='#8B835F', width=2)
		self.no_hl_brush = wx.Brush("#fff2d1")
		self.no_hl_text = '#8B835F'

		self.hl_pen = wx.Pen(colour='#4B4424', width=2)
		self.hl_brush = wx.Brush("#ffd976")
		self.hl_text = '#4B4424'

#		self.font = wx.Font(pointSize=self.size[1]/2.5, family=wx.FONTFAMILY_SWISS, style=wx.FONTWEIGHT_BOLD, weight=wx.FONTWEIGHT_BOLD)

		#add 'reset' button to index
		self.parent.catalog[str(self.parent.NextRGB()+(255,))] = 'reset'

	def update_ownUI(self):
		self.Draw()


##############################################

	def Draw(self):
		dc = wx.MemoryDC()
		self.bit1 = wx.EmptyBitmap(self.size[0], self.size[1])
		dc.SelectObject(self.bit1)
		dc.SetBackground(wx.Brush("White"))
		dc.Clear() # make sure you clear the bitmap!
		self.gcdc = wx.GCDC(dc) #make gcdc from the dc (for use of transparancy and antialiasing)
#		self.gcdc.SetFont(self.font)


		#make a hidden dc to which features can be drawn in uinique colors and later used for hittests
		self.hidden_dc = wx.MemoryDC()
		self.bit2 = wx.EmptyBitmap(self.size[0], self.size[1])
		self.hidden_dc.SelectObject(self.bit2)
		self.hidden_dc.SetBackground(wx.Brush("White"))
		self.hidden_dc.Clear() # make sure you clear the bitmap!

		if self.parent.highlighted == 'reset':
			self.gcdc.SetPen(self.hl_pen)
			self.gcdc.SetBrush(self.hl_brush)
			self.gcdc.SetTextForeground(self.hl_text)

		elif self.parent.highlighted != 'reset':
			self.gcdc.SetPen(self.no_hl_pen)
			self.gcdc.SetBrush(self.no_hl_brush)
			self.gcdc.SetTextForeground(self.no_hl_text)

		#place text in button. Make button bigger if text does not fit.
		text_size=self.gcdc.GetTextExtent(self.label)
		if text_size[0]>self.size[0]: 
			self.size[0] = text_size

		#draw button
		self.gcdc.DrawRoundedRectangle(0, 0, self.size[0], self.size[1], self.size[0]/5)

		#draw text
		btn_mid_x = self.size[0]/2 #middle of button in x
		btn_mid_y = self.size[1]/2 #middle of button in y
		self.gcdc.DrawText(self.label, btn_mid_x-text_size[0]/2, btn_mid_y-text_size[1]/2)


		#draw the hidden rectangle that is used for hittests
		self.hidden_dc.SetPen(wx.Pen(colour=self.parent.unique_color, width=0))
		self.hidden_dc.SetBrush(wx.Brush(colour=self.parent.unique_color))
		self.hidden_dc.DrawRoundedRectangle(0, 0, self.size[0], self.size[1], self.size[0]/5)

		dc.SelectObject(wx.NullBitmap) # need to get rid of the MemoryDC before Update() is called.
		self.hidden_dc.SelectObject(wx.NullBitmap) # need to get rid of the MemoryDC before Update() is called.



	def SetSize(self, size):
		self.size = size

	def GetSize(self, size):
		return self.size

	def SetLabel(self, label):
		self.label = label

	def GetLabel(self):
		return self.label

	def GetBitmaps(self):
		self.update_ownUI()
		return self.bit1, self.bit2


class TextCtrl:
	def __init__(self, parent, id=-1, size=(100,50), pointsize=12, label='', insertionpoint=False):
		'''A custom class that draws on the parent GCDC.
		The idea is to add some of the functionality and abstraction of a regular TextCtrl field.
		This approach seems easer than adding a bona fide TextCtrl on top of a gcdc.'''

		self.parent = parent
		self.id = id
		self.size=size #size as (x,y)
		self.pointsize=pointsize #size of text
		self.label=label #text to print
		self.insertionpoint = insertionpoint #where in the text the caret is

		self.underline_no_hl_pen = wx.Pen(colour='#FFFFFF', width=0) 
		self.no_hl_pen = wx.Pen(colour='#FFFFFF', width=0)
		self.no_hl_brush = wx.Brush("#FFFFFF")
		self.no_hl_text = '#8B835F'

		self.underline_hl_pen = wx.Pen(colour='#C5C4F5', width=3) 
		self.hl_pen = wx.Pen(colour='#C5C4F5', width=0)
		self.hl_brush = wx.Brush("#FFFFFF")
		self.hl_text = '#4B4424'

		self.caret_pen = wx.Pen(colour='#666666', width=2) 

#		self.font = wx.Font(pointSize=self.pointsize, family=wx.FONTFAMILY_SWISS, style=wx.FONTWEIGHT_BOLD, weight=wx.FONTWEIGHT_BOLD)


		#add 'reset' button to index
		self.parent.catalog[str(self.parent.NextRGB()+(255,))] = 'text'


	def update_ownUI(self):
		self.Draw()
		
	def Draw(self):
		dc = wx.MemoryDC()
		self.bit1 = wx.EmptyBitmap(self.size[0], self.size[1])
		dc.SelectObject(self.bit1)
		dc.SetBackground(wx.Brush("White"))
		dc.Clear() # make sure you clear the bitmap!
		self.gcdc = wx.GCDC(dc) #make gcdc from the dc (for use of transparancy and antialiasing)
#		self.gcdc.SetFont(self.font)


		#make a hidden dc to which features can be drawn in uinique colors and later used for hittests
		self.hidden_dc = wx.MemoryDC()
		self.bit2 = wx.EmptyBitmap(self.size[0], self.size[1])
		self.hidden_dc.SelectObject(self.bit2)
		self.hidden_dc.SetBackground(wx.Brush("White"))
		self.hidden_dc.Clear() # make sure you clear the bitmap!


		if self.parent.highlighted == 'text' or self.parent.text_edit_active is True:
			self.gcdc.SetPen(self.hl_pen)
			self.gcdc.SetBrush(self.hl_brush)
			self.gcdc.SetTextForeground(self.hl_text)

		elif self.parent.highlighted != 'text':
			self.gcdc.SetPen(self.no_hl_pen)
			self.gcdc.SetBrush(self.no_hl_brush)
			self.gcdc.SetTextForeground(self.no_hl_text)

		#place text in field. Make field x bigger if text does not fit.
		text_size=self.gcdc.GetTextExtent(self.label)
		if text_size[0]>self.size[0]: 
			self.size = (text_size[0], self.size[1])

		#make y size of field match the text size
		self.size = (self.size[0], text_size[1])

		#draw text field (basically a rectangle)
		self.gcdc.DrawRectangle(0, 0, self.size[0], self.size[1])

		#draw the underline
		if self.parent.highlighted == 'text' or self.parent.text_edit_active is True:
			self.gcdc.SetPen(self.underline_hl_pen)
		elif self.parent.highlighted != 'text':
			self.gcdc.SetPen(self.underline_no_hl_pen)

		self.gcdc.DrawLine(x1=1, y1=self.size[1]-self.size[1]/10, x2=1, y2=self.size[1])
		self.gcdc.DrawLine(x1=0, y1=self.size[1], x2=self.size[0], y2=self.size[1])
		self.gcdc.DrawLine(x1=self.size[0]-1, y1=self.size[1]-self.size[1]/10, x2=self.size[0]-1, y2=self.size[1])
		

		#draw the hidden rectangle that is used for hittests
		self.hidden_dc.SetPen(wx.Pen(colour=self.parent.unique_color, width=0))
		self.hidden_dc.SetBrush(wx.Brush(colour=self.parent.unique_color))
		self.hidden_dc.DrawRectangle(0, 0, self.size[0], self.size[1])



		#draw text
		btn_mid_x = self.size[0]/2 #middle of button in x
		btn_mid_y = self.size[1]/2 #middle of button in y
		self.gcdc.DrawText(self.label, btn_mid_x-text_size[0]/2, btn_mid_y-text_size[1]/2)

		#draw caret
		self.gcdc.SetPen(self.caret_pen)
		if self.insertionpoint is False:	
			pass	
		elif self.insertionpoint == 0:
			self.gcdc.DrawLine(x1=1, y1=0, x2=1, y2=0+self.size[1])
		elif self.insertionpoint == 1:
			self.gcdc.DrawLine(x1=self.gcdc.GetTextExtent(self.label[1])[0], y1=0, x2=self.gcdc.GetTextExtent(self.label[1])[0], y2=0+self.size[1])
		elif self.insertionpoint == 2:
			self.gcdc.DrawLine(x1=self.gcdc.GetTextExtent(self.label[0:2])[0], y1=0, x2=self.gcdc.GetTextExtent(self.label[0:2])[0], y2=0+self.size[1])
		elif self.insertionpoint == 3:
			self.gcdc.DrawLine(x1=self.size[0]-1, y1=0, x2=self.size[0]-1, y2=0+self.size[1])

#	def XY_to_insertionpoint(self, coordinate):
#		x, y = coordinate

#		return in

	def SendXY(self, coordinate):
		x, y = coordinate
		if x <= self.gcdc.GetTextExtent(self.label[0:1])[0]/2:
			self.insertionpoint = 0
		elif self.gcdc.GetTextExtent(self.label[0:1])[0]/2 < x <= self.gcdc.GetTextExtent(self.label[0:2])[0]-self.gcdc.GetTextExtent(self.label[1:2])[0]/2:
			self.insertionpoint = 1
		elif self.gcdc.GetTextExtent(self.label[0:2])[0]-self.gcdc.GetTextExtent(self.label[1:2])[0]/2 < x <= self.gcdc.GetTextExtent(self.label[0:3])[0]-self.gcdc.GetTextExtent(self.label[2:3])[0]/2:
			self.insertionpoint = 2
		elif self.gcdc.GetTextExtent(self.label[0:3])[0]-self.gcdc.GetTextExtent(self.label[2:3])[0]/2 < x:
			self.insertionpoint = 3



	def SetSize(self, size):
		self.size = size

	def GetSize(self, size):
		return self.size

	def SetText(self, label):
		self.label = label

	def GetText(self):
		return self.label

	def GetBitmaps(self):
		self.update_ownUI()
		return self.bit1, self.bit2




class CodonView(DNApyBaseDrawingClass):
	def __init__(self, parent, id):

		self.highlighted = False #a variable for keeping track of whether any object is highlighted
		self.codon = False
		self.target = []
		self.possible = []
		self.offtarget = []
		self.text_edit_active = False #to keep track of whether text is being edited

		#set up a dictionary to keep track of which color belongs to what object
		self.catalog = {} #for matching features with the unique colors
		self.catalog['(255, 255, 255, 255)'] = False #the background is white, have to add that key
		self.unique_color = (0,0,0)

		self.xc = 0
		self.yc = 0

		#create text field for showing/editing the ambigous codon
		self.text_field = TextCtrl(self, size=(80,60), label='')

		#create reset button
		self.reset_btn = Button(self, size=(80,30), label='Reset')

		#initialize
		super(CodonView, self).__init__(parent, wx.ID_ANY)



#		self.Bind(wx.EVT_LEFT_DOWN, self.OnLeftDown)
		self.Bind(wx.EVT_LEFT_UP, self.OnLeftUp)
#		self.Bind(wx.EVT_RIGHT_UP, self.OnRightUp)
		self.Bind(wx.EVT_MOTION, self.OnMotion)
#		self.Bind(wx.EVT_LEFT_DCLICK, self.OnLeftDouble)






############ Setting required methods ####################

	def update_globalUI(self):
		'''Method should be modified as to update other panels in response to changes in own panel.
		Preferred use is through sending a message using the pub module.
		Example use is: pub.Publisher.sendMessage('feature_list_updateUI', '').
		The first string is the "listening group" and deterimines which listeners get the message. 
		The second string is the message and is unimportant for this implementation.
		The listening group assigned here (to identify recipients) must be different from the listening group assigned in __init__ (to subscribe to messages).'''
#		pub.Publisher.sendMessage('from_plasmid_view', '')
		pass

	
	def update_ownUI(self):
		"""
		This would get called if the drawing needed to change, for whatever reason.

		The idea here is that the drawing is based on some data generated
		elsewhere in the system. If that data changes, the drawing needs to
		be updated.

		This code re-draws the buffer, then calls Update, which forces a paint event.
		"""
		dc = wx.MemoryDC()
		dc.SelectObject(self._Buffer)
		self.Draw(dc)
		dc.SelectObject(wx.NullBitmap) # need to get rid of the MemoryDC before Update() is called.
		self.Refresh()
		self.Update()



	def set_dna_selection(self, selection):
		'''Recieves requests for DNA selection and then sends it.'''
		assert type(selection) == tuple, 'Error, dna selection must be a tuple'
		selection = (int(selection[0]), int(selection[1]))
		genbank.dna_selection = selection


############### Done setting required methods #######################




	def Draw(self, dc):
		'''Class for drawing stuff on self.gcdc.
		It is important to only use the base_class drawing tools and not the inbuilt gcdc tools, 
		otherwise the hittests will not work'''
		self.xc = self.size[0]/3 #centre of codon circle in x
		self.yc = self.size[1]/2 #centre of codon circle in y
		self.Radius = self.yc/1.2
		self.unique_color = (0,0,0)

		dc.SetBackground(wx.Brush("White"))
		dc.Clear() # make sure you clear the bitmap!
		self.gcdc = wx.GCDC(dc) #make gcdc from the dc (for use of transparancy and antialiasing)

		#make a hidden dc to which features can be drawn in uinique colors and later used for hittests
		self.hidden_dc = wx.MemoryDC()
		self.hidden_dc.SelectObject(wx.EmptyBitmap(self.ClientSize[0], self.ClientSize[1]))
		self.hidden_dc.SetBackground(wx.Brush("White"))
		self.hidden_dc.Clear() # make sure you clear the bitmap!


		target_color = '#CCFF66'
		possible_color = '#FFFF66'
		offtarget_color = '#FF9966'
		nucleotide_color = '#8B835F'
		coding_nucleotide_color = '#4B4424' #for coloring the nucleotides encoded by the ambigous codon
		line_color = '#8B835F' #for lines
		

		first_nucleotide_thickness = self.Radius/3
		second_nucleotide_thickness = 2*(self.Radius/3)/3
		third_nucleotide_thickness = 1*(self.Radius/3)/3
		amino_acid_thickness = self.Radius/2.25



		
		#draw first nucleotide
		radius = first_nucleotide_thickness
		thickness = first_nucleotide_thickness
#		font = wx.Font(pointSize=thickness/1.5, family=wx.FONTFAMILY_SWISS, style=wx.FONTWEIGHT_BOLD, weight=wx.FONTWEIGHT_BOLD)
#		self.gcdc.SetFont(font)
		self.gcdc.SetPen(wx.Pen(colour='#666666', width=0))
		self.gcdc.SetBrush(wx.Brush("#ffe7ab"))
		nucleotides = ['U', 'C', 'A', 'G']
		for i in range(len(nucleotides)):
			start_angle = 0 + 90*i
			finish_angle = 90+90*i
			pointlist = self.make_arc(self.xc, self.yc, start_angle, finish_angle, radius, thickness, step=5)
			self.gcdc.DrawPolygon(pointlist)
			x1, y1 = self.AngleToPoints(self.xc, self.yc, radius/2, finish_angle-(finish_angle-start_angle)/2)

			#if nucleotide is part of degenerate codon it should have a different color
			self.gcdc.SetTextForeground((nucleotide_color))
			if self.codon is not False:
				if nucleotides[i].replace('U','T') in dna.UnAmb(self.codon[0]):
					self.gcdc.SetTextForeground((coding_nucleotide_color))
				
			self.gcdc.DrawText(nucleotides[i], x1-radius/4, y1-radius/2.8)


		#draw second nucleotide
		radius = first_nucleotide_thickness+second_nucleotide_thickness
		thickness = second_nucleotide_thickness
#		font = wx.Font(pointSize=thickness/2, family=wx.FONTFAMILY_SWISS, style=wx.FONTWEIGHT_BOLD, weight=wx.FONTWEIGHT_BOLD)
#		self.gcdc.SetFont(font)
		self.gcdc.SetBrush(wx.Brush("#ffd976"))
		nucleotides = ['UU', 'UC', 'UA', 'UG','CU', 'CC', 'CA', 'CG','AU', 'AC', 'AA', 'AG', 'GU', 'GC', 'GA', 'GG']
		for i in range(len(nucleotides)):
			start_angle = 0 + 22.5*i
			finish_angle = 22.5+22.5*i
			pointlist = self.make_arc(self.xc, self.yc, start_angle, finish_angle, radius, thickness, step=0.5)
			self.gcdc.DrawPolygon(pointlist)
			x1, y1 = self.AngleToPoints(self.xc, self.yc, radius/1.2, finish_angle-(finish_angle-start_angle)/2)

			#if nucleotide is part of degenerate codon it should have a different color
			self.gcdc.SetTextForeground((nucleotide_color))
			if self.codon is not False:
#				print('first')
#				print('nuc', nucleotides[i].replace('U','T'))
#				print('cod', dna.UnAmb(self.codon[0:2]))
				if nucleotides[i].replace('U','T') in dna.UnAmb(self.codon[0:2]):
					self.gcdc.SetTextForeground((coding_nucleotide_color))
			self.gcdc.DrawText(nucleotides[i][1], x1-radius/14, y1-radius/10)


		#draw third nucleotide
		radius = first_nucleotide_thickness+second_nucleotide_thickness+third_nucleotide_thickness
		thickness = third_nucleotide_thickness
#		font = wx.Font(pointSize=thickness/2, family=wx.FONTFAMILY_SWISS, style=wx.FONTWEIGHT_BOLD, weight=wx.FONTWEIGHT_BOLD)
#		self.gcdc.SetFont(font)
		self.gcdc.SetBrush(wx.Brush("#ffc700"))
		nucleotides = ['UUU', 'UUC', 'UUA', 'UUG','UCU', 'UCC', 'UCA', 'UCG','UAU', 'UAC', 'UAA', 'UAG', 'UGU', 'UGC', 'UGA', 'UGG',\
					'CUU', 'CUC', 'CUA', 'CUG','CCU', 'CCC', 'CCA', 'CCG','CAU', 'CAC', 'CAA', 'CAG', 'CGU', 'CGC', 'CGA', 'CGG',\
					'AUU', 'AUC', 'AUA', 'AUG','ACU', 'ACC', 'ACA', 'ACG','AAU', 'AAC', 'AAA', 'AAG', 'AGU', 'AGC', 'AGA', 'AGG',\
					'GUU', 'GUC', 'GUA', 'GUG','GCU', 'GCC', 'GCA', 'GCG','GAU', 'GAC', 'GAA', 'GAG', 'GGU', 'GGC', 'GGA', 'GGG']
		for i in range(len(nucleotides)):
			start_angle = 0 + 5.625*i
			finish_angle = 5.625+5.625*i
			pointlist = self.make_arc(self.xc, self.yc, start_angle, finish_angle, radius, thickness, step=0.1)
			self.gcdc.DrawPolygon(pointlist)
			x1, y1 = self.AngleToPoints(self.xc, self.yc, radius/1.1, finish_angle-(finish_angle-start_angle)/2)

			#if nucleotide is part of degenerate codon it should have a different color
			self.gcdc.SetTextForeground((nucleotide_color))
			if self.codon is not False:
				if nucleotides[i].replace('U','T') in dna.UnAmb(self.codon):
					self.gcdc.SetTextForeground((coding_nucleotide_color))
			self.gcdc.DrawText(nucleotides[i][2], x1-radius/30, y1-radius/24)

		#draw amino acids
		radius = first_nucleotide_thickness+second_nucleotide_thickness+third_nucleotide_thickness+amino_acid_thickness
		thickness = amino_acid_thickness
#		font = wx.Font(pointSize=third_nucleotide_thickness/2, family=wx.FONTFAMILY_SWISS, style=wx.FONTWEIGHT_NORMAL, weight=wx.FONTWEIGHT_NORMAL)
#		self.gcdc.SetFont(font)
		self.gcdc.SetTextForeground(('#000000'))

		AA = ['F', 'L', 'S', 'Y', 'stop', 'C', 'stop2', 'W', 'L2', 'P', 'H', 'Q', 'R', 'I', 'M', 'T', 'N', 'K', 'S2', 'R2', 'V', 'A', 'D', 'E', 'G']
		AA_width = {'F':2, 'L':2, 'S':4, 'Y':2, 'stop':2, 'C':2, 'stop2':1, 'W':1, 'L2':4, 'P':4, 'H':2, 'Q':2, 'R':4, 'I':3, 'M':1, 'T':4, 'N':2, 'K':2, 'S2':2, 'R2':2, 'V':4, 'A':4, 'D':2, 'E':2, 'G':4}
		AA_full = {'F':'Phenylalanine', 'L':'Leucine', 'S':'Serine', 'Y':'Tyrosine', 'stop':'Stop', 'C':'Cysteine', 'stop2':'Stop', 'W':'Tryptophan', 'L2':'Leucine', 'P':'Proline', 'H':'Histidine', 'Q':'Glutamine', 'R':'Arginine', 'I':'Isoleucine', 'M':'Methionine', 'T':'Threonine', 'N':'Asparagine', 'K':'Lysine', 'S2':'Serine', 'R2':'Arginine', 'V':'Valine', 'A':'Alanine', 'D':'Aspartic acid', 'E':'Glutamic acid', 'G':'Glycine'}
		finish_angle = 0

		#basic draw of the amino acid segments and names
		for i in range(len(AA)):
			#draw the amino acid segments
			self.gcdc.SetPen(wx.Pen(colour='#666666', width=0))
			if AA[i].replace('2','') in self.target: #if current AA is a selected one
				self.gcdc.SetBrush(wx.Brush(target_color))
			elif AA[i].replace('2','') in self.offtarget: #if it is in the off-targets list
				self.gcdc.SetBrush(wx.Brush(offtarget_color))
			elif AA[i].replace('2','') in self.possible: #if current AA is among the ones that may be selected without further off-targets
				self.gcdc.SetBrush(wx.Brush(possible_color))
			else:									#otherwise use standard color
				self.gcdc.SetBrush(wx.Brush("#fff2d1"))
			start_angle = finish_angle
			finish_angle = start_angle+5.625*AA_width[AA[i]]
			pointlist = self.make_arc(self.xc, self.yc, start_angle, finish_angle, radius, thickness, step=0.1)
			self.gcdc.DrawPolygon(pointlist)

			#draw hidden color which is used for hittests
			self.catalog[str(self.NextRGB()+(255,))] = AA[i]

			self.hidden_dc.SetPen(wx.Pen(colour=self.unique_color, width=0))
			self.hidden_dc.SetBrush(wx.Brush(colour=self.unique_color))
			self.hidden_dc.DrawPolygon(pointlist)			

			#draw lines
			angle = start_angle
			self.gcdc.SetPen(wx.Pen(colour=line_color, width=1))
			if angle in [0,90,180,270]:
				radius = 0
			elif angle % 22.5 == 0:
				radius = first_nucleotide_thickness
			elif angle % 5.625 ==0:
				radius = first_nucleotide_thickness+second_nucleotide_thickness
			x1, y1 = self.AngleToPoints(self.xc, self.yc, radius, angle)
			radius = radius = first_nucleotide_thickness+second_nucleotide_thickness+third_nucleotide_thickness+amino_acid_thickness
			x2, y2 = self.AngleToPoints(self.xc, self.yc, radius, angle)
			self.gcdc.DrawLine(x1, y1, x2, y2)

			#draw text
			text_angle = finish_angle-(finish_angle-start_angle)/2

#			if finish_angle <= 180:
#				text_extent = self.gcdc.GetTextExtent(AA_full[AA[i]])
#				text_radius = (first_nucleotide_thickness+second_nucleotide_thickness+third_nucleotide_thickness)*1.05
#
#				#need to adjust for text height. Imagine right angled triangle. Adjecent is radius. Opposite is half of the text height. Calculate tan angle.
#				tanangle = (0.5*text_extent[1])/text_radius #calculate the Tan(angle)
#				radians = math.atan(tanangle) #negate the Tan part and get radians
#				degrees = radians*(180/math.pi)	#convert radians to degrees
#				text_position_angle = text_angle-degrees			
#
#				tx, ty = self.AngleToPoints(self.xc, self.yc, text_radius, text_position_angle)
#				self.gcdc.DrawRotatedText(AA_full[AA[i]], tx, ty, -text_angle+90)
#			else:
#				text_extent = self.gcdc.GetTextExtent(AA_full[AA[i]])
#				text_radius = (first_nucleotide_thickness+second_nucleotide_thickness+third_nucleotide_thickness)*1.05 + text_extent[0]
#
#				#need to adjust for text height. Imagine right angled triangle. Adjecent is radius. Opposite is half of the text height. Calculate tan angle.
#				tanangle = (0.5*text_extent[1])/text_radius #calculate the Tan(angle)
#				radians = math.atan(tanangle) #negate the Tin part and get radians
#				degrees = radians*(180/math.pi)	#convert radians to degrees
#				text_position_angle = text_angle+degrees			
#
#				tx, ty = self.AngleToPoints(self.xc, self.yc, text_radius, text_position_angle)
#				self.gcdc.DrawRotatedText(AA_full[AA[i]], tx, ty, -text_angle-90)



		#now draw the highlighted ones
		finish_angle = 0
		self.gcdc.SetPen(wx.Pen(colour='#FF0000', width=1))
		self.gcdc.SetBrush(wx.Brush(colour=(0,0,0,0)))
		for i in range(len(AA)):
			#if current AA is highlighted, redraw that segment with a different pen
			start_angle = finish_angle
			finish_angle = start_angle+5.625*AA_width[AA[i]]
			if AA[i].replace('2', '') == self.highlighted: #if highlighted AA is the current one
				pointlist = self.make_arc(self.xc, self.yc, start_angle, finish_angle, radius, thickness, step=0.1)
				self.gcdc.DrawPolygon(pointlist)


		#write what the ambigous codon is 
		point_size = int(self.Radius/8)
#		font = wx.Font(pointSize=point_size, family=wx.FONTFAMILY_SWISS, style=wx.ITALIC, weight=wx.FONTWEIGHT_NORMAL)
#		self.gcdc.SetFont(font)
		self.gcdc.SetTextForeground((line_color))
		x = self.size[0]*0.62
		y = self.size[1]*0.1
		text = 'Codon:'
		self.gcdc.DrawText(text, x, y)

		x = self.size[0]*0.75
		point_size = int(self.Radius/6)
#		font = wx.Font(pointSize=point_size, family=wx.FONTFAMILY_SWISS, style=wx.FONTWEIGHT_BOLD, weight=wx.FONTWEIGHT_BOLD)
#		self.gcdc.SetFont(font)
		self.gcdc.SetTextForeground((coding_nucleotide_color))
		
		if self.codon is False:
			text = ''
		else:
			text = self.codon
		self.gcdc.DrawText(text, x, y)

		#draw custom text field to display codon
		self.text_field.SetText(text)
		bit1, bit2 = self.text_field.GetBitmaps()
#		self.gcdc.DrawBitmap(bit1, x, y)
#		self.hidden_dc.DrawBitmap(bit2, x, y)

		#draw reset button
#		bit1, bit2 = self.reset_btn.GetBitmaps()
#		self.gcdc.DrawBitmap(bit1, 15, self.size[1]-45)
#		self.hidden_dc.DrawBitmap(bit2, 15, self.size[1]-45)


		
		#and the bases they code for
		if self.codon is not False:
			#get text position based on the ambigous codon
			first_x = x + self.gcdc.GetTextExtent(text[0:1])[0] - self.gcdc.GetTextExtent(text[0])[0]/2
			second_x = x + self.gcdc.GetTextExtent(text[0:2])[0] - self.gcdc.GetTextExtent(text[1])[0]/2
			third_x = x + self.gcdc.GetTextExtent(text[0:3])[0] - self.gcdc.GetTextExtent(text[2])[0]/2
			start_y = y + self.gcdc.GetTextExtent(text[0])[1]

			#set new text size
			point_size = int(self.Radius/18)
#			font = wx.Font(pointSize=point_size, family=wx.FONTFAMILY_SWISS, style=wx.FONTSTYLE_ITALIC, weight=wx.FONTWEIGHT_NORMAL)
#			self.gcdc.SetFont(font)
			self.gcdc.SetTextForeground((coding_nucleotide_color))

			first = dna.UnAmb(self.codon[0])
			second = dna.UnAmb(self.codon[1])
			third = dna.UnAmb(self.codon[2])

			first_y = start_y
			for i in range(0, len(first)):
				text = first[i]
				#adjust for the size of that text
				pos_x = first_x - self.gcdc.GetTextExtent(text)[0]/2
				self.gcdc.DrawText(text, pos_x, first_y)
				first_y += point_size

			second_y = start_y
			for i in range(0, len(second)):
				text = second[i]
				#adjust for the size of that text
				pos_x = second_x - self.gcdc.GetTextExtent(text)[0]/2
				self.gcdc.DrawText(text, pos_x, second_y)
				second_y += point_size

			third_y = start_y
			for i in range(0, len(third)):
				text = third[i]
				#adjust for the size of that text
				pos_x = third_x - self.gcdc.GetTextExtent(text)[0]/2
				self.gcdc.DrawText(text, pos_x, third_y)
				third_y += point_size

		#draw key	
		width = self.Radius/16
		height = self.Radius/16
		
		point_size = int(self.Radius/20)
#		font = wx.Font(pointSize=point_size, family=wx.FONTFAMILY_SWISS, style=wx.FONTWEIGHT_NORMAL, weight=wx.FONTWEIGHT_NORMAL)
#		self.gcdc.SetFont(font)
		self.gcdc.SetTextForeground(('#666666'))

		#target key
		text = 'Target AA'
		x = 10
		y = 10
		self.gcdc.SetBrush(wx.Brush(target_color))
		self.gcdc.SetPen(wx.Pen(colour='#666666', width=0))
		self.gcdc.DrawRectangle(x, y, width, height)
		self.gcdc.DrawText(text, x+width*1.2, y)

		#possible key
		text = 'Possible AA'
		x = 10
		y += point_size*1.5
		self.gcdc.SetBrush(wx.Brush(possible_color))
		self.gcdc.SetPen(wx.Pen(colour='#E6E65C', width=1))
		self.gcdc.DrawRectangle(x, y, width, height)
		self.gcdc.DrawText(text, x+width*1.2, y)

		#possible key
		text = 'Off-target AA'
		x = 10
		y += point_size*1.5
		self.gcdc.SetBrush(wx.Brush(offtarget_color))
		self.gcdc.SetPen(wx.Pen(colour='#666666', width=0))
		self.gcdc.DrawRectangle(x, y, width, height)
		self.gcdc.DrawText(text, x+width*1.2, y)


		## draw graph ##
		AA_order = ('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', 'stop')
		AA_full = {'F':'Phenylalanine', 'L':'Leucine', 'S':'Serine', 'Y':'Tyrosine', 'stop':'Stop', 'C':'Cysteine', 'stop2':'Stop', 'W':'Tryptophan', 'L2':'Leucine', 'P':'Proline', 'H':'Histidine', 'Q':'Glutamine', 'R':'Arginine', 'I':'Isoleucine', 'M':'Methionine', 'T':'Threonine', 'N':'Asparagine', 'K':'Lysine', 'S2':'Serine', 'R2':'Arginine', 'V':'Valine', 'A':'Alanine', 'D':'Aspartic acid', 'E':'Glutamic acid', 'G':'Glycine'}

		originx = self.size[0]*0.75 #centre of plot in x
		originy = self.size[1]*0.45 #centre of plot in y
		sizex = self.size[0]*0.2
		sizey = self.size[1]*0.5

		xspacing = sizex/7
		yspacing = float(sizey)/float(21)
		tick_size = sizex/30
		

		#draw background rectangle
		self.gcdc.SetBrush(wx.Brush("#fff2d1"))
		self.gcdc.SetPen(wx.Pen(colour=line_color, width=0))
		self.gcdc.DrawRectangle(originx, originy, sizex, sizey)

		#frame
		self.gcdc.SetPen(wx.Pen(line_color))
		self.gcdc.DrawLine(originx, originy, originx+sizex, originy) #top line... #DrawLine(self, x1, y1, x2, y2)
		self.gcdc.DrawLine(originx+sizex, originy, originx+sizex, originy+sizey) #right line
		self.gcdc.DrawLine(originx+sizex, originy+sizey, originx, originy+sizey) #bottom line
		self.gcdc.DrawLine(originx, originy+sizey, originx, originy) #left line


		#title
		point_size = int(sizex/12)
#		font = wx.Font(pointSize=point_size, family=wx.FONTFAMILY_SWISS, style=wx.FONTWEIGHT_NORMAL, weight=wx.FONTWEIGHT_NORMAL)
#		self.gcdc.SetFont(font)
		self.gcdc.SetTextForeground((line_color))
		title = 'Codon count for each AA'
		self.gcdc.DrawText(title, originx, originy-self.gcdc.GetTextExtent(text)[1]*2)



		#y labels (amino acids)		
		point_size = int(self.Radius/18)
#		font = wx.Font(pointSize=point_size, family=wx.FONTFAMILY_SWISS, style=wx.FONTWEIGHT_NORMAL, weight=wx.FONTWEIGHT_NORMAL)
#		self.gcdc.SetFont(font)
		self.gcdc.SetTextForeground((line_color))
		for i in range(0, 21):	
			amino_acid = AA_full[str(AA_order[i])]
			self.gcdc.DrawText(amino_acid, originx-self.gcdc.GetTextExtent(amino_acid)[0]-tick_size, originy+(yspacing*i)-self.gcdc.GetTextExtent(amino_acid)[1]/2+yspacing/2)

		#x labels (count)
		for i in range(1, 7):	
			self.gcdc.DrawText(str(i), originx+xspacing*i-self.gcdc.GetTextExtent('6')[0]/2, originy-self.gcdc.GetTextExtent('6')[1]-tick_size/2)	



		#x ticks
#		for i in range(1, 7):
#			self.gcdc.DrawLine(originx+(xspacing*i), originy, originx+(xspacing*i), originy+tick_size)
#			self.gcdc.DrawLine(originx+(xspacing*i), originy+sizey, originx+(xspacing*i), originy+sizey-tick_size)

#		#y ticks
#		for i in range(1, 22):
#			self.gcdc.DrawLine(originx, originy+(yspacing*i), originx+tick_size, originy+(yspacing*i))
#			self.gcdc.DrawLine(originx+sizex, originy+(yspacing*i), originx+sizex-tick_size, originy+(yspacing*i))

		#draw bars
		if self.codon is not False:
			self.AA_count = mbc.count_codon_list(dna.UnAmb(self.codon))
			for i in range(0, 21):	
				AA = AA_order[i]
				if AA in self.target: #if current AA is a selected one
					self.gcdc.SetBrush(wx.Brush(target_color))
				elif AA in self.offtarget: #if it is in the off-targets list
					self.gcdc.SetBrush(wx.Brush(offtarget_color))
				else:	
					self.gcdc.SetBrush(wx.Brush('#666666'))

				count = self.AA_count[AA]
				self.gcdc.DrawRectangle(originx, originy+yspacing*i+yspacing*0.1, count*xspacing, yspacing*0.8) #(x, y, w, h)




##############################################################

	def HitTest(self):
		'''Tests whether the mouse is over any amino acid'''
		dc = wx.ClientDC(self) #get the client dc
		x, y = self.ScreenToClient(wx.GetMousePosition()) #get coordinate of mouse event
		pixel_color = self.hidden_dc.GetPixel(x,y) #use that coordinate to find pixel on the hidden dc
		print('color', pixel_color)
		print('item', self.catalog[str(pixel_color)])
		return self.catalog[str(pixel_color)] #return the amino acid


	def OnLeftUp(self, event):
		'''When left mouse button is lifted up, determine the DNA selection from angles generated at down an up events.'''
		amino_acid = self.HitTest()
		if amino_acid is not False:
			amino_acid = amino_acid.replace('2','')

			if amino_acid == 'reset': #special case to catch reset button
				self.codon = False
				self.target = []
				self.offtarget = []
				self.possible = []
				self.text_edit_active = False #exit text editing mode
			elif amino_acid == 'text': #special case to catch the codon edit region
				self.text_edit_active = True #switch to text editing mode
			elif amino_acid not in self.target:
				self.target.append(amino_acid)
				self.text_edit_active = False #exit text editing mode
			elif amino_acid in self.target:
				self.target.remove(amino_acid)
				self.text_edit_active = False #exit text editing mode
		else:
			self.text_edit_active = False #exit text editing mode
			
		if len(self.target)>0:
			self.codon, self.target, self.offtarget, self.possible = mbc.run(self.target)
		else:
			self.codon = False
			self.offtarget = []
			self.possible = []
		self.update_ownUI()



	def OnMotion(self, event):
		'''When mouse is moved with the left button down determine the DNA selection from angle generated at mouse down and mouse move event.'''
		amino_acid = self.HitTest()
		if amino_acid is not False:
			amino_acid = amino_acid.replace('2', '')
		if amino_acid is not self.highlighted: #if the index did not change
			self.highlighted = amino_acid
			self.update_ownUI()







##### main loop
class MyApp(wx.App):
	def OnInit(self):
		frame = wx.Frame(None, -1, title="Codon View", size=(900,500))
		panel =	CodonView(frame, -1)
		frame.Centre()
		frame.Show(True)
		self.SetTopWindow(frame)
		return True


if __name__ == '__main__': #if script is run by itself and not loaded	
	app = MyApp(0)
	app.MainLoop()
