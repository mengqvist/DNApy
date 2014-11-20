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
import math
from base_class import DNApyBaseDrawingClass
from base_class import DNApyBaseClass
import mixed_base_codons as mbc
import dna
import protein


		
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
		self.table = 1 #codon table

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
		'''Receives requests for DNA selection and then sends it.'''
		assert type(selection) == tuple, 'Error, dna selection must be a tuple'
		selection = (int(selection[0]), int(selection[1]))
		genbank.dna_selection = selection


############### Done setting required methods #######################




	def Draw(self, dc):
		'''Class for drawing stuff on gcdc.
		It is important to only use the base_class drawing tools and not the inbuilt gcdc tools, 
		otherwise the hittests will not work'''
		self.xc = self.size[0]/3 #centre of codon circle in x
		self.yc = self.size[1]/2 #centre of codon circle in y
		self.Radius = self.yc/1.2
		self.unique_color = (0,0,0)

		dc.SetBackground(wx.Brush("White"))
		dc.Clear() # make sure you clear the bitmap!
		gcdc = wx.GCDC(dc) #make gcdc from the dc (for use of transparency and antialiasing)

		#make a hidden dc to which features can be drawn in unique colors and later used for hittests
		self.hidden_dc = wx.MemoryDC()
		self.hidden_dc.SelectObject(wx.EmptyBitmap(self.ClientSize[0], self.ClientSize[1]))
		self.hidden_dc.SetBackground(wx.Brush("White"))
		self.hidden_dc.Clear() # make sure you clear the bitmap!

		#set what colors the different fields should have
		target_color = '#CCFF66'
		possible_color = '#FFFF66'
		offtarget_color = '#FF9966'
		nucleotide_color = '#8B835F'
		coding_nucleotide_color = '#4B4424' #for coloring the nucleotides encoded by the ambigous codon
		line_color = '#8B835F' #for lines
		
		#These parametrs determine the "thickness" of the amino acid sections
		first_nucleotide_thickness = self.Radius/3.0
		second_nucleotide_thickness = (self.Radius/3)/1.5
		third_nucleotide_thickness = (self.Radius/3)/4.0
		amino_acid_thickness = self.Radius/2.25



		
		#draw first nucleotide
		radius = first_nucleotide_thickness
		thickness = first_nucleotide_thickness
		font = wx.Font(pointSize=thickness/1.5, family=wx.FONTFAMILY_SWISS, style=wx.FONTWEIGHT_NORMAL, weight=wx.FONTWEIGHT_BOLD)
		gcdc.SetFont(font)
		gcdc.SetPen(wx.Pen(colour="#ffe7ab", width=0))
		gcdc.SetBrush(wx.Brush("#ffe7ab"))
		nucleotides = ['T', 'C', 'A', 'G']
		for i in range(len(nucleotides)):
			start_angle = 0 + 90*i
			finish_angle = 90+90*i
			pointlist = self.make_arc(self.xc, self.yc, start_angle, finish_angle, radius, thickness, step=5)
			gcdc.DrawPolygon(pointlist)
			x1, y1 = self.AngleToPoints(self.xc, self.yc, radius/2, finish_angle-(finish_angle-start_angle)/2)

			#if nucleotide is part of degenerate codon it should have a different color
			gcdc.SetTextForeground((nucleotide_color))
			if self.codon is not False:
				if nucleotides[i].replace('U','T') in dna.UnAmb(self.codon[0]):
					gcdc.SetTextForeground((coding_nucleotide_color))
				
			gcdc.DrawText(nucleotides[i], x1-radius/3.1, y1-radius/2.1)


		#draw second nucleotide
		radius = first_nucleotide_thickness+second_nucleotide_thickness
		thickness = second_nucleotide_thickness
		font = wx.Font(pointSize=thickness/2, family=wx.FONTFAMILY_SWISS, style=wx.FONTWEIGHT_NORMAL, weight=wx.FONTWEIGHT_BOLD)
		gcdc.SetFont(font)
		gcdc.SetPen(wx.Pen(colour="#ffd976", width=0))
		gcdc.SetBrush(wx.Brush("#ffd976"))
		nucleotides = ['TT', 'TC', 'TA', 'TG','CT', 'CC', 'CA', 'CG','AT', 'AC', 'AA', 'AG', 'GT', 'GC', 'GA', 'GG']
		for i in range(len(nucleotides)):
			start_angle = 0 + 22.5*i
			finish_angle = 22.5+22.5*i
			pointlist = self.make_arc(self.xc, self.yc, start_angle, finish_angle, radius, thickness, step=0.5)
			gcdc.DrawPolygon(pointlist)
			x1, y1 = self.AngleToPoints(self.xc, self.yc, radius/1.2, finish_angle-(finish_angle-start_angle)/2)

			#if nucleotide is part of degenerate codon it should have a different color
			gcdc.SetTextForeground((nucleotide_color))
			if self.codon is not False:
#				print('first')
#				print('nuc', nucleotides[i].replace('U','T'))
#				print('cod', dna.UnAmb(self.codon[0:2]))
				if nucleotides[i].replace('U','T') in dna.UnAmb(self.codon[0:2]):
					gcdc.SetTextForeground((coding_nucleotide_color))
			gcdc.DrawText(nucleotides[i][1], x1-radius/14.0, y1-radius/8.0)


		#draw third nucleotide
		radius = first_nucleotide_thickness+second_nucleotide_thickness+third_nucleotide_thickness
		thickness = third_nucleotide_thickness
		font = wx.Font(pointSize=thickness/2, family=wx.FONTFAMILY_SWISS, style=wx.FONTWEIGHT_NORMAL, weight=wx.FONTWEIGHT_BOLD)
		gcdc.SetFont(font)
		gcdc.SetPen(wx.Pen(colour="#ffc700", width=0))
		gcdc.SetBrush(wx.Brush("#ffc700"))
		codons = ['TTT', 'TTC', 'TTA', 'TTG','TCT', 'TCC', 'TCA', 'TCG','TAT', 'TAC', 'TAA', 'TAG', 'TGT', 'TGC', 'TGA', 'TGG',\
					'CTT', 'CTC', 'CTA', 'CTG','CCT', 'CCC', 'CCA', 'CCG','CAT', 'CAC', 'CAA', 'CAG', 'CGT', 'CGC', 'CGA', 'CGG',\
					'ATT', 'ATC', 'ATA', 'ATG','ACT', 'ACC', 'ACA', 'ACG','AAT', 'AAC', 'AAA', 'AAG', 'AGT', 'AGC', 'AGA', 'AGG',\
					'GTT', 'GTC', 'GTA', 'GTG','GCT', 'GCC', 'GCA', 'GCG','GAT', 'GAC', 'GAA', 'GAG', 'GGT', 'GGC', 'GGA', 'GGG']
		for i in range(len(codons)):
			start_angle = 0 + 5.625*i
			finish_angle = 5.625+5.625*i
			pointlist = self.make_arc(self.xc, self.yc, start_angle, finish_angle, radius, thickness, step=0.1)
			gcdc.DrawPolygon(pointlist)
			x1, y1 = self.AngleToPoints(self.xc, self.yc, radius/1.05, finish_angle-(finish_angle-start_angle)/2)

			#if nucleotide is part of degenerate codon it should have a different color
			gcdc.SetTextForeground((nucleotide_color))
			if self.codon is not False:
				if codons[i].replace('U','T') in dna.UnAmb(self.codon):
					gcdc.SetTextForeground((coding_nucleotide_color))
			gcdc.DrawText(codons[i][2], x1-radius/30, y1-radius/24)

		#set parameters for drawing amino acids
		radius = first_nucleotide_thickness+second_nucleotide_thickness+third_nucleotide_thickness+amino_acid_thickness
		thickness = amino_acid_thickness
		font = wx.Font(pointSize=third_nucleotide_thickness/2, family=wx.FONTFAMILY_SWISS, style=wx.FONTWEIGHT_NORMAL, weight=wx.FONTWEIGHT_BOLD)
		gcdc.SetFont(font)
		gcdc.SetTextForeground(('#000000'))
		finish_angle = 0
		
		#basic draw of the amino acid segments and names
		AA_width = 0
		current_AA = dna.Translate(codons[0], self.table)
		for codon in codons:
			AA = dna.Translate(codon, self.table)
			if codon == 'GGG': #catch the last codon
				AA_width += 1
				AA = None
				
			if current_AA == AA:
				AA_width += 1
			else:
				#draw the amino acid segments
				gcdc.SetPen(wx.Pen(colour='#666666', width=0))
				if current_AA in self.target: #if current AA is a selected one
					gcdc.SetBrush(wx.Brush(target_color))
				elif current_AA in self.offtarget: #if it is in the off-targets list
					gcdc.SetBrush(wx.Brush(offtarget_color))
				elif current_AA in self.possible: #if current AA is among the ones that may be selected without further off-targets
					gcdc.SetBrush(wx.Brush(possible_color))
				else:									#otherwise use standard color
					gcdc.SetBrush(wx.Brush("#fff2d1"))
				start_angle = finish_angle
				finish_angle = start_angle+5.625*AA_width
				pointlist = self.make_arc(self.xc, self.yc, start_angle, finish_angle, radius, thickness, step=0.1)
				gcdc.DrawPolygon(pointlist)

				#draw hidden color which is used for hittests
				self.catalog[str(self.NextRGB()+(255,))] = current_AA

				self.hidden_dc.SetPen(wx.Pen(colour=self.unique_color, width=0))
				self.hidden_dc.SetBrush(wx.Brush(colour=self.unique_color))
				self.hidden_dc.DrawPolygon(pointlist)			

				#draw lines
				angle = start_angle
				gcdc.SetPen(wx.Pen(colour=line_color, width=1))
				if angle in [0,90,180,270]:
					radius = 0
				elif angle % 22.5 == 0:
					radius = first_nucleotide_thickness
				elif angle % 5.625 ==0:
					radius = first_nucleotide_thickness+second_nucleotide_thickness
				x1, y1 = self.AngleToPoints(self.xc, self.yc, radius, angle)
				radius = radius = first_nucleotide_thickness+second_nucleotide_thickness+third_nucleotide_thickness+amino_acid_thickness
				x2, y2 = self.AngleToPoints(self.xc, self.yc, radius, angle)
				gcdc.DrawLine(x1, y1, x2, y2)

				#draw amino acid text
				text_angle = finish_angle-(finish_angle-start_angle)/2

				if finish_angle <= 180:
					text_extent = gcdc.GetTextExtent(protein.one_to_full(current_AA))
					text_radius = (first_nucleotide_thickness+second_nucleotide_thickness+third_nucleotide_thickness)*1.05

					#need to adjust for text height. Imagine right angled triangle. Adjecent is radius. Opposite is half of the text height. Calculate tan angle.
					tanangle = (0.5*text_extent[1])/text_radius #calculate the Tan(angle)
					radians = math.atan(tanangle) #negate the Tan part and get radians
					degrees = radians*(180/math.pi)	#convert radians to degrees
					text_position_angle = text_angle-degrees			

					tx, ty = self.AngleToPoints(self.xc, self.yc, text_radius, text_position_angle)
					gcdc.DrawRotatedText(protein.one_to_full(current_AA), tx, ty, -text_angle+90)
				else:
					text_extent = gcdc.GetTextExtent(protein.one_to_full(current_AA))
					text_radius = (first_nucleotide_thickness+second_nucleotide_thickness+third_nucleotide_thickness)*1.05 + text_extent[0]

					#need to adjust for text height. Imagine right angled triangle. Adjecent is radius. Opposite is half of the text height. Calculate tan angle.
					tanangle = (0.5*text_extent[1])/text_radius #calculate the Tan(angle)
					radians = math.atan(tanangle) #negate the Tin part and get radians
					degrees = radians*(180/math.pi)	#convert radians to degrees
					text_position_angle = text_angle+degrees			

					tx, ty = self.AngleToPoints(self.xc, self.yc, text_radius, text_position_angle)
					gcdc.DrawRotatedText(protein.one_to_full(current_AA), tx, ty, -text_angle-90)

				#now re-set the parameters for the next round
				current_AA = AA
				AA_width = 1


		#now draw the highlighted amino acids
#		finish_angle = 0
#		gcdc.SetPen(wx.Pen(colour='#FF0000', width=1))
#		gcdc.SetBrush(wx.Brush(colour=(0,0,0,0)))
#		for i in range(len(AA)):
#			#if current AA is highlighted, redraw that segment with a different pen
#			start_angle = finish_angle
#			finish_angle = start_angle+5.625*AA_width[AA[i]]
#			if AA[i] == self.highlighted: #if highlighted AA is the current one
#				pointlist = self.make_arc(self.xc, self.yc, start_angle, finish_angle, radius, thickness, step=0.1)
#				gcdc.DrawPolygon(pointlist)


		#write what the ambiguous codon is 
		point_size = int(self.Radius/8)
		font = wx.Font(pointSize=point_size, family=wx.FONTFAMILY_SWISS, style=wx.ITALIC, weight=wx.FONTWEIGHT_NORMAL)
		gcdc.SetFont(font)
		gcdc.SetTextForeground((line_color))
		x = self.size[0]*0.62
		y = self.size[1]*0.1
		text = 'Codon:'
		gcdc.DrawText(text, x, y)

		x = self.size[0]*0.75
		point_size = int(self.Radius/6)
		font = wx.Font(pointSize=point_size, family=wx.FONTFAMILY_SWISS, style=wx.FONTWEIGHT_NORMAL, weight=wx.FONTWEIGHT_NORMAL)
		gcdc.SetFont(font)
		gcdc.SetTextForeground((coding_nucleotide_color))
		
		if self.codon is False:
			text = ''
		else:
			text = self.codon
		gcdc.DrawText(text, x, y)

		
		#below the ambiguous codon, list the bases it codes for
		if self.codon is not False:
			#get text position based on the ambigous codon
			first_x = x + gcdc.GetTextExtent(text[0:1])[0] - gcdc.GetTextExtent(text[0])[0]/2
			second_x = x + gcdc.GetTextExtent(text[0:2])[0] - gcdc.GetTextExtent(text[1])[0]/2
			third_x = x + gcdc.GetTextExtent(text[0:3])[0] - gcdc.GetTextExtent(text[2])[0]/2
			start_y = y + gcdc.GetTextExtent(text[0])[1]

			#set new text size
			point_size = int(self.Radius/18)
			font = wx.Font(pointSize=point_size, family=wx.FONTFAMILY_SWISS, style=wx.FONTSTYLE_ITALIC, weight=wx.FONTWEIGHT_NORMAL)
			gcdc.SetFont(font)
			gcdc.SetTextForeground((coding_nucleotide_color))

			first = dna.UnAmb(self.codon[0])
			second = dna.UnAmb(self.codon[1])
			third = dna.UnAmb(self.codon[2])

			first_y = start_y
			for i in range(0, len(first)):
				text = first[i]
				#adjust for the size of that text
				pos_x = first_x - gcdc.GetTextExtent(text)[0]/2
				gcdc.DrawText(text, pos_x, first_y)
				first_y += point_size*1.2

			second_y = start_y
			for i in range(0, len(second)):
				text = second[i]
				#adjust for the size of that text
				pos_x = second_x - gcdc.GetTextExtent(text)[0]/2
				gcdc.DrawText(text, pos_x, second_y)
				second_y += point_size*1.2

			third_y = start_y
			for i in range(0, len(third)):
				text = third[i]
				#adjust for the size of that text
				pos_x = third_x - gcdc.GetTextExtent(text)[0]/2
				gcdc.DrawText(text, pos_x, third_y)
				third_y += point_size*1.2

		#draw key	
		width = self.Radius/16
		height = self.Radius/16
		
		point_size = int(self.Radius/20)
		font = wx.Font(pointSize=point_size, family=wx.FONTFAMILY_SWISS, style=wx.FONTWEIGHT_NORMAL, weight=wx.FONTWEIGHT_NORMAL)
		gcdc.SetFont(font)
		gcdc.SetTextForeground(('#666666'))

		#target AA key
		text = 'Target AA'
		x = 10
		y = 10
		gcdc.SetBrush(wx.Brush(target_color))
		gcdc.SetPen(wx.Pen(colour='#666666', width=0))
		gcdc.DrawRectangle(x, y, width, height)
		gcdc.DrawText(text, x+width*1.2, y)

		#possible AA key
		text = 'Possible AA'
		x = 10
		y += point_size*1.5
		gcdc.SetBrush(wx.Brush(possible_color))
		gcdc.SetPen(wx.Pen(colour='#E6E65C', width=1))
		gcdc.DrawRectangle(x, y, width, height)
		gcdc.DrawText(text, x+width*1.2, y)

		#off-target AA key
		text = 'Off-target AA'
		x = 10
		y += point_size*1.5
		gcdc.SetBrush(wx.Brush(offtarget_color))
		gcdc.SetPen(wx.Pen(colour='#666666', width=0))
		gcdc.DrawRectangle(x, y, width, height)
		gcdc.DrawText(text, x+width*1.2, y)


		## draw graph ##
		AA_order = ('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', '*')
		AA_full = {'F':'Phenylalanine', 'L':'Leucine', 'S':'Serine', 'Y':'Tyrosine', '*':'Stop', 'C':'Cysteine', 'stop2':'Stop', 'W':'Tryptophan', 'L2':'Leucine', 'P':'Proline', 'H':'Histidine', 'Q':'Glutamine', 'R':'Arginine', 'I':'Isoleucine', 'M':'Methionine', 'T':'Threonine', 'N':'Asparagine', 'K':'Lysine', 'S2':'Serine', 'R2':'Arginine', 'V':'Valine', 'A':'Alanine', 'D':'Aspartic acid', 'E':'Glutamic acid', 'G':'Glycine'}

		originx = self.size[0]*0.75 #centre of plot in x
		originy = self.size[1]*0.45 #centre of plot in y
		sizex = self.size[0]*0.2
		sizey = self.size[1]*0.5

		xspacing = sizex/7
		yspacing = float(sizey)/float(21)
		tick_size = sizex/30
		

		#draw background rectangle
		gcdc.SetBrush(wx.Brush("#fff2d1"))
		gcdc.SetPen(wx.Pen(colour=line_color, width=0))
		gcdc.DrawRectangle(originx, originy, sizex, sizey)

		#frame
		gcdc.SetPen(wx.Pen(line_color))
		gcdc.DrawLine(originx, originy, originx+sizex, originy) #top line... #DrawLine(self, x1, y1, x2, y2)
		gcdc.DrawLine(originx+sizex, originy, originx+sizex, originy+sizey) #right line
		gcdc.DrawLine(originx+sizex, originy+sizey, originx, originy+sizey) #bottom line
		gcdc.DrawLine(originx, originy+sizey, originx, originy) #left line


		#title
		point_size = int(sizex/15)
		font = wx.Font(pointSize=point_size, family=wx.FONTFAMILY_SWISS, style=wx.FONTWEIGHT_NORMAL, weight=wx.FONTWEIGHT_NORMAL)
		gcdc.SetFont(font)
		gcdc.SetTextForeground((line_color))
		title = 'Codon count for each AA'
		gcdc.DrawText(title, originx, originy-gcdc.GetTextExtent(text)[1]*2)



		#y labels (amino acids)		
		point_size = int(self.Radius/18)
		font = wx.Font(pointSize=point_size, family=wx.FONTFAMILY_SWISS, style=wx.FONTWEIGHT_NORMAL, weight=wx.FONTWEIGHT_NORMAL)
		gcdc.SetFont(font)
		gcdc.SetTextForeground((line_color))
		for i in range(0, 21):	
			amino_acid = AA_full[str(AA_order[i])]
			gcdc.DrawText(amino_acid, originx-gcdc.GetTextExtent(amino_acid)[0]-tick_size, originy+(yspacing*i)-gcdc.GetTextExtent(amino_acid)[1]/2+yspacing/2)

		#x labels (count)
		for i in range(1, 7):	
			gcdc.DrawText(str(i), originx+xspacing*i-gcdc.GetTextExtent('6')[0]/2, originy-gcdc.GetTextExtent('6')[1]-tick_size/2)	



		#x ticks
		for i in range(1, 7):
			gcdc.DrawLine(originx+(xspacing*i), originy, originx+(xspacing*i), originy+tick_size)
			gcdc.DrawLine(originx+(xspacing*i), originy+sizey, originx+(xspacing*i), originy+sizey-tick_size)

#		#y ticks
#		for i in range(1, 22):
#			gcdc.DrawLine(originx, originy+(yspacing*i), originx+tick_size, originy+(yspacing*i))
#			gcdc.DrawLine(originx+sizex, originy+(yspacing*i), originx+sizex-tick_size, originy+(yspacing*i))

		#draw bars according to how many times each AA is encoded
		if self.codon is not False:
			self.AA_count = protein.count_aa(dna.Translate(''.join(dna.UnAmb(self.codon)), self.table))
			for i in range(0, 21):	
				AA = AA_order[i]
				if AA in self.target: #if current AA is a selected one
					gcdc.SetBrush(wx.Brush(target_color))
				elif AA in self.offtarget: #if it is in the off-targets list
					gcdc.SetBrush(wx.Brush(offtarget_color))
				else:	
					gcdc.SetBrush(wx.Brush('#666666'))

				count = self.AA_count[AA]
				gcdc.DrawRectangle(originx, originy+yspacing*i+yspacing*0.1, count*xspacing, yspacing*0.8) #(x, y, w, h)




##############################################################

	def HitTest(self):
		'''Tests whether the mouse is over any amino acid'''
		dc = wx.ClientDC(self) #get the client dc
		x, y = self.ScreenToClient(wx.GetMousePosition()) #get coordinate of mouse event
		pixel_color = self.hidden_dc.GetPixel(x,y) #use that coordinate to find pixel on the hidden dc
#		print('color', pixel_color)
#		print('item', self.catalog[str(pixel_color)])
		return self.catalog[str(pixel_color)] #return the amino acid


	def OnLeftUp(self, event):
		'''When left mouse button is lifted up, determine the DNA selection from angles generated at down an up events.'''
		amino_acid = self.HitTest()
		if amino_acid is not False:

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
			codon_object = mbc.AmbigousCodon(self.target, self.table)
			self.codon = codon_object.getTriplet()
			self.target = codon_object.getTarget()
			self.offtarget = codon_object.getOfftarget()
			self.possible = codon_object.getPossible()
		else:
			self.codon = False
			self.offtarget = []
			self.possible = []
		self.update_ownUI()



	def OnMotion(self, event):
		pass
#		'''When mouse is moved with the left button down determine the DNA selection from angle generated at mouse down and mouse move event.'''
#		amino_acid = self.HitTest()
#		if amino_acid is not False:
#			amino_acid = amino_acid.replace('2', '') #Fix!
#		if amino_acid is not self.highlighted: #if the index did not change
#			self.highlighted = amino_acid
#			self.update_ownUI()





#make new class and bind in buttons and such			
class CodonButtonWrapper(DNApyBaseClass):
	'''
	This class is intended to glue together the plasmid drawing with control buttons.
	'''
	def __init__(self, parent, id):
		DNApyBaseClass.__init__(self, parent, id)
		self.codon_view = CodonView(self, -1)	


		

		
		##########  Add buttons and methods to control their behaviour ###########	
		#buttons
		reset = wx.Button(self, 1, 'Reset')		
		save = wx.Button(self, 3, 'Save Image')
		
		#the combobox
		options = ["1: Standard Code",
					"2: Vertebrate Mitochondrial Code",
					"3: Yeast Mitochondrial Code",
					"4: Mold, Protozoan, Coelenterate Mitochondrial Code",
					"5: Invertebrate Mitochondrial Code",
					"6: Ciliate, Dasycladacean and Hexamita Nuclear Code",
					"9: Echinoderm and Flatworm Mitochondrial Code",
					"10: Euplotid Nuclear Code",
					"11: Bacterial, Archaeal and Plant Plastid Code",
					"12: Alternative Yeast Nuclear Code",
					"13: Ascidian Mitochondrial Code",
					"14: Alternative Flatworm Mitochondrial Code",
					"15: Blepharisma Nuclear Code",
					"16: Chlorophycean Mitochondrial Code",
					"21: Trematode Mitochondrial Code",
					"22: Scenedesmus obliquus mitochondrial Code",
					"23: Thraustochytrium Mitochondrial Code",
					"24: Pterobranchia mitochondrial Code",
					"25: Candidate Division SR1 and Gracilibacteria Code"]

		self.combobox = wx.ComboBox(self, id=2, size=(-1, -1), choices=options, style=wx.CB_READONLY)
		self.combobox.Select(0)
		
		#bind feature list buttons
		self.Bind(wx.EVT_BUTTON, self.OnReset, id=1)
		self.Bind(wx.EVT_COMBOBOX, self.OnComboboxSelect, id=2)
		self.Bind(wx.EVT_BUTTON, self.OnSave, id=3)


		#arrange buttons vertically		
		sizer = wx.BoxSizer(wx.HORIZONTAL)
		sizer.Add(item=reset)
		sizer.Add(item=self.combobox)
		sizer.Add(item=save)




		#add feature list and buttons horizontally	
		sizer2 = wx.BoxSizer(wx.VERTICAL)
		sizer2.Add(item=sizer, proportion=0, flag=wx.EXPAND)
		sizer2.Add(item=self.codon_view, proportion=-1, flag=wx.EXPAND)

		self.SetSizer(sizer2)
		
	def update_ownUI(self):
		'''
		User interface updates.
		'''
		self.codon_view.update_ownUI()

	def OnReset(self, evt):
		self.codon_view.codon = False
		self.codon_view.target = []
		self.codon_view.offtarget = []
		self.codon_view.possible = []
		self.update_ownUI()
		
	def OnComboboxSelect(self, evt):
		self.codon_view.table = self.combobox.GetValue().split(':')[0]

		#compute result with new table
		if len(self.codon_view.target)>0:
			codon_object = mbc.AmbigousCodon(self.codon_view.target, self.codon_view.table)
			self.codon_view.codon = codon_object.getTriplet()
			self.codon_view.target = codon_object.getTarget()
			self.codon_view.offtarget = codon_object.getOfftarget()
			self.codon_view.possible = codon_object.getPossible()
		
		#update drawing
		self.update_ownUI()

		
	def OnSave(self, evt):	
		pass
	


##### main loop
class MyApp(wx.App):
	def OnInit(self):
		frame = wx.Frame(None, -1, title="Codon View", size=(900,500))
		panel =	CodonButtonWrapper(frame, -1)
		frame.Centre()
		frame.Show(True)
		self.SetTopWindow(frame)
		return True


if __name__ == '__main__': #if script is run by itself and not loaded	
	app = MyApp(0)
	app.MainLoop()
