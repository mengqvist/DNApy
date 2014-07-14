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



import wx
import wx.lib.graphics
import math
from base_class import DNApyBaseDrawingClass
import mixed_base_codons as mbc
import dna

class CodonView(DNApyBaseDrawingClass):
	def __init__(self, parent, id):
		
		self.highlighted = False #a variable for keeping track of whether any object is highlighted
		self.codon = False
		self.target = []
		self.possible = []
		self.offtarget = []
		DNApyBaseDrawingClass.__init__(self, parent, wx.ID_ANY)


		self.centre_x = 0
		self.centre_y = 0
		self.highlighted = False #variable for keeping track of whether an amino acid is highlighted


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
		pub.Publisher.sendMessage('from_plasmid_view', '')

	
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
		self.centre_x = self.size[0]/2 #centre of window in x
		self.centre_y = self.size[1]/2 #centro of window in y
		self.window_length = min(self.centre_x, self.centre_y)
		self.Radius = self.window_length/1.2
		x = self.centre_x
		y = self.centre_y
		target_color = '#CCFF66'
		possible_color = '#FFFF66'
		offtarget_color = '#FF9966'
		

		first_nucleotide_thickness = self.Radius/3
		second_nucleotide_thickness = 2*(self.Radius/3)/3
		third_nucleotide_thickness = 1*(self.Radius/3)/3
		amino_acid_thickness = self.Radius/2.25

		dc.SetBackground(wx.Brush("White"))
		dc.Clear() # make sure you clear the bitmap!
		gcdc = wx.GCDC(dc)
		gcdc.SetPen(wx.Pen(colour='#666666', width=0))

		#make a hidden dc to which features can be drawn in uinique colors and later used for hittests
		self.hidden_dc = wx.MemoryDC()
		self.hidden_dc.SelectObject(wx.EmptyBitmap(self.ClientSize[0], self.ClientSize[1]))
		self.hidden_dc.SetBackground(wx.Brush("White"))
		self.hidden_dc.Clear() # make sure you clear the bitmap!

		#set up a dictionary to keep track of which color belongs to what object
		self.catalog = {} #for matching features with the unique colors
		self.catalog['(255, 255, 255, 255)'] = False #the background is white, have to add that key
		unique_color = (0,0,0)

		#draw reset button
		## draw a reset button here ##

		print('codon', self.codon)
		#draw first nucleotide
		radius = first_nucleotide_thickness
		thickness = first_nucleotide_thickness
		font = wx.Font(pointSize=thickness/1.5, family=wx.FONTFAMILY_SWISS, style=wx.FONTWEIGHT_BOLD, weight=wx.FONTWEIGHT_BOLD)
		gcdc.SetFont(font)
		gcdc.SetBrush(wx.Brush("#ffe7ab"))
		nucleotides = ['U', 'C', 'A', 'G']
		for i in range(len(nucleotides)):
			start_angle = 0 + 90*i
			finish_angle = 90+90*i
			pointlist = self.make_arc(x, y, start_angle, finish_angle, radius, thickness, step=5)
			gcdc.DrawPolygon(pointlist)
			x1, y1 = self.AngleToPoints(x, y, radius/2, finish_angle-(finish_angle-start_angle)/2)

			#if nucleotide is part of degenerate codon it should have a different color
			gcdc.SetTextForeground(('#000000'))
			if self.codon is not False:
				if nucleotides[i].replace('U','T') in dna.UnAmb(self.codon[0]):
					gcdc.SetTextForeground(('#993366'))
				
			gcdc.DrawText(nucleotides[i], x1-radius/4, y1-radius/2.8)


		#draw second nucleotide
		radius = first_nucleotide_thickness+second_nucleotide_thickness
		thickness = second_nucleotide_thickness
		font = wx.Font(pointSize=thickness/2, family=wx.FONTFAMILY_SWISS, style=wx.FONTWEIGHT_BOLD, weight=wx.FONTWEIGHT_BOLD)
		gcdc.SetFont(font)
		gcdc.SetBrush(wx.Brush("#ffd976"))
		nucleotides = ['UU', 'UC', 'UA', 'UG','CU', 'CC', 'CA', 'CG','AU', 'AC', 'AA', 'AG', 'GU', 'GC', 'GA', 'GG']
		for i in range(len(nucleotides)):
			start_angle = 0 + 22.5*i
			finish_angle = 22.5+22.5*i
			pointlist = self.make_arc(x, y, start_angle, finish_angle, radius, thickness, step=0.5)
			gcdc.DrawPolygon(pointlist)
			x1, y1 = self.AngleToPoints(x, y, radius/1.2, finish_angle-(finish_angle-start_angle)/2)

			#if nucleotide is part of degenerate codon it should have a different color
			gcdc.SetTextForeground(('#000000'))
			if self.codon is not False:
#				print('first')
#				print('nuc', nucleotides[i].replace('U','T'))
#				print('cod', dna.UnAmb(self.codon[0:2]))
				if nucleotides[i].replace('U','T') in dna.UnAmb(self.codon[0:2]):
					gcdc.SetTextForeground(('#993366'))
			gcdc.DrawText(nucleotides[i][1], x1-radius/14, y1-radius/10)


		#draw third nucleotide
		radius = first_nucleotide_thickness+second_nucleotide_thickness+third_nucleotide_thickness
		thickness = third_nucleotide_thickness
		font = wx.Font(pointSize=thickness/2, family=wx.FONTFAMILY_SWISS, style=wx.FONTWEIGHT_BOLD, weight=wx.FONTWEIGHT_BOLD)
		gcdc.SetFont(font)
		gcdc.SetBrush(wx.Brush("#ffc700"))
		nucleotides = ['UUU', 'UUC', 'UUA', 'UUG','UCU', 'UCC', 'UCA', 'UCG','UAU', 'UAC', 'UAA', 'UAG', 'UGU', 'UGC', 'UGA', 'UGG',\
					'CUU', 'CUC', 'CUA', 'CUG','CCU', 'CCC', 'CCA', 'CCG','CAU', 'CAC', 'CAA', 'CAG', 'CGU', 'CGC', 'CGA', 'CGG',\
					'AUU', 'AUC', 'AUA', 'AUG','ACU', 'ACC', 'ACA', 'ACG','AAU', 'AAC', 'AAA', 'AAG', 'AGU', 'AGC', 'AGA', 'AGG',\
					'GUU', 'GUC', 'GUA', 'GUG','GCU', 'GCC', 'GCA', 'GCG','GAU', 'GAC', 'GAA', 'GAG', 'GGU', 'GGC', 'GGA', 'GGG']
		for i in range(len(nucleotides)):
			start_angle = 0 + 5.625*i
			finish_angle = 5.625+5.625*i
			pointlist = self.make_arc(x, y, start_angle, finish_angle, radius, thickness, step=0.1)
			gcdc.DrawPolygon(pointlist)
			x1, y1 = self.AngleToPoints(x, y, radius/1.1, finish_angle-(finish_angle-start_angle)/2)

			#if nucleotide is part of degenerate codon it should have a different color
			gcdc.SetTextForeground(('#000000'))
			if self.codon is not False:
				if nucleotides[i].replace('U','T') in dna.UnAmb(self.codon):
					gcdc.SetTextForeground(('#993366'))
			gcdc.DrawText(nucleotides[i][2], x1-radius/30, y1-radius/24)

		#draw amino acids
		radius = first_nucleotide_thickness+second_nucleotide_thickness+third_nucleotide_thickness+amino_acid_thickness
		thickness = amino_acid_thickness
		font = wx.Font(pointSize=third_nucleotide_thickness/2, family=wx.FONTFAMILY_SWISS, style=wx.FONTWEIGHT_NORMAL, weight=wx.FONTWEIGHT_NORMAL)
		gcdc.SetFont(font)
		gcdc.SetTextForeground(('#000000'))

		AA = ['F', 'L', 'S', 'Y', 'stop', 'C', 'stop2', 'W', 'L2', 'P', 'H', 'Q', 'R', 'I', 'M', 'T', 'N', 'K', 'S2', 'R2', 'V', 'A', 'D', 'E', 'G']
		AA_width = {'F':2, 'L':2, 'S':4, 'Y':2, 'stop':2, 'C':2, 'stop2':1, 'W':1, 'L2':4, 'P':4, 'H':2, 'Q':2, 'R':4, 'I':3, 'M':1, 'T':4, 'N':2, 'K':2, 'S2':2, 'R2':2, 'V':4, 'A':4, 'D':2, 'E':2, 'G':4}
		AA_full = {'F':'Phenylalanine', 'L':'Leucine', 'S':'Serine', 'Y':'Tyrosine', 'stop':'Stop', 'C':'Cysteine', 'stop2':'Stop', 'W':'Tryptophan', 'L2':'Leucine', 'P':'Proline', 'H':'Histidine', 'Q':'Glutamine', 'R':'Arginine', 'I':'Isoleucine', 'M':'Methionine', 'T':'Threonine', 'N':'Asparagine', 'K':'Lysine', 'S2':'Serine', 'R2':'Arginine', 'V':'Valine', 'A':'Alanine', 'D':'Aspartic acid', 'E':'Glutamic acid', 'G':'Glycine'}
		finish_angle = 0

		#basic draw of the amino acid segments and names
		for i in range(len(AA)):
			#draw the amino acid segments
			gcdc.SetPen(wx.Pen(colour='#666666', width=0))
			if AA[i].replace('2','') in self.target: #if current AA is a selected one
				gcdc.SetBrush(wx.Brush(target_color))
			elif AA[i].replace('2','') in self.offtarget: #if it is in the off-targets list
				gcdc.SetBrush(wx.Brush(offtarget_color))
			elif AA[i].replace('2','') in self.possible: #if current AA is among the ones that may be selected without further off-targets
				gcdc.SetBrush(wx.Brush(possible_color))
			else:									#otherwise use standard color
				gcdc.SetBrush(wx.Brush("#fff2d1"))
			start_angle = finish_angle
			finish_angle = start_angle+5.625*AA_width[AA[i]]
			pointlist = self.make_arc(x, y, start_angle, finish_angle, radius, thickness, step=0.1)
			gcdc.DrawPolygon(pointlist)

			#draw hidden color which is used for hittests
			unique_color = self.GetNextRGB(unique_color) #get color for drawing unique colors on the hidden dc
			self.catalog[str(unique_color+(255,))] = AA[i]

			self.hidden_dc.SetPen(wx.Pen(colour=unique_color, width=0))
			self.hidden_dc.SetBrush(wx.Brush(colour=unique_color))
			self.hidden_dc.DrawPolygon(pointlist)			

			#draw lines
			angle = start_angle
			gcdc.SetPen(wx.Pen(colour='#666666', width=1))
			if angle in [0,90,180,270]:
				radius = 0
			elif angle % 22.5 == 0:
				radius = first_nucleotide_thickness
			elif angle % 5.625 ==0:
				radius = first_nucleotide_thickness+second_nucleotide_thickness
			x1, y1 = self.AngleToPoints(x, y, radius, angle)
			radius = radius = first_nucleotide_thickness+second_nucleotide_thickness+third_nucleotide_thickness+amino_acid_thickness
			x2, y2 = self.AngleToPoints(x, y, radius, angle)
			gcdc.DrawLine(x1, y1, x2, y2)

			#draw text
			text_angle = finish_angle-(finish_angle-start_angle)/2
			if finish_angle <= 180:
				text_extent = gcdc.GetTextExtent(AA_full[AA[i]])
				text_radius = (first_nucleotide_thickness+second_nucleotide_thickness+third_nucleotide_thickness)*1.05

				#need to adjust for text height. Imagine right angled triangle. Adjecent is radius. Opposite is half of the text height. Calculate tan angle.
				tanangle = (0.5*text_extent[1])/text_radius #calculate the Tan(angle)
				radians = math.atan(tanangle) #negate the Tin part and get radians
				degrees = radians*(180/math.pi)	#convert radians to degrees
				text_position_angle = text_angle-degrees			

				tx, ty = self.AngleToPoints(x, y, text_radius, text_position_angle)
				gcdc.DrawRotatedText(AA_full[AA[i]], tx, ty, -text_angle+90)
			else:
				text_extent = gcdc.GetTextExtent(AA_full[AA[i]])
				text_radius = (first_nucleotide_thickness+second_nucleotide_thickness+third_nucleotide_thickness)*1.05 + text_extent[0]

				#need to adjust for text height. Imagine right angled triangle. Adjecent is radius. Opposite is half of the text height. Calculate tan angle.
				tanangle = (0.5*text_extent[1])/text_radius #calculate the Tan(angle)
				radians = math.atan(tanangle) #negate the Tin part and get radians
				degrees = radians*(180/math.pi)	#convert radians to degrees
				text_position_angle = text_angle+degrees			

				tx, ty = self.AngleToPoints(x, y, text_radius, text_position_angle)
				gcdc.DrawRotatedText(AA_full[AA[i]], tx, ty, -text_angle-90)



		#now draw the highlighted ones
		finish_angle = 0
		gcdc.SetPen(wx.Pen(colour='#FF0000', width=1))
		gcdc.SetBrush(wx.Brush(colour=(0,0,0,0)))
		for i in range(len(AA)):
			#if current AA is highlighted, redraw that segment with a different pen
			start_angle = finish_angle
			finish_angle = start_angle+5.625*AA_width[AA[i]]
			if AA[i].replace('2', '') == self.highlighted: #if highlighted AA is the current one
				pointlist = self.make_arc(x, y, start_angle, finish_angle, radius, thickness, step=0.1)
				gcdc.DrawPolygon(pointlist)


		#draw key	
		width = self.Radius/16
		height = self.Radius/16
		
		point_size = int(self.Radius/16)
		font = wx.Font(pointSize=point_size, family=wx.FONTFAMILY_SWISS, style=wx.FONTWEIGHT_NORMAL, weight=wx.FONTWEIGHT_NORMAL)
		gcdc.SetFont(font)
		gcdc.SetTextForeground(('#666666'))

		#target key
		text = 'Target AA'
		x = self.centre_x-self.Radius
		y = self.Radius*0.1
		gcdc.SetBrush(wx.Brush(target_color))
		gcdc.SetPen(wx.Pen(colour='#666666', width=0))
		gcdc.DrawRectangle(x, y, width, height)
		gcdc.DrawText(text, x+width*1.2, y-height*0.2)

		#possible key
		text = 'Possible AA'
		x = self.centre_x-self.Radius
		y = self.Radius*0.1+point_size*1.5
		gcdc.SetBrush(wx.Brush(possible_color))
		gcdc.SetPen(wx.Pen(colour='#E6E65C', width=1))
		gcdc.DrawRectangle(x, y, width, height)
		gcdc.DrawText(text, x+width*1.2, y-height*0.2)

		#possible key
		text = 'Off-target AA'
		x = self.centre_x-self.Radius
		y = self.Radius*0.1+point_size*3
		gcdc.SetBrush(wx.Brush(offtarget_color))
		gcdc.SetPen(wx.Pen(colour='#666666', width=0))
		gcdc.DrawRectangle(x, y, width, height)
		gcdc.DrawText(text, x+width*1.2, y-height*0.2)


##############################################################

	def HitTest(self):
		'''Tests whether the mouse is over any amino acid'''
		dc = wx.ClientDC(self) #get the client dc
		x, y = self.ScreenToClient(wx.GetMousePosition()) #get coordinate of mouse event
		pixel_color = self.hidden_dc.GetPixel(x,y) #use that coordinate to find pixel on the hidden d
		return self.catalog[str(pixel_color)] #return the amino acid


	def OnLeftUp(self, event):
		'''When left mouse button is lifted up, determine the DNA selection from angles generated at down an up events.'''
		amino_acid = self.HitTest().replace('2','')
		if amino_acid not in self.target:
			self.target.append(amino_acid)
		elif amino_acid in self.target:
			self.target.remove(amino_acid)
		elif amino_acid == 'reset':
			self.codon = False
			self.target = []
			self.offtarget = []
			self.possible = []
			
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
		frame = wx.Frame(None, -1, title="Codon View", size=(500,600), style = wx.NO_FULL_REPAINT_ON_RESIZE)
		panel =	CodonView(frame, -1)
		frame.Centre()
		frame.Show(True)
		self.SetTopWindow(frame)
		return True


if __name__ == '__main__': #if script is run by itself and not loaded	
	app = MyApp(0)
	app.MainLoop()
