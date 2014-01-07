#!/usr/bin/python


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

import wx.richtext as rt
import wx


class create(rt.RichTextCtrl):
	'''A class to print colored output to a rich textctrl'''
	def __init__(self, parent, style):
		rt.RichTextCtrl.__init__(self, parent, style)
		self.SetEditable(False) #make it not editable
		font = wx.Font(pointSize=10, family=wx.FONTFAMILY_DEFAULT, style=wx.FONTSTYLE_NORMAL, weight=wx.FONTWEIGHT_NORMAL, underline=False, faceName='Source Code Pro', encoding=wx.FONTENCODING_DEFAULT) #could also use Inconsolata
		self.SetFont(font)



	def write(self, string, stringtype):
		'''General method for printing to the Output Panel'''
		if stringtype == 'DNA':
			self.BeginTextColour('#009999')
		elif stringtype == 'DNAcolor':
			#this is too slow!
			self.WriteText(string)
			
			color = '#333333'
			self.attr.SetTextColour(color)
			self.SetStyleEx(rt.RichTextRange(insertionpoint, insertionpoint+30), self.attr)
			
			previousbase = ''
			
			string=string.upper()
			i = 0

			[(self.attr.SetTextColour('#33CC00'), self.SetStyleEx(rt.RichTextRange(insertionpoint + i, insertionpoint + i+1), self.attr)) for base in string for i in xrange(len(string)) if base =='A'] 
#			for base in string:
#				start = insertionpoint + i
#				end = start + 1
#				if base == '-': color = '#000000'
#				elif base == 'A': color = '#33CC00'
#				elif base == 'T': color = '#CC0000'
#				elif base == 'C': color = '#0066CC'	
#				elif base == 'G': color = '#000000'	
#				elif base == 'N': color = '#FF00CC'	
#				else: color = '#FF6600'
#
#
#				self.attr.SetTextColour(color)
#				self.SetStyleEx(rt.RichTextRange(start, end), self.attr)
#				i += 1
		elif stringtype == 'Protein':
			self.BeginTextColour('#CC6600')
		elif stringtype == 'Text':
			self.BeginTextColour('#333333')
		elif stringtype == 'File':
			self.BeginTextColour('#330099')
		elif stringtype == 'Barcode':	
			self.BeginTextColour('#FF00FF')
			
		self.WriteText(string)
		self.EndTextColour()
		
		
		if stringtype == 'Replace':
			self.BeginTextColour('#333333')
			self.SetValue(string)
			self.EndTextColour()
			
	def write_image(self, image):
		'''General method for printing images to the Output Panel'''
		pass
#		self.WriteImage(images._rt_smiley.GetImage())	
#		
#		bool 	WriteBitmap(self, bitmap, bitmapType)
#				Write a bitmap at the current insertion point.
#		bool 	WriteImage(self, image, bitmapType)
#				Write an image at the current insertion point.
#		bool 	WriteImageBlock(self, imageBlock)
#				Write an image block at the current insertion point.
#		bool 	WriteImageFile(self, filename, bitmapType)
#				Load an image from file and write at the current insertion point.



if __name__ == '__main__': #if script is run by itself and not loaded	
	app = wx.App() # creation of the wx.App object (initialisation of the wxpython toolkit)
	frame = wx.Frame(None, title="Output Panel") # creation of a Frame with a title
	frame.output = create(frame, style=wx.VSCROLL|wx.HSCROLL) # creation of a richtextctrl in the frame
	frame.output.write('CACC', 'DNA') #testing..
	frame.Show() # frames are invisible by default so we use Show() to make them visible
	app.MainLoop() # here the app enters a loop waiting for user input


