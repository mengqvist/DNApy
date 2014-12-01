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
from base_class import DNApyBaseClass
#from wx.lib.pubsub import pub

class create(DNApyBaseClass):
	'''A class to print colored output to a rich textctrl'''
	def __init__(self, parent, style):
		super(create, self).__init__(parent, style) 

		self.rtc = rt.RichTextCtrl(self)
		
		
		self.rtc.SetEditable(False) #make it not editable
		font = wx.Font(10, wx.MODERN, wx.NORMAL, wx.NORMAL, False, u'Consolas')
		self.rtc.SetFont(font)
#		self.rtc.Bind(wx.EVT_KEY_DOWN, self.OnKeyPress)	

		#determing which listening group from which to recieve messages about UI updates
#		self.listening_group = 'placeholder' 		
#		pub.Publisher.subscribe(self.listen_to_updateUI, self.listening_group)

		sizer = wx.BoxSizer(wx.HORIZONTAL)
		sizer.Add(item=self.rtc, proportion=-1, flag=wx.EXPAND)
		self.SetSizer(sizer)



####### Modify methods from base class to fit current needs #########

	def update_globalUI(self):
		'''Method should be modified as to update other panels in response to changes in own panel.
		Preferred use is through sending a message using the pub module.
		Example use is: pub.Publisher.sendMessage('feature_list_updateUI', '').
		The first string is the "listening group" and deterimines which listeners get the message. 
		The second string is the message and is unimportant for this implementation.
		The listening group assigned here (to identify recipients) must be different from the listening group assigned in __init__ (to subscribe to messages).'''
		pass


	def update_ownUI(self):
		'''Updates all fields depending on which feature is chosen'''
		pass
		
#####################################################################

	def OnKeyPress(self, evt):
		print('keypress')

	def clear(self):
		'''Remove any text already in the Output Panel'''
		self.rtc.Clear()

	def write(self, string, stringtype):
		'''General method for printing to the Output Panel'''
		if stringtype == 'DNA':
			self.rtc.BeginTextColour('#009999')
		elif stringtype == 'DNAcolor':
			#this is too slow!
			self.rtc.WriteText(string)
			
			color = '#333333'
			self.rtc.attr.SetTextColour(color)
			self.rtc.SetStyleEx(rt.RichTextRange(insertionpoint, insertionpoint+30), self.attr)
			
			previousbase = ''
			
			string=string.upper()
			i = 0

			[(self.rtc.attr.SetTextColour('#33CC00'), self.rtc.SetStyleEx(rt.RichTextRange(insertionpoint + i, insertionpoint + i+1), self.rtc.attr)) for base in string for i in xrange(len(string)) if base =='A'] 
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
			self.rtc.BeginTextColour('#CC6600')
		elif stringtype == 'Text':
			self.rtc.BeginTextColour('#333333')
		elif stringtype == 'File':
			self.rtc.BeginTextColour('#330099')
		elif stringtype == 'Barcode':	
			self.rtc.BeginTextColour('#FF00FF')
			
		self.rtc.WriteText(string)
		self.rtc.EndTextColour()
		
		
		if stringtype == 'Replace':
			self.rtc.BeginTextColour('#333333')
			self.rtc.SetValue(string)
			self.rtc.EndTextColour()
			
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


