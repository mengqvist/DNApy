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
from base_class import DNApyBaseDrawingClass

import featurelist_GUI
import plasmid_GUI

import cairo
from wx.lib.wxcairo import ContextFromDC
import pango # for text
import pangocairo # for glyphs and text
import math

import colcol

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




########### class for text ######################
class TextEdit(DNApyBaseDrawingClass):
	def __init__(self, parent, id):
		

		# initiate parent
		self.parent = parent
		
		# get data
		self.dna 			= genbank.gb.GetDNA()
		self.cdna 			= genbank.gb.GetDNA()

		# text styler
		self.textSpacing 	= 50
		self.fontsize 		= 9
		self.lineGap 		= 5
		self.lineHeight		= self.textSpacing +  self.fontsize 
		
		# interactive cursor
		self.PositionPointer = 1
		
		# some helper varuiable #TODO remove this
		self.minHeight = 1000
		
		# storage dict, to prevent recalculating to much
		self.cairoStorage = {
				'features'	: None, # feature Object
				'fPaths' 	: [], # path of alls features as List
				'cursorPath': None,
				'cursor'	: None
				
		}
		
		
		#self.drawinPanel = wx.Panel(self.parent)
		DNApyBaseDrawingClass.__init__(self, parent, wx.ID_ANY)
		
		# set a sizer
		#sizer = wx.BoxSizer(wx.HORIZONTAL)
		#sizer.Add(item=self.drawinPanel)
		#self.SetSizer(sizer)
		
		self.Bind(wx.EVT_LEFT_DOWN, self.OnLeftDown)
		self.Bind(wx.EVT_LEFT_UP, self.OnLeftUp)
		self.Bind(wx.EVT_MOTION, self.OnMotion) # drag and drop, selection etc.
		self.Bind(wx.EVT_KEY_DOWN, self.OnKeyPress) #This is important for controlling the input into the editor
		#wx.EVT_KEY_DOWN(self, self.OnKeyPress)
		
		#self.SetAutoLayout(True)
		#wx.Panel.__init__(self, parent)
		

####### Modify methods from base calss to fit current needs #########

	def update_globalUI(self):
		'''
		Method should be modified as to update other panels in response to changes in own panel.
		'''
		MSG_CHANGE_TEXT = "change.text"
		pub.sendMessage(MSG_CHANGE_TEXT, text="DNA view says update!")



	def update_ownUI(self):
		"""
		This would get called if the drawing needed to change, for whatever reason.

		The idea here is that the drawing is based on some data generated
		elsewhere in the system. If that data changes, the drawing needs to
		be updated.

		This code re-draws the buffer, then calls Update, which forces a paint event.
		"""
		print self.get_selection()


		# start dc
		dc = wx.MemoryDC()
		dc.SelectObject(self._Buffer)
		dc.SetBackground(wx.Brush("White"))	# start the DC and clear it
		dc.Clear() 							# make sure you clear the bitmap!
		# make a cairo context
		self.ctx = ContextFromDC(dc)		# load it into cairo
		
		# start pango as a font backend
		self.pango = pangocairo.CairoContext(self.ctx)
		self.pango.set_antialias(cairo.ANTIALIAS_SUBPIXEL)
		
		self.dna = genbank.gb.GetDNA()
		if self.dna != None:
			self.cdna = dna.C(self.dna)
		else:
			self.dna = ''	
			self.cdna = ''
			
		# render the text
		self.displayText()
		# draw a cursor
		self.drawCursor()
		# create colorful features here
		self.drawFeatures()

		dc.SelectObject(wx.NullBitmap) # need to get rid of the MemoryDC before Update() is called.
		self.Refresh()
		self.Update()


		return None

		#sequence = genbank.gb.GetDNA()
		#if sequence == None:
		#	sequence = ''
		#self.stc.SetText(sequence) #put the DNA in the editor
		#start, finish, zero = genbank.dna_selection
		#if finish == -1: #a caret insertion (and no selection). Will actually result in a selection where start is one larger than finish
		#	self.stc.GotoPos(start-1) #set caret insertion
		#else:
		#	self.stc.SetSelection(start-1, finish) #update selection

	def displayText(self):
		
		# get extend to wrap text correct
		width, height 	= self.parent.GetVirtualSize()
		if width > 30:
			w = width - 8
		else:
			w = width
		self.w = w # width for other drawing operators

		# start layout
		layout 			= self.pango.create_layout()
		layout.set_wrap(pango.WRAP_WORD_CHAR)
		layout.set_width(pango.SCALE * w)
		layout.set_spacing(pango.SCALE *self.textSpacing) 
		layout.set_alignment(pango.ALIGN_LEFT)

		
		# new attr list for the dna
		dnaAttr = pango.AttrList()
		# color
		gray 		= 56 * 65535/255
		fg_color	= pango.AttrForeground(gray,gray,gray, 0, len(self.dna))
		spacing 	= pango.AttrLetterSpacing(1500, 0, len(self.dna))
		# add attr
		dnaAttr.insert(fg_color)
		dnaAttr.insert(spacing)
		# set attr
		layout.set_attributes(dnaAttr)
		
		#set font
		font = pango.FontDescription("monospace normal "+ "1")
		font.set_size(pango.SCALE * self.fontsize)
		layout.set_font_description(font)
		


	
	

		# output sense strain
		starty = 20
		self.ctx.move_to(8,starty)
		#layout.set_text(self.dna)
		
		
		
		# make some markup
		# all markup should be performed here
		start, end, zero = self.get_selection()
		
	
		# correct direktion of selection
		if start > end:
			s 	= end
			e 	= start
		else:
			s 	= start
			e	= end

		# make markup selection
		if abs(end-start) > 0 and zero == -1:
			s = s-1
			e = e-1
			senseText = self.dna[0:s] + '<span background="#0082ED">'+self.dna[s:e] + '</span>' + self.dna[e:]
			antiText = self.cdna[0:s] + '<span background="#0082ED">'+self.cdna[s:e] + '</span>' + self.cdna[e:]	
			
		elif abs(end-start) > 0 and zero == 1:
			senseText = '<span background="#0082ED">'+self.dna[0:s] + '</span>' + self.dna[s:e] + '<span background="#0082ED">'+self.dna[e:] + '</span>'
			antiText = '<span background="#0082ED">'+self.cdna[0:s] + '</span>' + self.cdna[s:e] + '<span background="#0082ED">'+self.cdna[e:] + '</span>'
		else:
			senseText = self.dna
			antiText  = self.cdna
		
		layout.set_markup(senseText)
		self.pango.update_layout(layout)
		self.pango.show_layout(layout)
		
		self.senseLayout = layout.copy() # copy sense layout
		
		# get the sense lines
		self.senseLines = layout.get_iter()
		
		#while self.senseLines.next_line():
		#	a,b = self.senseLines.get_line_extents()
		#	line = self.senseLines.get_line()
		
		# determine the average width of a char
		n = 0
		w = 0
		#self.lineLength =  self.senseLines.get_line().length
		self.lineLength =  layout.get_line(0).length
		
		x1, y1, w1, h1 	= self.senseLayout.index_to_pos(1)
		x2, y2, w2, h2 	= self.senseLayout.index_to_pos(1+self.lineLength)
		self.lineHeight = abs(y1 - y2) / pango.SCALE
		
		#print line.start_index, line.length
		#self.lineLength = line.length
		#lineWith = line.get_extents()[0][2]/pango.SCALE
		#print lineWith
		#while self.senseLines.next_char() and n < 20:
		#	w = w + self.senseLines.get_char_extents()[3]
		#	n = n + 1
		#self.averageWidth = w/(n*pango.SCALE)
		#self.nCharLine = lineWith/self.averageWidth
		#print  n, self.averageWidth, lineWith, self.nCharLine
		#print a, b
		
		
		
		
		# make them in usefull cairo coordinates like (start, end) each as (x,y)
		self.sensLineCoord = []
		
		#while self.senseLines.next_line():
		#	a,b = self.senseLines.get_line_extents()
		#	xa,ya,wa,ha = a
		#	xb,yb,wb,hb = b
		#	
		#	# scale to pixels
		#	xb = xb / pango.SCALE
		#	yb = yb / pango.SCALE
		#	wb = wb / pango.SCALE
		#	hb = hb / pango.SCALE				
		#	self.sensLineCoord.append(((xb, yb),(wb, hb)))#/pango.SCALE )

		# print anti sense strain
		self.ctx.move_to(8,starty + self.fontsize + self.lineGap)
		layout.set_markup(antiText)
		self.pango.update_layout(layout)
		self.pango.show_layout(layout)
		# get the lines
		#self.antiLines = layout.get_lines_readonly()
		
		# update the height of the context, based on the dna stuff
		w, h = layout.get_pixel_size()
		self.minHeight = h + self.textSpacing + 50

		# resize window
		w,h = self.GetSize()					# prevents missfomration on resize
		self.SetSize(wx.Size(w,self.minHeight)) # prevents missfomration on resize
		
		#cursor1s, cursor1w = layout.get_cursor_pos(4)
		#cursor2s, cursor2w = layout.get_cursor_pos(8)
		
		#layout.move_cursor_visually(True, 4, 0, 10)

		
		
		# make ticks:
		#set font
		font = pango.FontDescription("sans normal "+ "1")
		font.set_size(pango.SCALE * 8)
		layout.set_font_description(font)
		
		# new attr list for the ticks
		ticksAttr = pango.AttrList()
		# color
		gray 		= 110 * 65535/255
		fg_color	= pango.AttrForeground(gray,gray,gray, 0, len(self.dna))
		spacing 	= pango.AttrLetterSpacing(1500, 0, len(self.dna))
		# add attr
		ticksAttr.insert(fg_color)
		ticksAttr.insert(spacing)
		# set attr
		layout.set_attributes(ticksAttr)
		
		# number of lines:
		if len(self.dna) != 0:
			nLines = math.ceil((len(self.dna)-1) / self.lineLength)
		
			n = 0 # counter
			y = 8
			while n <= nLines:
				
				self.ctx.move_to(8, y)
				text = "%d" % (n *  self.lineLength + 1)
				layout.set_text(text)
				
				y = y + self.lineHeight
				
				n = n + 1
				self.pango.update_layout(layout)
				self.pango.show_layout(layout)
		
		
		
		
		

		return None
	
	def drawFeatures(self):
		features 	= genbank.gb.get_all_feature_positions()
		paths 		= {} # (color, path)
		
		if self.cairoStorage['features'] != features :
			# something changed, we have to draw the paths new
			for feature in features:
				 
				featuretype, complement, start, finish, name, index = feature
				
				if featuretype != "source":
					# hittets name
					hname 	= self.hitName(name, index)
					name	= self.nameBeautiful(name)
				
					xs, ys, ws, hs = self.senseLayout.index_to_pos(start)
					xf, yf, wf, hf = self.senseLayout.index_to_pos(finish)
				
					# if ys and yf are not the same, we have feature over multiple lines!
					dy 		= abs(ys - yf) / pango.SCALE
					ySteps 	= dy / self.lineHeight
					xs 		= xs / pango.SCALE
					xf 		= xf / pango.SCALE
					ys 		= ys / pango.SCALE + self.lineHeight + self.fontsize - (0.5 * self.textSpacing)
					yf		= yf / pango.SCALE
					
					# draw the name somewhere at the start
					nameheight = 9 # px
					self.ctx.select_font_face('Arial', cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
					self.ctx.set_font_size(nameheight)
					xbearing, ybearing, TextWidth, TextHeight, xadvance, yadvance = self.ctx.text_extents(name)
					self.ctx.move_to(xs + 10 , ys + nameheight + 2)
					self.ctx.text_path(name)
					textpath = self.ctx.copy_path()
					self.ctx.new_path()
					# end name
					
					# box for name
					self.ctx.move_to(xs, ys+10)
					self.ctx.rectangle(xs, ys, TextWidth + 20 , nameheight + 5)
					
					# lines below the dna
					while ySteps >= 0:
						if ySteps == 0:
							xE = xf
						else:
							xE = self.w

						# make the path
						self.ctx.move_to(xs, ys)
						self.ctx.line_to(xE, ys)
						
						xs = 8
						ys = ys + self.lineHeight
						ySteps = ySteps - 1
					
					# get color
					color = eval(featuretype)['fw'] #get the color of feature (as string)

					
					paths[hname] = (color, self.ctx.copy_path(),textpath)
					self.ctx.new_path() # clear canvas
			# save the paths
			self.cairoStorage['fPaths'] =paths
			# update storage status
			self.cairoStorage['features'] = features
			
		# draw the beautiful features
		
		for a in self.cairoStorage['fPaths']:
			color, path, tpath = self.cairoStorage['fPaths'][a]
			
			# get color
			assert type(color) == str
			r,g,b = colcol.hex_to_rgb(color)
			r = float(r)/255
			g = float(g)/255
			b = float(b)/255
			
			# append and draw path
			self.ctx.append_path(path)
			self.ctx.set_line_width(1)
			self.ctx.set_source_rgba(r ,g ,b , 1) # Solid color
			self.ctx.stroke_preserve()
			self.ctx.fill()
			# write text
			self.ctx.append_path(tpath)
			self.ctx.set_line_width(1)
			self.ctx.set_source_rgba(1 ,1 ,1 , 1) # Solid color
			self.ctx.fill()
			
			
			
		
			
			
		return None
		
	
	def drawCursor(self):
		''' draw and display cursor in the textfield'''
		indexO 		= self.cairoStorage['cursor']
		index 		= self.PositionPointer
		if indexO != index:	
			# turn index to x,y
			x, y,w,h 	= self.senseLayout.index_to_pos(index)
			x 			= x /pango.SCALE
			y 			= y / pango.SCALE
			y1 			= y + 0.25 * self.lineHeight
			y2 			= y + 0.75 * self.lineHeight
			self.ctx.move_to(x,y1)
			self.ctx.line_to(x,y2)
		
			self.cairoStorage['cursorPath'] = self.ctx.copy_path()
			self.ctx.new_path()
		
		# draw cursor:
		self.ctx.append_path(self.cairoStorage['cursorPath'])
		self.ctx.set_line_width(0.3)
		self.ctx.set_source_rgba(0,0,0, 1) # Solid color
		self.ctx.stroke()
	
######################################################

	def OnSize(self, event):
		# remove buffered features, to force redraw:
		self.cairoStorage['features'] = None
		
		# The Buffer init is done here, to make sure the buffer is always
		# the same size as the Window
		#Size  = self.GetClientSizeTuple()
		w, h = self.GetSize()
		Size = wx.Size(w,self.minHeight)
		self.SetSizeWH(w,self.minHeight) # 
		#self.SetSize(Size)
		self.parent.SetVirtualSize(Size)
		#self.SetSize(Size)



		
		self.size  = self.parent.GetVirtualSize()

		# Make new offscreen bitmap: this bitmap will always have the
		# current drawing in it, so it can be used to save the image to
		# a file, or whatever.
		w, h = self.size

		#self.SetSize(Size)		
		self._Buffer = wx.EmptyBitmap(w, h)
		#self.SetSize(Size)
		self.update_ownUI()	
		self.SetSize(Size)  # must be here, so after ui update the size adjudust again correctly
		
	
	def OnLeftDown(self, event):
		# remove selection:
		self.set_dna_selection()
		
		# get position
		x, y 				= self.ScreenToClient(wx.GetMousePosition())
		xp 					= x * pango.SCALE
		yp 					= int(y - 0.5 * self.lineHeight) * pango.SCALE 
		self.LeftDownIndex 	= None
		index 				= self.senseLayout.xy_to_index(xp,yp)
		
		if index != False:
			self.LeftDownIndex = index[0] 
		
		# update cursor
		self.PositionPointer = index[0] 
			
		
		event.Skip() #very important to make the event propagate and fulfill its original function

	def OnLeftUp(self, event):
		#self.set_dna_selection() #update the varable keeping track of DNA selection
		#self.set_cursor_position() # update cursor position
		self.update_ownUI()
		#self.update_globalUI()
		event.Skip() #very important to make the event propagate and fulfill its original function

	def OnMotion(self, event):
		if event.Dragging() and event.LeftIsDown():
			oldSel = genbank.dna_selection
			
			# get position
			x, y 				= self.ScreenToClient(wx.GetMousePosition())
			xp 					= x * pango.SCALE
			yp 					= (y - 20) * pango.SCALE  # TODO replace 20 with something calculated
			self.LeftMotionIndex 	= None
			index 				=  self.senseLayout.xy_to_index(xp,yp)
			if index != False:
				self.LeftMotionIndex = index[0]
				self.PositionPointer = index[0] 
			
			if self.LeftDownIndex < self.LeftMotionIndex:
				selection = (self.LeftDownIndex, self.LeftMotionIndex, -1)
			else:
				selection = (self.LeftMotionIndex, self.LeftDownIndex, -1)
			
			if oldSel != selection:
				self.set_dna_selection(selection) #update the varable keeping track of DNA selection
				
			

			
		event.Skip() #very important to make the event propagate and fulfill its original function


	def OnRightUp(self, event):
		event.Skip() #very important to make the event propagate and fulfill its original function


	#### rendering


####################################################################
	def hitName(self, name, index):
		# unique and reproducible id for each feature
		name 	= self.nameBeautiful(name)
		seq 	= "%s%s" % (name, index)
		name 	= ''.join(seq.split())
		return name
	
	def nameBeautiful(self, name):
		# remove any '"' in front or at the end
		if name[0:1] == '"':
			name = name[1:]
		if name[-1:] == '"':
			name = name[:-1]
		return name
	
	def set_cursor_position(self):
		# set position of cursor for status bar:
		genbank.cursor_position = self.PositionPointer+1 # +1 because we have no 0 position, A|CGA is pos 2


	def set_dna_selection(self, selection = (0,0, -1)):
		'''
		Updates DNA selection.
		'''
		##selection = self.get_selection()
		a, b, none = selection			# new to enable selections over 0
		selection = (a, b , -1)
		genbank.dna_selection = selection
		self.update_globalUI()
		self.update_ownUI()

	def moveCursor(self, offset, linejump = False):
		### handle cursor movement:
		# offset +-int
		# linejump False/up/down
		
		index = self.PositionPointer
		topology = genbank.gb.gbfile['locus']['topology']
		
		length = len(genbank.gb.GetDNA())

		if linejump == False:
			index = index + offset
		elif linejump == "up":
			# jump to the last row, but keep index to the left: #TODO, not correctly behaving if jumping between start and end
			index = index - self.lineLength
		elif linejump == "down":
			index = index + self.lineLength
		if topology == 'circular':
			if index < 1:
				index = index + length
			elif index > length:
				index = index - length
		else:
			if index < 1:
				index = 1
			elif index > length - 1:
				index = length - 1

		
		self.PositionPointer = index
		self.set_cursor_position()
		
	def OnKeyPress(self, evt):
		key = evt.GetUniChar()

		index = self.PositionPointer
		
		shift = evt.ShiftDown() # is shift down?
		
		if key in [97, 65, 116, 84, 99, 67, 103, 71]: # [a, A, t, T, c, C, g, G']
			start, finish, null = self.get_selection()
			if start != finish: # if a selection, delete it, then paste
				genbank.gb.Delete(start+1, finish, visible=False)
				
			if shift == True:
				genbank.gb.Paste(index, chr(key).upper())
			elif shift == False:
				genbank.gb.Paste(index, chr(key).lower())
			
			# move pointer
			self.PositionPointer = self.PositionPointer + 1
			self.update_ownUI()
			self.update_globalUI()
			#self.stc.SetSelection(start+1, start+1)

		
		elif key == 8: #backspace
			start, finish, null = self.get_selection()
			if start != finish: # if a selection, delete it
				genbank.gb.Delete(start+1, finish)
				self.set_dna_selection() # remove selection
				self.update_ownUI()
				self.update_globalUI()

			else:
				genbank.gb.Delete(index-1, index-1)
				self.PositionPointer = self.PositionPointer - 1 # set pointer back
				self.update_ownUI()
				self.update_globalUI()

				#self.stc.SetSelection(start-1, start-1)
		
		elif key == 127: #delete
			start, finish, null = self.get_selection()
			if start != finish: # if a selection, delete it
				genbank.gb.Delete(start+1, finish)
			else:
				genbank.gb.Delete(index, index)
			self.set_dna_selection() # remove selection
			self.update_ownUI()
			self.update_globalUI()

		elif key == 314 and shift == False: #left
			self.moveCursor(-1)		
			self.set_dna_selection() # remove selection
			self.update_ownUI()
			self.update_globalUI()

		elif key == 314 and shift == True: #left + select
			start, finish, null = self.get_selection()
			if start == 0:
				start = self.PositionPointer
			if finish == self.PositionPointer:
				self.moveCursor(-1)
				finish = self.PositionPointer
			else:
				self.moveCursor(-1)
				start = self.PositionPointer
			self.set_dna_selection((start, finish, null)) # remove selection
			self.update_ownUI()
			self.update_globalUI()


		elif key == 316 and shift == False: #right
			self.moveCursor(1)#self.PositionPointer = self.PositionPointer + 1
			self.set_dna_selection() # remove selection
			self.update_ownUI()
			self.update_globalUI()

		elif key == 316 and shift == True: #right + select
			start, finish, null = self.get_selection()
			if start == 0:
				start = self.PositionPointer
			
			if finish == self.PositionPointer:
				self.moveCursor(1)
				finish = self.PositionPointer
			else:
				self.moveCursor(1)
				start = self.PositionPointer
			self.set_dna_selection((start, finish, null)) # remove selection
			self.update_ownUI()
			self.update_globalUI()


		elif key == 315 and shift == False: #up
			self.moveCursor(1, "up")
			self.update_ownUI()
			self.update_globalUI()

		elif key == 315 and shift == True: #up + select
			print "select with keyboard"

		elif key == 317 and shift == False: #down
			self.moveCursor(1, "down")
			self.update_ownUI()
			self.update_globalUI()

		elif key == 317 and shift == True: #down + select
			print "select with keyboard"

		#self.set_dna_selection() #update the varable keeping track of DNA selection

		self.set_cursor_position() # udpate cursor position, just in case

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
		'''start, finish = self.stc.GetSelection()
		if start == finish: #not a selection
			finish = 0
		elif start > finish: #the selection was made backwards
			start, finish = finish, start

		selection = (start+1, finish)'''
		
		return genbank.dna_selection

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

			mposition = self.stc.CharPositionFromPoint(xposition, yposition)


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

	genbank.dna_selection = (0, 0, -1)	 #variable for storing current DNA selection
	genbank.feature_selection = False #variable for storing current feature selection
	genbank.search_hits = []

	import sys
	assert len(sys.argv) == 2, 'Error, this script requires a path to a genbank file as an argument.'
	print('Opening %s' % str(sys.argv[1]))

	genbank.gb = genbank.gbobject(str(sys.argv[1])) #make a genbank object and read file


	app = MyApp(0)
	app.MainLoop()
