#!/usr/bin/env python


#DNApy is a DNA editor written purely in python.
#The program is intended to be an intuitive and fully featured 
#editor for molecular and synthetic biology.
#Enjoy!
#
#copyright (C) 2014-2015  Martin Engqvist |
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#LICENSE:
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
#Get source code at: https://github.com/mengqvist/DNApy
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

import featureedit_GUI

# drawing
import cairo
from wx.lib.wxcairo import ContextFromDC
import pango # for text
import pangocairo # for glyphs and text
import math

# for color function
import colcol

# for hashing the enzymes to detect change in enzymes
import hashlib


files	=	{}   #list with all configuration files
files['default_dir'] 	= os.path.abspath(os.path.dirname(sys.argv[0]))+"/"
files['default_dir']	= string.replace(files['default_dir'], "\\", "/")
files['default_dir']	= string.replace(files['default_dir'], "library.zip", "")
settings=files['default_dir']+"settings"   ##path to the file of the global settings
exec(compile(open(settings).read(), settings, 'exec')) #gets all the pre-assigned settings


########### class for text ######################
class TextEdit(DNApyBaseDrawingClass):
	'''
		new Text editor that allows to display dna and
		complementary dna in a usefull manner
		
		Also it can and should show:
			- features (done, buggy)
			- restriction enzymes (partly)
			- dna modifications (todo)
			- primer (todo)
	'''
	def __init__(self, parent, id):
		
		#parent
		self.parent = parent
		
		# get data
		self.dna 			= genbank.gb.GetDNA()
		self.cdna 			= genbank.gb.GetDNA()

		# settings for the editor 
		self.sY				= 60	# offset top
		self.dX				= 8 	# border left and right
		self.textSpacing 	= 80	# seperation of lines
		self.fontsize 		= 9		# line size
		self.lineGap 		= 6
		
		# interactive cursor
		self.PositionPointer 		= 1
		self.PositionPointerEnd 	= False
		
		# initial minimal height of the editor 
		self.minHeight = 200
		
		# storage dict, to prevent recalculating to much
		self.cairoStorage = {
				'features'	: None, # feature Object
				'fPaths' 	: [], 	# path of alls features as List
				'cursorPath': None,
				'cursor'	: None,
				'enzymes'	: {},	# enzymes object
				'ePaths'	: [],	# enzymes drawn
				'ticks'		: None,	# ticks object --> self.dna
				'tPaths'	: []	# ticks drawn		
		}
		
		
		
		# call window
		DNApyBaseDrawingClass.__init__(self, parent, wx.ID_ANY)

		# bind events
		self.Bind(wx.EVT_LEFT_DCLICK, self.OnLeftDouble)
		self.Bind(wx.EVT_LEFT_DOWN, self.OnLeftDown)
		self.Bind(wx.EVT_LEFT_UP, self.OnLeftUp)
		self.Bind(wx.EVT_MOTION, self.OnMotion) # drag and drop, selection etc.
		self.Bind(wx.EVT_KEY_DOWN, self.OnKeyPress) #This is important for controlling the input into the editor
		

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

		#####################
		# remove buffered features, to force redraw:
		self.cairoStorage['features'] = None
		
		# The Buffer init is done here, to make sure the buffer is always
		# the same size as the Window
		wC, hC = self.parent.GetClientSize()

		wP, hP = self.parent.GetSize()
		self.SetSize(wx.Size(wC, self.minHeight))
		self.parent.SetVirtualSize(wx.Size(wC-30, self.minHeight))
		

		# Make new offscreen bitmap: this bitmap will always have the
		# current drawing in it, so it can be used to save the image to
		# a file, or whatever.
#		if wC != 0 and hC !=0:
#			self._Buffer = wx.EmptyBitmap(wC, self.minHeight)

#		else:
#			self._Buffer = wx.EmptyBitmap(wP, hP)
				
		##########################

		# start dc
		dc = wx.MemoryDC()
		dc.SelectObject(self._Buffer)
		dc.SetBackground(wx.Brush("White"))
		dc.Clear() 
		# make a cairo context
		self.ctx = ContextFromDC(dc)		
		# start pango as a font backend
		self.pango = pangocairo.CairoContext(self.ctx)
		self.pango.set_antialias(cairo.ANTIALIAS_SUBPIXEL)
		
		# reload dna
		self.dna = genbank.gb.GetDNA()
		if self.dna != None:
			self.cdna 	= dna.C(self.dna)
		else:
			self.dna 	= ''	# no dna, if there is nothing
			self.cdna 	= ''	# empty cdna
			
		# render the text
		self.displayText()

		# create colorful features here
		self.drawFeatures()
		
		# draw enzymes
		self.drawEnzymes()
		
		# draw a cursor
		self.drawCursor()
		
		# draw Tics above all!
		self.drawTicks()
	
		
		# end of canvas handling
		dc.SelectObject(wx.NullBitmap) # need to get rid of the MemoryDC before Update() is called.
		self.Refresh()
		self.Update()

		return None
	
	
		
	
	def displayText(self):
		''' main function for the editors text capabilities
			here the text gets positioned and all its color.
			Also the selection highlight is done here'''
		
		# get window size to wrap text
		width, height 	= self.parent.GetVirtualSize()
		if width > 30:
			w = width - self.dX
		else:
			w = width
		# save with for other functions
		self.w = w

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
		
		# font and color are prepared for the DNA

		
		
		# get selection, if any
		start, end, zero = self.get_selection()
		if end != -1:
			# correct direktion of selection
			if start > end:
				s 	= end - 1
				e 	= start
			else:
				s 	= start - 1
				e	= end
		# make markup selection
		if abs(end-start) >= 0 and zero == -1 and end != -1 :
			senseText = self.dna[0:s] + '<span background="#0082ED">'+self.dna[s:e] + '</span>' + self.dna[e:]
			antiText = self.cdna[0:s] + '<span background="#0082ED">'+self.cdna[s:e] + '</span>' + self.cdna[e:]	
			
		elif abs(end-start) >= 0 and zero == 1 and end != -1:
			senseText = '<span background="#0082ED">'+self.dna[0:s] + '</span>' + self.dna[s:e] + '<span background="#0082ED">'+self.dna[e:] + '</span>'
			antiText = '<span background="#0082ED">'+self.cdna[0:s] + '</span>' + self.cdna[s:e] + '<span background="#0082ED">'+self.cdna[e:] + '</span>'
		else:
			senseText = self.dna
			antiText  = self.cdna
		
		# output dna
		self.ctx.move_to(self.dX,self.sY)
		layout.set_markup(senseText)
		self.pango.update_layout(layout)
		self.pango.show_layout(layout)
		
		# copy sense layout
		self.senseLayout = layout.copy() 

		
		# determine the length and height of a line
		self.lineLength =  layout.get_line(0).length
		x1, y1, w1, h1 	= self.senseLayout.index_to_pos(1)
		x2, y2, w2, h2 	= self.senseLayout.index_to_pos(1+self.lineLength)
		self.lineHeight = abs(y1 - y2) / pango.SCALE
		
		# print anti sense strain
		self.ctx.move_to(self.dX,self.sY + self.fontsize + self.lineGap)
		layout.set_markup(antiText)
		self.pango.update_layout(layout)
		self.pango.show_layout(layout)
		
		# update the height of the context, based on the inserted dna
		w, h = layout.get_pixel_size()
		self.minHeight = h + 2*self.sY
		# resize window
#		w,h = self.GetSize()					# prevents missfomration on resize
#		self.SetSize(wx.Size(w,self.minHeight)) # prevents missfomration on resize
		
		return None
	
	
	def drawFeatures(self):
		features 	= genbank.gb.get_all_feature_positions()
		paths 		= {} # (color, path)
		
		featureheight = 10
		
		if self.cairoStorage['features'] != features :
			# something changed, we have to draw the paths new
			self.drawn_fw_locations = []
			self.drawn_rv_locations = []
			
			# font settings
			# draw the name somewhere at the start
			nameheight = 9 # px
			self.ctx.select_font_face('Arial', cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
			self.ctx.set_font_size(nameheight)
			
			for feature in features:
				 
				featuretype, complement, start, finish, name, index = feature
				
				if featuretype != "source":
				
					# get names
					hname 	= self.hitName(name, index)
					name	= self.nameBeautiful(name)
					#get position
					xs, ys, ws, hs = self.senseLayout.index_to_pos(start)
					xf, yf, wf, hf = self.senseLayout.index_to_pos(finish+1)
					# linecount:
					# if ys and yf are not the same, we have feature over multiple lines!
					dy 		= abs(ys - yf) / pango.SCALE
					ySteps 	= dy / self.lineHeight # number of lines!
					xs 		= xs / pango.SCALE
					xf 		= xf / pango.SCALE
					ys		= (ys / pango.SCALE) + self.sY 
					yf		= yf / pango.SCALE
					
					# first sense strain over the text
					if complement == False:
						self.drawn_fw_locations, l = self.find_overlap(self.drawn_fw_locations, (start, finish))				
						level = l * (featureheight+6)
					

						xbearing, ybearing, TextWidth, TextHeight, xadvance, yadvance = self.ctx.text_extents(name)				
						ys = ys  - level - 5
						self.ctx.move_to(xs + 10, ys )
						self.ctx.text_path(name)
						textpath = self.ctx.copy_path()
						self.ctx.new_path()
						# end name
					
						# box for name
						self.ctx.move_to(xs, ys)
						self.ctx.rectangle(xs, ys-featureheight, TextWidth + 20 , nameheight + 5)
					
						# lines above the dna
						while ySteps >= 0:
							if ySteps == 0:
								xE = xf
							else:
								xE = self.w + 5 # TODO: figure out, why we need to add 5 px, to make it look correct
 
							# make a box-path
							self.ctx.move_to(xs, ys)
							self.ctx.rectangle(xs, ys-5, xE-xs , 5)
							
							# the second x start is always 8px, or what we will opt for
							xs = self.dX
							ys = ys + self.lineHeight
							ySteps = ySteps - 1
						
						# path
						fillpath = self.ctx.copy_path()
						self.ctx.new_path()
						# get color
						color = eval(featuretype)['fw'] #get the color of feature (as string)
					else:
						self.drawn_rv_locations, l = self.find_overlap(self.drawn_rv_locations, (start, finish))				
						level = l * (featureheight+6)
						# set y
						ys = ys + 2*self.fontsize + self.lineGap + 9 + level 

						# lines below the dna
						while ySteps >= 0:
							if ySteps == 0:
								xE = xf
							else:
								xE = self.w

							# make the path
							#self.ctx.move_to(xs, ys)
							#self.ctx.line_to(xE, ys)
							self.ctx.move_to(xs, ys)
							self.ctx.rectangle(xs, ys+5, xE-xs , 5)
						
							xs = self.dX
							ys = ys + self.lineHeight
							ySteps = ySteps - 1
						
						# change y:
						ys = ys - self.lineHeight
						# buffer path
						bufferPath = self.ctx.copy_path()
						self.ctx.new_path()
						
						# name at the end:
						xbearing, ybearing, TextWidth, TextHeight, xadvance, yadvance = self.ctx.text_extents(name)	
						self.ctx.move_to(xf - TextWidth - 10, ys + 10 )
						self.ctx.text_path(name)
						textpath = self.ctx.copy_path()
						self.ctx.new_path()
						# end name

						# box for name
						self.ctx.append_path(bufferPath)
						self.ctx.move_to(xf - TextWidth, ys + featureheight)
						self.ctx.rectangle(xf - TextWidth - 20, ys, TextWidth + 20 , nameheight + 5)

						# get path
						fillpath = self.ctx.copy_path()
						self.ctx.new_path()
						# get color
						color = eval(featuretype)['fw'] #get the color of feature (as string)

					
					paths[hname] = (color,fillpath ,textpath)
					self.ctx.new_path() # clear canvas
			# save the paths
			self.cairoStorage['fPaths'] =	paths
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
			#self.ctx.stroke_preserve()
			self.ctx.fill()
			# write text
			self.ctx.append_path(tpath)
			self.ctx.set_line_width(1)
			self.ctx.set_source_rgba(1 ,1 ,1 , 1) # Solid color
			self.ctx.fill()

		return None
	
	
	def drawEnzymes(self):
		# load enzymes
		enzymes 	=  genbank.gb.restrictionEnzymes.selection
		# create hash of selection and DNA
		hashable 	= "%s%s" % (frozenset(enzymes), genbank.gb.gbfile['dna'])
		e1 			= hashlib.md5(hashable).hexdigest()
		
		if e1 != self.cairoStorage['enzymes']:
			nameheight = 9 # px
			self.ctx.select_font_face('Arial', cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
			self.ctx.set_font_size(nameheight)
			dna 			= self.dna
			length			= len(dna)
			index			= 0	# same as in plasmid_GUI.py
			enzymesPaths 	= {}
			# loop enzymes ordered dict
			for e in enzymes:
				for site in enzymes[e].restrictionSites:
					name, start, end, cut51, cut52, dnaMatch = site
					startRef 	= start - 1
					hitname 	= self.hitName(name, index)
					
					# get position to draw
					xs, ys, ws, hs = self.senseLayout.index_to_pos(startRef)
					xs = xs / pango.SCALE + self.dX
					ys = ys / pango.SCALE + self.sY
					hs = hs / pango.SCALE
					# move there, but above the line
					y = ys# - hs
					self.ctx.move_to(xs,y)
					# print text as path
					self.ctx.text_path(name)
					textpath = self.ctx.copy_path()
					self.ctx.new_path()
					# make a line, representing the cut
					xC51, yC51, wC51, hC51 	= self.senseLayout.index_to_pos(startRef + enzymes[e].c51)
					xC31, yC31, wC31, hC31 	= self.senseLayout.index_to_pos(startRef + enzymes[e].c31)
					xC51 					= xC51 / pango.SCALE + self.dX
					yC51 					= yC51 / pango.SCALE + self.sY
					tHeight 				= hC51 / pango.SCALE
					xC31 					= xC31 / pango.SCALE +  self.dX
					yC31 					= yC31 / pango.SCALE + self.sY
					
					if yC51 == yC31:
						self.ctx.move_to(xC51, yC51)
						self.ctx.line_to(xC51, yC51 + tHeight + 0.25 * self.lineGap)
						self.ctx.line_to(xC31, yC31 + tHeight + 0.25 * self.lineGap)
						self.ctx.line_to(xC31, yC31 + 2 * tHeight + 0.5* self.lineGap)
					else:
						if yC51 <= yC31:
							lineend 	= self.dX + self.w
							linestart 	= self.dX
						else:
							# switch
							linestart	= self.dX + self.w
							lineend 	= self.dX
							
						# split the drawing in two parts
						self.ctx.move_to(xC51, yC51)
						self.ctx.line_to(xC51, yC51 + tHeight + 0.25 * self.lineGap)
						self.ctx.line_to(lineend, yC51 + tHeight + 0.25 * self.lineGap)
						# part 2
						self.ctx.move_to(linestart, yC31 + tHeight + 0.25 * self.lineGap)
						self.ctx.line_to(xC31, yC31 + tHeight + 0.25 * self.lineGap)
						self.ctx.line_to(xC31, yC31 + tHeight + 0.25 * self.lineGap)
						self.ctx.line_to(xC31, yC31 + 2 * tHeight + 0.5* self.lineGap)
						
					linepath 				= self.ctx.copy_path()
					self.ctx.new_path()
					enzymesPaths[hitname] 	= (textpath, linepath)
					
					# raise index
					index = index + 1
			# update storage
			self.cairoStorage['enzymes'] 	= e1				
			self.cairoStorage['ePaths'] 	= enzymesPaths
		
		# draw enzymes:
		self.ctx.set_line_width(1)
		self.ctx.set_source_rgba(0.6,0.3,0.3 , 1) # Solid color
		for e in self.cairoStorage['ePaths']:
			path, pathc = self.cairoStorage['ePaths'][e]
			self.ctx.append_path(path)
			self.ctx.fill()
			self.ctx.append_path(pathc)
			self.ctx.stroke()
	
	def drawTicks(self):
		''' draw the nice infos about the dna position in each line'''
		nameheight = 9 # px
		self.ctx.set_font_size(nameheight)
		self.ctx.select_font_face('Arial', cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
		# number of lines:
		if len(self.dna) != 0 and self.cairoStorage['ticks'] != self.dna:
			# count the lines:
			nLines = math.ceil((len(self.dna)-1) / self.lineLength)
		
			# empty paths:
			self.cairoStorage['tPaths'] = []
			
			n = 0 # counter
			y =  self.sY  -2 #- self.fontsize - 2
			
			# clear paths
			self.ctx.new_path()
			while n <= nLines:
				
				self.ctx.move_to(self.dX + 3, y)
				# make the tick
				text = "%d" % (n *  self.lineLength + 1)
				xbearing, ybearing, TextWidth, TextHeight, xadvance, yadvance = self.ctx.text_extents(text)	
				self.ctx.text_path(text)
				textpath = self.ctx.copy_path()
				self.ctx.new_path()
				# make a box
				boxPath = self.carioRoundBox(self.dX,y  - TextHeight - 2, TextWidth + 6, TextHeight + 4, 0.6)
				self.ctx.new_path()
				
				#save paths
				self.cairoStorage['tPaths'].append((textpath, boxPath))
				
				# raise counter
				y = y + self.lineHeight
				
				n = n + 1
				
			# set storage:
			self.cairoStorage['ticks'] = self.dna
			

		# print the ticks:
		self.ctx.set_line_width(0)

		for i in self.cairoStorage['tPaths']:
			textpath, boxPath = i
			# draw box
			self.ctx.set_source_rgba(1,1,1, 0.6) # white, transparent
			self.ctx.append_path(boxPath)
			self.ctx.fill()
			# draw position
			self.ctx.set_source_rgba(0,0,0,1) # dark black
			self.ctx.append_path(textpath)
			self.ctx.fill()

		return None
	
	def carioRoundBox(self, x, y, width=10, height=10, aspect=0.8):
		''' draw a box with rounded corners '''
		# from http://cairographics.org/samples/rounded_rectangle/
		corner_radius = height / 10.0   

		radius = corner_radius / aspect
		degrees = math.pi / 180.0

		self.ctx.new_sub_path()
		self.ctx.arc( x + width - radius, y + radius, radius, -90 * degrees, 0 * degrees)
		self.ctx.arc( x + width - radius, y + height - radius, radius, 0 * degrees, 90 * degrees)
		self.ctx.arc(x + radius, y + height - radius, radius, 90 * degrees, 180 * degrees)
		self.ctx.arc(x + radius, y + radius, radius, 180 * degrees, 270 * degrees)
		self.ctx.close_path()

		path = self.ctx.copy_path()
		self.ctx.new_path()
		return path

	def drawCursor(self):
		''' draw and display cursor in the textfield'''
		indexO 		= self.cairoStorage['cursor']
		index 		= self.PositionPointer
		pos1, pos2, zero = self.get_selection()
		if indexO != index and pos2 == -1:
			index = index - 1 # we save the genetic curors 1-len(), but mean the cursor 0-len()-1
			
			# turn index to x,y
			x, y, w, h 	= self.senseLayout.index_to_pos(index)
			x 			= x / pango.SCALE + self.dX
			y 			= (y / pango.SCALE) + self.sY
			y1 			= y 
			y2 			= y + 0.5 * self.lineHeight
			
			# move cursor if its the line end:
			if self.PositionPointerEnd == True:
				x 		= x + self.dX
			self.ctx.move_to(x,y1)
			self.ctx.line_to(x,y2)
		
			self.cairoStorage['cursorPath'] = self.ctx.copy_path()
			self.ctx.new_path()
		elif pos2 != -1:
			# there is a selection, so no cursor:
			self.cairoStorage['cursorPath'] = None
		
		if self.cairoStorage['cursorPath'] != None:
			# draw cursor:
			self.ctx.append_path(self.cairoStorage['cursorPath'])
			self.ctx.set_line_width(0.3)
			self.ctx.set_source_rgba(0,0,0, 1) # Solid color
			self.ctx.stroke()
	
######################################################



		
	
	def OnLeftDown(self, event):
		# remove selection:
		self.set_dna_selection()
		
		
		# get position
		x, y 				= self.ScreenToClient(wx.GetMousePosition())
		xp 					= (x - self.dX) * pango.SCALE
		yp 					= int(y - self.sY - self.fontsize) * pango.SCALE 
		self.LeftDownIndex 	= None
		index 				= self.senseLayout.xy_to_index(xp,yp)
		
		if index != False:
			# click can be left or right of the center of the char
			if index[1] == 0:
				self.LeftDownIndex = index[0] + 1 
			else:
				self.LeftDownIndex = index[0] + 2
		
		# update cursor
		self.PositionPointer = self.LeftDownIndex
		self.set_cursor_position()
			
		
		event.Skip() #very important to make the event propagate and fulfill its original function

	def OnLeftUp(self, event):
		# get selection
		pos1, pos2, zero 	= self.get_selection()
		if pos2 == -1:
			# reset selection
			#self.set_dna_selection()

			# maybe we cliked on a feature, label or a line?
			hit = self.HitTest()
			if hit != None:
				self.hit = hit
				# we hit one, so we need to draw the selection!
				# wich feature was hit?
				featurelist             = genbank.gb.get_all_feature_positions()
				for i in range(0,len(featurelist)):
					# load all the feature infos
					featuretype, complement, start, finish, name, index = featurelist[i]
					hName = self.hitName(name, index)
					if hit == hName:
						# a feature with start smaller then finish, is a featur laying over zero
						if start > finish:
							zero = 1 # 1 --> feature starts left and ends right of +1
						# set selection 
						self.set_dna_selection((start,finish, zero))
						self.set_cursor_position(start)
			else:
				self.hit = None
		
		# update the UI maybe?
		self.update_ownUI()
		self.update_globalUI()
		event.Skip() #very important to make the event propagate and fulfill its original function
	
	def OnLeftDouble(self, event):
		'''doubleclick aktion'''
		# remove selection:
		self.set_dna_selection()
		# if we hit a feature, wen can open a dialog
		hit = self.HitTest()
		if hit != None:
			# wich feature was hit?
			featurelist             = genbank.gb.get_all_feature_positions()
			for i in range(0,len(featurelist)):
				# load all the feature infos
				featuretype, complement, start, finish, name, index = featurelist[i]
				hName = self.hitName(name, index)
				if hit == hName:
					genbank.feature_selection = copy.copy(index)
					dlg = featureedit_GUI.FeatureEditDialog(None, 'Edit Feature') # creation of a dialog with a title'''
					dlg.ShowModal()
					dlg.Center()
					# update the UI maybe?
					self.update_ownUI()
					self.update_globalUI()

	
	def OnMotion(self, event):
		if event.Dragging() and event.LeftIsDown():
			oldSel = genbank.dna_selection
			
			# get position
			x, y 				= self.ScreenToClient(wx.GetMousePosition())
			xp 					= (x - self.dX) * pango.SCALE
			yp 					= int(y - self.sY - self.fontsize) * pango.SCALE 
			self.LeftMotionIndex 	= None
			index 				=  self.senseLayout.xy_to_index(xp,yp)
			# get motion cursor
			if index != False:
				if index[1] == 0:
					self.LeftMotionIndex = index[0] + 1
				else:
					self.LeftMotionIndex = index[0] + 2 


			# update cursor
			if self.LeftMotionIndex >= self.LeftDownIndex:
				self.PositionPointer = self.LeftMotionIndex + 1
			else:
				self.PositionPointer = self.LeftMotionIndex 
			
			# end of line excepions
			if (self.PositionPointer - 1)  % self.lineLength == 0 and self.PositionPointer >= self.lineLength :
				self.PositionPointerEnd = True
				self.PositionPointer = self.PositionPointer -1
			else:
				self.PositionPointerEnd = False
				
			# get selection
			if self.LeftDownIndex <= self.LeftMotionIndex:
				selection = (self.LeftDownIndex, self.LeftMotionIndex, -1)
			else:
				selection = (self.LeftMotionIndex , self.LeftDownIndex , -1)
			
			if oldSel != selection:
				self.set_dna_selection(selection) #update the varable keeping track of DNA selection
				self.update_ownUI()
				self.update_globalUI()
				
			

			
		event.Skip() #very important to make the event propagate and fulfill its original function


#	def OnRightUp(self, event):
#		event.Skip() #very important to make the event propagate and fulfill its original function


	#### rendering


####################################################################
	def find_overlap(self, drawn_locations, new_range):
		'''
		Takes two ranges and determines whether the new range has overlaps with the old one.
		If there are overlaps the overlap locations are returned.
		This is used when drawing features. If two features overlap I want them drawn on different levels.
		'''
		assert type(drawn_locations) == list
		assert type(new_range) == tuple

		if drawn_locations == []:
			drawn_locations.append([new_range])
			return drawn_locations, 0
		else:
			i = 0
			while i < len(drawn_locations):
				overlap_found = False
				for n in range(0,len(drawn_locations[i])):
					if drawn_locations[i][n][0]<=new_range[0]<=drawn_locations[i][n][1] or drawn_locations[i][n][0]<=new_range[1]<=drawn_locations[i][n][1]: #if they overlap
						overlap_found = True
					elif new_range[0]<=drawn_locations[i][n][0]<=new_range[1] or new_range[0]<=drawn_locations[i][n][1]<=new_range[1]: #if they overlap
						overlap_found = True
				if overlap_found == False:
					drawn_locations[i].append(new_range)
					return drawn_locations, i
					break
				elif i+1==len(drawn_locations):
					drawn_locations.append([new_range])
					return drawn_locations, i+1
					break
				i += 1
				
	def hitName(self, name, index):
		# unique and reproducible id for each feature
		name 	= self.nameBeautiful(name)
		seq 	= "%s%s" % (name, index)
		name 	= ''.join(seq.split())
		return name
		
		
	def HitTest(self):
		'''Tests whether the mouse is over any feature or a enzyme'''
		hit 	= None


		# get the mouse positions
		x, y 	= self.ScreenToClient(wx.GetMousePosition())
		x2, y2 	= self.ctx.device_to_user(x,y)

		
		# list of all the paths:
		loop = [self.cairoStorage['fPaths'], 	# features
				self.cairoStorage['ePaths']]	# enzymes
		# loop over all possible paths to find the one we have under the mouse
		for paths in loop:
			for i in paths:
				path = paths[i]
				if type(path) == list or type(path) == tuple:
					path = path[1]
				
				# load the path
				self.ctx.append_path(path)
				self.ctx.set_line_width(0.1) # reduce linewith, so we can have nice clicking detection
				inFill 		= self.ctx.in_fill(x2,y2)
				inStroke 	= self.ctx.in_stroke(x2,y2)
				# check if this path is hit
				if inFill == True or inStroke == True:
					hit = i
				self.ctx.new_path()

		return hit
	
	def nameBeautiful(self, name):
		# remove any '"' in front or at the end
		if name[0:1] == '"':
			name = name[1:]
		if name[-1:] == '"':
			name = name[:-1]
		return name
	
	def set_cursor_position(self, pointer=False):
		if pointer == False:
			pointer = self.PositionPointer
		else:
			self.PositionPointer = pointer
		# set position of cursor for status bar:
		genbank.cursor_position = pointer


	def set_dna_selection(self, selection = (1,-1, -1)): # standart selection
		'''
		Updates DNA selection.
		'''
		##selection = self.get_selection()
		a, b, none = selection			# new to enable selections over 0
		# if start after end, swap them
		if a > b and b != -1:
			a,b = b,a
		selection = (a, b , -1)
		genbank.dna_selection = selection
		#self.update_globalUI()
		#self.update_ownUI()

	def moveCursor(self, offset, linejump = False):
		### handle cursor movement:
		# offset +-int
		# linejump False/up/down
		
		index = self.PositionPointer
		topology = genbank.gb.gbfile['locus']['topology']
		
		length = len(genbank.gb.GetDNA())

		if linejump == False:
				# check if we moved to a Line end
			# if so, we first do not increase the counter:
			if self.PositionPointerEnd == False and index % self.lineLength == 0:
				self.PositionPointerEnd = True
			elif self.PositionPointerEnd == True and index % self.lineLength == 0:
				self.PositionPointerEnd = False
				index = index + offset
			else:
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
			if finish != -1: # if a selection, delete it, then paste
				genbank.gb.Delete(start+1, finish, visible=False)
				self.PositionPointer = start
			if shift == True:
				genbank.gb.Paste(index, chr(key).upper())
			elif shift == False:
				genbank.gb.Paste(index, chr(key).lower())
			
			# move pointer
			self.PositionPointer = self.PositionPointer + 1
			
			if finish != -1:
				self.set_dna_selection()
			self.update_ownUI()
			self.update_globalUI()
			
			


		
		elif key == 8: #backspace
			start, finish, null = self.get_selection()
			if finish != -1: # if a selection, delete it
				genbank.gb.Delete(start+1, finish)
				self.set_dna_selection() # remove selection
				self.PositionPointer = start
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
			if finish != -1: # if a selection, delete it
				genbank.gb.Delete(start, finish)
				self.PositionPointer = start
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
			# start a selection
			if finish == -1:
				start = self.PositionPointer - 1
				finish = self.PositionPointer
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
			# start a selection
			if finish == -1:
				start = self.PositionPointer
				finish = self.PositionPointer
			
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
			print("select with keyboard")

		elif key == 317 and shift == False: #down
			self.moveCursor(1, "down")
			self.update_ownUI()
			self.update_globalUI()

		elif key == 317 and shift == True: #down + select
			print("select with keyboard")

		#self.set_dna_selection() #update the varable keeping track of DNA selection

		self.set_cursor_position() # udpate cursor position, just in case

##########################
	#def make_outputpopup(self):
	#	'''Creates a popup window in which output can be printed'''
	#	self.outputframe = wx.Frame(None, title="Output Panel") # creation of a Frame with a title
	#	self.output = output.create(self.outputframe, style=wx.VSCROLL|wx.HSCROLL) # creation of a richtextctrl in the frame




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
		start 	= 1
		finish 	= len(genbank.gb.GetDNA())+1 # TODO figure the coorect way to select. Write it in some dokumentation
		self.set_dna_selection((start, finish, -1))
		#self.set_dna_selection()
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
		start, finish, zero = self.get_selection()
		if finish == -1:
			raise ValueError('Cannot modify an empty selection')
		else:
			genbank.gb.Upper(start, finish)
			self.update_ownUI()
			self.set_dna_selection()

	def lowercase(self):
		'''Change selection to lowercase'''
		start, finish, zero = self.get_selection()
		if finish == -1:
			raise ValueError('Cannot modify an empty selection')
		else:
			genbank.gb.Lower(start, finish)
			self.update_ownUI()
			self.set_dna_selection()

	def reverse_complement_selection(self):
		'''Reverse-complement current selection'''
		start, finish, zero = self.get_selection()
		if finish == -1:
			raise ValueError('Cannot modify an empty selection')
		else:
			genbank.gb.RCselection(start, finish)
			self.update_ownUI()
			self.update_globalUI()
			self.set_dna_selection()

	def delete(self):
		'''Deletes a selection and updates dna and features'''
		start, finish, zero = self.get_selection()
		if finish == -1:
			raise ValueError('Cannot delete an empty selection')
		else:
			genbank.gb.Delete(start, finish)
			self.update_ownUI()
			self.update_globalUI()
			self.set_dna_selection()

	def cut(self):
		'''Cut DNA and store it in clipboard together with any features present on that DNA'''
		start, finish, zero = self.get_selection()
		if finish == -1:
			raise ValueError('Cannot cut an empty selection')
		else:
			genbank.gb.Cut(start, finish)
			self.set_cursor_position(start)
			self.update_ownUI()
			self.update_globalUI()
			self.set_dna_selection()

	def cut_reverse_complement(self):
		'''Cut reverse complement of DNA and store it in clipboard together with any features present on that DNA'''
		start, finish, zero = self.get_selection()
		if finish == -1:
			raise ValueError('Cannot cut an empty selection')
		else:
			genbank.gb.CutRC(start, finish)
			self.update_ownUI()
			self.update_globalUI()
			self.set_dna_selection()
			
	def paste(self):
		'''Paste DNA and any features present on that DNA'''
		start, finish, zero = self.get_selection()
		if finish == -1:
			pass
		else: #If a selection, remove sequence, reposition pointer
			genbank.gb.Delete(start, finish, visible=False)
			self.set_cursor_position(start)
		genbank.gb.RichPaste(genbank.cursor_position)
		self.update_ownUI()
		self.update_globalUI()
		self.set_dna_selection() # reset selection

	def paste_reverse_complement(self):
		'''Paste reverse complement of DNA and any features present on that DNA'''
		start, finish, zero = self.get_selection()
		if finish == -1:
			pass
		else: #If a selection, remove sequence, reposition pointer
			genbank.gb.Delete(start, finish, visible=False)
			self.set_cursor_position(start)
		genbank.gb.PasteRC(genbank.cursor_position)
		self.update_ownUI()
		self.update_globalUI()
		self.set_dna_selection() # reset selection

	def copy(self):
		'''Copy DNA and features into clipboard'''
		start, finish, zero = self.get_selection()

		if finish == -1:
			raise ValueError('Cannot copy an empty selection')
		else:
			#genbank.gb.Copy(start, finish)

			# try the new copy:
			genbank.gb.RichCopy(start, finish)

	def copy_reverse_complement(self):
		'''Copy reverse complement of DNA'''
		start, finish, zero = self.get_selection()
		if finish == -1:
			raise ValueError('Cannot copy an empty selection')
		else:
			genbank.gb.CopyRC(start, finish)


#######################################################

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
			raise ValueError('Cannot translate an empty selection')
		else:
			DNA = genbank.gb.GetDNA(start, finish)
			protein = dna.Translate(DNA)
			self.translate_output(protein, DNA, 'leading strand')

	def translate_selection_reverse_complement(self):
		'''Translate reverse-complement of selected DNA'''
		start, finish = self.get_selection()
		if finish == -1:
			raise ValueError('Cannot translate an empty selection')
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

