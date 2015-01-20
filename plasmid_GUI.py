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

#TODO
#fix long plasmid names
#add 'dna ruler'
#add rightclick menus

import wx

import cairo
from wx.lib.wxcairo import ContextFromDC

from wx.lib.pubsub import setupkwargs 		#this line not required in wxPython2.9.
 	                                  	#See documentation for more detail
from wx.lib.pubsub import pub


import genbank
import copy
import math

import os, sys
import string
from base_class import DNApyBaseDrawingClass
from base_class import DNApyBaseClass
import featureedit_GUI

import colcol

import options					# new option file to make options easy to change in one file (or settings.txt)

files={}   #list with all configuration files
files['default_dir'] = os.path.abspath(os.path.dirname(sys.argv[0]))+"/"
files['default_dir']=string.replace(files['default_dir'], "\\", "/")
files['default_dir']=string.replace(files['default_dir'], "library.zip", "")
settings=files['default_dir']+"settings"   ##path to the file of the global settings
execfile(settings) #gets all the pre-assigned settings




class PlasmidView(DNApyBaseDrawingClass):
	'''
		Draws a plasmid view using cairo on a wx.dc
		uses current file if any
	'''	
	def __init__(self, parent, id):
		
		# load options into self.opt
		self.opt = options.options()

		# variables to store current events like highlight etc.
		# for hover and click
		self.hittest = {}			# to store the hittest paths
		self.Highlight 		= None		# hightlight this on move
		self.HighlightClick 	= None		# hightlight constantly
		
		# for selection
		self.selectionDrawing  		= []
		self.selectionDrawingDirection 	= "right"
		self.selectionDrawingDiff	= 0
		
		# initialise the window
		DNApyBaseDrawingClass.__init__(self, parent, wx.ID_ANY)

		self.parent = parent
	
		genbank.dna_selection = (1,1)

		self.Bind(wx.EVT_LEFT_DOWN, self.OnLeftDown)
		self.Bind(wx.EVT_LEFT_UP, self.OnLeftUp)
		#self.Bind(wx.EVT_RIGHT_UP, self.OnRightUp)
		self.Bind(wx.EVT_MOTION, self.OnMotion)
		#self.Bind(wx.EVT_LEFT_DCLICK, self.OnLeftDouble)
		



		

############ Setting required methods ####################

	def update_globalUI(self):
		'''
		Method should be modified as to update other panels in response to changes in own panel.
		'''
		MSG_CHANGE_TEXT = "change.text"
		pub.sendMessage(MSG_CHANGE_TEXT, text="Plasmid view says update!")


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
		
		dc.SetBackground(wx.Brush("White"))	# start the DC and clear it
		dc.Clear() 				# make sure you clear the bitmap!
		ctx = ContextFromDC(dc)		# load it into cairo
		self.DrawCairo(ctx)
		
		dc.SelectObject(wx.NullBitmap) # need to get rid of the MemoryDC before Update() is called.
		self.Refresh()
		self.Update()


	def set_dna_selection(self, selection):
		'''Receives requests for DNA selection and then sends it.'''
		assert type(selection) == tuple, 'Error, dna selection must be a tuple'
		selection = (int(selection[0]-1), int(selection[1]))
		genbank.dna_selection = selection
		self.update_globalUI()


############### Done setting required methods #######################
	
	def exportMapAs(self, filepath, fileType=".svg"):
		'''	This method is similar to update_ownUI only that it 
			exports a file
			it makes a surface as svg instead of wx.dc
		'''
		
		# create new surface:
		fo = file(filepath, 'w')
		width, height = 1000,1000 # set to 1000x1000, can be anything
		
		if fileType==".svg":
			surface = cairo.SVGSurface (fo, width, height)
		elif fileType==".png":
			# does not work yet!
			surface = cairo.Surface (fo, width, height) 
		else:
			return False

		ctx = cairo.Context (surface)
		self.DrawCairo(ctx, (width,height) )
		surface.finish()
		
		self.update_ownUI()
		return True
	
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





	def DrawCairo(self, ctx, exportSize=False):
		''' Function that draws the whole plasmid in every aspect'''

		self.ctx = ctx

		
		if exportSize == False:
			width, height = self.GetVirtualSize()	# set hight and width. Its 100, 100 
			canvasWidth  = self.size[0]
			canvasHeight = self.size[1]
			ratio = float(canvasHeight)/ float(canvasWidth)
		else:
			width, height = exportSize	# set hight and width. Its 2000, 1000 maybe 
			canvasWidth  = width
			canvasHeight = height
			ratio = float(canvasHeight)/ float(canvasWidth)
		
		self.labelHelper		= {}	# a helper variable for labels drawn into arrows!
		
		# only resize the canvas if software is loaded. This prevents error messages
		if canvasWidth != 0 and canvasHeight != 0 and canvasWidth > 100 and canvasHeight > 100:
			self.ctx.scale(canvasWidth*ratio, canvasHeight) # Normalizing the canvas
			self.ctx.scale(0.01, 0.01) 			# make it 100*ratiox100
			self.ctx.translate (50/ratio,50) 		# set center to 0,0		
		
		
		# get radius rom options:
		radius 		= self.opt.pRadius
		self.radius 	= radius
		
		# load label positions:
		self.labelPos = self.labelPositions()
		
		
		# draw the helper plain white circle to make selection possible:
		self.ctx.arc(0,0,radius+2.5,0,2*math.pi)		# outer circle
		self.ctx.set_source_rgb (1, 1, 1) 			# white
		self.ctx.set_line_width (0)				# no line
		self.markerArea1 = self.ctx.copy_path()			# copy for later
		self.ctx.fill()						# fill
		
		self.ctx.arc(0,0,radius-2.5,0,2*math.pi)		# inner circle
		self.ctx.set_source_rgb (1,1,1)				# Solid white		
		self.ctx.set_line_width (0)				# no line
		self.markerArea2 = self.ctx.copy_path()			# copy for later
		self.ctx.fill()						# fill



		# draw the plasmid (two circles with different radius)
		self.ctx.arc(0,0,radius+0.4,0, 2*math.pi)		# circle
		self.ctx.set_source_rgb (0, 0, 0) 			# Solid color
		self.ctx.set_line_width (0.4)				# line width
		self.ctx.stroke()					# stroke only no fill!

		# inner plasmid
		self.ctx.arc(0,0,radius-0.4,0,2*math.pi)		# inner circle
		self.ctx.stroke()					# stroke the same color and settings
		
		# if a file is loaded draw the rest
		if genbank.gb.GetDNA():



			#draw plasmid name
			self.drawCairoPlasmidName(self.ctx)
		

			
			# draw features and labels in arrow
			self.drawCairoFeatures()
			
			# should we draw selection?
			if len(self.selectionDrawing) > 1:
				start_rad 	= self.selectionDrawing[0]
				end_rad 	= self.selectionDrawing[1]
				

				# now draw blue selection:
				self.ctx.set_source_rgba(0,0.67,1,0.75) # Solid color
				self.ctx.set_line_width(3)
				x,y = self.radial2cartesian(radius, start_rad)
				
				
				if start_rad - self.selectionDrawingDirection > 0:
					# start is bigger than first finish, so 5'-3'
					direction = "left"
				else:
					direction = "right"
				
				
				self.ctx.move_to(x,y)
				if direction == "right": # 3' to 5'
					self.ctx.arc(0,0,radius, start_rad, end_rad )
				else:
					self.ctx.arc_negative(0,0,radius, start_rad, end_rad )
				self.ctx.stroke()
				
				
				## LINES ###################
				# set style
				self.ctx.set_source_rgb(0,0,0) # Solid color
				self.ctx.set_line_width(0.2) # or 0.1
				
				
				# load the angles and draw the line!
				
				sx1, sy1 = self.radial2cartesian(radius+2, start_rad)
				sx2, sy2 = self.radial2cartesian(radius-2, start_rad)
				self.ctx.move_to(sx1,sy1)
				self.ctx.line_to(sx2,sy2)
				self.ctx.stroke()
				
				# and the endline, same thing
				
				ex1, ey1 = self.radial2cartesian(radius+2, end_rad)
				ex2, ey2 = self.radial2cartesian(radius-2, end_rad)
				self.ctx.move_to(ex1,ey1)
				self.ctx.line_to(ex2,ey2)
				self.ctx.stroke()
				
				
			elif len(self.selectionDrawing) == 1:
				# just show one line for hover here
					
				## LINES ###################
				# set style
				self.ctx.set_source_rgb(0,0,0) # Solid color
				self.ctx.set_line_width(0.2) # or 0.1
				
				
				# load the angles and draw the line!
				start_rad = self.selectionDrawing[0]
				sx1, sy1 = self.radial2cartesian(radius+2, start_rad)
				sx2, sy2 = self.radial2cartesian(radius-2, start_rad)
				self.ctx.move_to(sx1,sy1)
				self.ctx.line_to(sx2,sy2)
				self.ctx.stroke()
		
		

		return True





	
	
	
	
	
	def drawCairoFeatures(self):
		'''Loops trough features and calls other functions to draw the arrows and
		the labels'''
		
		# set some variables
		featurelist             = genbank.gb.get_all_feature_positions()
		self.drawn_fw_locations = [] 		# for location tracing foreward
		self.drawn_rv_locations = [] 		# for location tracing reverse
		drawHightlight 			= False		
		highPath       			= []
		highColor				= []	
		
		
						
		
		
		
		#####################
		# Settings loaded from self.opt
		radius				= self.opt.pRadius
		arrow_thickness   	= self.opt.pArrowThick		# like everything in percent
		radius_change     	= self.opt.pRadChange		# initial space between arrow and plasmid
		lev               	= self.opt.pLevAdd		# level distance
		arrow_head_length 	= self.opt.pArrowHeadLength	
		
		length           	= float(len(genbank.gb.GetDNA()))
		degreeMult       	= float(360/length)
		
		bordercolor      	= self.opt.pFeatureBC
		borderwidth       	= self.opt.pBW
		
		bordercolorHigh   	= self.opt.pFeatureBChigh
		borderwidthHigh   	= self.opt.pBWhigh
		#####################
		
		# loop through all features
		for i in range(0,len(featurelist)): 

			# load all the feature infos
			featuretype, complement, start, finish, name, index = featurelist[i]
			feature     = featurelist[i]
			
			# variable to save as a hittest
			hittestName 					= "%s%s" % (name, index)	# name and index as indicator!
			self.hittest[hittestName] 		= []						# initialise list to store arrow and label of a feature as paths
			self.labelHelper[hittestName] 	= [True]					# save radius of arrow to draw the label maybe

			
			# highlitght this element?
			if self.Highlight == hittestName or self.HighlightClick == hittestName:
				drawHightlight 	= True
				passBC    		= bordercolorHigh 	# border Color
			else:
				passBC    		= bordercolor		# border Color
				drawHightlight	= False
		
			# draw the arrow:
			arrowResult, arrowResult2 = self.drawCairoArrow(feature, hittestName, radius, radius_change, lev,arrow_head_length, length, degreeMult, borderwidth, passBC, drawHightlight )
			
			# if we want to highlight this
			if arrowResult != False:
				# we get color and the path from the function to save!
				highPath.append(arrowResult)
				highColor.append(arrowResult2)

		
			
		# now draw the highlitgh arrow, so its on top:
		# it is then drawn twice at the same position!
		if len(highPath) >0:
			i = 0
			for  path, color in zip(highPath, highColor):
				self.ctx.append_path(path)
				r,g,b = colcol.hex_to_rgb(color)
				r = float(r)/255
				g = float(g)/255
				b = float(b)/255
				self.ctx.set_source_rgb (r,g,b,) 		# Solid color
				self.ctx.set_line_width (borderwidthHigh) 	# or 0.1
				self.ctx.fill_preserve ()
				self.ctx.set_source_rgb (bordercolorHigh[0],bordercolorHigh[1],bordercolorHigh[2])	
				self.ctx.stroke()
			
	
		
		# loop through all features to draw the labels in arrows!
		for i in range(0,len(featurelist)):	
			# load all the feature infos
			featuretype, complement, start, finish, name, index = featurelist[i]
			feature     = featurelist[i]
		
			# variable to save as a hittest
			hittestName = "%s%s" % (name, index)		# random pattern
			
			# draw the label and the line connecting them
			labelresult = self.drawArrowLabels(feature, hittestName)
			if labelresult != False:
				# this label was not drawn into arrow
				# save it hittest name as missing
				self.labelHelper[hittestName] = [False,feature]
			
			
		
		# here draw the other Labels:
		self.drawCairoLabels()
					
	
		return True
	
	
	
	
	# only draw arrow of feature
	def drawCairoArrow(self, feature, hittestName, radius,  radius_change, lev,arrow_head_length, length, degreeMult, bw, bc, high ):
		'''Function to draw the arrow onto the canvas and save 
		its path for a hittest'''
		# load infos:
		featuretype, complement, start, finish, name, index = feature
		
		# draw
		if complement == False:
			self.drawn_fw_locations, levelMult = self.find_overlap(self.drawn_fw_locations, (start, finish))
			if levelMult>3 or finish-start<=5: #Only allow for tree levels. Also, for very short features, draw them at bottom level
				levelMult = 0
				radius_change = radius_change * 1.25

			level 				= lev * (1+levelMult)
	
			s_deg				= float(start)  * degreeMult			# in degree
			f_deg 				= float(finish) * degreeMult
			s_rad 				= (s_deg * math.pi/180) - math.pi/2 
			f_rad   			= (f_deg * math.pi/180) - math.pi/2 
		
			f_rad_wo_head  		= ((f_deg-arrow_head_length) * math.pi/180) - math.pi/2 
			f_rad_wo_head_help 	= ((f_deg-2*arrow_head_length) * math.pi/180) - math.pi/2 
			

			
			# weather to draw a head or not, only if feature is large enough
			drawHead = True
			if f_rad_wo_head_help < s_rad:
				f_rad_wo_head 	= f_rad
				drawHead 		= False


			arrow_head_x, arrow_head_y  	= self.radial2cartesian(radius+level, f_rad)
			arrow_head_x1, arrow_head_y1  	= self.radial2cartesian(radius+level+radius_change+0.5, f_rad_wo_head)
			arrow_head_x2, arrow_head_y2  	= self.radial2cartesian(radius+level-radius_change-0.5, f_rad_wo_head)


			# get the color:
			color = eval(featuretype)['fw'] #get the color of feature (as string)
			assert type(color) == str

			# move to start
			x_start, y_start = self.radial2cartesian(radius+level+radius_change, s_rad)
			self.ctx.move_to(x_start, y_start)
			
			self.labelHelper[hittestName].append(radius+level-radius_change)	# inner radius
			self.labelHelper[hittestName].append(radius+level+radius_change)	# outer radius
			
			# the actual drawing:
			self.ctx.arc(0, 0, radius+level+radius_change,s_rad, f_rad_wo_head);

			if drawHead:
				self.ctx.line_to(arrow_head_x1,arrow_head_y1)
				self.ctx.line_to(arrow_head_x,arrow_head_y)
				self.ctx.line_to(arrow_head_x2,arrow_head_y2)
			self.ctx.arc_negative(0, 0, radius+level-radius_change, f_rad_wo_head, s_rad);
			self.ctx.close_path ()
			
		else:
			self.drawn_rv_locations, levelMult = self.find_overlap(self.drawn_rv_locations, (start, finish))
			if levelMult>3 or finish-start<=5: #Only allow for tree levels. Also, for very short features, draw them at bottom level
				levelMult = 0
				radius_change = radius_change * 1.25 # for bigger items!

			level 				= lev * (1+levelMult)
	
			s_deg				= float(finish)  * degreeMult			# in degree
			f_deg 				= float(start) * degreeMult
			s_rad 				= (s_deg * math.pi/180) - math.pi/2 
			f_rad   			= (f_deg * math.pi/180) - math.pi/2 
		
			f_rad_wo_head  		= ((f_deg + arrow_head_length) * math.pi/180) - math.pi/2 
			f_rad_wo_head_help 	= ((f_deg + 1.5 * arrow_head_length) * math.pi/180) - math.pi/2 
			
			# weather to draw a head or not, only if feature is large enough
			drawHead = True
			if f_rad_wo_head_help > s_rad:
				f_rad_wo_head 	= f_rad
				drawHead 		= False

				
			arrow_head_x, arrow_head_y  = self.radial2cartesian(radius-level, f_rad)
			arrow_head_x1, arrow_head_y1  = self.radial2cartesian(radius-level-radius_change-0.5, f_rad_wo_head)
			arrow_head_x2, arrow_head_y2  = self.radial2cartesian(radius-level+radius_change+0.5, f_rad_wo_head)

			# get the color:
			color = eval(featuretype)['rv'] #get the color of feature (as string)
			assert type(color) == str

			# move to start
			x_start, y_start = self.radial2cartesian(radius-level-radius_change, s_rad)
			self.ctx.move_to(x_start, y_start)
			
			self.labelHelper[hittestName].append(radius-level-radius_change)	# smallest radius
			self.labelHelper[hittestName].append(radius-level+radius_change)	# outer radius
			
			# the actual drawing:
			self.ctx.arc_negative(0, 0, radius-level-radius_change,s_rad, f_rad_wo_head);
			if drawHead:
				self.ctx.line_to(arrow_head_x1,arrow_head_y1)
				self.ctx.line_to(arrow_head_x,arrow_head_y)
				self.ctx.line_to(arrow_head_x2,arrow_head_y2)
			self.ctx.arc(0, 0, radius-level+radius_change, f_rad_wo_head, s_rad);
			self.ctx.close_path ()
	

			

		# now draw the arrows with color and lines
		r,g,b = colcol.hex_to_rgb(color)
		r = float(r)/255
		g = float(g)/255
		b = float(b)/255

		self.ctx.set_source_rgba (r,g,b,1.0) # Solid color
		self.ctx.set_line_width (bw) # or 0.1
		
		self.ctx.fill_preserve ()
		self.ctx.set_source_rgb (bc[0],bc[1],bc[2])
		
		# save feature to hittest:
		self.hitpath = self.ctx.copy_path()
		self.hittest[hittestName].append(self.hitpath)

		self.ctx.stroke()
		
		# return the path if it should be highlighted!
		if high == True:
			return self.hitpath, color
		else:
			return False, False
	
	
	def drawArrowLabels(self, f, hittestName):
		'''
			Draws labels into arrows
		'''
		

		featuretype, complement, start, finish, name, index = f
		
		arrowRadius = float(self.labelHelper[hittestName][1]) + 0.5
		back_rad = 0
		returnValue	= False
	
		# calculate length of coresponding arrow, ignoring the level!
		length         		= float(len(genbank.gb.GetDNA()))
		degreeMult     	  	= float(360/length)
		s_deg				= float(start)  * degreeMult			# in degree
		f_deg 				= float(finish) * degreeMult
		middle_deg			= f_deg - ((f_deg-s_deg)/2)				# the center of the arrow
		s_rad 				= (s_deg * math.pi/180)  - math.pi/2 
		f_rad   			= (f_deg * math.pi/180)  - math.pi/2 
	
		ar_l = 2 * math.pi * (arrowRadius-3.2) * (f_deg - s_deg)/360 # 3.2 is arrow head
	
		dist2rad			= (f_rad - s_rad)/ar_l
	
	
	
		# calculate length of written text
		# remove any '"' in front or at the end 
		if name[0:1] == '"':
			name = name[1:]
			lengthName = len(name)
			if name[-1:] == '"':
				name = name[:-1]
			
	
		self.ctx.select_font_face('Arial', cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD)
		self.ctx.set_font_size(2);
		self.ctx.set_source_rgb (1,1,1)
		xbearing, ybearing, TextWidth, TextHeight, xadvance, yadvance = self.ctx.text_extents(name)
	
		# we want some spacer:
		spacer 		= 0.1 * TextWidth
		ar_l_spacer = ar_l - spacer
	
		# check if text would fix into arrow
		if TextWidth < ar_l_spacer:
		
			
			# if yes we have to draw every symbole one for one
			startspacer = (ar_l - TextWidth)/2
			
			#
			# we have to draw the labbels different fore foreward and reverse features
			# also upper and lower features are different
			# so we need for time a similar code
			#
			
			# draw foreward and upper labels here:
			if complement == False and (middle_deg <= 90 or middle_deg >= 270)  :

				start_rad	= s_rad + (startspacer * dist2rad)
		
				back_rad 	=  0.5* math.pi + start_rad  # sum of rotations to later rotate back
		
				# rotate context so the drawing position is up:
				self.ctx.rotate(back_rad)
			
		
				startposition_x, startposition_y = self.radial2cartesian(arrowRadius, -0.5* math.pi)
				self.ctx.move_to(startposition_x, startposition_y); 
				for n in name:
					
					# to remove variations in width we take 20 time the same char
					helptext = "%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s" %(n,n,n,n,n,n,n,n,n,n,n,n,n,n,n,n,n,n,n,n)
					xbearing, ybearing, TextWidth, TextHeight, xadvance, yadvance = self.ctx.text_extents(helptext)
					self.ctx.show_text(n)

			
					new_rad	= (float(TextWidth) / 20.0) * dist2rad
					back_rad = back_rad + new_rad
					self.ctx.rotate(new_rad)
			elif complement == False and (middle_deg > 90 or middle_deg < 270) :

				start_rad	= f_rad - (startspacer * dist2rad) # f_rad because we write from head to start
		
				back_rad 	=  - 0.5* math.pi + start_rad  # sum of rotations to later rotate back

				# rotate context so the drawing position is down:
				self.ctx.rotate(back_rad)
			
		
				startposition_x, startposition_y = self.radial2cartesian(arrowRadius+1.3, 0.5* math.pi )
				self.ctx.move_to(startposition_x, startposition_y); 
				for n in name:
					# to remove variations in width we take 20 time the same char
					helptext = "%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s" %(n,n,n,n,n,n,n,n,n,n,n,n,n,n,n,n,n,n,n,n)
					xbearing, ybearing, TextWidth, TextHeight, xadvance, yadvance = self.ctx.text_extents(helptext)
					self.ctx.show_text(n)

			
					new_rad	= (float(TextWidth) / 20.0) * dist2rad
					back_rad = back_rad - new_rad
					self.ctx.rotate(-new_rad)
			elif complement == True and (middle_deg < 90 or middle_deg > 270) :

				start_rad	= s_rad + (startspacer * dist2rad)
		
				back_rad 	=  0.5* math.pi + start_rad  # sum of rotations to later rotate back
		
				# rotate context so the drawing position is up:
				self.ctx.rotate(back_rad)
			
		
				startposition_x, startposition_y = self.radial2cartesian(arrowRadius, -0.5* math.pi)
				self.ctx.move_to(startposition_x, startposition_y); 
				for n in name:
					
					# to remove variations in width we take 20 time the same char
					helptext = "%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s" %(n,n,n,n,n,n,n,n,n,n,n,n,n,n,n,n,n,n,n,n)
					xbearing, ybearing, TextWidth, TextHeight, xadvance, yadvance = self.ctx.text_extents(helptext)
					self.ctx.show_text(n)

			
					new_rad	= (float(TextWidth) / 20.0) * dist2rad
					back_rad = back_rad + new_rad
					self.ctx.rotate(+new_rad)
			elif complement == True and (middle_deg > 90 or middle_deg < 270) :

				start_rad	= f_rad - (startspacer * dist2rad)
		
				back_rad 	=  -0.5* math.pi + start_rad  # sum of rotations to later rotate back
		
				# rotate context so the drawing position is up:
				self.ctx.rotate(back_rad)
			
		
				startposition_x, startposition_y = self.radial2cartesian(arrowRadius+1.3, 0.5* math.pi)
				self.ctx.move_to(startposition_x, startposition_y); 
				for n in name:
					
					# to remove variations in width we take 20 time the same char
					helptext = "%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s" %(n,n,n,n,n,n,n,n,n,n,n,n,n,n,n,n,n,n,n,n)
					xbearing, ybearing, TextWidth, TextHeight, xadvance, yadvance = self.ctx.text_extents(helptext)
					self.ctx.show_text(n)

			
					new_rad	= (float(TextWidth) / 20.0) * dist2rad
					back_rad = back_rad - new_rad
					self.ctx.rotate(-new_rad)
		else:
			# we could not draw this label
			# we should draw this label later!
			returnValue = hittestName
		self.ctx.rotate(-back_rad)

		return returnValue
	
	
	
	
	
	def drawCairoLabels(self):
		'''Function to draw labels of:
			- features wich are to short for drawing inside
			- restriction enzymes
			- any other
		'''

		# cairo context in -> self.ctx
		
		# all positions should be in --> self.labelPos
		
		# all missing labels are saved as self.labelHelper[hittestName] = [False,feature]
		# so you can loop trough them
		
		for name in self.labelHelper:
			if self.labelHelper[name][0] == False:
				featuretype, complement, start, finish, name, index = self.labelHelper[name][1]
				
				#radius = self.labelHelper[name][3]

				
				# middle of feature:
				length 		= float(len(genbank.gb.GetDNA()))
				middle 		= start + (finish-start)/2 % length
				middleAngle = self.pos_to_angle_rad(middle)
				x,y 	= self.radial2cartesian(self.opt.pRadius+0.4, middleAngle) 	# middle of item, outer plasmid radius
				x1,y1 	= self.radial2cartesian(self.opt.labelRadius, middleAngle)	# endpoint of line
		
				# draw a line from label to feature				
				self.ctx.set_source_rgba(0.0,0.0,0.0,1) # Solid color
				self.ctx.set_line_width(0.2) # or 0.1
				self.ctx.move_to(x,y)
				self.ctx.line_to(x1,y1)
				self.ctx.stroke()
				
				# draw a box:
				
				#self.hitpath = self.ctx.copy_path()
				#self.ctx.fill()
				
				# draw a label
				
				# add box-path to hittest list of this feature or enzyme		
				#self.hittest[name].append(self.hitpath) 
				
				

		return True
	
	
	
	
	
	def drawCairoPlasmidName(self, ctx):
		'''Draw the plasmid name and basepairs in the middle'''
		
		#width  = 100
		#height = 100
		self.hittest['plasmidname'] = []
		
		name = genbank.gb.fileName.split('.')[0]
		basepairs = '0 bp'
		if genbank.gb.GetDNA() != None:
			basepairs = str(len(genbank.gb.GetDNA())) + ' bp'
		
		# name
		self.ctx.select_font_face('Arial', cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
		self.ctx.set_font_size(2.5);
		xbearing, ybearing, TextWidth, TextHeight, xadvance, yadvance = self.ctx.text_extents(name)
		self.ctx.move_to(-TextWidth/2, 0);  
		
		# save feature to hittest:
		self.hitpath = self.ctx.copy_path()
		self.hittest['plasmidname'].append(self.hitpath)
	
		
		self.ctx.show_text(name)
		
		# bp
		self.ctx.select_font_face('Arial', cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
		self.ctx.set_font_size(2);
		xbearing, ybearing, TextWidth, TextHeight, xadvance, yadvance = self.ctx.text_extents(basepairs)
		self.ctx.move_to(-TextWidth/2, 2.7);  
		self.ctx.show_text(basepairs);
		
		
		return True
	
	
	
	
	
	def drawCairoSelection(self, start=0, finish=0, feature=False):
			
		# we need:
		radius			  = self.radius

	
		self.selectionDrawing = []
		
		# feature to highlight or drag and drop?
		if feature != False:
			featurelist             = genbank.gb.get_all_feature_positions()
			for i in range(0,len(featurelist)):	
				# load all the feature infos
				featuretype, complement, start, finish, name, index = featurelist[i]
		
				temHit = "%s%s" % (name, index)
				if feature == temHit:
					length         		= float(len(genbank.gb.GetDNA()))
					degreeMult     	  	= float(2*math.pi/length)
					start_rad			= self.pos_to_angle_rad(start)			# in degree
					finish_rad 			= self.pos_to_angle_rad(finish)	
					
					self.selectionDrawingDirection = start_rad + 0.01
					
					# set selection for real!
					self.set_dna_selection((start,finish))
					
		
		else:
		
			# we have a selection:
			start_rad		= start#   - math.pi/2 
			finish_rad   	= finish # + math.pi/2 
		
		# set Drawing angles
		self.selectionDrawing.append(start_rad)	
		self.selectionDrawing.append(finish_rad)	
		
		
		# save a third value for direction.
		# it is the first value of the end
		# so it is the second start point!
		# if its greater than start is is a 3'-5' direction
		# else it should be 5'-3' selection!
		
		# we have to reset this as soon as we hit the start again!:

		# check if the movement just hit the start and maybe changed direction of selection!
		# this part is pretty awesome and cost me some thinking. Hope it works well
		oldDrawingDiff 				= self.selectionDrawingDiff
		self.selectionDrawingDiff 	= start_rad - finish_rad
		
		oldSign = math.copysign(1, oldDrawingDiff)
		newSign = math.copysign(1, self.selectionDrawingDiff)
		
		if math.fabs(round(self.selectionDrawingDiff)) < 2 and oldSign != newSign:
			self.selectionDrawingDirection = False
	
		
		if self.selectionDrawingDirection == False:
			self.selectionDrawingDirection = finish_rad

		
		self.update_ownUI()
		return True
	
	
	
	# function to make a List of allowed label positions
	# the later can be used to position as much labels on the screen as possible
	def labelPositions(self):
		allowedPositions = []
		
		# how high should a label be?
		generalHight = 15 # px on screen
		generalHight = self.ctx.device_to_user_distance(0,generalHight)[1]

		# get plasmid view width and height
		width, height = self.GetVirtualSize()
		

		
		# so labels are arranged around the plasmid.
		# at least the radius of the plasmid and 3 layers of arrows have to fit below:		
		labelRadius = self.opt.labelRadius


		
		
		
		
		return allowedPositions


	
############### Setting methods for converting stuff ##############

	def radial2cartesian(self,radius, angle, cx=0,cy=0): # angle in radial!
		radius = float(radius)
		angle  = float(angle)
		x = radius * math.cos(angle) + cx
		y = radius * math.sin(angle) + cy
		
		return x,y
		
	def cartesian2radial(self,x, y):

		x=-x # to make it correct
		
		add =  math.radians(90)
		a = math.atan2(x,y) + 1 * math.pi - add
		
		return a

	def angle_to_pos(self, angle):
		'''Convert an angle of a circle to a DNA position'''
		len_dna = float(len(genbank.gb.GetDNA()))
		dna_pos = int(self.AngleToFraction(angle)*len_dna)
		return dna_pos

	def pos_to_angle(self, pos):
		'''Calculate angles from DNA positions'''
		assert type(pos) == int, 'Error, position needs to be an integer'
		len_dna = float(len(genbank.gb.GetDNA()))
		if len_dna == 0:
			angle = 0
		else:
			angle = self.FractionToAngle(pos/float(len_dna))
		return angle	
		
	def angle_to_pos_rad(self, angle):
		len_dna = float(len(genbank.gb.GetDNA()))
		radMult        	= float(2*math.pi/len_dna)
		position 		= (angle-math.pi/2)/radMult + 0.5 * len_dna
		
		return position
	
	def pos_to_angle_rad(self, position):
		len_dna = float(len(genbank.gb.GetDNA()))
		radMult        	= float(2*math.pi/len_dna)
		angle = (position - 0.5 * len_dna) * radMult + math.pi/2
		
		return angle

########## Done with angle to dna methods ####################



######### Mouse methods #####################


	def HitTest(self):
		'''Tests whether the mouse is over any feature or label'''
		hit = None
		x, y = self.ScreenToClient(wx.GetMousePosition())	
		
		# get the mouse positions
		x2, y2 = self.ctx.device_to_user(x,y)

		
		for path in self.hittest:
			cairoCtx = self.ctx
			for i in self.hittest[path]:
				# load the path
				self.ctx.append_path(i)
				inFill = self.ctx.in_fill(x2,y2) 
				inStroke = self.ctx.in_stroke(x2,y2)
				# check if this path is hit
				if inFill == True or inStroke == True:
					hit = path
			
					
				self.ctx.fill()
		
		return hit
	
	def HitTestMarkerArea(self):
		'''Check if the marker area (area where selection is possible) was hit'''
		
		# get the mouse position
		x, y = self.ScreenToClient(wx.GetMousePosition())	
		x2, y2 = self.ctx.device_to_user(x,y)

		self.ctx.append_path(self.markerArea1)
		inFill1 = self.ctx.in_fill(x2,y2) 
		self.ctx.stroke()
		
		self.ctx.append_path(self.markerArea2)
		inFill2 = self.ctx.in_fill(x2,y2) 
		self.ctx.stroke()
		
		# check if the path is hit for selection
		if inFill1 == True and inFill2 == False:
			return True
		else:
			return False

		
			

	def OnLeftDown(self, event):
		'''When left mouse button is pressed down, store angle at which this happened.'''
		
		# get the mouse position
		x, y = self.ScreenToClient(wx.GetMousePosition())	
		x2, y2 = self.ctx.device_to_user(x,y)
		
		self.selectionDrawing 			= [] # empty selection
		self.selectionDrawingDirection  = False
		
		inArea = self.HitTestMarkerArea()
		
		if inArea == True:		
			a  = self.cartesian2radial(x2,y2)
			self.left_down_angle = a
		else:
			self.left_down_angle = False		
				
			
		
		#return hit


	def OnLeftUp(self, event):
		# remove start point, so we know selection is done:
		self.left_down_angle == False
		
		# maybe he hit a feature?
		hit = self.HitTest()
		if hit != None:
			self.HighlightClick = hit
			# we hit one, so we need to draw the selection!
			# wich feature was hit?
			self.drawCairoSelection(0, 0, hit)
			
		else:
			self.HighlightClick = None
		#self.set_dna_selection((start, finish))
		self.update_ownUI()


	def OnMotion(self, event):
		'''When mouse is moved with the left button down determine the DNA selection from angle generated at mouse down and mouse move event.'''
		
		if event.Dragging() and event.LeftIsDown() and self.left_down_angle != False:
			# we have to make a selection now!
			start_angle = self.left_down_angle
			
			# get the mouse positions
			x, y = self.ScreenToClient(wx.GetMousePosition())	
			x2, y2 = self.ctx.device_to_user(x,y)

			a  = self.cartesian2radial(x2,y2)
			finish_angle = a
			self.drawCairoSelection(start_angle, finish_angle)
			

			
			#set genebank selection
			start 	= self.angle_to_pos_rad(start_angle)
			finish 	= self.angle_to_pos_rad(finish_angle)
			
			self.set_dna_selection((start,finish))

		
		
		else:
			oldhit 	= self.Highlight
			hit 	= self.HitTest()
			if hit:
				self.Highlight = hit
				# only update after change
				if oldhit != hit:
					self.update_ownUI()
				
			elif len(self.selectionDrawing) < 2: # only hover position, if none is selected --> maybe use two variables to make both possible?
				# no feature to highlight, but maybe selection hover?
				self.Highlight = None
				
				# check if we are in the "hover area"
				inArea = self.HitTestMarkerArea()
		
				if inArea == True:
					# get the mouse position
					x, y = self.ScreenToClient(wx.GetMousePosition())	
					x2, y2 = self.ctx.device_to_user(x,y)
					self.selectionDrawing = [] # empty it
					self.selectionDrawing.append(self.cartesian2radial(x2,y2))	
					self.update_ownUI()
				elif oldhit != hit:
					self.update_ownUI()
				elif len(self.selectionDrawing) == 1: # TODO 
					self.selectionDrawing = [] # empty it so the hover line wont be seen all the time
					self.update_ownUI()
			

		

		
		
		
		#if event.Dragging() and event.LeftIsDown():
		
			
			#up_angle = self.PointsToAngle(self.centre_x, self.centre_y, x, y)
			#down_angle = self.left_down_angle

			#if down_angle <= up_angle:
			#	start = self.angle_to_pos(down_angle)
			#	finish = self.angle_to_pos(up_angle)
			#elif down_angle > up_angle:
			#	start = self.angle_to_pos(up_angle)
			#	finish = self.angle_to_pos(down_angle)			

			#self.set_dna_selection((start, finish))
			#self.update_ownUI()
		'''else:
			new_index = self.HitTest()
			if new_index is self.highlighted_feature: #if the index did not change
				pass
			else:
				self.highlighted_feature = new_index
				self.update_ownUI()'''


	def OnLeftDouble(self, event):
		'''When left button is double clicked, launch the feature edit dialog.'''
		'''new_index = self.HitTest() #this does not get the "true" feature index. Some featues are split and this is an index that accounts for that.
		if new_index is not False: #False is returned for the background
			featurelist = genbank.gb.get_all_feature_positions()
			featuretype, complement, start, finish, name, index = featurelist[new_index]
			genbank.feature_selection = copy.copy(index)

			dlg = featureedit_GUI.FeatureEditDialog(None, 'Edit Feature') # creation of a dialog with a title'''
		dlg.ShowModal()
		dlg.Center()


	def OnRightUp(self, event):
		print('plasmid right')


############ Done with mouse methods ####################

class PlasmidView2(PlasmidView):
	'''
	This class is intended to glue together the plasmid drawing with control buttons.
	'''
	def __init__(self, parent, id):
		DNApyBaseClass.__init__(self, parent, id)
		panel1 = wx.Panel(self)
		panel2 = wx.Panel(self)

		##########  Add buttons and methods to control their behaviour ###########	
		#buttons
		padding = 10 #how much to add around the picture

		imageFile = files['default_dir']+"/icon/circle.png"
		image1 = wx.Image(imageFile, wx.BITMAP_TYPE_ANY).ConvertToBitmap()
		circle = wx.BitmapButton(panel1, id=10, bitmap=image1, size = (image1.GetWidth()+padding, image1.GetHeight()+padding), name = "share")

		imageFile = files['default_dir']+"/icon/group.png"
		image1 = wx.Image(imageFile, wx.BITMAP_TYPE_ANY).ConvertToBitmap()
		group = wx.BitmapButton(panel1, id=11, bitmap=image1, size = (image1.GetWidth()+padding, image1.GetHeight()+padding), name = "share")

		imageFile = files['default_dir']+"/icon/radiating.png"
		image1 = wx.Image(imageFile, wx.BITMAP_TYPE_ANY).ConvertToBitmap()
		radiating = wx.BitmapButton(panel1, id=12, bitmap=image1, size = (image1.GetWidth()+padding, image1.GetHeight()+padding), name = "share")

		
		imageFile = files['default_dir']+"/icon/new_small.png"
		image1 = wx.Image(imageFile, wx.BITMAP_TYPE_ANY).ConvertToBitmap()
		newfeature = wx.BitmapButton(panel1, id=1, bitmap=image1, size = (image1.GetWidth()+padding, image1.GetHeight()+padding), name = "share")

		imageFile = files['default_dir']+"/icon/remove_small.png"
		image1 = wx.Image(imageFile, wx.BITMAP_TYPE_ANY).ConvertToBitmap()
		deletefeature = wx.BitmapButton(panel1, id=2, bitmap=image1, size = (image1.GetWidth()+padding, image1.GetHeight()+padding), name = "share")

		imageFile = files['default_dir']+"/icon/edit.png"
		image1 = wx.Image(imageFile, wx.BITMAP_TYPE_ANY).ConvertToBitmap()
		edit = wx.BitmapButton(panel1, id=6, bitmap=image1, size = (image1.GetWidth()+padding, image1.GetHeight()+padding), name = "edit")
		
		#bind feature list buttons
		self.Bind(wx.EVT_BUTTON, self.OnCircularLabels, id=10)
		self.Bind(wx.EVT_BUTTON, self.OnGroupLabels, id=11)
		self.Bind(wx.EVT_BUTTON, self.OnRadiatingLabels, id=12)
		
		self.Bind(wx.EVT_BUTTON, self.OnNew, id=1)
		self.Bind(wx.EVT_BUTTON, self.OnDelete, id=2)
		self.Bind(wx.EVT_BUTTON, self.OnEditFeature, id=6)

		#arrange buttons vertically		
		sizer = wx.BoxSizer(wx.VERTICAL)
		sizer.Add(item=circle)
		sizer.Add(item=group)
		sizer.Add(item=radiating)
		sizer.Add(item=newfeature)
		sizer.Add(item=deletefeature)
		sizer.Add(item=edit)

		panel1.SetSizer(sizer)



		#add the actual plasmid view
		self.plasmid_view = PlasmidView(panel2, -1)	

		sizer1 = wx.BoxSizer(wx.HORIZONTAL)
		sizer1.Add(item=self.plasmid_view, proportion=-1, flag=wx.EXPAND)
		panel2.SetSizer(sizer1)



		#add feature list and buttons horizontally	
		sizer2 = wx.BoxSizer(wx.HORIZONTAL)
		sizer2.Add(item=panel2, proportion=-1, flag=wx.EXPAND)
		sizer2.Add(item=panel1)

		self.SetSizer(sizer2)
	
	
	
	def update_ownUI(self):
		'''
		User interface updates.
		'''
		self.plasmid_view.update_ownUI()

	def OnCircularLabels(self, evt):
		self.plasmid_view.label_type = 'circular'
		self.update_ownUI()
		
	def OnGroupLabels(self, evt):
		self.plasmid_view.label_type = 'group'
		self.update_ownUI()
		
	def OnRadiatingLabels(self, evt):
		self.plasmid_view.label_type = 'radiating'
		self.update_ownUI()
		
	def OnNew(self, evt):
		pass
		
	def	OnDelete(self, evt):
		pass
		
	def	OnEditFeature(self, evt):
		pass
		
	def OnLeftUp(self, event):
		event.Skip()
	#######################################################################
	
##### main loop
class MyApp(wx.App):
	def OnInit(self):

		frame = wx.Frame(None, -1, title="Plasmid View", size=(700,600), style = wx.NO_FULL_REPAINT_ON_RESIZE)


		self.plasmid_view = PlasmidView2(frame, -1)	

		
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

	import sys
	assert len(sys.argv) == 2, 'Error, this script requires a path to a genbank file as an argument.'
	print('Opening %s' % str(sys.argv[1]))

	genbank.gb = genbank.gbobject(str(sys.argv[1])) #make a genbank object and read file


	app = MyApp(0)
	app.MainLoop()
