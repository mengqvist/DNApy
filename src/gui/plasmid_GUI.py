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

#TODO
#fix long plasmid names
#add 'dna ruler'
#add rightclick menus

import wx

import cairo
from wx.lib.wxcairo import ContextFromDC
import pango # for text
import pangocairo # for glyphs and text

from wx.lib.pubsub import setupkwargs 		#this line not required in wxPython2.9.
 	                                  	#See documentation for more detail
from wx.lib.pubsub import pub



import genbank
import copy
import math
import random # testweise

import os, sys
import string
from base_class import DNApyBaseDrawingClass
from base_class import DNApyBaseClass
import featureedit_GUI

import colcol

import hashlib

import options					# new option file to make options easy to change in one file (or settings.txt)

files					={}   #list with all configuration files
files['default_dir'] 	= os.path.abspath(os.path.dirname(sys.argv[0]))+"/"
files['default_dir']	= string.replace(files['default_dir'], "\\", "/")
files['default_dir']	= string.replace(files['default_dir'], "library.zip", "")
settings				= files['default_dir'] + "settings"   ##path to the file of the global settings
exec(compile(open(settings).read(), settings, 'exec')) #gets all the pre-assigned settings





class plasmidstore():
	''' class that defines a storing object for the plasmid drawing process'''
	def __init__(self):
		# BASIC INFO
		self.features 		= None 		# hold the feature info
		self.enzymes		= None 		# hold the selected enzyme info
		
		# FEATURES
		self.drawnfeatures 	= {} 		# hold all cairo objects to draw Features
		self.drawnlabels 	= {} 		# hold all written stuff
		self.drawnLLines 	= {} 		# labellines
		self.labelBoxes		= {} 		# make selecting labels more easy
		self.labelULine		= {}
		
		# ENZYMES
		self.drawnlabelsE 	= {} 		# hold all written stuff
		self.drawnLLinesE 	= {} 		# labellines
		self.labelBoxesE	= {} 		# make selecting labels more easy
		
		# SELECTION
		self.drawnSelection = None		# one path for selection only
		self.drawnPosition 	= None 		# hold the line, indicating the cursor position
			
		# LABELING HELPER
		self.LiF 			= [] # labels to be written in feature arrows
		self.LoF 			= [] # labels to be written as textlabels of features
		self.LoE 			= [] # labels to be writtem as restrictionlabels
		self.LiFold 		= [] # labels formerly to be written in feature arrows
		self.LoFold 		= [] # labels formerly to be written as textlabels of features
		self.LoEold 		= [] # labels formerly to be writtem as restrictionlabels
		
		# outer radius of all features for drawing lines:
		self.radiusOuter	= {}
		
		# interaction:
		self.interaction	= {	'hit'				: None,		# save the hitTest name for highlight
								'leftDown'			: None,
								'markerArea1' 		: None,
								'markerArea2' 		: None,
								'selection'			: (0, 0, -1),
								'selectionOld'		: (0, 0, -1),
								'position'			: 1,
								'positionOld'		: 0
							  } 		
		

		return None


class drawPlasmid(DNApyBaseDrawingClass):
	''' This class handle the drawing of a plasmid'''
	'''
		it should draw:
		- Plasmid 
		- Features
		- feature names
		- RestriktionSites
		- Names of Restriction sites
		- Selection
		- higlights

		It has to Calculate in this order:
		- The Size of the Plasmid
		- The angular positon of the features
		- If the Label of a feature can be drawn inside
		- the angular positions of all Restrictionsites
		- The arrangements of the labels
			- featurelabels
			- Restriktionsite Labels

		Possible Interactions:
		- select dna
			- drag and drop select
			- click on feature
		- zoom in and out
		- double click feature

		Each cycle must be something like this:

		1. Check if we have to calcaulate the length of the features
			self.calc.fl = True
			if true: check if this changes anything for the other labels (do labels have to be drawn outside or inside compared to before)
		2. Check if we have to calcaulate Label Positions again:
			self.calc.lp = true
		3. If self.calc.lp and fl == False:
			Just redraw the previous plasmid
		   else:
			recalculate the desired stuff

		x. Highlight and interaction'''
	def __init__(self, parent, id):
	
		# object to store our cacluations over time
		self.plasmidstore 		= plasmidstore()		
		self.radius 			= 25								# cairo unit
		self.radiusI 			= self.radius - 0.013 * self.radius	# cairo unit
		self.radiusO 			= self.radius + 0.013 * self.radius	# cairo unit
		self.radiusLabels 		= self.radius + 17 					# cairo unit
		self.arrowWidth 		= 1.8 								# cairo unit
		self.labelHeight 		= 12 								# pixel
		

		self.Highlight 			= None								# store who to highlight
		self.enzymeGap			= "-|-"								# seperator to make unique hitname less unique
		
		#self.EnzymeMonitor = Monitor()
		
		# start the window
		DNApyBaseDrawingClass.__init__(self, parent, wx.ID_ANY)
		
		# bind events
		self.Bind(wx.EVT_LEFT_UP, self.OnLeftUp)
		self.Bind(wx.EVT_LEFT_DOWN, self.OnLeftDown)
		self.Bind(wx.EVT_MOTION, self.OnMotion)
		self.Bind(wx.EVT_LEFT_DCLICK, self.OnLeftDouble)


		
	
		
		return None
	
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

		# update selection from editor
		self.updateSelection() 
		self.updatePosition()
		
		# start the calculations if any.
		self.checkFeatureCalculations() 		# changed features may change labels also
		self.checkEnzymeCalculations() 			# check for restriction enzymes
		self.checkLabelCalculations() 			# changed features may change labels also
		
		# check if we even have to redraw anything:

		dc = wx.MemoryDC()
		dc.SelectObject(self._Buffer)

		dc.SetBackground(wx.Brush("White"))	# start the DC and clear it
		dc.Clear() 				# make sure you clear the bitmap!
		self.ctx = ContextFromDC(dc)		# load it into cairo
		
		# seletion starts here
		self.checkSelectionCalculations()		# check and make a selection arc

		
		# glyphs
		#self.pango = pangocairo.CairoContext(self.ctx)
		#self.pango.set_antialias(cairo.ANTIALIAS_SUBPIXEL)

		self.draw()

		dc.SelectObject(wx.NullBitmap) # need to get rid of the MemoryDC before Update() is called.
		self.Refresh()
		self.Update()

		return None
	
	def set_dna_selection(self, selection):
		'''Receives requests for DNA selection and then sends it.'''
		assert type(selection) == tuple, 'Error, dna selection must be a tuple'
		selection = (int(selection[0]), int(selection[1]), int(selection[2]))
		genbank.dna_selection = selection
		self.plasmidstore.interaction["selection"] = selection
		self.update_globalUI()
	


	def set_cursor_position(self, position):
		# set position of cursor for status bar:
		genbank.cursor_position 					= position 
		self.plasmidstore.interaction["position"] 	= position

	def updateSelection(self):
		''' get selection from editor and update own store '''
		selection =  genbank.dna_selection
		assert type(selection) == tuple, 'Error, dna selection must be a tuple'
		selection = (int(selection[0]), int(selection[1]), int(selection[2]))
		self.plasmidstore.interaction["selection"] = selection
	
	def updatePosition(self):
		''' get position from editor and update own store '''
		position 									=  genbank.cursor_position
		self.plasmidstore.interaction["position"] 	= position

	def hitName(self, name, index):
		# unique and reproducible id for each feature
		name 	= self.nameBeautiful(name)
		seq 	= "%s%s" % (name, index)
		name 	= ''.join(seq.split())
		return name
	
	def HitTest(self):
		'''Tests whether the mouse is over any feature or label'''
		hit 	= None


		# get the mouse positions
		x, y 	= self.ScreenToClient(wx.GetMousePosition())
		x2, y2 	= self.ctx.device_to_user(x,y)

		# list of all the paths:
		loop = [self.plasmidstore.drawnfeatures, 	# features
				self.plasmidstore.drawnlabels, 		# featurelabes
				self.plasmidstore.labelBoxes,		# featurelabels
				self.plasmidstore.drawnlabelsE,  	# enzymes
				self.plasmidstore.labelBoxesE]		# enzymes
		# loop over all possible paths to find the one we have under the mouse
		for paths in loop:
			for i in paths:
				path = paths[i]
				if type(path) == list:
					path = path[0]
				
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
	
	def HitTestSelection(self):
		'''Check if the selection area (area where selection is possible) was hit'''

		# get the mouse position
		x, y 	= self.ScreenToClient(wx.GetMousePosition())
		x2, y2 	= self.ctx.device_to_user(x,y)

		self.ctx.append_path(self.plasmidstore.interaction["markerArea1"])
		inFill1 = self.ctx.in_fill(x2,y2)
		self.ctx.stroke()

		self.ctx.append_path(self.plasmidstore.interaction["markerArea2"])
		inFill2 = self.ctx.in_fill(x2,y2)
		self.ctx.stroke()

		# check if the path is hit for selection
		if inFill1 == True and inFill2 == False:
			return True
		else:
			return False
	
	def saveSelection(self, start=0, finish=0, zero=-1, featureHit=False):

		# feature to highlight or drag and drop?
		if featureHit != False:
			featurelist             = genbank.gb.get_all_feature_positions()
			for i in range(0,len(featurelist)):
				# load all the feature infos
				featuretype, complement, start, finish, name, index = featurelist[i]

				hName = self.hitName(name, index)
				if featureHit == hName:
					# a feature with start smaller then finish, is a featur laying over zero
					if start > finish:
						zero = 1 # 1 --> feature starts left and ends right of +1
					# set selection for real!
					self.set_dna_selection((start,finish, zero))
		else:
			self.set_dna_selection((start,finish, zero))
		
		return None
	
	def nameBeautiful(self, name):
		# remove any '"' in front or at the end
		if name[0:1] == '"':
			name = name[1:]
		if name[-1:] == '"':
			name = name[:-1]
		return name

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
	
	
	############### Inteaction with mosue and keyboard ################
	def OnLeftUp(self, event):
		''' handle left ouseclick up'''
		# get selection
		pos1, pos2, zero 	= self.plasmidstore.interaction["selection"]
		if self.plasmidstore.interaction["leftDown"] == True and pos2 != -1:
			# selection is finish
			self.plasmidstore.interaction["leftDown"] = False
			# save the mouse position
			x, y 		= self.ScreenToClient(wx.GetMousePosition())
			x2, y2 		= self.ctx.device_to_user(x,y)
			pos 		= self.cartesian2position(x2,y2)
			
			self.plasmidstore.interaction["selection"] = (pos1, pos, zero)

		elif pos2 == -1:
			# reset selection
			self.plasmidstore.interaction["selection"] = (1, -1, -1)

			# maybe we cliked on a feature, label or a line?
			hit = self.HitTest()
			if hit != None:
				self.plasmidstore.interaction["hit"] = hit
				# we hit one, so we need to draw the selection!
				# wich feature was hit?
				self.saveSelection(1, -1, -1, hit)
			else:
				self.plasmidstore.interaction["hit"] = None
		
		# update the UI maybe?
		self.update_ownUI()

	
	def OnLeftDown(self, event):
		# remove the position cursor:
		self.set_cursor_position(0)
		
		selection = self.HitTestSelection()
		if selection:
			# save the mouse position
			x, y = self.ScreenToClient(wx.GetMousePosition())
			x2, y2 = self.ctx.device_to_user(x,y)
			pos = self.cartesian2position(x2,y2)
			self.plasmidstore.interaction["leftDown"] = True
			self.plasmidstore.interaction["selection"] = (pos, -1, -1)
		else:
			# remove selection
			self.plasmidstore.interaction["hit"] = None
			self.saveSelection(1, -1, -1, False)
		
		return None
		
		
	def OnLeftDouble(self, event):
		'''When left button is double clicked, launch the feature edit dialog.'''
		hit = self.HitTest() #this does not get the "true" feature index. Some featues are split and this is an index that accounts for that.
		if hit is not False:
			# get the index of this feature
			featurelist             = genbank.gb.get_all_feature_positions()
			for i in range(0,len(featurelist)):
				# load all the feature infos
				featuretype, complement, start, finish, name, index = featurelist[i]
				hName = self.hitName(name, index)
				if hit == hName:
					genbank.feature_selection = copy.copy(index)

			dlg = featureedit_GUI.FeatureEditDialog(None, 'Edit Feature') # creation of a dialog with a title
			dlg.ShowModal()
			dlg.Center()
		
		
	def OnMotion(self, event):
		'''When mouse is moved with the left button down determine the DNA selection from angle generated at mouse down and mouse move event.'''
		
		
		if event.Dragging() and event.LeftIsDown() and self.plasmidstore.interaction["leftDown"] == True:
			# save the mouse position
			x, y 				= self.ScreenToClient(wx.GetMousePosition())
			x2, y2 				= self.ctx.device_to_user(x,y)
			pos 				= self.cartesian2position(x2,y2)
			pos1, pos2, zero 	= self.plasmidstore.interaction["selection"]
			self.plasmidstore.interaction["selection"] = (pos1, pos, zero)
			
			# if we have a sudden increase, we might have moven over zero
			if abs(pos2 - pos) > self.dnaLength/2 and pos2 != -1:
				zero = zero * -1 # turn zero around
				self.plasmidstore.interaction["selection"] = (pos1, pos, zero)

			#set genebank selection
			self.saveSelection(pos1, pos, zero)
			self.update_ownUI()



		

		# check if something is beneath the mousecursor
		oldhit 	= self.Highlight
		hit 	= self.HitTest()
		# save the state
		self.Highlight = hit
		# only update after change
		if oldhit != hit:
			self.update_ownUI()

				
		
		
		#elif len(self.selectionDrawing) < 2: # only hover position, if none is selected --> maybe use two variables to make both possible?
			# no feature to highlight, but maybe selection hover?
		#	self.Highlight = None


			# check if we are in the "hover area"
		#	inArea = self.HitTestMarkerArea()

		#	if inArea == True:
		#		# get the mouse position
		#		x, y = self.ScreenToClient(wx.GetMousePosition())
		#		x2, y2 = self.ctx.device_to_user(x,y)
		#		self.selectionDrawing = [] # empty it
		#		self.selectionDrawing.append(self.cartesian2radial(x2,y2))
		#		self.update_ownUI()
		#	elif oldhit != hit:
		#		self.update_ownUI()
		#	elif len(self.selectionDrawing) == 1: # TODO
		#		self.selectionDrawing = [] # empty it so the hover line wont be seen all the time
		#		self.update_ownUI()




	
	
	############### Setting methods for converting stuff ##############

	def radial2cartesian(self,radius, angle, cx=0,cy=0): # angle in radial!
		radius = float(radius)
		angle  = float(angle)
		x = radius * math.cos(angle) + cx
		y = radius * math.sin(angle) + cy
		return x,y

	def cartesian2radial(self,x, y):
		x 		= -x # cairo somehow swaps the x axis
		diff 	=  -math.radians(90)	# remove some radians, because 0 is on top
		r 		= math.atan2(x,y) + 1 * math.pi + diff
		return r

	def position2angle(self, position):
		radMult        	= float(2*math.pi/self.dnaLength)
		angle = (position - 0.5 * self.dnaLength) * radMult + math.pi/2
		return angle
		
	def position2cartesian(self, position, radius):
		angle = self.position2angle(position)
		x,y = self.radial2cartesian(radius, angle)
		return x, y
	
	def cartesian2position(self, x,y):
		angle 	= self.cartesian2radial(x,y)
		radMult	= float(2*math.pi/self.dnaLength)
		pos 	= 0.5 * self.dnaLength  + (angle - math.pi/2)/radMult 
		pos 	= int(round(pos, 0))
		return pos
	


	def merge_dicts(self, *dict_args):
		'''
		Given any number of dicts, shallow copy and merge into a new dict,
		precedence goes to key value pairs in latter dicts.
		'''
		result = {}
		for dictionary in dict_args:
		    result.update(dictionary)
		return result


	########################################
	# caclcualting functions
	def checkFeatureCalculations(self):
		''' this functions takes the current positions and calcualtions
			of the features and checks fast if it has to recalculated anything'''
		redrawfeatures = False

		
		
		
		# compare the features to the previously calculated ones
		featuresOld = self.plasmidstore.features
		featuresNew = genbank.gb.get_all_feature_positions()
		
		if featuresNew != None:
			# get the length once
			self.dnaLength 			= float(len(genbank.gb.GetDNA()))
		else:
			self.dnaLength 			= 1 # or 0 but this might cause problems if there is a division x/length

		# compare old and new
		if featuresOld != featuresNew:

			# print ("DO recalc features")
			# we only have to redraw if the length, the color or direction changed
			# calculate and save all the features
			allFeatures = {} # storing list
			# some basic info we have to update:
			
			self.drawn_fw_locations = []
			self.drawn_rv_locations = []
			

			labelFeature = [] 	# labels for features outside of the arrow
								# save as [[middle, name, index, hitname, y],[...]]

			
			# loop the features and get the path
			for feature in 	featuresNew:
				featuretype, complement, start, finish, name, index = feature
				hitname 	= self.hitName(name, index)
				
				# do not display the annoying source
				if featuretype != "source": 
					# path creation
					allFeatures[hitname] = [self.cairoFeature(feature), featuretype]
				
					# label handling
					length 		= self.dnaLength
					middle 		= start + (finish-start)/2 % length
				
					name		= self.nameBeautiful(name)
					y 			= None # to be populated
					labelFeature.append([middle, name, index, hitname, y])
				
			# save the new paths in the fancy storing class for further processing
			self.plasmidstore.drawnfeatures = allFeatures
			
			# save label drawn outside as text
			self.plasmidstore.LoF = labelFeature

		
		


		self.plasmidstore.features = featuresNew # everything is over and the new are the current (old)
		return redrawfeatures
		
		
	def cairoFeature(self, feature):
		featurepath = None
		featuretype, complement, start, finish, name, index = feature
		name 	= self.nameBeautiful(name)
		hname 	= self.hitName(name, index)
		

		
		# create the correct path for foreward and reverse
		if complement == False:
			self.drawn_fw_locations, levelMult = self.find_overlap(self.drawn_fw_locations, (start, finish))
			arrowMult 	= 1 # to account for the different direction of addition
			radiusAdd	= levelMult * (self.arrowWidth + 0.5) + 0.4
			radiusO 	= self.radiusO + radiusAdd + self.arrowWidth
			radiusI 	= self.radiusO + radiusAdd
			
			
		else:
			self.drawn_rv_locations, levelMult = self.find_overlap(self.drawn_rv_locations, (start, finish))		
			
			arrowMult 	= -1 # to account for the different direction of addition
			radiusAdd	= levelMult * (self.arrowWidth + 0.5) + 0.4
			radiusO 	= self.radiusI - radiusAdd - self.arrowWidth
			radiusI 	= self.radiusI - radiusAdd
			
		# calculate the radians for end and start
		s = self.position2angle(start)
		e = self.position2angle(finish)

		
		# large features get a arrow head:
		arrowHead = False
		if abs(e - s) > math.radians(5): # more then x degrees
			arrowHead = True
		
		
		if arrowHead == True and complement == False:
			# arrow head has x radians
			r 		= radiusI + arrowMult * self.arrowWidth/2
			xA, yA 	= self.radial2cartesian(r, e)
			e 		= e - math.radians(2) # new end
			
			self.ctx.new_path()
			self.ctx.arc(0,0,radiusI, s,e)
		
			# draw arrow head sometimes:
			self.ctx.line_to(xA, yA)
			#self.ctx.move_to(cx1, cy1)
			self.ctx.arc_negative(0,0,radiusO, e,s)
			self.ctx.close_path()
			featurepath  = self.ctx.copy_path()
			self.ctx.new_path() # clear canvas
		elif arrowHead == True and complement == True:
			# arrow head has x radians
			r 		= radiusI + arrowMult * self.arrowWidth/2
			xA, yA 	= self.radial2cartesian(r, s)
			s 		= s + math.radians(2) # new "start"
			
			self.ctx.new_path()
			self.ctx.arc_negative(0,0,radiusI, e,s)
			self.ctx.line_to(xA, yA)
			self.ctx.arc(0,0,radiusO, s,e)
			self.ctx.close_path()
			featurepath  = self.ctx.copy_path()
			self.ctx.new_path() # clear canvas
		else:
			self.ctx.new_path()
			self.ctx.arc(0,0,radiusI, s,e)
			self.ctx.arc_negative(0,0,radiusO, e,s)
			self.ctx.close_path()
			featurepath  = self.ctx.copy_path()
			self.ctx.new_path() # clear canvas
		
		# save the outer radius to allow draing of the lines:
		if complement == False:
			self.plasmidstore.radiusOuter[hname] = radiusO
		else:
			self.plasmidstore.radiusOuter[hname] = radiusI
		
		return featurepath
	
	
	
	def caironameLabels(self, labels, site, kind):
		''' create paths for labels given'''
		
		if kind == "features":
			self.ctx.select_font_face('Arial', cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD)
		else:
			self.ctx.select_font_face('Arial', cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
		self.ctx.set_font_size(self.lH)
		textPaths			= {}
		linePaths 			= {}
		underlinePaths 		= {}
		boxPaths			= {}
		
		shift = 4 # shift the labels closer to the plasmid, as the have a bigger radius
		
		
		for l in labels:
			name 	= l[1]
			y 		= l[4]
			if kind == "enzymes":
				name2		= "%s%s" %(name, self.enzymeGap)
			else:
				name2		= name
			hname 	= self.hitName(name2, l[2])

					
			xbearing, ybearing, TextWidth, TextHeight, xadvance, yadvance = self.ctx.text_extents(name)

			if site == "left":
				x = - math.cos(math.asin(y/self.radiusLabels)) * self.radiusLabels - TextWidth + shift
			else:
				x = math.cos(math.asin(y/self.radiusLabels)) * self.radiusLabels - shift
			self.ctx.move_to(x, y);
			self.ctx.text_path(name)
						
			# save label:
			path		 = self.ctx.copy_path()
			textPaths[hname] = path
			self.ctx.new_path() # clear canvas
			
			# make a box behind the text, so we can select more easys:
			self.ctx.rectangle(x, y - self.lH, TextWidth, self.lH)
			boxPaths[hname] = self.ctx.copy_path()
			self.ctx.new_path()
			
			# featurs have a underline:
			if kind == "features":
				self.ctx.move_to(x, y)
				self.ctx.line_to(x + TextWidth, y)
				path		 = self.ctx.copy_path()
				underlinePaths[hname] = path
				self.ctx.new_path() # clear canvas
			
			# make a line from text to plasmid:
			if site == "left":
				x = x + TextWidth
			else:
				x = x
			
			# get the radius to where the line should lead:
			try:
				r = self.plasmidstore.radiusOuter[hname]
			except:
				r = self.radiusI
					
			y1 = y - self.lH/2 + 0.3
			self.ctx.move_to(x, y1)
			x2, y2 = self.position2cartesian(l[0], self.radiusO + 7)
			x3, y3 = self.position2cartesian(l[0], r)
			self.ctx.line_to(x2,y2)
			self.ctx.line_to(x3,y3)
			linePaths[hname] = self.ctx.copy_path() # save line path
			self.ctx.new_path() # clear canvas
			
			
		return textPaths, linePaths, boxPaths, underlinePaths

	def checkLabelCalculations(self):
		''' this functions takes the current positions and calcualtions
			of the labels and checks fast if it has to recalculated anything'''
		# labels to handle:
		LoF = self.plasmidstore.LoF
		LoE = self.plasmidstore.LoE
		
		# check if the labels changed, onyl recalculate everything, if they did
		if LoF != self.plasmidstore.LoFold or LoE != self.plasmidstore.LoEold:
			# print("DO recalculate label")
			# sort labels for their middle position:
			LoF = sorted(LoF, key=lambda x: x[0]) 
			LoE = sorted(LoE, key=lambda x: x[0])
		
			# group the labels
			# TODO, should be interactive anyway
		
			# now split in two sites:
			half = self.dnaLength/2
			# Features
			LoF1 = [] # first half
			LoF2 = [] # second half
			LoF1I = [] # indexing
			LoF2I = [] # indexing
			for l in LoF:
				if l[0] < half:
					LoF1.append(l)
					LoF1I.append(l[3])
				else:
					LoF2.append(l)
					LoF2I.append(l[3])
			# Enzymes
			LoE1 = [] # first half
			LoE2 = [] # second half
			LoE1I = [] # indexing
			LoE2I = [] # indexing
			for l in LoE:
				if l[0] < half:
					LoE1.append(l)
					LoE1I.append(l[3])
				else:
					LoE2.append(l)
					LoE2I.append(l[3])
			

		
			# then check the numbers
			LoF1, LoE1 = self.checkNumberOfLabels(LoF1, LoE1)
			
			# arrange and draw labels:
			Labels1 = self.arrangeLabels(LoF1, LoE1, -1)
			Labels2 = self.arrangeLabels(LoF2, LoE2, +1)
			
			# sperate into features and enzymes... 
			# this must be, so we can order them together, but color them
			# different
			LabelsFeatures1 	= []
			LabelsFeatures2 	= []
			LabelsEnzymes1	 	= []
			LabelsEnzymes2 		= []
			for label in Labels1:
				hit = label[3]
				if hit in LoF1I:
					LabelsFeatures1.append(label)
				if hit in LoE1I:
					LabelsEnzymes1.append(label)
			for label in Labels2:
				hit = label[3]
				if hit in LoF2I:
					LabelsFeatures2.append(label)
				if hit in LoE2I:
					LabelsEnzymes2.append(label)
			# done serperating the labels agai				

			
			# FEATURES: make the paths for features
			laP1, liP1, bP1, UlP1	= self.caironameLabels(LabelsFeatures1, "right", "features")
			laP2, liP2, bP2, UlP2  	= self.caironameLabels(LabelsFeatures2, "left", "features")
			# FEATURES: save paths for displaying
			self.plasmidstore.drawnlabels 			= self.merge_dicts(laP1, laP2) 
			self.plasmidstore.drawnLLines 			= self.merge_dicts(liP1, liP2) 
			self.plasmidstore.labelBoxes 			= self.merge_dicts(bP1, bP2) 
			self.plasmidstore.labelULine 			= self.merge_dicts(UlP1, UlP2) 
			
			
			
			# ENZYMES: make the paths for enzymes
			laP1, liP1, bP1, UlP1	= self.caironameLabels(LabelsEnzymes1, "right", "enzymes")
			laP2, liP2, bP2, UlP2 	= self.caironameLabels(LabelsEnzymes2, "left", "enzymes")
			# ENZYMES: save paths for displaying
			self.plasmidstore.drawnlabelsE 			= self.merge_dicts(laP1, laP2) 
			self.plasmidstore.drawnLLinesE 			= self.merge_dicts(liP1, liP2) 
			self.plasmidstore.labelBoxesE 			= self.merge_dicts(bP1, bP2) 

			
			# save the new status quo of the labels 
			self.plasmidstore.LoFold = self.plasmidstore.LoF
			self.plasmidstore.LoEold = self.plasmidstore.LoE
			
			return None
		else:
			# print("NOT recalculate label")
			return None

	
	def checkNumberOfLabels(self, LoF=None, LoE=None):
		''' this function is called before label sorting'''
		
	
		# determine height of label displaying area
		height = self.radius * 2.2 
		# divide by the height of each label
		lH = self.ctx.device_to_user_distance(self.labelHeight,0)[0] # label height
		lH = round(lH, 3)
		self.lH = lH # for other functions to use
		# round down, to make it work
		nL = math.trunc(height/lH)
	
		# number of labels to handle
		nFL = len(LoF)
		nEL = len(LoE)
	
		if nFL + nEL <= nL:
			# no problem here
			return LoF, LoE # return feature Label and Text Label
		elif nFL <= nL:
			# reduce the nTL so, that it may work.
			# reduce most populated areas
			# print("group some Restrictionsite labels agressivly")
			# call function and check again if we can fit everything in
			# if it did work--> fine
			# else show it anyway --> or figure something out
			return LoF, LoE # return feature Label and Text Label
		elif nFL > nL:
			# we have to many feature labels.
			# so we hide all TL and remove the labels of the smallest features!
			# print("delete some Restrictionsite labels")
			return LoF, LoE # return feature Label and Text Label
		else:
			# print("check Number of labels: undefined")
			return LoF, LoE # return feature Label and Text Label
	

	#def groupLabels(self, labels):
	#	''' restrictionenzymelabels can be grouped 
	#		as they somtimes occure at the very same place '''
	#	lastPos = None
	#	i		= 0
	#	for l in labels:
	#		pos = l[0]
	#
	#		if pos == lastPos:
	#			# two labels at the very same position
	#			# we merge the names
	#			cname 	= l[1]
	#			oname 	= labels[i][1]
	#			name 	= "%s, %s"
	#		i = i + 1
	#	# we can group labels of same type and same position
	#	return None

	def groupAgressivly(self):
		# we can group some close restrictionlabels 
		# and only show the position indicating by a number
		return None

	def arrangeLabels(self, F, E, index):
		''' this function does the actual positioning of all labels'''
		
		# now they are equal
		# enzymes, features etc: merge
		Labels = F + E
		Labels = sorted(Labels, key=lambda x: x[0]) 
		
		# find the most horiziontal one
		horzI = None
		horzN = self.dnaLength
		if index == -1:
			target = self.dnaLength / 4 * 1
		else:
			target = self.dnaLength / 4 * 3

		i = 0
		# start by putting the initial position
		for l in Labels:
			# append a new entry with the current y coordinate:
			x,y = self.position2cartesian(l[0], self.radiusLabels - 2)
			l[4] =  round(y, 3)
			
			# most horizontal one:
			if abs(l[0] - target) < horzN:
				horzN = abs(l[0] - target) 	# save new differenze
				horzI = i					# save index
			
			i = i + 1
		
		
		
		# create series for loop for better positioning
		# this way we start positioning in the hotizontal position and 
		# arrange the labels in both directions more accurate
		h 		= 1
		loop 	= [horzI]
		newI 	= horzI
		switch 	= 1
		both	= True
		while h < len(Labels):
			if newI > 0 and newI < len(Labels)- 1 and both == True:
				newI = newI + switch * h
				switch = -1 * switch
			elif newI == 0:
				switch = 1
				newI = newI + switch * h
				both = False
			elif newI == len(Labels) - 1:
				switch = -1
				newI = newI + switch * h
				both = False
			else:
				newI = newI + switch * 1

			loop.append(newI)
			h = h + 1

		# now iterate n times and check if x[i] overlaps with either x[i+1] or x[i-1]
		# this is the key element of the label positioning algorythm
		# n is critical for performance, higher n = more steps
		# dy is the step size in px
		dy 			= self.ctx.device_to_user_distance(5,0)[0] # amount of movement in pixel
		i 			= 0 		# just a counter
		n 			= 100		# n of steps
		overlapp 	= 1			# initial overlap, so it is not 0
		while i < n and overlapp > 0:

			overlapp = 0
			#for l in Labels:
			for a in loop:
				l = Labels[a]
				# for those labels with one upper and lower label
				y = l[4]			# y
				U = y + self.lH		# upper border
				L = y				# lower border
				
				# get the correct upper and lower bound of the adjacend features or the end of the radius
				if index > 0:
					if a == 0:
						# next label
						nL = +self.radiusLabels
					else:
						# next label
						nextLabel = Labels[a - 1]
						nL = nextLabel[4]		# next lower
						
					if a == len(Labels)-1:
						lU = -self.radiusLabels
					else: 
						# previous label
						prevLabel = Labels[a + 1]
						lU = prevLabel[4] + self.lH # last upper
				else:
					if a == len(Labels)-1:
						# next label
						nL = self.radiusLabels
					else:
						# next label
						nextLabel = Labels[a + 1]
						nL = nextLabel[4]		# next lower
						
					if a == 0:
						lU = -self.radiusLabels
					else: 
						# previous label
						prevLabel = Labels[a - 1]
						lU = prevLabel[4] + self.lH # last upper
				
				
				
				# check upper bound
				dU = U - nL
				# check lower bound 
				dL = lU - L
	

				# move the item up or down, according to its overlapping
				if dU > 0 or dL > 0:
					
					# we need to move somewhere
					if dU > 0 and dL <= 0:
						# move down
						y = y - dy

						# raise the overlap
						overlapp = overlapp + dU
					elif dU <= 0 and dL > 0:
						# we need to move up
						y = y + dy

						# raise the overlap
						overlapp = overlapp + dL
					elif dU > 0 and dL > 0:
						# we need to make them even
						y = y + (dU - (dU + dL)/2)

						# raise the overlap
						overlapp = overlapp + dU + dL
					
					# prevent to large y
					if y > self.radiusLabels:
						y = self.radiusLabels
					elif y < -self.radiusLabels:
						y = - self.radiusLabels
					

					# save y:
					l[4] = y
	
			# print("DO reposition labels, Round: ", i)
			i = i + 1
	
		# after n cycles stop calculations and return the positions
		return Labels

	
	
	def checkEnzymeCalculations(self):
		
		# load enzymes
		enzymes 	=  genbank.gb.restrictionEnzymes.selection
		# create hash of selection and DNA
		hashable 	= "%s%s" % (frozenset(enzymes), genbank.gb.gbfile['dna'])
		e1 			= hashlib.md5(hashable).hexdigest()

		labels = []

		if e1 !=self.plasmidstore.enzymes:
			# we have new enzymes
			index = 0
			for e in enzymes:
				for site in enzymes[e].restrictionSites:
					name, start, end, cut51, cut52, dnaMatch = site
					length 		= self.dnaLength
					middle 		= start + (end-start)/2 % length
					name		= self.nameBeautiful(name)
					name2		= "%s%s" %(name, self.enzymeGap)
					hitname 	= self.hitName(name2, index)
					y 			= None # to be populated
					labels.append([middle, name, index, hitname, y])
					
					index = index + 1
			self.plasmidstore.enzymes= e1
			self.plasmidstore.LoE = labels
		
		return None
	
	
	def checkSelectionCalculations(self):
		''' only chnage the arc if we changed something '''
		# selection
		selOld 		= self.plasmidstore.interaction["selectionOld"]
		sel 		= self.plasmidstore.interaction["selection"]

		
		if sel != selOld:
			if sel[1] == -1:
				# remove selection arc
				self.plasmidstore.drawnSelection = None
			else:
				# draw selection arc
				pos1, pos2, zero = sel

				a1 = self.position2angle(pos1)
				a2 = self.position2angle(pos2)
				
				if (pos1 < pos2 and zero == -1) or (pos2 < pos1 and zero == 1):
					self.ctx.arc(0,0,self.radius,a1,a2)
				else :
					self.ctx.arc_negative(0,0,self.radius,a1,a2)
				path = self.ctx.copy_path() # save line path
				self.ctx.new_path() 		# clear canvas
				self.plasmidstore.drawnSelection = path
			self.plasmidstore.interaction["selectionOld"] = sel
		
		
		# draw the line for indicating the position
		position 	= self.plasmidstore.interaction["position"]
		positionOld	= self.plasmidstore.interaction["positionOld"]

		if position != positionOld:
			if position < self.dnaLength and position > 0:
				# get the coordinates for the line
				x1, y1 = self.position2cartesian(position, self.radiusI - 1)
				x2, y2 = self.position2cartesian(position, self.radiusO + 1)
				# draw the actual line
				self.ctx.move_to(x1,y1)
				self.ctx.line_to(x2,y2)
				# copy, save and delete
				path = self.ctx.copy_path()
				self.plasmidstore.drawnPosition = path
				self.ctx.new_path() 		
				# save settings
				
			else:
				self.plasmidstore.drawnPosition = None
		
		self.plasmidstore.interaction["positionOld"] = position
		return None
	
	############################################
	# drawing functions

	def draw(self):
		''' # print the stored elements on th canvas'''
		# prepare the canvas
		width, height = self.GetVirtualSize()
		if width < height:
			ratio = float(width)/float(height)
			# only resize the canvas if software is loaded. This prevents error messages
			if width != 0 and height != 0 and width > 100 and height > 100:
				self.ctx.scale(width, height*ratio) # Normalizing the canvas
				self.ctx.scale(0.01, 0.01) 			# make it 100*ratiox100
				self.ctx.translate (50,50/ratio) 		# set center to 0,0
		if width >= height:
			ratio = float(height)/float(width)
			# only resize the canvas if software is loaded. This prevents error messages
			if width != 0 and height != 0 and width > 100 and height > 100:
				self.ctx.scale(width*ratio, height) # Normalizing the canvas
				self.ctx.scale(0.01, 0.01) 			# make it 100*ratiox100
				self.ctx.translate (50/ratio,50) 		# set center to 0,0
		self.drawFeatures()
		self.drawLabels()
		self.drawSelection()
		return None

	def drawFeatures(self):
		# print("DO output path feature")
		# only resize the canvas if software is loaded. This prevents error messages

		# draw the helper plain white circle to make selection possible:
		self.ctx.arc(0,0,self.radiusO+1,0,2*math.pi)		# outer circle
		self.ctx.set_source_rgb (1, 1, 1) 				# white
		self.ctx.set_line_width (0)						# no line
		self.plasmidstore.interaction["markerArea1"] = self.ctx.copy_path()	# copy for later
		self.ctx.new_path()								# remove

		self.ctx.arc(0,0,self.radiusI-1,0,2*math.pi)		# inner circle
		self.ctx.set_source_rgb (1,1,1)					# Solid white
		self.ctx.set_line_width (0)						# no line
		self.plasmidstore.interaction["markerArea2"] = self.ctx.copy_path()	# copy for later
		self.ctx.new_path()								# remove

		# draw the plasmid (two circles with different radius)
		radius=10
		self.ctx.arc(0,0,self.radiusO,0, 2*math.pi)		# circle
		self.ctx.set_source_rgb (0, 0, 0) 			# Solid color
		self.ctx.set_line_width (0.24)				# line width
		self.ctx.stroke()					# stroke only no fill!

		# inner plasmid
		self.ctx.arc(0,0,self.radiusI,0,2*math.pi)		# inner circle
		self.ctx.stroke()					# stroke the same color and settings


		# draw the buffered featrues:
		i  = 0
		for a in self.plasmidstore.drawnfeatures:
			path 		= self.plasmidstore.drawnfeatures[a][0]
			# get color
			featuretype = self.plasmidstore.drawnfeatures[a][1]
			color = eval(featuretype)['fw'] #get the color of feature (as string)
			assert type(color) == str
			r,g,b = colcol.hex_to_rgb(color)
			r = float(r)/255
			g = float(g)/255
			b = float(b)/255

			# put the path on the canvas and fill
			#self.ctx.new_path()
			self.ctx.append_path(path)
			self.ctx.set_source_rgba (r,g,b,1.0) # Solid color
			self.ctx.fill_preserve()
			
			# outline for the feature that is beneeth the mouse:
			if a == self.Highlight:
				#self.ctx.set_source_rgb (1, 0, 0)
				self.ctx.set_line_width(0.3)
				self.ctx.stroke()
			else:
				# remove path
				self.ctx.new_path()
			
			i = i + 1
		
		return None



	def drawLabels(self):
		''' function to put the labels on the canvas'''
		
		
		# FEATURES
		self.ctx.set_source_rgb (0, 0, 0) 			# Solid color
		self.ctx.set_line_width(0.1)
		# draw feature label names
		for a in self.plasmidstore.drawnlabels:
			path = self.plasmidstore.drawnlabels[a]
			self.ctx.append_path(path)
			# darker color for the feature that is beneeth the mouse:
			if a == self.Highlight:
				self.ctx.set_source_rgba (0, 0, 0, 1)
			else:
				self.ctx.set_source_rgba (0, 0, 0, 0.7)
			self.ctx.fill()
		
		# draw the lines
		for a in self.plasmidstore.drawnLLines:
			path = self.plasmidstore.drawnLLines[a]
			self.ctx.append_path(path)
			
			# darker color for the feature that is beneeth the mouse:
			if a == self.Highlight:
				self.ctx.set_source_rgba (0, 0, 0, 1)
			else:
				self.ctx.set_source_rgba (0, 0, 0, 0.3)
			self.ctx.stroke()
			
		# draw the underlines
		#for a in self.plasmidstore.labelULine:
		#	path = self.plasmidstore.labelULine[a]
		#	self.ctx.append_path(path)
		#	self.ctx.set_line_width(2)
		#	
		#	# darker color for the feature that is beneeth the mouse:
		#	if a == self.Highlight:
		##		self.ctx.set_source_rgba (0, 0, 0, 1)
		#	else:
		#		self.ctx.set_source_rgba (0, 0, 0, 0.3)
		#	self.ctx.stroke()
		
		
		# ENZYMES
		self.ctx.set_source_rgb (0, 0, 0) 			# Solid color
		self.ctx.set_line_width(0.1)
		# draw feature label names
		for a in self.plasmidstore.drawnlabelsE:
			path = self.plasmidstore.drawnlabelsE[a]
			self.ctx.append_path(path)
			# darker color for the feature that is beneeth the mouse:
			if self.Highlight != None and a.split(self.enzymeGap)[0] == self.Highlight.split(self.enzymeGap)[0]:
				self.ctx.set_source_rgba (0, 0, 0, 1)
			else:
				self.ctx.set_source_rgba (0, 0, 0, 0.7)
			self.ctx.fill()
		
		# draw the lines
		for a in self.plasmidstore.drawnLLinesE:
			path = self.plasmidstore.drawnLLinesE[a]
			self.ctx.append_path(path)
			
			# darker color for the feature that is beneeth the mouse:
			if self.Highlight != None and a.split(self.enzymeGap)[0] == self.Highlight.split(self.enzymeGap)[0]:
				self.ctx.set_source_rgba (0, 0, 0, 1)
			else:
				self.ctx.set_source_rgba (0, 0, 0, 0.3)
			self.ctx.stroke()
		
	def drawSelection(self):
		''' draw arc for selection '''
		if self.plasmidstore.drawnSelection != None:
			#0082ED
			r,g,b = colcol.hex_to_rgb('#0082ED') 

			r = float(r)/255
			g = float(g) /255
			b = float(b) /255

			self.ctx.set_source_rgba (r,g,b, 0.6) 			# Solid color
			self.ctx.set_line_width(3)
			path = self.plasmidstore.drawnSelection
			self.ctx.append_path(path)
			self.ctx.stroke()
		
		if self.plasmidstore.drawnSelection == None and self.plasmidstore.drawnPosition != None:
			self.ctx.set_source_rgba (0,0,0) 			# Solid color
			self.ctx.set_line_width(0.2)
			path = self.plasmidstore.drawnPosition
			self.ctx.append_path(path)
			self.ctx.stroke()



class PlasmidView2(DNApyBaseDrawingClass):
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

		imageFile 	= files['default_dir']+"/icon/circle.png"
		image1 		= wx.Image(imageFile, wx.BITMAP_TYPE_ANY).ConvertToBitmap()
		circle 		= wx.BitmapButton(panel1, id=10, bitmap=image1, size = (image1.GetWidth()+padding, image1.GetHeight()+padding), name = "share")

		imageFile 	= files['default_dir']+"/icon/group.png"
		image1 		= wx.Image(imageFile, wx.BITMAP_TYPE_ANY).ConvertToBitmap()
		group 		= wx.BitmapButton(panel1, id=11, bitmap=image1, size = (image1.GetWidth()+padding, image1.GetHeight()+padding), name = "share")

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
		#self.plasmid_view = PlasmidView(panel2, -1)
		self.plasmid_view = drawPlasmid(panel2, -1)

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
	exec(compile(open(settings).read(), settings, 'exec')) #gets all the pre-assigned settings

	genbank.dna_selection = (1, -1)	 #variable for storing current DNA selection
	genbank.feature_selection = False #variable for storing current feature selection

	import sys
	assert len(sys.argv) == 2, 'Error, this script requires a path to a genbank file as an argument.'
	# print('Opening %s' % str(sys.argv[1]))

	genbank.gb = genbank.gbobject(str(sys.argv[1])) #make a genbank object and read file


	app = MyApp(0)
	app.MainLoop()
