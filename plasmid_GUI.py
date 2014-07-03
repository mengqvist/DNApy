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
#fix labels so they don't overlap
#fix labels so they don't run offscreen
#fix long labels
#fix long plasmid names
#add 'dna ruler'

import wx
import wx.lib.graphics
from wx.lib.pubsub import pub


import math
import genbank
import copy

import os, sys
import string
from base_class import DNApyBaseClass
import featureedit_GUI

files={}   #list with all configuration files
files['default_dir'] = os.path.abspath(os.path.dirname(sys.argv[0]))+"/"
files['default_dir']=string.replace(files['default_dir'], "\\", "/")
files['default_dir']=string.replace(files['default_dir'], "library.zip", "")
settings=files['default_dir']+"settings"   ##path to the file of the global settings
execfile(settings) #gets all the pre-assigned settings

class Base(DNApyBaseClass):
	def __init__(self, parent, id):
		wx.Panel.__init__(self, parent)

		#determing which listening group from which to recieve messages about UI updates
		self.listening_group = 'from_feature_list' #needs to be assigned or will raise an error		
		pub.Publisher.subscribe(self.listen_to_updateUI, self.listening_group)

		self.listening_group2 = 'from_feature_edit'		
		pub.Publisher.subscribe(self.listen_to_updateUI, self.listening_group2)		

		self.listening_group3 = 'from_dna_edit'		
		pub.Publisher.subscribe(self.listen_to_updateUI, self.listening_group3)	

		self.listening_group4 = 'from_main'
		pub.Publisher.subscribe(self.listen_to_updateUI, self.listening_group4)

		self.listening_group5 = 'private_group_for_those_that_affect_DNA_selection_from_DNA_editor'
		pub.Publisher.subscribe(self.listen_to_updateUI, self.listening_group5)


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


class BufferedWindow(Base):

    """

    A Buffered window class.

    To use it, subclass it and define a Draw(DC) method that takes a DC
    to draw to. In that method, put the code needed to draw the picture
    you want. The window will automatically be double buffered, and the
    screen will be automatically updated when a Paint event is received.

    When the drawing needs to change, you app needs to call the
    UpdateDrawing() method. Since the drawing is stored in a bitmap, you
    can also save the drawing to file by calling the
    SaveToFile(self, file_name, file_type) method.

    """

    def __init__(self, *args, **kwargs):
        # make sure the NO_FULL_REPAINT_ON_RESIZE style flag is set.
 
#        kwargs['style'] = kwargs.setdefault('style', wx.NO_FULL_REPAINT_ON_RESIZE) | wx.NO_FULL_REPAINT_ON_RESIZE
        Base.__init__(self, *args, **kwargs)

        self.Bind(wx.EVT_PAINT, self.OnPaint)
        self.Bind(wx.EVT_SIZE, self.OnSize)

        # OnSize called to make sure the buffer is initialized.
        # This might result in OnSize getting called twice on some
        # platforms at initialization, but little harm done.
        self.OnSize(None)
        self.paint_count = 0


    def Draw(self, dc):
        ## just here as a place holder.
        ## This method should be over-ridden when subclassed
        raise NotImplementedError

    def OnPaint(self, event):
        # All that is needed here is to draw the buffer to screen
        dc = wx.BufferedPaintDC(self, self._Buffer)

    def OnSize(self,event):
        # The Buffer init is done here, to make sure the buffer is always
        # the same size as the Window
        #Size  = self.GetClientSizeTuple()
        Size  = self.ClientSize
        self.size = Size

        # Make new offscreen bitmap: this bitmap will always have the
        # current drawing in it, so it can be used to save the image to
        # a file, or whatever.
        self._Buffer = wx.EmptyBitmap(*Size)
        self.update_ownUI()

    def SaveToFile(self, FileName, FileType=wx.BITMAP_TYPE_PNG):
        ## This will save the contents of the buffer
        ## to the specified file. See the wxWindows docs for 
        ## wx.Bitmap::SaveFile for the details
        self._Buffer.SaveFile(FileName, FileType)




class PlasmidView(BufferedWindow):
	def __init__(self, *args, **kwargs):	
		genbank.dna_selection = (1,1)

		self.centre_x = 0
		self.centre_y = 0
		self.highlighted_feature = False
		BufferedWindow.__init__(self, *args, **kwargs)
	

		self.Bind(wx.EVT_LEFT_DOWN, self.OnLeftDown)
		self.Bind(wx.EVT_LEFT_UP, self.OnLeftUp)
		self.Bind(wx.EVT_RIGHT_UP, self.OnRightUp)
		self.Bind(wx.EVT_MOTION, self.OnMotion)
		self.Bind(wx.EVT_LEFT_DCLICK, self.OnLeftDouble)

	def find_overlap(self, drawn_locations, new_range):
		'''Takes two ranges and determines whether the new range has overlaps with the old one.
		If there are overlaps the overlap locations are returned.'''
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



	def Draw(self, dc):
		# make a path that contains a circle and some lines
		self.centre_x = self.size[0]/2 #centre of window in x
		self.centre_y = self.size[1]/2 #centro of window in y
		self.window_length = min(self.centre_x, self.centre_y)
		self.Radius = self.window_length/1.5

#		dc.SetDeviceOrigin(size_x/2, size_y/2)

		dc.SetBackground(wx.Brush("White"))
		dc.Clear() # make sure you clear the bitmap!
		gcdc = wx.GCDC(dc)

		#make a hidden dc to which features can be drawn in uinique colors and later used for hittests
		self.hidden_dc = wx.MemoryDC()
		self.hidden_dc.SelectObject(wx.EmptyBitmap(self.ClientSize[0], self.ClientSize[1]))
		self.hidden_dc.SetBackground(wx.Brush("White"))
		self.hidden_dc.Clear() # make sure you clear the bitmap!

		#draw DNA circles
		gcdc.SetPen(wx.Pen(colour='#444444', width=3))
		gcdc.SetBrush(wx.Brush("White"))
		gcdc.DrawCircle(x=self.centre_x, y=self.centre_y, radius=self.Radius) #outer DNA circle

		#draw plasmid name
		self.Draw_plasmid_name(gcdc)

		#draw features
		self.Draw_features(gcdc)

		#draw enzymes
		self.Draw_enzymes(gcdc)

		#draw selection
		self.Draw_selection(gcdc)
	
			
#		self.hidden_dc.SelectObject(wx.NullBitmap) # need to get rid of the MemoryDC before Update() is called.


	def Draw_plasmid_name(self, gcdc):
		'''Draw the plasmid name and basepairs in the middle'''
		name = genbank.gb.fileName.split('.')[0]
		basepairs = str(len(genbank.gb.GetDNA())) + ' bp'

		font = wx.Font(pointSize=self.window_length/16, family=wx.FONTFAMILY_SWISS, style=wx.FONTWEIGHT_BOLD, weight=wx.FONTWEIGHT_BOLD)
		gcdc.SetFont(font)
		gcdc.SetTextForeground(('#666666'))
		name_length = gcdc.GetTextExtent(name) #length of text in pixels
		gcdc.DrawText(name, self.centre_x-name_length[0]/2, self.centre_y-name_length[1])

		font = wx.Font(pointSize=self.window_length/20, family=wx.FONTFAMILY_SWISS, style=wx.FONTWEIGHT_NORMAL, weight=wx.FONTWEIGHT_NORMAL)
		gcdc.SetFont(font)
		gcdc.SetTextForeground(('#666666'))
		basepairs_length = gcdc.GetTextExtent(basepairs)
		gcdc.DrawText(basepairs, self.centre_x-basepairs_length[0]/2, self.centre_y+basepairs_length[1]/2)



	def Draw_selection(self, gcdc):
		'''Draws the current selection'''
		gcdc.SetBrush(wx.Brush(colour=wx.Colour(0,75,255,128))) #blue
		gcdc.SetPen(wx.Pen(colour='#444444', width=1))

		start, finish = copy.copy(genbank.dna_selection)
		start_angle = self.pos_to_angle(start-1)
		finish_angle = self.pos_to_angle(finish)

#		print('plasmid start finsh', start, finish)
#		print('plasmid angles', start_angle, finish_angle)

		if start == finish+1 or finish_angle-start_angle<0.3:
			xc=self.centre_x
			yc=self.centre_y
			x1 = xc + self.Radius * math.cos((finish_angle-90)*(math.pi/180))
			y1 = yc + self.Radius * math.sin((finish_angle-90)*(math.pi/180))
			gcdc.DrawLine(xc, yc, x1, y1)
		
		else:	
			xc=self.centre_x
			yc=self.centre_y
			x1 = xc + self.Radius * math.cos((finish_angle-90)*(math.pi/180)) #the latter needs to be first as the arc draws backwards
			y1 = yc + self.Radius * math.sin((finish_angle-90)*(math.pi/180))
			x2 = xc + self.Radius * math.cos((start_angle-90)*(math.pi/180))
			y2 = yc + self.Radius * math.sin((start_angle-90)*(math.pi/180))
			gcdc.DrawArc(x1, y1, x2, y2, xc, yc);


	def Draw_enzymes(self, gcdc):
		pass

	def Draw_features(self, gcdc):
		'''Function dedicated to drawing feature arrows. The highlighted variable is used to pass an integer in case one wishes to highlight features.'''
		self.window_length = min(self.centre_x, self.centre_y) #length of shortest part of window
		feature_thickness = self.window_length/12 #thickness of feature arrows and is used for a bunch of derived measurements
		outside_space = feature_thickness/4 #space between the outermost feature and the DNA circle
		arrowhead_length = 5 #length of arrowhead
		step = 0.25 #degree interval at which polygon point should be drawn
		spacer = feature_thickness/4 #for in-between features


		drawn_locations = [] #for keeping track of how many times a certain region has been painted on

		#features
		featurelist = genbank.gb.get_all_feature_positions()
		self.feature_catalog = {} #for matching features with the unique colors
		self.feature_catalog['(255, 255, 255, 255)'] = False #the background is white, have to add that key
		R = 0
		G = 0
		B = 0
		for i in range(0,len(featurelist)): 
			if R == 255 and G == 255 and B == 255:
				raise ValueError
			elif  R == 255 and G == 255:
				R = 0
				G = 0
				B += 1
			elif R == 255:
				R = 0
				G += 1
			else:
				R += 1
			unique_color = (R,G,B) #for drawing unique colors on the hidden dc
			self.feature_catalog[str(unique_color+(255,))] = i	

			featuretype, complement, start, finish, name, index = featurelist[i]
			drawn_locations, level = self.find_overlap(drawn_locations, (start, finish))

			#check so that stuff is not drawn too close to center
			#if it is too close, draw at top level
			if outside_space+feature_thickness+(feature_thickness+spacer)*level > self.Radius:
				level = 0

			featuretype = featuretype.replace('-', 'a') #for -10 and -35 region
			featuretype = featuretype.replace("5'", "a5") #for 5' features
			featuretype = featuretype.replace("3'", "a3") #for 5' features


			start_angle = self.pos_to_angle(start)
			finish_angle = self.pos_to_angle(finish)
			pointlist = [] #for storing drawing points for polygon
			xc=self.centre_x #centre of circle
			yc=self.centre_y #centre of circle

			#set color surrounding feature. Normally black, red if feature is highlighted
			if i is self.highlighted_feature: #if the current feature corresponds to that which should be highlighted
				gcdc.SetPen(wx.Pen(colour='#FF0000', width=2))
			else:
				gcdc.SetPen(wx.Pen(colour='#444444', width=1))

			#draw feature
			if complement == False:
				color = eval(featuretype)['fw'] #get the color of feature (as string)
				assert type(color) == str
				gcdc.SetBrush(wx.Brush(color))

				if arrowhead_length > int(finish_angle-start_angle): #if feature is too short to make arrow, make box
					#near side of box
					i = 0
					while i <= int(finish_angle-start_angle):
						x1 = xc + (self.Radius-outside_space-((feature_thickness+spacer)*level)) * math.cos((finish_angle-i-90)*(math.pi/180)) #the latter needs to be first as the arc draws backwards
						y1 = yc + (self.Radius-outside_space-((feature_thickness+spacer)*level)) * math.sin((finish_angle-i-90)*(math.pi/180))
						pointlist.append((x1,y1))
						i += step

					#far side of box
					i = int(finish_angle-start_angle)
					while i >= 0:
						x1 = xc + (self.Radius-outside_space-feature_thickness-((feature_thickness+spacer)*level)) * math.cos((finish_angle-i-90)*(math.pi/180)) #the latter needs to be first as the arc draws backwards
						y1 = yc + (self.Radius-outside_space-feature_thickness-((feature_thickness+spacer)*level)) * math.sin((finish_angle-i-90)*(math.pi/180))
						pointlist.append((x1,y1))
						i -= step

				else: #if not too short, draw arrow
					#near side of arrow
					i = 0
					while i <= int(finish_angle-start_angle):
						if i == 0: #the point of the arrow
							x1 = xc + (self.Radius-outside_space-feature_thickness/2-((feature_thickness+spacer)*level)) * math.cos((finish_angle-i-90)*(math.pi/180)) #the latter needs to be first as the arc draws backwards
							y1 = yc + (self.Radius-outside_space-feature_thickness/2-((feature_thickness+spacer)*level)) * math.sin((finish_angle-i-90)*(math.pi/180))
						elif i < arrowhead_length: #don't draw in-between the arrowhead point and the arrowhead body
							pass
#						elif i == arrowhead_length: #uncomment this if the arrowhead should have "wings on the side"
#							x1 = xc + (self.Radius-outside_space+feature_thickness/2-((feature_thickness+spacer)*level)) * math.cos((finish_angle-i-90)*(math.pi/180)) #the latter needs to be first as the arc draws backwards
#							y1 = yc + (self.Radius-outside_space+feature_thickness/2-((feature_thickness+spacer)*level)) * math.sin((finish_angle-i-90)*(math.pi/180))
						else:
							x1 = xc + (self.Radius-outside_space-((feature_thickness+spacer)*level)) * math.cos((finish_angle-i-90)*(math.pi/180)) #the latter needs to be first as the arc draws backwards
							y1 = yc + (self.Radius-outside_space-((feature_thickness+spacer)*level)) * math.sin((finish_angle-i-90)*(math.pi/180))
						pointlist.append((x1,y1))
						i += step

					#far side of arrow
					i = int(finish_angle-start_angle)
					while i >= 0:
						if i < arrowhead_length:
							pass
#						elif i == arrowhead_length:  #uncomment this if the arrowhead should have "wings on the side"
#							x1 = xc + (self.Radius-outside_space-feature_thickness-feature_thickness/2-((feature_thickness+spacer)*level)) * math.cos((finish_angle-i-90)*(math.pi/180)) #the latter needs to be first as the arc draws backwards
#							y1 = yc + (self.Radius-outside_space-feature_thickness-feature_thickness/2-((feature_thickness+spacer)*level)) * math.sin((finish_angle-i-90)*(math.pi/180))
						else:
							x1 = xc + (self.Radius-outside_space-feature_thickness-((feature_thickness+spacer)*level)) * math.cos((finish_angle-i-90)*(math.pi/180)) #the latter needs to be first as the arc draws backwards
							y1 = yc + (self.Radius-outside_space-feature_thickness-((feature_thickness+spacer)*level)) * math.sin((finish_angle-i-90)*(math.pi/180))
						pointlist.append((x1,y1))
						i -= step


			elif complement == True:
				color = eval(featuretype)['rv'] #get the color of feature (as string)
				assert type(color) == str
				gcdc.SetBrush(wx.Brush(color))

				if arrowhead_length > int(finish_angle-start_angle): #if feature is too short to make arrow, make box
					#near side of box
					i = 0
					while i <= int(finish_angle-start_angle):
						x1 = xc + (self.Radius-outside_space-((feature_thickness+spacer)*level)) * math.cos((finish_angle-i-90)*(math.pi/180)) #the latter needs to be first as the arc draws backwards
						y1 = yc + (self.Radius-outside_space-((feature_thickness+spacer)*level)) * math.sin((finish_angle-i-90)*(math.pi/180))
						pointlist.append((x1,y1))
						i += step

					#far side of box
					i = int(finish_angle-start_angle)
					while i >= 0:
						x1 = xc + (self.Radius-outside_space-feature_thickness-((feature_thickness+spacer)*level)) * math.cos((finish_angle-i-90)*(math.pi/180)) #the latter needs to be first as the arc draws backwards
						y1 = yc + (self.Radius-outside_space-feature_thickness-((feature_thickness+spacer)*level)) * math.sin((finish_angle-i-90)*(math.pi/180))
						pointlist.append((x1,y1))
						i -= step		
	
				else: #otherwise make arrow
					#near side of arrow
					i = 0
					while i <= int(finish_angle-start_angle):
						if i == int(finish_angle-start_angle):
							x1 = xc + (self.Radius-outside_space-feature_thickness/2-((feature_thickness+spacer)*level)) * math.cos((finish_angle-i-90)*(math.pi/180)) #the latter needs to be first as the arc draws backwards
							y1 = yc + (self.Radius-outside_space-feature_thickness/2-((feature_thickness+spacer)*level)) * math.sin((finish_angle-i-90)*(math.pi/180))
						elif i > int(finish_angle-start_angle)-arrowhead_length:
							pass
#						elif i == int(finish_angle-start_angle)-arrowhead_length: #uncomment this if the arrowhead should have "wings on the side"
#							x1 = xc + (self.Radius-outside_space+feature_thickness/2-((feature_thickness+spacer)*level)) * math.cos((finish_angle-i-90)*(math.pi/180)) #the latter needs to be first as the arc draws backwards
#							y1 = yc + (self.Radius-outside_space+feature_thickness/2-((feature_thickness+spacer)*level)) * math.sin((finish_angle-i-90)*(math.pi/180))
						else:
							x1 = xc + (self.Radius-outside_space-((feature_thickness+spacer)*level)) * math.cos((finish_angle-i-90)*(math.pi/180)) #the latter needs to be first as the arc draws backwards
							y1 = yc + (self.Radius-outside_space-((feature_thickness+spacer)*level)) * math.sin((finish_angle-i-90)*(math.pi/180))
						pointlist.append((x1,y1))
						i += step
				
					#far side of arrow
					i = int(finish_angle-start_angle)
					while i >= 0:
						if i > int(finish_angle-start_angle)-arrowhead_length:
							pass
#						elif i == int(finish_angle-start_angle)-arrowhead_length:  #uncomment this if the arrowhead should have "wings on the side"
#							x1 = xc + (self.Radius-outside_space-feature_thickness-feature_thickness/2-((feature_thickness+spacer)*level)) * math.cos((finish_angle-i-90)*(math.pi/180)) #the latter needs to be first as the arc draws backwards
#							y1 = yc + (self.Radius-outside_space-feature_thickness-feature_thickness/2-((feature_thickness+spacer)*level)) * math.sin((finish_angle-i-90)*(math.pi/180))
						else:
							x1 = xc + (self.Radius-outside_space-feature_thickness-((feature_thickness+spacer)*level)) * math.cos((finish_angle-i-90)*(math.pi/180)) #the latter needs to be first as the arc draws backwards
							y1 = yc + (self.Radius-outside_space-feature_thickness-((feature_thickness+spacer)*level)) * math.sin((finish_angle-i-90)*(math.pi/180))
						pointlist.append((x1,y1))
						i -= step

			#first draw the hidden features which are used for hittests on click
			self.hidden_dc.SetPen(wx.Pen(colour=unique_color, width=0))
			self.hidden_dc.SetBrush(wx.Brush(colour=unique_color))
			self.hidden_dc.DrawPolygon(pointlist)

			#now draw the real features
			gcdc.DrawPolygon(pointlist)



			###############
			# Draw labels #
			###############

			#draw label for feature
			#text parameters
			font_size = int(feature_thickness*0.6)
			font = wx.Font(pointSize=font_size, family=wx.FONTFAMILY_SWISS, style=wx.FONTWEIGHT_BOLD, weight=wx.FONTWEIGHT_BOLD)
			gcdc.SetFont(font)
			gcdc.SetTextForeground((0,0,0))
			label_line_length = self.window_length/8
			max_label_length = 120 #max length of label in pixels

			xc=self.centre_x #centre of circle
			yc=self.centre_y #centre of circle
			feature_name = name
			name_length = gcdc.GetTextExtent(feature_name) #length of text in pixels
			feature_radius = self.Radius-outside_space-((feature_thickness+spacer)*level) #the feature radius depends on where on the plasmid it is drawn
			feature_length = feature_radius*((finish_angle-start_angle)*math.pi)/float(180) #length of feature in pixels  len=radius*((theta*pi)/180)
		
			if name_length[0] < feature_length*0.9: #if feature is long enough to put text inside, do it
				mid_text = start_angle+(finish_angle-start_angle)/2	#middle of text should be here (angle)
				for i in range(0,len(feature_name)):
					if i == 0:
						plasmid_angle = mid_text-((1+gcdc.GetTextExtent(feature_name)[0]/2)*180)/(math.pi*feature_radius) #determine beginning angle
						text_angle = plasmid_angle + ((gcdc.GetTextExtent(feature_name[i])[0]/2)*180)/(math.pi*feature_radius) #determine beginning angle

					elif feature_name[i-1] == ' ': #GetTextExtent does not work on spaces
						plasmid_angle += ((font_size/3)*180)/(math.pi*feature_radius) #add length of space if one is present
#						text_angle = plasmid_angle + 
					else: 
						plasmid_angle += (1+gcdc.GetTextExtent(feature_name[i-1])[0]*180)/(math.pi*feature_radius) #add length of previous letter
						text_angle = plasmid_angle + gcdc.GetTextExtent(feature_name[i])[0]/2
					
					x1 = xc + feature_radius * math.cos((plasmid_angle-90)*(math.pi/180)) #the latter needs to be first as the arc draws backwards
					y1 = yc + feature_radius * math.sin((plasmid_angle-90)*(math.pi/180))
					gcdc.DrawRotatedText(feature_name[i], x1, y1, -(text_angle))
			
			else: #if feature is too short to put text inside,  put text on the outside
				while name_length[0] > max_label_length: #shorten text if it is too long 
					feature_name = feature_name[:-3]+'..'
					name_length = gcdc.GetTextExtent(feature_name) #length of text in pixels
		
				#draw the lines to the label and the label itself		
				angle = start_angle+(finish_angle-start_angle)/2
				x1 = xc + (self.Radius-outside_space-((feature_thickness+spacer)*level)) * math.cos((angle-90)*(math.pi/180)) #the latter needs to be first as the arc draws backwards
				y1 = yc + (self.Radius-outside_space-((feature_thickness+spacer)*level)) * math.sin((angle-90)*(math.pi/180))
				x2 = xc + (self.Radius+label_line_length) * math.cos((angle-90)*(math.pi/180)) #the latter needs to be first as the arc draws backwards
				y2 = yc + (self.Radius+label_line_length) * math.sin((angle-90)*(math.pi/180))
				gcdc.DrawLine(x1,y1,x2,y2)
				if angle <= 180:
					gcdc.DrawLine(x2,y2,x2+name_length[0]+3,y2)
					gcdc.DrawText(feature_name,x2+3,y2-gcdc.GetTextExtent(feature_name)[1])
				elif angle > 180:
					gcdc.DrawLine(x2,y2,x2-name_length[0],y2)
					gcdc.DrawText(feature_name,x2-name_length[0],y2-gcdc.GetTextExtent(feature_name)[1])

				#draw hidden box at text positon, used for hittests
				self.hidden_dc.SetPen(wx.Pen(colour=unique_color, width=0))
				self.hidden_dc.SetBrush(wx.Brush(colour=unique_color))
				self.hidden_dc.DrawRectangle(x2, y2, gcdc.GetTextExtent(feature_name)[0], -gcdc.GetTextExtent(feature_name)[1])

					

	def angle_to_pos(self, angle):
		'''Calculate DNA position from a start and end angle'''
		len_dna = float(len(genbank.gb.GetDNA()))
		dna_pos = int((angle/float(360))*len_dna)
		return dna_pos

	def pos_to_angle(self, pos):
		'''Calculate angles from DNA positions'''
		assert type(pos) == int, 'Error, position needs to be an integer'
		len_dna = float(len(genbank.gb.GetDNA()))
		if len_dna == 0:
			angle = 0
		else:
			angle = 360*(pos/float(len_dna))
		return angle

	def calc_angle(self, x, y):
		'''Calculates the angle corresponding to mouse clicks (as a fraction float)'''
		z = math.sqrt(math.pow((x-self.centre_x),2) + math.pow((y-self.centre_y),2)) #determine z
#		cosangle = (x-centre_x)/z
#		cos_radians = math.acos(cosangle)
#		degrees = cos_radians*(180/math.pi)
#		print('cos', degrees)

		sinangle = (y-self.centre_y)/z 	
		sin_radians = math.asin(sinangle)
		degrees = sin_radians*(180/math.pi)	
		#now I have the angle of the triangle. Based on where it is placed I have to calculate the 'real' angle

		#for these calculations it is important to remember that y increases as you go down...
		x_difference = x-self.centre_x
		y_difference = y-self.centre_y
		if x_difference >= 0 and y_difference >= 0: #triangle is in bottom right of circle
			angle = 90+degrees
#			print('bottom right')
		elif  x_difference <= 0 and y_difference <= 0: #triangle is in top left of circle
			angle = 180+90-degrees
#			print('top left')
		elif x_difference >= 0 and y_difference <= 0: #triangle is in in top right of circle
			angle = 90+degrees
#			print('top right')
		elif x_difference <= 0 and y_difference >= 0: #triangle is in bottom left of circle
			angle = 180+90-degrees
		else:
			ValueError
		return angle

	def HitTest(self):
		'''Tests whether the mouse is over any feature or label'''
		dc = wx.ClientDC(self) #get the client dc
		x, y = self.ScreenToClient(wx.GetMousePosition()) #get coordinate of mouse event
		pixel_color = self.hidden_dc.GetPixel(x,y) #use that coordinate to find pixel on the hidden d
		return self.feature_catalog[str(pixel_color)] #return the index
			

	def OnLeftDown(self, event):
		'''When left mouse button is pressed down, store angle at which this happened.'''
		self.centre_x = self.size[0]/2 #centre of window in x
		self.centre_y = self.size[1]/2 #centro of window in y
		x, y = self.ScreenToClient(wx.GetMousePosition())	
		angle = self.calc_angle(x, y)
		self.left_down_angle = angle #save the angle at which left button was clicked for later use


	def OnLeftUp(self, event):
		'''When left mouse button is lifted up, determine the DNA selection from angles generated at down an up events.'''
		self.centre_x = self.size[0]/2 #centre of window in x
		self.centre_y = self.size[1]/2 #centro of window in y
		x, y = self.ScreenToClient(wx.GetMousePosition())	

		up_angle = self.calc_angle(x, y)
		down_angle = self.left_down_angle

		if abs(down_angle-up_angle) <= 0.2: # want to do 'down == up' but I need some tolerance
			self.highlighted_feature = self.HitTest()
			if self.highlighted_feature is False: #if there is no feature, draw the line
				start = self.angle_to_pos(down_angle) 
				finish = self.angle_to_pos(up_angle)
			else:
				featuretype, complement, start, finish, name, index = genbank.gb.get_all_feature_positions()[self.highlighted_feature] #get info for the feature that was 'hit'
				start += 1 #need to adjust for some reason
		elif down_angle < up_angle:
			start = self.angle_to_pos(down_angle)
			finish = self.angle_to_pos(up_angle)
		elif down_angle > up_angle:
			start = self.angle_to_pos(up_angle)
			finish = self.angle_to_pos(down_angle)

		self.set_dna_selection((start, finish))
		self.update_ownUI()
		pub.Publisher.sendMessage('private_group_for_those_that_affect_DNA_selection_from_plasmid_view', '') #tell others that DNA selection changed


	def OnMotion(self, event):
		'''When mouse is moved with the left button down determine the DNA selection from angle generated at mouse down and mouse move event.'''
		if event.Dragging() and event.LeftIsDown():
			x, y = self.ScreenToClient(wx.GetMousePosition())	

			up_angle = self.calc_angle(x, y)
			down_angle = self.left_down_angle

			if down_angle <= up_angle:
				start = self.angle_to_pos(down_angle)
				finish = self.angle_to_pos(up_angle)
			elif down_angle > up_angle:
				start = self.angle_to_pos(up_angle)
				finish = self.angle_to_pos(down_angle)			

			self.set_dna_selection((start, finish))
			self.update_ownUI()
			pub.Publisher.sendMessage('private_group_for_those_that_affect_DNA_selection_from_plasmid_view', '') #tell others that DNA selection changed
		else:
			new_index = self.HitTest()
			if new_index is self.highlighted_feature: #if the index did not change
				pass
			else:
				self.highlighted_feature = new_index
				self.update_ownUI()

	def OnLeftDouble(self, event):
		'''When left button is duble clicked'''
		new_index = self.HitTest() #this does not get the "true" feature index. Some featues are split and this is an index that accounts for that.
		if new_index is not False:
			print('not false')
		
			featurelist = genbank.gb.get_all_feature_positions()
			featuretype, complement, start, finish, name, index = featurelist[new_index]
			genbank.feature_selection = copy.copy(index)

			dlg = featureedit_GUI.FeatureEditDialog(None, 'Edit Feature') # creation of a dialog with a title
			dlg.ShowModal()
			dlg.Center()

	def OnRightUp(self, event):
		print('plasmid right')





##### main loop
class MyApp(wx.App):
	def OnInit(self):
		frame = wx.Frame(None, -1, title="Plasmid View", size=(700,600), style = wx.NO_FULL_REPAINT_ON_RESIZE)
		panel =	PlasmidView(frame, -1)
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
