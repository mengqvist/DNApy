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
from wx.lib.pubsub import pub

import math
import genbank

import os, sys
import string
from base_class import DNApyBaseClass


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

#		self.listening_group3 = 'dna_selection_request'		
#		pub.Publisher.subscribe(self.set_dna_selection, self.listening_group3)	


	def update_globalUI(self):
		'''Method should be modified as to update other panels in response to changes in own panel.
		Preferred use is through sending a message using the pub module.
		Example use is: pub.Publisher.sendMessage('feature_list_updateUI', '').
		The first string is the "listening group" and deterimines which listeners get the message. 
		The second string is the message and is unimportant for this implementation.
		The listening group assigned here (to identify recipients) must be different from the listening group assigned in __init__ (to subscribe to messages).'''
		#pub.Publisher.sendMessage('from_dna_edit', '')
		pass

	
	def update_ownUI(self):
		'''For changing background color of text ranges'''
		pass

	def set_dna_selection(self):
		pass


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

        wx.EVT_PAINT(self, self.OnPaint)
        wx.EVT_SIZE(self, self.OnSize)

        # OnSize called to make sure the buffer is initialized.
        # This might result in OnSize getting called twice on some
        # platforms at initialization, but little harm done.
        self.OnSize(None)
        self.paint_count = 0


    def Draw(self, dc):
        ## just here as a place holder.
        ## This method should be over-ridden when subclassed
        pass

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
        self.UpdateDrawing()

    def SaveToFile(self, FileName, FileType=wx.BITMAP_TYPE_PNG):
        ## This will save the contents of the buffer
        ## to the specified file. See the wxWindows docs for 
        ## wx.Bitmap::SaveFile for the details
        self._Buffer.SaveFile(FileName, FileType)

    def UpdateDrawing(self):
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
        del dc # need to get rid of the MemoryDC before Update() is called.
        self.Refresh()
        self.Update()




class PlasmidView(BufferedWindow):
	def __init__(self, *args, **kwargs):	
		self.x = 0
		self.y = 0

		self.selection_down = None
		self.selection_up = None

		self.centre_x = 0
		self.centre_y = 0
		BufferedWindow.__init__(self, *args, **kwargs)
	

		wx.EVT_LEFT_DOWN(self, self.OnLeftDown)
		wx.EVT_LEFT_UP(self, self.OnLeftUp)
		wx.EVT_RIGHT_UP(self, self.OnRightUp)
		wx.EVT_MOTION(self, self.OnMotion)



	def Draw(self, dc):
		# make a path that contains a circle and some lines
		self.centre_x = self.size[0]/2 #centre of window in x
		self.centre_y = self.size[1]/2 #centro of window in y

#		dc.SetDeviceOrigin(size_x/2, size_y/2)

		dc.SetBackground(wx.Brush("White"))
		dc.Clear() # make sure you clear the bitmap!

		gcdc = wx.GCDC(dc)

		if self.centre_x > self.centre_y:
			self.Radius = self.centre_y/1.5
		else:
			self.Radius = self.centre_x/1.5

		#DNA circles
		gcdc.SetPen(wx.Pen(colour='#444444', width=3))
		gcdc.SetBrush(wx.Brush("White"))
		gcdc.DrawCircle(x=self.centre_x, y=self.centre_y, radius=self.Radius) #outer DNA circle
		gcdc.DrawCircle(x=self.centre_x, y=self.centre_y, radius=self.Radius-5) #inner DNA circle


		#features
		featurelist = genbank.gb.get_all_feature_positions()
		for entry in featurelist: 
			featuretype, complement, start, finish = entry
			featuretype = featuretype.replace('-', 'a') #for -10 and -35 region
			featuretype = featuretype.replace("5'", "a5") #for 5' features
			featuretype = featuretype.replace("3'", "a3") #for 5' features

			start_angle, finish_angle = self.pos_to_angle(start, finish)


			if complement == False:
				color = eval(featuretype)['fw'] #get the color of feature (as string)
				assert type(color) == str

				gcdc.SetPen(wx.Pen(colour='#444444', width=1))
				gcdc.SetBrush(wx.Brush(color))

				xc=self.centre_x
				yc=self.centre_y
			
				pointlist = []

				step = 0.25			
				i = 0
				while i <= int(finish_angle-start_angle):
					x1 = xc + (self.Radius+25) * math.cos((finish_angle-i-90)*(math.pi/180)) #the latter needs to be first as the arc draws backwards
					y1 = yc + (self.Radius+25) * math.sin((finish_angle-i-90)*(math.pi/180))
					pointlist.append((x1,y1))
					i += step

				i = int(finish_angle-start_angle)
				while i >= 0:
					x1 = xc + (self.Radius+0) * math.cos((finish_angle-i-90)*(math.pi/180)) #the latter needs to be first as the arc draws backwards
					y1 = yc + (self.Radius+0) * math.sin((finish_angle-i-90)*(math.pi/180))
					pointlist.append((x1,y1))
					i -= step

			elif complement == True:
				color = eval(featuretype)['rv'] #get the color of feature (as string)
				assert type(color) == str

				gcdc.SetPen(wx.Pen(colour='#444444', width=1))
				gcdc.SetBrush(wx.Brush(color))

				xc=self.centre_x
				yc=self.centre_y
			
				pointlist = []

				step = 0.5			
				i = 0
				while i <= int(finish_angle-start_angle):
					x1 = xc + (self.Radius-6) * math.cos((finish_angle-i-90)*(math.pi/180)) #the latter needs to be first as the arc draws backwards
					y1 = yc + (self.Radius-6) * math.sin((finish_angle-i-90)*(math.pi/180))
					pointlist.append((x1,y1))
					i += step

				i = int(finish_angle-start_angle)
				while i >= 0:
					x1 = xc + (self.Radius-31) * math.cos((finish_angle-i-90)*(math.pi/180)) #the latter needs to be first as the arc draws backwards
					y1 = yc + (self.Radius-31) * math.sin((finish_angle-i-90)*(math.pi/180))
					pointlist.append((x1,y1))
					i -= step

			gcdc.DrawPolygon(pointlist)



		#selection
		gcdc.SetPen(wx.Pen(colour='#444444', width=1))
		gcdc.SetBrush(wx.Brush(colour=wx.Colour(0,75,255,128))) #blue
		if self.selection_down != None and self.selection_up != None:
			xc=self.centre_x
			yc=self.centre_y
			x1 = xc + self.Radius * math.cos((self.selection_up-90)*(math.pi/180)) #the latter needs to be first as the arc draws backwards
			y1 = yc + self.Radius * math.sin((self.selection_up-90)*(math.pi/180))
			x2 = xc + self.Radius * math.cos((self.selection_down-90)*(math.pi/180))
			y2 = yc + self.Radius * math.sin((self.selection_down-90)*(math.pi/180))
			gcdc.DrawArc(x1, y1, x2, y2, xc, yc)
			
		#enzymes
		#names

	def pos_to_angle(self, start, finish):
		'''Calculate angles from DNA positions'''
		len_dna = float(len(genbank.gb.GetDNA()))
		start_angle = 360*(start/len_dna)
		finish_angle = 360*(finish/len_dna)

		return start_angle, finish_angle

	def calc_location(self):
		'''Calculates the DNA location corresponding to mouse clicks (as a fraction float)'''
		z = math.sqrt(math.pow((self.x-self.centre_x),2) + math.pow((self.y-self.centre_y),2)) #determine z
#		cosangle = (self.x-centre_x)/z
#		cos_radians = math.acos(cosangle)
#		degrees = cos_radians*(180/math.pi)
#		print('cos', degrees)

		sinangle = (self.y-self.centre_y)/z 	
		sin_radians = math.asin(sinangle)
		degrees = sin_radians*(180/math.pi)	
		#now I have the angle of the triangle. Based on where it is placed I have to calculate the 'real' angle

		#for these calculations it is important to remember that y increases as you go down...
		x_difference = self.x-self.centre_x
		y_difference = self.y-self.centre_y
		if x_difference > 0 and y_difference > 0: #triangle is in bottom right of circle
			angle = 90+degrees
#			print('bottom right')
		elif  x_difference < 0 and y_difference < 0: #triangle is in top left of circle
			angle = 180+90-degrees
#			print('top left')
		elif x_difference > 0 and y_difference < 0: #triangle is in in top right of circle
			angle = 90+degrees
#			print('top right')
		elif x_difference < 0 and y_difference > 0: #triangle is in bottom left of circle
			angle = 180+90-degrees
#			print('bottom left')
#		print('angle', angle)

		#now calculate what percentage of the circle that is
		#length = (diameter*pi*angle)/360
#		length = (self.Radius*2*math.pi*angle)/360
#		percent_length = length/(self.Radius*2*math.pi)
#		print('percent', percent_length)
		return angle


	def OnLeftDown(self, event):
		self.centre_x = self.size[0]/2 #centre of window in x
		self.centre_y = self.size[1]/2 #centro of window in y
		x, y = self.ScreenToClient(wx.GetMousePosition())	
		self.x = x
		self.y = y

		angle = self.calc_location()
		self.selection_down = angle
		self.selection_up = angle +1
		self.UpdateDrawing()

	def OnLeftUp(self, event):
		self.centre_x = self.size[0]/2 #centre of window in x
		self.centre_y = self.size[1]/2 #centro of window in y
		x, y = self.ScreenToClient(wx.GetMousePosition())	
		self.x = x
		self.y = y

		angle = self.calc_location()
		if (self.selection_down-1)<=angle<=(self.selection_down+1):
			self.selection_up = self.selection_down +1
		else:
			self.selection_up = angle
		self.UpdateDrawing()




	def OnMotion(self, event):
	#	time.sleep(0.01)
		if event.Dragging() and event.LeftIsDown():
			x, y = self.ScreenToClient(wx.GetMousePosition())	
			self.x = x
			self.y = y

			angle = self.calc_location()
			self.selection_up = angle

			self.UpdateDrawing()

	def OnRightUp(self, event):
		self.Close()





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
