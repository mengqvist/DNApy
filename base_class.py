#!/usr/bin/env python


#This file is part of DNApy. DNApy is a DNA editor written purely in python. 
#The program is intended to be an intuitive, fully featured, 
#extendible, editor for molecular and synthetic biology.  
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


import wx
import os, sys
import math
import colcol



class DNApyBaseClass(wx.Panel):
	'''
	This is a base class which should serve as the base for any panel used in DNApy.
	'''
	def __init__(self, parent, id):
		wx.Panel.__init__(self, parent)


#### own UI updates ####
	def update_ownUI(self):
		'''
		Method should be modified as to update graphical elements of own panel only. 
		Preferred use is through a self.update_ownUI() call.
		'''
		raise NotImplementedError

	
class DNApyBaseDrawingClass(DNApyBaseClass):
############# This first piece of code comes from a drawing example by ChrisBarker-NOAA  #####################
	'''
	A Buffered window class which should serve as the base for drawing panels used in DNApy.

	To use it, subclass it and define a Draw(DC) method that takes a DC
	to draw to. In that method, put the code needed to draw the picture
	you want. The window will automatically be double buffered, and the
	screen will be automatically updated when a Paint event is received.

	When the drawing needs to change, you app needs to call the
	UpdateDrawing() method. Since the drawing is stored in a bitmap, you
	can also save the drawing to file by calling the
	SaveToFile(self, file_name, file_type) method.
	
	Draw, OnPaint, OnSize and SaveToFile methods from ChrisBarker-NOAA at GitHub
	https://github.com/PythonCHB/wxPythonDemos/blob/master/ClientDCTest.py
	
	Other code by Martin K. M. Engqvist for the DNApy project.
	'''
	def __init__(self, parent, id):
		wx.Panel.__init__(self, parent)
		self.parent = parent
		
		self.Bind(wx.EVT_PAINT, self.OnPaint)
		self.Bind(wx.EVT_SIZE, self.OnSize)

		# OnSize called to make sure the buffer is initialized.
		# This might result in OnSize getting called twice on some
		# platforms at initialization, but little harm done.
		self.OnSize(None)
		self.paint_count = 0

		self.unique_color = None

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
		
######### ChrisBarker-NOAA code ends ##################################################################
#######################################################################################################
		
	def NextRGB(self):
		'''
		Get next RGB color.
		'''
		if self.unique_color == None:
			self.unique_color = (0,0,0)
		else:
			self.unique_color = colcol.NextRGB(self.unique_color)
		return self.unique_color

############### Setting methods for inter-converting angles to percent from angles ##############

	def AngleToFraction(self, angle):
		'''
		Calculate where on a circle (as a fraction of 1) a certain angle is.
		The input is an angle 0 to 360 degrees.
		The output is a fraction from 0 to 1.
		'''
		assert type(angle) is int or type(angle) is float, 'Error, the input needs to be an integer or a float'''
		assert 0<=angle<=360, 'Error, "%s" is not a valid angle. It needs to be 0 to 360.' % str(angle)
		fraction = angle/float(360)
		assert type(fraction) is float and 0<=fraction<=1, 'Error, fraction needs to be a float between 0 and 1.'
		return fraction


	def FractionToAngle(self, fraction):
		'''
		Calculate which angle (as degrees) a certain fraction of a circle corresponds to.
		'''
		assert type(fraction) is float and 0<=fraction<=1, 'Error, fraction needs to be a float between 0 and 1.'
		angle = fraction*360
		assert type(angle) is float, 'Error, the angle needs to be a float'''
		assert 0<=angle<=360, 'Error, "%s" is not a valid angle. It needs to be 0 to 360.' % str(angle)
		return angle	


	def PointsToAngle(self, centre_x, centre_y, x, y):
		'''
		Calculates the angle (in degrees as a float) corresponding to mouse clicks.
		'''
		assert type(centre_x) is int and type(centre_y) is int and type(x) is int and type(y) is int, 'Error, the input coordinates need to be integers.'

		#the click is one of the corners of a right-angled triangle
		z = math.sqrt(math.pow((x-centre_x),2) + math.pow((y-centre_y),2)) #determine z (the hypotenuse)
		sinangle = (y-centre_y)/z #calculate the Sin(angle)
		radians = math.asin(sinangle) #negate the Sin part and get radians
		degrees = radians*(180/math.pi)	#convert radians to degrees
		#now I have the angle of the triangle. Based on where it is placed I have to calculate the 'real' angle of a circle where 0 is at the top.

		#for these calculations it is important to remember that y increases as you go down...
		x_difference = x-centre_x
		y_difference = y-centre_y
		if x_difference >= 0 and y_difference >= 0: #triangle is in bottom right of circle
			degrees = 90+degrees
#			print('bottom right')
		elif  x_difference <= 0 and y_difference <= 0: #triangle is in top left of circle
			degrees = 180+90-degrees
#			print('top left')
		elif x_difference >= 0 and y_difference <= 0: #triangle is in in top right of circle
			degrees = 90+degrees
#			print('top right')
		elif x_difference <= 0 and y_difference >= 0: #triangle is in bottom left of circle
			degrees = 180+90-degrees
#			print('bottom left')
		else:
			ValueError
		assert type(degrees) is float and 0<=degrees<=360 
		return degrees


	def AngleToPoints(self, centre_x, centre_y, radius, angle):	
		'''
		Takes the centre of a circle, an angle (in degrees) and a radius and calculates the correspoinding XY coordinate on the circle
		'''
		assert type(centre_x) is int and type(centre_y) is int, 'Error, the input coordinates need to be integers.'
 		assert (type(angle) is int or type(angle) is float) and (type(radius) is int or type(radius) is float), 'Error, the input angle and radius need to be integers or floats.'
		assert 0<=angle<=360, 'Error, the input angle is %s, but needs to be between 0 and 360.' % angle
		x = int(centre_x + radius * math.cos((angle-90)*(math.pi/180)))
		y = int(centre_y + radius * math.sin((angle-90)*(math.pi/180)))
		return (x, y)

################### Done with angle methods ############################

	def make_arc(self, xc, yc, start_angle, finish_angle, radius, thickness, step=0.25, arrowhead_length=5, arrow=False):
		'''
		Compile a list of points that describe an arc of defined range and thickness.
		The arc can have the shape of an arrow.
		'''
		pointlist = [] #for storing drawing points for polygon
		if arrow==False:
			#far side of box
			i = 0
			while i <= finish_angle-start_angle:
				angle = finish_angle-i
				x, y = self.AngleToPoints(xc, yc, radius, angle)
				pointlist.append((x,y))
				i += step
			pointlist.append(self.AngleToPoints(xc, yc, radius, start_angle))

			#near side of box
			i = finish_angle-start_angle
			radius = radius-thickness
			while i >= 0:	
				angle = finish_angle-i
				x, y = self.AngleToPoints(xc, yc, radius, angle)						
				pointlist.append((x,y))
				i -= step
			pointlist.append(self.AngleToPoints(xc, yc, radius, finish_angle))

		elif arrow=='rv':
			#far side of arrow
			i = 0
			while i <= int(finish_angle-start_angle):
				if i >int(finish_angle-start_angle) - arrowhead_length:
					pass
				else:
					angle = finish_angle-i
					x, y = self.AngleToPoints(xc, yc, radius, angle)
					pointlist.append((x,y))
				i += step

			#point of arrow
			arrowpoint = radius - thickness/2
			angle = finish_angle-int(finish_angle-start_angle)
			x, y = self.AngleToPoints(xc, yc, arrowpoint, angle)
			pointlist.append((x,y))
		
			#near side of arrow
			i = int(finish_angle-start_angle)
			radius = radius-thickness
			while i >= 0:
				if i > int(finish_angle-start_angle) - arrowhead_length:
					pass
				else:
					angle = finish_angle-i
					x, y = self.AngleToPoints(xc, yc, radius, angle)
					pointlist.append((x,y))
				i -= step

		elif arrow=='fw':
			#far side of arrow
			i = 0
			while i <= int(finish_angle-start_angle):
				if i == 0: #the point of the arrow
					arrowpoint = radius - thickness/2
					angle = finish_angle-i
					x, y = self.AngleToPoints(xc, yc, arrowpoint, angle)
					pointlist.append((x,y))
				elif i < arrowhead_length:
					pass
				else:
					angle = finish_angle-i
					x, y = self.AngleToPoints(xc, yc, radius, angle)						
					pointlist.append((x,y))
				i += step

			#near side of arrow
			i = int(finish_angle-start_angle)
			radius = radius - thickness
			while i >= 0:
				if i < arrowhead_length:
					pass
				else:
					angle = finish_angle-i
					x, y = self.AngleToPoints(xc, yc, radius, angle)
					pointlist.append((x,y))
				i -= step
		return pointlist
	
