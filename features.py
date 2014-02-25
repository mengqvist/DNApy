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

import ast
import wx
import genbank
import sys, os
from string import *

import featureviewer
import featureeditor

## to do ##
#fix buttons
#make pretty
#add right-click menu
#add search function
#how to modify the location?
#make the xref and protein id and such copyable


files={}   #list with all configuration files

files['default_dir'] = os.path.abspath(os.path.dirname(sys.argv[0]))+"/"
files['default_dir']=replace(files['default_dir'], "\\", "/")
files['default_dir']=replace(files['default_dir'], "library.zip", "")

settings=files['default_dir']+"settings"   ##path to the file of the global settings
execfile(settings) #gets all the pre-assigned settings

	

		

class MyPanel(wx.Panel):
	def __init__(self, parent):
		wx.Panel.__init__(self, parent)
		
		self.current_highlight_color = '#FFFFFF'
		
		splitter1 = wx.SplitterWindow(self, -1, style=wx.SP_3D)
		

		##
		#first panel, for showing feature overview
		self.feature_list = featureviewer.FeatureView(splitter1, id=wx.ID_ANY)

		#second panel, for editing
		self.feature_edit = featureeditor.EditFeatureView(splitter1, id=wx.ID_ANY)

		splitter1.SplitHorizontally(self.feature_list, self.feature_edit)

		#global sizer		
		globsizer = wx.BoxSizer(wx.HORIZONTAL)
		globsizer.Add(splitter1, -1, wx.EXPAND)
		self.SetSizer(globsizer)
		self.Centre()
	
	def updateUI(self):
		"""Update feature list content"""
		self.feature_list.updateUI()



##### main loop
class MyApp(wx.App):
	def OnInit(self):
		frame = wx.Frame(None, -1, "featureedit.py")
		frame.featureeditor = MyPanel(frame)
		frame.Show(True)
		self.SetTopWindow(frame)
		return True
        
if __name__ == '__main__': #if script is run by itself and not loaded	
	app = MyApp(0)
	app.MainLoop()



	
