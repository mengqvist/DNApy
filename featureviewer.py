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

import featureeditor

files={}   #list with all configuration files

files['default_dir'] = os.path.abspath(os.path.dirname(sys.argv[0]))+"/"
files['default_dir']=replace(files['default_dir'], "\\", "/")
files['default_dir']=replace(files['default_dir'], "library.zip", "")

settings=files['default_dir']+"settings"   ##path to the file of the global settings
execfile(settings) #gets all the pre-assigned settings


class FeatureView(wx.Panel):
	"""Class for viewing features as a list"""
	def __init__(self, parent, id):
		wx.Panel.__init__(self, parent)
		self.listview = wx.ListCtrl(self, id=3001, style=wx.LC_REPORT|wx.LC_SINGLE_SEL)
		self.listview.InsertColumn(0, "Feature", format=wx.LIST_FORMAT_LEFT, width=200)
		self.listview.InsertColumn(1, "Type", format=wx.LIST_FORMAT_LEFT, width=100)
		self.listview.InsertColumn(2, "Location on DNA", format=wx.LIST_FORMAT_LEFT, width=200)
		self.listview.InsertColumn(3, "Strand", format=wx.LIST_FORMAT_LEFT, width=120)

#		font = wx.Font(pointSize=10, family=wx.FONTFAMILY_DEFAULT, style=wx.FONTSTYLE_NORMAL, weight=wx.FONTWEIGHT_NORMAL, underline=False, faceName='Monospace', encoding=wx.FONTENCODING_DEFAULT)
#		self.listview.SetFont(font)
		self.listview.Bind(wx.EVT_LIST_ITEM_FOCUSED, self.ListOnSelect)
		self.listview.Bind(wx.EVT_LIST_ITEM_ACTIVATED, self.ListOnActivate)

		
#		newfeature = wx.Button(self, 1, 'New Feature')
#		deletefeature = wx.Button(self, 2, 'Delete Feature')
#		moveup = wx.Button(self, 4, 'Move Up')
#		movedown = wx.Button(self, 5, 'Move Down')
#		copytranslation = wx.Button(self, 5, 'Copy Translation')

		#buttons
		imageFile = files['default_dir']+"/icon/new_small.png"
		image1 = wx.Image(imageFile, wx.BITMAP_TYPE_ANY).ConvertToBitmap()
		newfeature = wx.BitmapButton(self, id=1, bitmap=image1, size = (image1.GetWidth()+15, image1.GetHeight()+15), name = "share")

		imageFile = files['default_dir']+"/icon/remove_small.png"
		image1 = wx.Image(imageFile, wx.BITMAP_TYPE_ANY).ConvertToBitmap()
		deletefeature = wx.BitmapButton(self, id=2, bitmap=image1, size = (image1.GetWidth()+15, image1.GetHeight()+15), name = "share")

		imageFile = files['default_dir']+"/icon/move_up.png"
		image1 = wx.Image(imageFile, wx.BITMAP_TYPE_ANY).ConvertToBitmap()
		moveup = wx.BitmapButton(self, id=4, bitmap=image1, size = (image1.GetWidth()+15, image1.GetHeight()+15), name = "share")

		imageFile = files['default_dir']+"/icon/move_down.png"
		image1 = wx.Image(imageFile, wx.BITMAP_TYPE_ANY).ConvertToBitmap()
		movedown = wx.BitmapButton(self, id=5, bitmap=image1, size = (image1.GetWidth()+15, image1.GetHeight()+15), name = "share")


		#bind feature list buttons
		self.Bind(wx.EVT_BUTTON, self.OnNew, id=1)
		self.Bind(wx.EVT_BUTTON, self.OnDelete, id=2)
		self.Bind(wx.EVT_BUTTON, self.OnMoveFeatureUp, id=4)
		self.Bind(wx.EVT_BUTTON, self.OnMoveFeatureDown, id=5)
		
		sizer = wx.BoxSizer(wx.VERTICAL)
		sizer.Add(newfeature)
		sizer.Add(deletefeature)
		sizer.Add(moveup)
		sizer.Add(movedown)
		
		sizer2 = wx.BoxSizer(wx.HORIZONTAL)
		sizer2.Add(self.listview, 3, wx.EXPAND)
		sizer2.Add(sizer, 0, wx.EXPAND)

		self.SetSizer(sizer2)


	def get_selection(self):
		"""Get currently selected feature"""
		return self.listview.GetFocusedItem()

	def ListOnSelect(self, event):	
		'''Updates selection depending on which feature is chosen'''
		index = self.get_selection()
		genbank.gb.set_feature_selection(index)

	def ListOnActivate(self, event):
		'''Updates feature AND DNA selection based on which feature is chosen'''
		index = self.get_selection()
		genbank.gb.set_feature_selection(index)

		locations = genbank.gb.get_feature_location(index)
		start = genbank.gb.get_location(locations[0])[0]
		finish = genbank.gb.get_location(locations[-1])[1]
		genbank.gb.set_dna_selection((start, finish))  #how to I propagate this to the DNA view???

	def OnNew(self, event):
		'''Make new feature'''
		#make feature and update interface


		dlg = wx.Dialog(self, style=wx.YES_NO|wx.CANCEL)
		self.feature_edit = featureeditor.EditFeatureView(dlg, id=wx.ID_ANY)


		result = dlg.ShowModal()
		dlg.Destroy()
		if result == wx.ID_YES: #if yes, remove clutter
			genbank.gb.clean_clutter()
		
		


			genbank.gb.add_feature() #add arguments here!!!!!!!!!
			self.updateUI()

			number = self.listview.GetItemCount()-1
			self.update_feature_selection(number)
		

	def OnDelete(self, event):
		'''Delete selected feature'''
		#identify feature, remove it and update interface
		index = self.get_selection()
		feature = genbank.gb.get_feature(index)
		number = self.listview.GetFirstSelected()
		genbank.gb.remove_feature(feature)
		self.updateUI()
	
		#set highlight, focus and selection
		if number == self.listview.GetItemCount():
			number = number-1
		self.update_feature_selection(number)


	def OnMoveFeatureUp(self, event):
		'''Move feature up one step'''
		index = self.get_selection()
		feature = genbank.gb.get_feature(index)
		number = self.listview.GetFirstSelected()
		genbank.gb.move_feature(feature, 'u')	
		self.updateUI()
		if number != 0:
			number = number-1
		self.update_feature_selection(number)
		

	def OnMoveFeatureDown(self, event):
		'''Move feature up down step'''
		index = self.get_selection()
		feature = genbank.gb.get_feature(index)
		number = self.listview.GetFirstSelected()
		genbank.gb.move_feature(feature, 'd')
		self.updateUI()
	
		if number != self.listview.GetItemCount()-1:
			number = number+1
		self.update_feature_selection(number)


	def update_feature_selection(self, number):
		'''Updates which feature is selected'''
		genbank.gb.set_feature_selection(number)
		self.listview.SetItemState(number, wx.LIST_STATE_SELECTED,wx.LIST_STATE_SELECTED) #for the highlight
		self.listview.Select(number, True) #to select it
		self.listview.Focus(number) #to focus it


	def OnCopyFASTA(self, event):
		pass

	def OnCopyDNA(self, event):
		pass
		
	def OnCopyTranslation(self, event):
		pass

	def get_feature_color(self, feature):
		'''Takes single feature and finds its color, if any'''
		#piece of code to check if there are assigned colors to the features
		
		try: 
			Key = eval(feature['key'])
							
		except:
			Key = grey	
		
		complement = feature['complement']
		if complement == True: self.current_highlight_color = Key['rv']
		else: self.current_highlight_color = Key['fw']

	def updateUI(self):
		'''Refreshes table from features stored in the genbank object'''
		
		#need to figure out how to do this without changing the selection....

		self.listview.DeleteAllItems()
		n = 0 #for feautrecolor
		for entry in genbank.gb.gbfile['features']:
			col0 = entry['qualifiers'][0].split('=')[1]
	#		col0 = 'T7\terminator'
			col1 = entry['key']
			col2 = str(entry['location'])[1:-1]
			if entry['complement'] == True:
				col3 = 'complement'
			else:
				col3 = 'leading'
#			col4 = entry['qualifiers']
			
			self.listview.Append([col0, col1, col2, col3])	
		
			#coloring
			self.get_feature_color(entry)
			color = self.current_highlight_color
			item = n
			self.listview.SetItemBackgroundColour(item, color)	
			n += 1
		#self.autosize()

