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
#nothing at the moment


import ast
import wx
from wx.lib.agw import ultimatelistctrl as ULC
from wx.lib.pubsub import Publisher as pub

import genbank
import sys, os
import string
import copy

#for asserts and tests
import types
import unittest

import featureedit_GUI

class NewFeatureEdit(featureedit_GUI.FeatureEdit):
	'''Class that is bacially the same as the FeatureList, but with updates to how update_globalUI works'''
	def update_globalUI(self):
		pub.sendMessage('feature_list_updateUI', 'Update!!!') #send a message to update UI, a trick to make panals from different frames talk


class NewFeatureDialog(wx.Dialog):
	def __init__(self, parent, title):
		super(NewFeatureDialog, self).__init__(parent=parent,id=wx.ID_ANY, title=title, size=(700, 300)) 		


		self.feature_edit = NewFeatureEdit(self, id=wx.ID_ANY)	#get the feature edit panel
		self.NewFeatureButtonPanel = wx.Panel(self) #make new panel to hold the buttons

		self.OK = wx.Button(self.NewFeatureButtonPanel, 7, 'OK')
		self.OK.Bind(wx.EVT_BUTTON, self.OnOK, id=7)
#		self.Cancel = wx.Button(self.NewFeatureButtonPanel, 8, 'Cancel')
#		self.Cancel.Bind(wx.EVT_BUTTON, self.OnCancel, id=8)

		#organize buttons
		buttonsizer = wx.BoxSizer(wx.HORIZONTAL)
		buttonsizer.Add(self.OK)
#		buttonsizer.Add(self.Cancel)
		self.NewFeatureButtonPanel.SetSizer(buttonsizer)
		
		#organize feature edit panel and button panel
		sizer = wx.BoxSizer(wx.VERTICAL)
		sizer.Add(item=self.feature_edit, proportion=-1, flag=wx.EXPAND)
		sizer.Add(item=self.NewFeatureButtonPanel, proportion=0)
		self.SetSizer(sizer)		

#		wx.EVT_CLOSE(self, self.OnCancel) #when window is closed on x, cancel



	def OnOK(self, event):
		'''Accept new feature from the "new feature" popup"'''

#		key = self.feature_edit.
#		qualifiers = self.feature_edit.
#		location = self.feature_edit.
#		complement = self.feature_edit
#		join = self.feature_edit
#		order = self.feature_edit

#		self.update_globalUI()  #update DNA view
		self.Destroy()

#	def OnCancel(self, event):
#		'''Reject new feeature from the "new feature" popup'''
#		genbank.gb.remove_feature(genbank.gb.get_feature(index=-1))
#		genbank.feature_selection = 0
##		pub.sendMessage('feature_index', 0)
#		pub.sendMessage('feature_list_updateUI', 'update UI pls')
##		pub.sendMessage('dna_edit_updateUI', 'update UI pls')
#		self.Destroy()


class FeatureList(wx.Panel):
	"""Class for viewing features as a list"""
	def __init__(self, parent, id):
		wx.Panel.__init__(self, parent)
		self.feature_list = ULC.UltimateListCtrl(self, id=3001, agwStyle=ULC.ULC_REPORT|ULC.ULC_SINGLE_SEL)

		self.feature_list.InsertColumn(0, "Feature", format=wx.LIST_FORMAT_LEFT, width=200)
		self.feature_list.InsertColumn(1, "Type", format=wx.LIST_FORMAT_LEFT, width=100)
		self.feature_list.InsertColumn(2, "Location on DNA", format=wx.LIST_FORMAT_LEFT, width=200)
		self.feature_list.InsertColumn(3, "Strand", format=wx.LIST_FORMAT_LEFT, width=120)

#		font = wx.Font(pointSize=10, family=wx.FONTFAMILY_DEFAULT, style=wx.FONTSTYLE_NORMAL, weight=wx.FONTWEIGHT_NORMAL, underline=False, faceName='Monospace', encoding=wx.FONTENCODING_DEFAULT)
#		self.feature_list.SetFont(font)
		self.feature_list.Bind(wx.EVT_LIST_ITEM_SELECTED, self.ListOnSelect)
		self.feature_list.Bind(wx.EVT_LIST_ITEM_ACTIVATED, self.ListOnActivate)

		#subscribe to messages sent from feature edit module
#		pub.subscribe(self.listen_to_index_msg, 'feature_index')
		pub.subscribe(self.listen_to_updateUI, 'feature_list_updateUI')

		
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

		imageFile = files['default_dir']+"/icon/edit.png"
		image1 = wx.Image(imageFile, wx.BITMAP_TYPE_ANY).ConvertToBitmap()
		edit = wx.BitmapButton(self, id=6, bitmap=image1, size = (image1.GetWidth()+15, image1.GetHeight()+15), name = "edit")
		
		#bind feature list buttons
		self.Bind(wx.EVT_BUTTON, self.OnNew, id=1)
		self.Bind(wx.EVT_BUTTON, self.OnDelete, id=2)
		self.Bind(wx.EVT_BUTTON, self.OnMoveFeatureUp, id=4)
		self.Bind(wx.EVT_BUTTON, self.OnMoveFeatureDown, id=5)
		self.Bind(wx.EVT_BUTTON, self.OnEditFeature, id=6)


		#arrange buttons vertically		
		sizer = wx.BoxSizer(wx.VERTICAL)
		sizer.Add(item=newfeature)
		sizer.Add(item=deletefeature)
		sizer.Add(item=moveup)
		sizer.Add(item=movedown)
		sizer.Add(item=edit)

		#add feature list and buttons horisontally	
		sizer2 = wx.BoxSizer(wx.HORIZONTAL)
		sizer2.Add(item=self.feature_list, proportion=3, flag=wx.EXPAND)
		sizer2.Add(item=sizer, proportion=0, flag=wx.EXPAND)

		self.SetSizer(sizer2)
		self.update_ownUI()

	def listen_to_index_msg(self, index):
		'''For recieving index update request from the edit feature class'''
		print('type', type(index))
		print(index)
		newindex = copy.copy(index)
		print('type', type(newindex))
		print(newindex)
		self.update_feature_selection(newindex)
		print('yay', newindex)

	def listen_to_updateUI(self, msg):
		'''For recieving requests for UI updates.
			Sending messages is a roundabout way of doing this, 
			but I have not figured out a different way of getting the frames to talk.'''
		self.update_ownUI()

	def ListOnSelect(self, event):	
		'''Updates selection depending on which feature is chosen'''
		index = self.feature_list.GetFirstSelected()
		genbank.feature_selection = copy.copy(index)
		self.update_globalUI() #update feature editor

	def ListOnActivate(self, event):
		'''Updates feature AND DNA selection based on which feature is chosen'''
		index = self.feature_list.GetFirstSelected()
		genbank.feature_selection = copy.copy(index)

		locations = genbank.gb.get_feature_location(index)
		start = genbank.gb.get_location(locations[0])[0]
		finish = genbank.gb.get_location(locations[-1])[1]
		genbank.dna_selection = copy.copy((start-1, finish))
		self.GetTopLevelParent().dnaview.stc.SetSelection(start-1, finish) #update DNA selection

######
######
######
######
	def OnNew(self, event):
		'''Make new feature'''
		#make feature and update interface
		start, finish = self.get_dna_selection()
		
		genbank.gb.add_feature(key='misc_feature', qualifiers=['/note=New feature'], location=['%s..%s' % (start, finish)], complement=False, join=False, order=False)
		genbank.feature_selection = copy.copy(len(genbank.gb.get_all_features())-1)
		pub.sendMessage('feature_list_updateUI', 'update UI pls')

		dlg = NewFeatureDialog(None, 'New Feature') # creation of a dialog with a title
		dlg.ShowModal()
		dlg.Center()
		

	def OnDelete(self, event):
		'''Delete selected feature'''
		index = self.feature_list.GetFirstSelected()
		feature = genbank.gb.get_feature(index)
		genbank.gb.remove_feature(feature)
		self.update_ownUI() #update the feature view
		self.update_globalUI() 
	
		#set highlight, focus and selection
		if index == self.feature_list.GetItemCount():
			index = index-1
		if index < 0:
			self.update_feature_selection(None)	
		else:
			self.update_feature_selection(index)


	def OnMoveFeatureUp(self, event):
		'''Move feature up one step'''
		index = self.feature_list.GetFirstSelected()
		feature = genbank.gb.get_feature(index)
		genbank.gb.move_feature(feature, 'u')	
		self.update_ownUI()

		#set highlight, focus and selection
		if index != 0:
			index = index-1
		self.update_feature_selection(index)
		

	def OnMoveFeatureDown(self, event):
		'''Move feature up down step'''
		index = self.feature_list.GetFirstSelected()
		feature = genbank.gb.get_feature(index)
		genbank.gb.move_feature(feature, 'd')
		self.update_ownUI()

		#set highlight, focus and selection	
		if index != self.feature_list.GetItemCount()-1:
			index = index+1
		self.update_feature_selection(index)

	def OnEditFeature(self, event):
		'''Edit a feature that is already present'''
#		genbank.feature_selection = 		
#		pub.sendMessage('feature_list_updateUI', 'update UI pls')

		dlg = NewFeatureDialog(None, 'Edit Feature') # creation of a dialog with a title
		dlg.ShowModal()
		dlg.Center()		

	def focus_feature_selection(self):
		index = copy.copy(genbank.feature_selection)
#		print('type', type(index))
#		print('index', index)

		self.feature_list.SetItemState(item=index, state=ULC.ULC_STATE_SELECTED, stateMask=wx.LIST_STATE_SELECTED) #for the highlight
		self.feature_list.Select(index, True) #to select it
		self.feature_list.Focus(index) #to focus it

	def update_feature_selection(self, index):
		'''Updates which feature is selected'''
		genbank.feature_selection = copy.copy(index)
		if index == None:
			pass
		else:
			self.focus_feature_selection()



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


	def get_dna_selection(self):
		'''This method is intended for later modification by children.
			It is needed to get the dna selection for creating new features.'''
		return (1,1)

	def update_globalUI(self):
		'''This method is intended for later modification by children.
			It is needed to be flexible on which other components are updated in a full GUI application.'''		
		pass

	def update_ownUI(self):
		'''Refreshes only the UI of this panel by re-filling the table from features stored in the genbank object'''
		#need to figure out how to do this without changing the selection....
		self.feature_list.DeleteAllItems()
		item = 0 #for feautrecolor
		features = genbank.gb.get_all_features()
		for entry in features:
#			print(entry)
			if len(entry['qualifiers']) == 0:
				col0 = entry['key']	
			else:
				col0 = entry['qualifiers'][0].split('=')[1]
			col1 = entry['key']
			locationstring = ''
			for location in entry['location']:
				if locationstring != '':
					locationstring += ', '
				locationstring += str(location)
			col2 = locationstring
			if entry['complement'] == True:
				col3 = 'complement'
			else:
				col3 = 'leading'
#			col4 = entry['qualifiers']
		
			self.feature_list.Append([col0, col1, col2, col3])	
	
			#coloring
			self.get_feature_color(entry)
			hexcolor = self.current_highlight_color #get hex color
			r, g, b = self.hex_to_rgb(hexcolor) #convert to RGB
			color = wx.Colour(r, g, b) #make color object
			self.feature_list.SetItemBackgroundColour(item, color)	
			item += 1
		try:
			if genbank.feature_selection != None: #focus on the selected feature
				self.focus_feature_selection()
		except:
			pass
		

	def hex_to_rgb(self, value):
		value = value.lstrip('#')
		lv = len(value)
		return tuple(int(value[i:i+lv/3], 16) for i in range(0, lv, lv/3))

	def rgb_to_hex(self, rgb):
		return '#%02x%02x%02x' % rgb



######################################
######################################



##### main loop
class MyApp(wx.App):
	def OnInit(self):
		frame = wx.Frame(None, -1, title="Feature List", size=(700,500))
		panel =	FeatureList(frame, -1)
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

	genbank.dna_selection = (0, 0)	 #variable for storing current DNA selection
	genbank.feature_selection = False #variable for storing current feature selection

	import sys
	assert len(sys.argv) == 2, 'Error, this script requires a path to a genbank file as an argument.'
	print('Opening %s' % str(sys.argv[1]))

	genbank.gb = genbank.gbobject(str(sys.argv[1])) #make a genbank object and read file


	app = MyApp(0)
	app.MainLoop()
