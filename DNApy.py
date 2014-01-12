#!/usr/bin/python


#DNApy is a DNA editor written purely in python. 
#The program is intended to be an intuitive, fully featured, 
#extendable, editor for molecular and synthetic biology.  
#Enjoy!
#
#Copyright (C) 2014  Martin Engqvist | 
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
#Get source code at: https://github.com/0b0bby0/DNApy
#
from copy import deepcopy
from string import *
from os import access,listdir
import sys, os
import wx
import wx.richtext as rt

from wxPython.lib.buttons import *
from wxPython.lib.colourselect import *

import pyperclip

import output
import dna
import genbank

import dnaeditor
import featureeditor


files={}   #dictionary with all configuration files

files['default_dir'] = os.path.abspath(os.path.dirname(sys.argv[0]))+"/"
files['default_dir']=replace(files['default_dir'], "\\", "/")
files['default_dir']=replace(files['default_dir'], "library.zip", "")

variables=files['default_dir']+"variables"   ##path to the file of the global variables
settings=files['default_dir']+"settings"   ##path to the file of the global settings
execfile(variables) #gets all the pre-assigned variables
execfile(settings) #gets all the pre-assigned settings




class MyFrame(wx.Frame):
	tab_list=[] #list of tabs 
	current_tab=0 #contains the current tab
	panel=[] #list of panels for the textbox

	def __init__(self, parent, id, title):
		wx.Frame.__init__(self, parent, id, title)
		ID=wx.NewId()
		self.DNApy = wx.Notebook(self, ID, style=0) ######create blank notebook
#		wx.EVT_NOTEBOOK_PAGE_CHANGED(self, ID, self.page_change)
#		wx.EVT_CLOSE(self, self.OnCloseWindow)

		self.generate_dnaview_tab("")
		self.generate_featureview_tab("")
		self.generate_genbankview_tab("")
		
		#create toolbars
		self.__generate_toolbar()
		
		
		#create Menu Bar
		self.create_menu()

		#create statusbar
		self.statusbar = self.CreateStatusBar(2)
		self.statusbar.SetStatusStyles(styles=[wx.SB_FLAT, wx.SB_FLAT])

		sizer = wx.BoxSizer(wx.VERTICAL)
		sizer.Add(self.frame_1_toolbar, 0, wx.EXPAND)
		sizer.Add(self.DNApy, -1, wx.EXPAND)
		self.SetSizer(sizer)	
		
		self.Centre()	


	def generate_dnaview_tab(self, evt):
		number=len(self.tab_list)

		self.panel.append(wx.Panel(self.DNApy, -1))

		self.dnaview = dnaeditor.MyPanel(self.panel[number])
	
		self.tab_list.append(self.dnaview)


		sizer_1=wx.BoxSizer(wx.HORIZONTAL)
		sizer_1.Add(self.tab_list[number], 1, wx.EXPAND, 0)
		self.DNApy.AddPage(self.panel[number], "DNA view")
		self.panel[number].SetSizer(sizer_1)


	def generate_featureview_tab(self, evt):
		number=len(self.tab_list)

		self.panel.append(wx.Panel(self.DNApy, -1))

		self.featureview = featureeditor.MyPanel(self.panel[number])
	
		self.tab_list.append(self.featureview)


		sizer_1=wx.BoxSizer(wx.HORIZONTAL)
		sizer_1.Add(self.tab_list[number], 1, wx.EXPAND, 0)
		self.DNApy.AddPage(self.panel[number], "Feature view")
		self.panel[number].SetSizer(sizer_1)


	def generate_vectorview_tab(self, evt):
		pass

	def generate_squencingview_tab(self, evt):
		pass

	def generate_genbankview_tab(self, evt):
		number=len(self.tab_list)

		self.panel.append(wx.Panel(self.DNApy, -1))

		self.genbankview = output.create(self.panel[number], style=wx.VSCROLL|wx.HSCROLL)
		self.genbankview.SetEditable(False)

		self.tab_list.append(self.genbankview)


		sizer_1=wx.BoxSizer(wx.HORIZONTAL)
		sizer_1.Add(self.tab_list[number], 1, wx.EXPAND, 0)
		self.DNApy.AddPage(self.panel[number], "GenBank view")
		self.panel[number].SetSizer(sizer_1)

################ file functions #################

	def new_document(self, evt):
		'''Create new gb file'''
		#self.gb=genbank.gbobject() #make new gb in panel	
		pass

	def open_file(self, evt):
		'''Function for opening file'''
		self.dir_to_open = '/home/martin/Python_files/DNApy/'
		dlg = wx.FileDialog( self, style=wx.OPEN|wx.FILE_MUST_EXIST,   defaultDir=self.dir_to_open ,wildcard='GenBank files (*.gb)|*|Any file (*)|*')
		dlg.ShowModal()
		fileName = dlg.GetFilename()
		all_path=dlg.GetPath()
		dire=dlg.GetDirectory()
		dlg.Destroy()
		if(fileName == None or fileName == "" ):
			return1
		
		name, extension = fileName.split('.')
		if extension == 'gb':
			#self.gb=genbank.gbobject() #make new gb in panel
			genbank.gb.readgb(all_path) #read the file

			#update DNA
			self.dnaview.gbviewer.SetValue(genbank.gb.get_dna()) #put dna in box
			self.dnaview.paint_features() #update features	

			#update features
			self.featureview.populate_feature_list()	

			#update genbank
			self.genbankview.SetValue(genbank.gb.make_gbstring())		

		else:
			print("error, not a gb file")		

#		self.Bind(wx.EVT_UPDATE_UI, self.update_statusbar)
#		wx.EVT_CLOSE(self, self.OnCloseWindow)
		
#		wx.EVT_KEY_DOWN(self, self.OnKeyPress)
	
	def save_all(self, evt):
		pass
	
	def save_file(self, evt):
		'''Function for saving file'''
#		try:


		genbank.gb.write_file(genbank.gb.get_filepath())
		
#			self.OnUpdateUI(1) #update statusbar
#		except:
#			self.save_as_file("")
		

	def save_as_file(self, evt):
		'''Function for saving file as'''
		filepath = genbank.gb.get_filepath()	
		for i in range(len(filepath)): #get directory for file
			if filepath[i] == '/':
				dire = filepath[0:i+1]
		#get save dialog
		dlg = wx.FileDialog( self, style=wx.SAVE | wx.OVERWRITE_PROMPT,defaultDir=dire,wildcard='TXT files (*)|*|Any file (*)|*')
		dlg.ShowModal()
		all_path=dlg.GetPath()
		fileName=dlg.GetFilename()
		dire=dlg.GetDirectory()
		dlg.Destroy()
		if (fileName == None or fileName == ""):
			return
		else:
			#try:
			genbank.gb.update_filepath(all_path)
			self.save_file("")
			#except:
			#	error_window(7, self)

	def quit(self, evt):
		'''Function for quiting program'''
		print("close")
#    		self.close_all("")
       		self.Close()
       		#self.Destroy()		


	def OnCloseWindow(self, evt):
		'''not currently used'''
#		self.close_all("")
#		foo=self.GetSize()  ###except for the window size of file 
#		if(self.IsMaximized()==0):
#			file=open(files['size'], "w")
#			file.write(str(foo[0])+"\n"+str(foo[1]))
#			file.close()
		self.Close()
		#self.Destroy()


##########################################################


	def __generate_toolbar(self):
		'''For generating toolbar with buttons'''
				
		self.frame_1_toolbar = wx.ToolBar(self, wx.ID_ANY, style=wx.TB_HORIZONTAL|wx.TB_FLAT|wx.TB_DOCKABLE)

   		#syntax for toolbar
   		#AddLabelTool(self, id, label, bitmap, bmpDisabled, kind, shortHelp, longHelp, clientData) 
   		

   		
   		#New Document
   		self.frame_1_toolbar.AddLabelTool(500, "New Document", wx.Bitmap(files['default_dir']+"/icon/new.png", wx.BITMAP_TYPE_ANY), wx.NullBitmap, wx.ITEM_NORMAL, 'New File', "New File") #last one is the one displayed in status bar
   		wx.EVT_TOOL(self, 500, self.new_document)
		#Open File
   		self.frame_1_toolbar.AddLabelTool(501, "Open File", wx.Bitmap(files['default_dir']+"/icon/open.png", wx.BITMAP_TYPE_ANY), wx.NullBitmap, wx.ITEM_NORMAL, 'Load File', 'Load File')
   		wx.EVT_TOOL(self, 501, self.open_file)
		#Save current file
   		self.frame_1_toolbar.AddLabelTool(502, "Save current file", wx.Bitmap(files['default_dir']+"/icon/save.png", wx.BITMAP_TYPE_ANY), wx.NullBitmap, wx.ITEM_NORMAL, 'Save File', 'Save File')
   		wx.EVT_TOOL(self, 502, self.save_file)
		#Save as
   		self.frame_1_toolbar.AddLabelTool(503, "Save as", wx.Bitmap(files['default_dir']+"/icon/saveas.png", wx.BITMAP_TYPE_ANY), wx.NullBitmap, wx.ITEM_NORMAL, 'Save File As', 'Save File As')
   		wx.EVT_TOOL(self, 503, self.save_as_file)
		#Cut
   		self.frame_1_toolbar.AddLabelTool(504, "Cut", wx.Bitmap(files['default_dir']+"/icon/cut.png", wx.BITMAP_TYPE_ANY), wx.NullBitmap, wx.ITEM_NORMAL, 'Cut', 'Cut')
   		wx.EVT_TOOL(self, 504, self.dnaview.Cut)
		#Copy
   		self.frame_1_toolbar.AddLabelTool(505, "Copy", wx.Bitmap(files['default_dir']+"/icon/copy.png", wx.BITMAP_TYPE_ANY), wx.NullBitmap, wx.ITEM_NORMAL, 'Copy', 'Copy')
   		wx.EVT_TOOL(self, 505, self.dnaview.Copy)
		#Paste
   		self.frame_1_toolbar.AddLabelTool(506, "Paste", wx.Bitmap(files['default_dir']+"/icon/paste.png", wx.BITMAP_TYPE_ANY), wx.NullBitmap, wx.ITEM_NORMAL, 'Paste', 'Paste')
   		wx.EVT_TOOL(self, 506, self.dnaview.Paste)
   		#Undo
   		self.frame_1_toolbar.AddLabelTool(513, "Undo", wx.Bitmap(files['default_dir']+"/icon/undo.png", wx.BITMAP_TYPE_ANY), wx.NullBitmap, wx.ITEM_NORMAL, 'Undo', 'Undo')
   		wx.EVT_TOOL(self, 513, self.dnaview.Undo)   
   		#Redo
   		self.frame_1_toolbar.AddLabelTool(514, "Redo", wx.Bitmap(files['default_dir']+"/icon/redo.png", wx.BITMAP_TYPE_ANY), wx.NullBitmap, wx.ITEM_NORMAL, 'Redo', 'Redo')
   		wx.EVT_TOOL(self, 514, self.dnaview.Redo) 
		#Search Upward
   		self.frame_1_toolbar.AddLabelTool(507, "Search Upward", wx.Bitmap(files['default_dir']+"/icon/up.png", wx.BITMAP_TYPE_ANY), wx.NullBitmap, wx.ITEM_NORMAL, 'Search Upward', 'Search Upward')
   		wx.EVT_TOOL(self, 507, self.dnaview.search_up)
		#Search window
		self.dnaview.search_tool()
		#Search Downward
   		self.frame_1_toolbar.AddLabelTool(509, "Search Downward", wx.Bitmap(files['default_dir']+"/icon/down.png", wx.BITMAP_TYPE_ANY), wx.NullBitmap, wx.ITEM_NORMAL, 'Search Downward', 'Search Downward')
   		wx.EVT_TOOL(self, 509, self.dnaview.search_down)
		#Print current window
#   		self.frame_1_toolbar.AddLabelTool(510, "Print current window", wx.Bitmap(files['default_dir']+"/icon/print.png", wx.BITMAP_TYPE_ANY), wx.NullBitmap, wx.ITEM_NORMAL, 'Print Current Window', 'Print Current Window')
 #  		wx.EVT_TOOL(self, 510, self.print_setup)
		
		self.frame_1_toolbar.Realize()
		
		
		##### Toolbar 2 #####
		
#		self.frame_2_toolbar = wx.ToolBar(self, wx.ID_ANY, style=wx.TB_HORIZONTAL|wx.TB_FLAT|wx.TB_DOCKABLE)
			
		#Select or Mutate
#		self.selormut = wx.ComboBox(self.frame_2_toolbar, id=wx.ID_ANY, size=(85, 28), choices=['Select', 'Mutate'], style=wx.CB_READONLY)
#		self.frame_2_toolbar.AddControl(self.selormut)		
#		self.selormut.SetSelection(0)
		
		#nucleotide or amino acid
#		self.nucleotideoraminoacid = wx.ComboBox(self.frame_2_toolbar, id=wx.ID_ANY, size=(120, 28), choices=['Nucleotide', 'Amino Acid', 'Feature'], style=wx.CB_READONLY)
#		self.frame_2_toolbar.AddControl(self.nucleotideoraminoacid)
#		self.nucleotideoraminoacid.SetSelection(0)
		
		#'position'
#		self.positionbox=wx.TextCtrl(self.frame_2_toolbar, id=wx.ID_ANY, size=(70, 25), value="position")
#		self.frame_2_toolbar.AddControl(self.positionbox)
#		self.positionbox.SetEditable(False)	
		
		
		#'pos #'
#		self.posnum=wx.TextCtrl(self.frame_2_toolbar, id=wx.ID_ANY, size=(90, 25), value="")
#		self.frame_2_toolbar.AddControl(self.posnum)
		
		
		#'in'
#		self.inbox=wx.TextCtrl(self.frame_2_toolbar, id=wx.ID_ANY, size=(25, 25), value="in")
#		self.frame_2_toolbar.AddControl(self.inbox)
#		self.inbox.SetEditable(False)
				
		#featurebox
#		self.featurebox = wx.ComboBox(self.frame_2_toolbar, id=wx.ID_ANY, size=wx.DefaultSize, choices=['File', 'feature1'], style=wx.CB_READONLY)
#		self.frame_2_toolbar.AddControl(self.featurebox)
		
		#'go'

#		self.frame_2_toolbar.Realize()






	def create_menu(self):     #method for creating menu
		self.menubar = wx.MenuBar()
		fileitem = wx.Menu()
			
		#new document
		fileitem.Append(1, "New\tCtrl+Q", "New Document")
		wx.EVT_MENU(self, 1, self.new_document)

		#open document
		fileitem.Append(2, "Open\tCtrl+O", "Open Document")
		wx.EVT_MENU(self, 2, self.open_file)
		fileitem.AppendSeparator()

		#save document
		fileitem.Append(3, "Save\tCtrl+S", "Save current document")
		wx.EVT_MENU(self, 3, self.save_file)

		#save document as
		fileitem.Append(4, "Save as", "Save a copy of current document")
		wx.EVT_MENU(self, 4, self.save_as_file)

		#save all
		fileitem.Append(5, "Save all", "Save all open tabs")
		wx.EVT_MENU(self, 5, self.save_all)
		fileitem.AppendSeparator()

		#close single
		fileitem.Append(5, "Close", "Close current document")
		wx.EVT_MENU(self, 5, self.dnaview.close_single)

		#close all
		fileitem.Append(6, "Close all", "Close all tabs")
		wx.EVT_MENU(self, 6, self.dnaview.close_all)
		fileitem.AppendSeparator()

		#quit
		fileitem.Append(7, "Exit", "Exit program")
		wx.EVT_MENU(self, 7, self.quit)

		self.menubar.Append(fileitem, "&File")

		######################### For 'Edit DNA' menu item #############################################
		self.edit = wx.Menu()
		#undo
		self.edit.Append(9, "Undo\tCtrl+Z", "Undo")
		wx.EVT_MENU(self, 9, self.dnaview.Undo)

		#redo
		self.edit.Append(10, "Redo\tCtrl+Y", "Redo")
		wx.EVT_MENU(self, 10, self.dnaview.Redo)
		self.edit.AppendSeparator() #________________________devider

		#cut
		self.edit.Append(11, "Cut\tCtrl+X", "Cut selected DNA")
		wx.EVT_MENU(self,11, self.dnaview.Cut)

		#copy
		self.edit.Append(12, "Copy\tCtrl+C", "Copy selected DNA")
		wx.EVT_MENU(self, 12, self.dnaview.Copy)

		#paste
		self.edit.Append(13, "Paste\tCtrl+V", "Paste DNA")
		wx.EVT_MENU(self, 13, self.dnaview.Paste)
		
		#cut reverse complement
		self.edit.Append(111, "Cut Rev-Comp\tCtrl+Shift+X", "Cut reverse-complement of selected DNA")
		wx.EVT_MENU(self,111, self.dnaview.Cut_RevComp)

		#copy reverse complement
		self.edit.Append(121, "Copy Rev-Comp\tCtrl+Shift+C", "Copy reverse-complement of selected DNA")
		wx.EVT_MENU(self, 121, self.dnaview.Copy_RevComp)		
		
		#paste reverse complement
		self.edit.Append(131, "Paste Rev-Comp\tCtrl+Shift+V", "Paste reverse-complement of DNA")
		wx.EVT_MENU(self, 131, self.dnaview.Paste_RevComp)

		#reverse-complement selection
		self.edit.Append(141, "Rev-Comp selection\tCtrl+Shift+R", "Reverse-complement seleected DNA")
		wx.EVT_MENU(self,141, self.dnaview.RevComp_sel)
		self.edit.AppendSeparator() #________________________devider
		
		#select all
		self.edit.Append(14, "Select all", "Select all text")
		wx.EVT_MENU(self, 14, self.dnaview.select_all)
		self.edit.AppendSeparator() #________________________devider

		#uppercase
		self.edit.Append(34, "Uppercase\tCtrl+W", "Convert selected text to uppercase")
		wx.EVT_MENU(self, 34, self.dnaview.Uppercase)

		#lowercase
		self.edit.Append(35, "Lowercase\tCtrl+E", "Convert selected text to lowercase")
		wx.EVT_MENU(self, 35, self.dnaview.Lowercase)
		self.edit.AppendSeparator() #________________________devider

		self.menubar.Append(self.edit, "Edit DNA")
		
	
	
		######## Features menu item ########
		self.features = wx.Menu()
		
		self.features.Append(40, "List Features", "List all features in file")
		wx.EVT_MENU(self,40, self.dnaview.list_features)
		
#		self.features.Append(41, "Edit Features", "Edit features in file")
#		wx.EVT_MENU(self,41, self.edit_features)		
		
		self.menubar.Append(self.features, "Features")
		


		######## Protein menu item #########
		self.protein = wx.Menu()
		
		#translate
		self.protein.Append(30, "Translate\tCtrl+T", "Translate DNA to protein")
		wx.EVT_MENU(self,30, self.dnaview.translate_selection)

		#translate reverse complement
		self.protein.Append(31, "Translate Rev-Comp\tCtrl+Shift+T", "Translate DNA to protein")
		wx.EVT_MENU(self,31, self.dnaview.translate_selection_reverse_complement)

		#translate feature
		self.protein.Append(32, "Translate feature", "Translate DNA feature to protein")
		wx.EVT_MENU(self,32, self.dnaview.translate_feature)

		self.menubar.Append(self.protein, "Protein")



		########## For 'Help' menu item #############
		self.help = wx.Menu()
		#about program
		self.help.Append(21, "About", "About this program")
		wx.EVT_MENU(self, 21, self.dnaview.info)

		#print IUPAC codes for dna and AA
		self.help.Append(22, "IUPAC codes", "IUPAC codes of DNA and amino acids")
		wx.EVT_MENU(self, 22, self.dnaview.IUPAC_codes)


		self.menubar.Append(self.help, "Help")

		self.SetMenuBar(self.menubar)










##### main loop
class MyApp(wx.App):
    def OnInit(self):
        frame = MyFrame(None, -1, "DNApy.py")
        frame.Show(True)
        self.SetTopWindow(frame)
        return True
        
if __name__ == '__main__': #if script is run by itself and not loaded	
	app = MyApp(0)
	app.MainLoop()
