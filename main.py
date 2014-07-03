#!/usr/bin/python


#DNApy is a DNA editor written purely in python. 
#The program is intended to be an intuitive, fully featured, 
#extendable, editor for molecular and synthetic biology.  
#Enjoy!
#
#copyright (C) 2014  Martin Engqvist | 
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
import string
from os import access,listdir
import sys, os
import wx
import wx.richtext as rt


import pyperclip
import subprocess



import dna
import genbank
import output

#GUI components
import dnaeditor_GUI
import featureedit_GUI
import featurelist_GUI
import plasmid_GUI
import genbank_GUI
import mixed_base_codons_GUI
from wx.lib.pubsub import pub

#TODO
#add pretty dna view
#make rightklick menus
#fix mutate
#fix find previus and find next

files={}   #dictionary with all configuration files

files['default_dir'] = os.path.abspath(os.path.dirname(sys.argv[0]))+"/"
files['default_dir']=string.replace(files['default_dir'], "\\", "/")
files['default_dir']=string.replace(files['default_dir'], "library.zip", "")

variables=files['default_dir']+"variables"   ##path to the file of the global variables
settings=files['default_dir']+"settings"   ##path to the file of the global settings
execfile(variables) #gets all the pre-assigned variables
execfile(settings) #gets all the pre-assigned settings


class MyFrame(wx.Frame):
	'''Main frame of DNApy'''
	tab_list=[] #list of tabs 
	current_tab=0 #contains the current tab
	panel=[] #list of panels for the textbox
	
	
	def __init__(self, parent, id, title):
		wx.Frame.__init__(self, parent, id, title, size=windowsize) #size of program
		ID=wx.NewId()
#		self.DNApy = wx.Notebook(self, ID, style=0) ######create blank notebook
#		wx.EVT_NOTEBOOK_PAGE_CHANGED(self, ID, self.page_change)

		self.Bind(wx.EVT_CLOSE, self.OnCloseWindow)

		self.listening_group = 'from_feature_list'	
		pub.Publisher.subscribe(self.listen_to_updateUI, self.listening_group)

		self.listening_group1 = 'from_dna_edit' #recieve updates from DNA editor
		pub.Publisher.subscribe(self.listen_to_updateUI, self.listening_group1)

		self.listening_group2 = 'from_feature_edit'		
		pub.Publisher.subscribe(self.listen_to_updateUI, self.listening_group2)		

		self.listening_group4 = 'from_plasmid_view'
		pub.Publisher.subscribe(self.listen_to_updateUI, self.listening_group4)	



		#create toolbars
		self.__generate_toolbar()
		self.__generate_searchandmutate_toolbar()
		
		
		#create Menu Bar
		self.create_menu()

		#create statusbar
		self.statusbar = self.CreateStatusBar(2)
		self.statusbar.SetStatusStyles(styles=[wx.SB_FLAT, wx.SB_FLAT])


		self.fileopen = False #used to see if a file is open
		self.new_file(None) #create new genbank file

		genbank.dna_selection = (1, 1)	 #variable for storing current DNA selection
		genbank.feature_selection = False #variable for storing current feature selection
		genbank.search_hits = []
		genbank.gb.fileName = ''

		#create splitter and panels
#		self.splitter2 = wx.SplitterWindow(self, 0, style=wx.SP_3D)

		self.splitter1 = wx.SplitterWindow(self, 0, style=wx.SP_3D)	
		self.feature_list = featurelist_GUI.FeatureList(self.splitter1, id=wx.ID_ANY)
		self.dnaview = dnaeditor_GUI.DNAedit(self.splitter1, id=wx.ID_ANY)
		self.splitter1.SplitHorizontally(self.feature_list, self.dnaview, sashPosition=-(windowsize[1]-240))


		#plasmid view
		self.plasmid_frame = wx.Frame(self, -1, title="Plasmid view", size=(500,500))
		self.plasmid_view = plasmid_GUI.PlasmidView(self.plasmid_frame, -1)		

#		self.splitter2.SplitVertically(self.splitter1, self.plasmid_view, sashPosition=-(windowsize[0]/2.2))
		
		self.do_layout()
		self.Centre()
	
		

	def do_layout(self):
		'''Pack toolbar and the tabs in their sizers'''
		#add to sizer
		sizer_1 = wx.BoxSizer(wx.VERTICAL)
		sizer_1.Add(item=self.frame_1_toolbar, proportion=0, flag=wx.EXPAND)
		sizer_1.Add(item=self.splitter1, proportion=-1, flag=wx.EXPAND)
		#if second toolbar is present, add that too.
		try:
			sizer_1.Add(self.frame_2_toolbar, 0, wx.EXPAND)
		except:
			pass
		self.SetSizer(sizer_1)
		



	


#	def generate_genbankview_tab(self, evt):
#		number=len(self.tab_list)
#
#		self.panel.append(wx.Panel(self.DNApy, -1))
#		self.genbankview = genbank_GUI.MyPanel(self.panel[number], style=wx.VSCROLL|wx.HSCROLL)
#		self.genbankview.rtc.SetEditable(False)
#
#		self.tab_list.append(self.genbankview)
#
#
#		sizer_1=wx.BoxSizer(wx.HORIZONTAL)
#		sizer_1.Add(self.tab_list[number], 1, wx.EXPAND, 0)
#		self.DNApy.AddPage(self.panel[number], "GenBank")
#		self.panel[number].SetSizer(sizer_1)

################ file functions #################

	def new_file(self, evt):
		'''Create new gb file'''
		if self.fileopen == False: #if no file is open, make blank gb file
			genbank.gb = genbank.gbobject() #make new gb in panel	


			self.SetTitle('NewFile - DNApy')
#			self.page_change("")

			self.frame_1_toolbar.EnableTool(502, 1)
			self.frame_1_toolbar.EnableTool(503, 1)
			self.frame_1_toolbar.EnableTool(504, 1)
			self.frame_1_toolbar.EnableTool(505, 1)
			self.frame_1_toolbar.EnableTool(506, 1)
			self.frame_1_toolbar.EnableTool(511, 1)
			self.frame_1_toolbar.EnableTool(512, 1)
			self.Bind(wx.EVT_UPDATE_UI, self.update_statusbar)
			self.fileopen = True
			

		elif self.fileopen == True: #if a file IS open, make a new window
			subprocess.Popen("python ~/Python_files/DNApy/main.py 1", shell=True)

	def open_file(self, evt):
		'''Function for opening file'''
		self.dir_to_open = default_filepath
		dlg = wx.FileDialog( self, style=wx.OPEN|wx.FILE_MUST_EXIST,   defaultDir=self.dir_to_open ,wildcard='GenBank files (*.gb)|*|Any file (*)|*')
		dlg.ShowModal()
		genbank.gb.fileName = dlg.GetFilename()
		all_path=dlg.GetPath()
		dire=dlg.GetDirectory()
		dlg.Destroy()
		if(genbank.gb.fileName == None or genbank.gb.fileName == "" ):
			return1
		
		name, extension = genbank.gb.fileName.split('.')
		if extension.lower() == 'gb':
			genbank.gb = genbank.gbobject(all_path) #make a genbank object and read file
			self.dnaview.update_ownUI()
			self.dnaview.update_globalUI()

			self.SetTitle(genbank.gb.fileName+' - DNApy')
			if genbank.gb.clutter == True: #if tags from ApE or Vector NTI is found in file
				dlg = wx.MessageDialog(self, style=wx.YES_NO|wx.CANCEL, message='This file contains tags from the Vector NTI or ApE programs. Keeping these tags may break compatibility with other software. Removing them will clean up the file, but may result in the loss of some personalized styling options when this file is viewed in Vector NTI or ApE. Do you wish to REMOVE these tags?')
				result = dlg.ShowModal()
				dlg.Destroy()
				if result == wx.ID_YES: #if yes, remove clutter
					genbank.gb.clean_clutter()
					self.dnaview.update_ownUI()
					self.dnaview.update_globalUI()				

			self.frame_1_toolbar.EnableTool(502, True)
			self.frame_1_toolbar.EnableTool(503, True)
			self.frame_1_toolbar.EnableTool(504, True)
			self.frame_1_toolbar.EnableTool(505, True)
			self.frame_1_toolbar.EnableTool(506, True)
			self.frame_1_toolbar.EnableTool(511, True)
			self.frame_1_toolbar.EnableTool(512, True)
			self.fileopen = True
			
		else:
			print("error, not a gb file")		

		self.Bind(wx.EVT_UPDATE_UI, self.update_statusbar)

	
	def save_file(self, evt):
		'''Function for saving file'''
		genbank.gb.Save()


	def save_as_file(self, evt):
		'''Function for saving file as'''
		filepath = genbank.gb.GetFilepath()
		print(filepath)	
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
			genbank.gb.fileName = fileName
			if genbank.gb.fileName[-3:].lower() != '.gb': #make sure it has gb file ending
				all_path += '.gb'
				genbank.gb.fileName += '.gb'
			#try:
			genbank.gb.SetFilepath(all_path)
			self.save_file("")
			self.SetTitle(genbank.gb.fileName+' - DNApy')
			#except:
			#	error_window(7, self)


	def quit(self, evt):
		'''Function for quiting program'''
		self.Close()
		self.Destroy()		


	def OnCloseWindow(self, evt):
		'''Function for when window is closed'''

#		foo=self.GetSize()  ###except for the window size of file 
#		if(self.IsMaximized()==0):
#			file=open(files['size'], "w")
#			file.write(str(foo[0])+"\n"+str(foo[1]))
#			file.close()
		self.Destroy()


##########################################################



	def update_statusbar(self, evt):
		'''Updates statusbar'''
		#this stuff is for the statusbar
#		if len(self.tab_list) == 0:
#			string = 'File unmodified'
#		elif self.tab_list[self.current_tab].modify==0:
#			string = 'File unmodified'
#		elif self.tab_list[self.current_tab].modify==1:
#			string = 'File not yet saved'
#		self.current_tab=self.DNApy.GetSelection()
#		if self.current_tab == 0: #if dna editor is active
			
		#mposition, Feature = self.dnaview.mouse_position("") #get mouse position
		mposition = 'None'
		Feature = "None"
	
		try:
			Position = str(mposition+1)
		except:
			Position = ""
	
		try:
			Feature = str(Feature)
		except:
			Feature = ""
	
		try:		
			SelectionFrom, SelectionTo = (str(self.dnaview.stc.GetSelection()[0]+1), str(self.dnaview.stc.GetSelection()[1]))
			if SelectionFrom == '-1' and SelectionTo == '-2': #no selection if true
				SelectionFrom, SelectionTo = ("0", "0")
		except:
			SelectionFrom, SelectionTo = ("0", "0")
		try:	
			Length = str(self.dnaview.stc.GetSelection()[1] - self.dnaview.stc.GetSelection()[0])
		except:
			Length = ""


		self.SetStatusText('Position: %s      Feature: %s' % (Position, Feature), 0) #text in first field
	
		if float(Length)/3 == 1: #if one triplet is selected, show the AA
			AA = ': %s' % dna.Translate(self.dnaview.stc.GetSelectedText())
		else:
			AA = ''
		
		self.SetStatusText('Selection: %s to %s,   %s bp,   %.1f AA%s' % (SelectionFrom, SelectionTo, Length, float(Length)/3, AA), 1) #text in second field

#		else:
#			self.SetStatusText('', 0)
#			self.SetStatusText('', 1)		


######### get and set methods #########

	#DNA#
	def get_dna_selection(self):
		'''Method for getting which DNA range is currently selected'''
		return genbank.dna_selection

		## I should probably modify this method to broadcast by pypub		

	def set_dna_selection(self, selection):
		'''Method for selecting a certain DNA range'''
		#input needs to be a touple of two values
		genbank.dna_selection = selection
		start = selection[0]
		finish = selection[1]
#		print('Selection from %s to %s') % (start, finish)


		## I should probably modify this method to broadcast by pypub

#	#Feature#
#	def get_feature_selection(self):
#		'''Get index of currently selected feature, if any'''
#		return genbank.feature_selection #returns integer


#	def set_feature_selection(self, index):
#		'''Set currently selected feature'''
#		assert type(index) == int, "Error, index must be an integer."
#		#this is currently independent from DNA selection
#		num_features = len(genbank.gb.get_all_features()) #number of features already present
#		if index == -1 and num_features == 0: #just so that it is easier to select the last feature
#			index = 0
#		elif index == -1 and num_features != 0: #just so that it is easier to select the last feature			
#			index = num_features-1

#		genbank.feature_selection = int(copy.copy(index))
#		#add logic to find first and last position for feature and make DNA selection match.
#		#print('Feature "%s" selected') % (self.get_feature_label(genbank.feature_selection))

######### cut, paste, copy methods ########

	def cut(self, evt):
		'''Check which panel is select and cut accordingly'''
		control = wx.Window.FindFocus() #which field is selected?
		if control == self.searchinput: #the searchbox
			start, finish = self.searchinput.GetSelection()
			if start != -2 and finish != -2: #must be a selection
				pyperclip.copy(self.searchinput.GetValue()[start:finish])
				control.SetValue(self.searchinput.GetValue()[:start]+self.searchinput.GetValue()[finish:])
		elif control == self.dnaview.stc: #the main dna window	
			self.dnaview.cut()

	def paste(self, evt):
		'''Check which panel is select and paste accordingly'''
		control = wx.Window.FindFocus() #which field is selected?
		if control == self.searchinput: #the searchbox
			control.SetValue(pyperclip.paste())
		elif control == self.dnaview.stc: #the main dna window
			self.dnaview.paste()		

	def copy(self, evt):
		'''Check which panel is select and copy accordingly'''
		control = wx.Window.FindFocus() #which field is selected?
		if control == self.searchinput: #the searchbox
			start, finish = self.searchinput.GetSelection()
			if start != -2 and finish != -2: #must be a selection
				pyperclip.copy(self.searchinput.GetValue()[start:finish])
		elif control == self.dnaview.stc: #the main dna window	
			self.dnaview.copy()

	def cut_reverse_complement(self, evt):
		self.dnaview.cut_reverse_complement()

	def paste_reverse_complement(self, evt):
		self.dnaview.paste_reverse_complement()

	def copy_reverse_complement(self, evt):
		self.dnaview.copy_reverse_complement()

#	def delete(self, evt):
#		self.dnaview.delete()

#		print('deleted')

	def reverse_complement_selection(self, evt):
		self.dnaview.reverse_complement_selection()

	def select_all(self, evt):
		self.dnaview.select_all()

##########################################

	def make_outputpopup(self):
		'''Creates a popup window in which output can be printed'''
		self.outputframe = wx.Frame(None, title="Output Panel") # creation of a Frame with a title
		self.output = output.create(self.outputframe, style=wx.VSCROLL|wx.HSCROLL) # creation of a richtextctrl in the frame
		

############## Info methods

	def info(self, event):
		string = '''
This file is part of DNApy. DNApy is a DNA editor written purely in python. 
The program is intended to be an intuitive, fully featured, 
extendable, editor for molecular and synthetic biology.  
Enjoy!

copyright (C) 2014  Martin Engqvist | 

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
LICENSE:
This file is part of DNApy.
DNApy is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.
 
DNApy is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Library General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software Foundation,
Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Get source code at: https://github.com/0b0bby0/DNApy

'''
		self.make_outputpopup()		
		self.output.write(string+'\n', 'Text')
		self.outputframe.Show()
		
	def IUPAC_codes(self, event):
		'''Prints info about the IUPAC codes for DNA and amino acids'''
		string = '''
Nucleotide code 	Base
A 				Adenine
C 				Cytosine
G 				Guanine
T (or U) 			Thymine (or Uracil)
R 				A or G
Y 				C or T
S 				G or C
W 				A or T
K 				G or T
M 				A or C
B 				C or G or T
D 				A or G or T
H 				A or C or T
V 				A or C or G
N 				any base
. or - 				gap

Amino acid 	Three letter 	Amino acid
A 			Ala 			Alanine
C 			Cys 			Cysteine
D 			Asp 			Aspartic Acid
E 			Glu 			Glutamic Acid
F 			Phe 			Phenylalanine
G 			Gly 			Glycine
H 			His 			Histidine
I 			Ile 			Isoleucine
K 			Lys 			Lysine
L 			Leu 			Leucine
M 			Met 			Methionine
N 			Asn 			Asparagine
P 			Pro 			Proline
Q 			Gln 			Glutamine
R 			Arg 			Arginine
S 			Ser 			Serine
T 			Thr 			Threonine
V 			Val 			Valine
W 			Trp 			Tryptophan
Y 			Tyr 			Tyrosine
'''

		self.make_outputpopup()		
		self.output.write(string+'\n', 'Text')
		self.outputframe.Show()


	
	def codon_table(self, event):
		'''Prints the standard codon table for translating DNA to protein'''
		string = '''
Put Table here
'''
		self.make_outputpopup()	
		self.output.write(string+'\n', 'Text')
		self.outputframe.Show()

	def mixed_base_codon_dlg(self, event):
		"""Make popup dialog where amino acids can be chosen and the resulting mixed base codon is computed"""
		self.dlg = wx.Frame(None, title="Find mixed base codons", size=(400, 500)) # creation of a Frame with a title
		self.mixed_base_codon = mixed_base_codons_GUI.MixedBaseCodon(self.dlg)	
		self.dlg.Show()
		self.dlg.Center()

	def uppercase(self, event):
		self.dnaview.uppercase()

	def lowercase(self, event):
		self.dnaview.lowercase()

	def translate_selection(self, event):
		self.dnaview.translate_selection()

	def translate_selection_reverse_complement(self, event):
		self.dnaview.translate_selection_reverse_complement()

	def translate_feature(self, event):
		self.dnaview.translate_feature()
#####################################

	def Undo(self, evt):
		genbank.gb.Undo()
		self.update_globalUI()


	def Redo(self, evt):
		genbank.gb.Redo()
		self.update_globalUI()

######################################
	def listen_to_updateUI(self, msg):
		self.update_ownUI()


	def update_ownUI(self):
		#and now actually update the UI
		print('filever index', genbank.gb.get_file_version_index())
		if genbank.gb.get_file_version_index() <= 0: #if there are no undos available
			self.frame_1_toolbar.EnableTool(513, False)
		else:
			self.frame_1_toolbar.EnableTool(513, True)

		if genbank.gb.get_file_version_index() >= len(genbank.gb.file_versions)-1: #if there are no redos available
			self.frame_1_toolbar.EnableTool(514, False)
		else:
			self.frame_1_toolbar.EnableTool(514, True)

	def update_globalUI(self):
		pub.Publisher.sendMessage('from_main', '')

######################################

	def list_features(self, evt):
		'''List all features in output panel'''
		self.make_outputpopup()
#		tabtext = str(genbank.gbviewer.GetPageText(genbank.gbviewer.GetSelection()))
		tabtext = 'Replace!'
		self.output.write('%s | List features\n' % tabtext, 'File')
		featurelist = genbank.gb.ListFeatures()
		for feature in featurelist:
			self.output.write('%s\n' % feature, 'Text')
		self.outputframe.Show()

	def view_genbank(self, evt):
		'''View the genbank file as text'''
		dlg = wx.Frame(self, id=wx.ID_ANY, title="Edit Qualifier", size=(600,600)) #make frame
		genbankview = genbank_GUI.MyPanel(dlg, style=wx.VSCROLL|wx.HSCROLL) #put the genbank view panel inside
		genbankview.rtc.SetEditable(False)	
		sizer = wx.BoxSizer(wx.VERTICAL) #make sizer
		sizer.Add(item=genbankview, proportion=-1, flag=wx.EXPAND)	#put panel in sizer
		dlg.SetSizer(sizer) #assign sizer to window
		dlg.Show()

	def view_output(self, evt):	
		'''Make an output window in which things can be printed'''
		self.outputframe = wx.Frame(self, title="Output Panel") # creation of a Frame with a title
		self.outputwindow = output.create(self.outputframe, style=wx.VSCROLL|wx.HSCROLL) # creation of a richtextctrl in the frame
		sizer = wx.BoxSizer(wx.VERTICAL) #make sizer
		sizer.Add(item=self.outputwindow, proportion=-1, flag=wx.EXPAND)	#put panel in sizer
		self.outputframe.SetSizer(sizer) #assign sizer to window
		self.outputframe.Show()		

######### Toolbar and Menu ############

	def __generate_toolbar(self):
		'''For generating toolbar with buttons'''
				
		self.frame_1_toolbar = wx.ToolBar(self, wx.ID_ANY, style=wx.TB_HORIZONTAL|wx.TB_FLAT|wx.TB_DOCKABLE)

   		#syntax for toolbar
   		#AddLabelTool(self, id, label, bitmap, bmpDisabled, kind, shortHelp, longHelp, clientData) 
   		

   		
   		#New Document
   		self.frame_1_toolbar.AddLabelTool(500, "New Document", wx.Bitmap(files['default_dir']+"/icon/new.png", wx.BITMAP_TYPE_ANY), wx.NullBitmap, wx.ITEM_NORMAL, 'New File', "New File") #last one is the one displayed in status bar
   		wx.EVT_TOOL(self, 500, self.new_file)
		#Open File
   		self.frame_1_toolbar.AddLabelTool(501, "Open File", wx.Bitmap(files['default_dir']+"/icon/open.png", wx.BITMAP_TYPE_ANY), wx.NullBitmap, wx.ITEM_NORMAL, 'Load File', 'Load File')
   		wx.EVT_TOOL(self, 501, self.open_file)
		#Save current file
   		self.frame_1_toolbar.AddLabelTool(502, "Save current file", wx.Bitmap(files['default_dir']+"/icon/save.png", wx.BITMAP_TYPE_ANY), wx.NullBitmap, wx.ITEM_NORMAL, 'Save File', 'Save File')
   		wx.EVT_TOOL(self, 502, self.save_file)
		#Save as
   		self.frame_1_toolbar.AddLabelTool(503, "Save as", wx.Bitmap(files['default_dir']+"/icon/saveas.png", wx.BITMAP_TYPE_ANY), wx.NullBitmap, wx.ITEM_NORMAL, 'Save File As', 'Save File As')
   		wx.EVT_TOOL(self, 503, self.save_as_file)
		#cut
   		self.frame_1_toolbar.AddLabelTool(504, "cut", wx.Bitmap(files['default_dir']+"/icon/cut.png", wx.BITMAP_TYPE_ANY), wx.NullBitmap, wx.ITEM_NORMAL, 'cut', 'cut')
   		wx.EVT_TOOL(self, 504, self.cut)
		#copy
   		self.frame_1_toolbar.AddLabelTool(505, "copy", wx.Bitmap(files['default_dir']+"/icon/copy.png", wx.BITMAP_TYPE_ANY), wx.NullBitmap, wx.ITEM_NORMAL, 'copy', 'copy')
   		wx.EVT_TOOL(self, 505, self.copy)
		#paste
   		self.frame_1_toolbar.AddLabelTool(506, "paste", wx.Bitmap(files['default_dir']+"/icon/paste.png", wx.BITMAP_TYPE_ANY), wx.NullBitmap, wx.ITEM_NORMAL, 'paste', 'paste')
   		wx.EVT_TOOL(self, 506, self.paste)
   		#Undo
   		self.frame_1_toolbar.AddLabelTool(513, "Undo", wx.Bitmap(files['default_dir']+"/icon/undo.png", wx.BITMAP_TYPE_ANY), wx.NullBitmap, wx.ITEM_NORMAL, 'Undo', 'Undo')
   		wx.EVT_TOOL(self, 513, self.Undo)   
   		#Redo
   		self.frame_1_toolbar.AddLabelTool(514, "Redo", wx.Bitmap(files['default_dir']+"/icon/redo.png", wx.BITMAP_TYPE_ANY), wx.NullBitmap, wx.ITEM_NORMAL, 'Redo', 'Redo')
   		wx.EVT_TOOL(self, 514, self.Redo) 
		#Print current window
#   		self.frame_1_toolbar.AddLabelTool(510, "Print current window", wx.Bitmap(files['default_dir']+"/icon/print.png", wx.BITMAP_TYPE_ANY), wx.NullBitmap, wx.ITEM_NORMAL, 'Print Current Window', 'Print Current Window')
 #  		wx.EVT_TOOL(self, 510, self.print_setup)

   		self.frame_1_toolbar.AddCheckTool(516, wx.Bitmap(files['default_dir']+"/icon/plasmid.png", wx.BITMAP_TYPE_ANY), wx.Bitmap(files['default_dir']+"/icon/plasmid.png", wx.BITMAP_TYPE_ANY), 'Plasmid view', 'Plasmid view')
   		wx.EVT_TOOL(self, 516, self.toggle_plasmid_view)

   		self.frame_1_toolbar.AddCheckTool(511, wx.Bitmap(files['default_dir']+"/icon/search.png", wx.BITMAP_TYPE_ANY), wx.Bitmap(files['default_dir']+"/icon/search.png", wx.BITMAP_TYPE_ANY), 'Find', 'Find')
   		wx.EVT_TOOL(self, 511, self.toggle_searchandmutate_toolbar)

   		self.frame_1_toolbar.AddCheckTool(512, wx.Bitmap(files['default_dir']+"/icon/mutate.png", wx.BITMAP_TYPE_ANY), wx.Bitmap(files['default_dir']+"/icon/search.png", wx.BITMAP_TYPE_ANY), 'Mutate', 'Mutate')
   		wx.EVT_TOOL(self, 512, self.toggle_searchandmutate_toolbar)

		
		self.frame_1_toolbar.Realize()
		
		self.frame_1_toolbar.EnableTool(502, 0) #save current file button
		self.frame_1_toolbar.EnableTool(503, 0) #save as button
		self.frame_1_toolbar.EnableTool(504, 0) #cut button
		self.frame_1_toolbar.EnableTool(505, 0) #copy button
		self.frame_1_toolbar.EnableTool(506, 0)	#paste button
		self.frame_1_toolbar.EnableTool(511, 0) #search button
		self.frame_1_toolbar.EnableTool(512, 0) #mutate button
		self.frame_1_toolbar.EnableTool(513, 0) #Undo button
		self.frame_1_toolbar.EnableTool(514, 0) #Redo button

	def __generate_searchandmutate_toolbar(self):
		##### Toolbar 2 #####
		self.nucleotideoraminoacidSelection = 0

		self.frame_2_toolbar = wx.ToolBar(self, wx.ID_ANY, style=wx.TB_HORIZONTAL|wx.TB_FLAT|wx.TB_DOCKABLE)
		
		self.findormutSelection = 'Find'
		self.add_search_tools(self.findormutSelection)
			

		self.frame_2_toolbar.Realize()
		self.frame_2_toolbar.Hide()

	def toggle_plasmid_view(self, event):
		'''When the plasmid view button is toggled, show/hide the plasmid veiw'''
		

		if self.frame_1_toolbar.GetToolState(516) == True:
			self.plasmid_frame.Show()

		elif self.frame_1_toolbar.GetToolState(516) == False:
			try:
				self.plasmid_frame.Hide()
			except:
				pass


	def add_search_tools(self, typeof):
		'''Adds tools to the find/mutate toolbar. A string "find" or "mutate" is passed to determine which toolbar to build'''
		featurelist = ['Molecule']
		try:
			#make a list of all feature names
			features = genbank.gb.get_all_features()
			for entry in features:
				featurelist.append(entry['qualifiers'][0].split('=')[1])
		except:
			pass

		if typeof == 'Find':
			self.nucleotideoraminoacid = wx.ComboBox(self.frame_2_toolbar, id=601, size=(120, 28), choices=['Nucleotide', 'Amino Acid', 'Feature'], style=wx.CB_READONLY)
		elif typeof == 'Mutate':
			self.nucleotideoraminoacid = wx.ComboBox(self.frame_2_toolbar, id=601, size=(120, 28), choices=['Nucleotide', 'Amino Acid'], style=wx.CB_READONLY)
		self.frame_2_toolbar.AddControl(self.nucleotideoraminoacid)
		self.nucleotideoraminoacid.SetSelection(self.nucleotideoraminoacidSelection)
		wx.EVT_COMBOBOX(self, 601, self.OnChangeSearchParams)

		#'input'
		self.searchinput=wx.TextCtrl(self.frame_2_toolbar, id=wx.ID_ANY, size=(100, 28), value="")
		self.frame_2_toolbar.AddControl(self.searchinput)
	

		#'in'
		nucleotideoraa = self.nucleotideoraminoacid.GetValue()
		if nucleotideoraa == 'Nucleotide' or nucleotideoraa == 'Amino Acid':
			self.inbox=wx.TextCtrl(self.frame_2_toolbar, id=wx.ID_ANY, size=(25, 28), value="in")
			self.frame_2_toolbar.AddControl(self.inbox)
			self.inbox.SetEditable(False)
			
		#featurebox
		if nucleotideoraa == 'Nucleotide' or nucleotideoraa == 'Amino Acid':
			try:
				#get features...
				self.featurebox = wx.ComboBox(self.frame_2_toolbar, id=603, size=(120, 28), choices=['Molecule' and features], style=wx.CB_READONLY)
			except:
				self.featurebox = wx.ComboBox(self.frame_2_toolbar, id=603, size=(120, 28), choices=featurelist, style=wx.CB_READONLY)
			self.frame_2_toolbar.AddControl(self.featurebox)
			self.featurebox.SetSelection(0)
		
		#'go'
		if typeof == 'Find':
			#find previous
	   		self.frame_2_toolbar.AddLabelTool(507, "Find previous", wx.Bitmap(files['default_dir']+"/icon/up.png", wx.BITMAP_TYPE_ANY), wx.NullBitmap, wx.ITEM_NORMAL, 'Find previous', 'Find previous')
	   		wx.EVT_TOOL(self, 507, self.find_previous)
			#find next
	   		self.frame_2_toolbar.AddLabelTool(509, "Find next", wx.Bitmap(files['default_dir']+"/icon/down.png", wx.BITMAP_TYPE_ANY), wx.NullBitmap, wx.ITEM_NORMAL, 'Find next', 'Find next')
	   		wx.EVT_TOOL(self, 509, self.find_next)
			#search glass
	   		self.frame_2_toolbar.AddLabelTool(604, "Find", wx.Bitmap(files['default_dir']+"/icon/search.png", wx.BITMAP_TYPE_ANY), wx.NullBitmap, wx.ITEM_NORMAL, 'Find', 'Find')
	   		wx.EVT_TOOL(self, 604, self.find)
		elif typeof == 'Mutate':
			self.frame_2_toolbar.AddLabelTool(604, "Mutate", wx.Bitmap(files['default_dir']+"/icon/up.png", wx.BITMAP_TYPE_ANY), wx.NullBitmap, wx.ITEM_NORMAL, 'Mutate', 'Mutate')
	   		wx.EVT_TOOL(self, 604, self.mutate)

	def find(self, evt):
		'''Find nucleotide in molecule'''
		searchtype = self.nucleotideoraminoacid.GetValue() #type of search
		searchframe = int(self.featurebox.GetSelection())-1 #where to search
		searchstring = self.searchinput.GetValue()
		if searchstring.isdigit():
			searchstring = int(searchstring)

		if searchtype == 'Nucleotide':
			genbank.search_hits = genbank.gb.FindNucleotide(searchstring, searchframe)
		elif searchtype == 'Amino Acid':
			genbank.search_hits = genbank.gb.FindAminoAcid(searchstring, searchframe)
		elif searchtype == 'Feature':
			genbank.search_hits = genbank.gb.FindFeature(searchstring)

		self.update_globalUI()
		start, finish = self.get_dna_selection()
		self.dnaview.stc.SetSelection(start, finish)

	
	def find_previous(self, evt):
		'''Select prevous search hit'''
		genbank.gb.find_previous()
		start, finish = self.get_dna_selection()
		self.dnaview.stc.SetSelection(start, finish)

	def find_next(self, evt):
		'''Select next search hit'''
		genbank.gb.find_next()
		start, finish = self.get_dna_selection()
		self.dnaview.stc.SetSelection(start, finish)

	def mutate(self, evt):
		'''Mutate DNA or protein'''
		mutationtype = self.nucleotideoraminoacid.GetValue() #type of mutation
		mutationframe = int(self.featurebox.GetSelection())-1 #where to mutate (feature or molecule)
		mutationinput = self.searchinput.GetValue() #which mutation to perform
		if ',' in mutationinput: #input allows for many comma-seperated mutations
			mutation_list = mutationinput.split(',')
		else:
			mutation_list = [mutationinput]

		for mutation in mutation_list:
			if mutationtype == 'Nucleotide':
				genbank.gb.mutate('D', mutationframe, mutation)
			elif mutationtype == 'Amino Acid':
				genbank.gb.mutate('A', mutationframe, mutation)
		
		self.update_globalUI()

	def OnChangeSearchParams(self, evt):
		'''When changes are made to options in the searchbar, update which downstream options are available'''
		#get selection for find or mutate
		if evt.GetId() == 601:
			self.nucleotideoraminoacidSelection = self.nucleotideoraminoacid.GetSelection()

		self.frame_2_toolbar.ClearTools()
		self.add_search_tools(self.findormutSelection)


	def toggle_searchandmutate_toolbar(self, evt):
		'''Toggles the visibility of the search toolbar'''
		if self.frame_1_toolbar.GetToolState(511) == True and self.frame_1_toolbar.GetToolState(512) == False:
			self.frame_2_toolbar.Show()
			self.frame_1_toolbar.EnableTool(512, 0)
			self.frame_2_toolbar.ClearTools()
			self.findormutSelection = 'Find'
			self.add_search_tools(self.findormutSelection)

		elif self.frame_1_toolbar.GetToolState(511) == False and self.frame_1_toolbar.GetToolState(512) == True:
			self.frame_2_toolbar.Show()
			self.frame_1_toolbar.EnableTool(511, 0)
			self.frame_2_toolbar.ClearTools()
			self.findormutSelection = 'Mutate'
			self.add_search_tools(self.findormutSelection)

		elif self.frame_1_toolbar.GetToolState(511) == False and self.frame_1_toolbar.GetToolState(512) == False:
			self.frame_1_toolbar.EnableTool(511, 1)
			self.frame_1_toolbar.EnableTool(512, 1)
			self.frame_2_toolbar.Hide()

		self.Layout()


	def create_menu(self):     #method for creating menu
		self.menubar = wx.MenuBar()
		fileitem = wx.Menu()
			
		#new document
		fileitem.Append(1, "New\tCtrl+Q", "New File")
		wx.EVT_MENU(self, 1, self.new_file)

		#open document
		fileitem.Append(2, "Open\tCtrl+O", "Open File")
		wx.EVT_MENU(self, 2, self.open_file)
		fileitem.AppendSeparator()

		#save document
		fileitem.Append(3, "Save\tCtrl+S", "Save current file")
		wx.EVT_MENU(self, 3, self.save_file)

		#save document as
		fileitem.Append(4, "Save as\tShift+Ctrl+S", "Save a copy of current file")
		wx.EVT_MENU(self, 4, self.save_as_file)

		#save all
#		fileitem.Append(5, "Save all", "Save all open tabs")
#		wx.EVT_MENU(self, 5, self.save_all)
#		fileitem.AppendSeparator()

		#close single
#		fileitem.Append(5, "Close", "Close current file")
#		wx.EVT_MENU(self, 5, self.dnaview.close_file)

		#close all
#		fileitem.Append(6, "Close all", "Close all tabs")
#		wx.EVT_MENU(self, 6, self.dnaview.close_all)
#		fileitem.AppendSeparator()

		#quit
		fileitem.Append(7, "Quit", "Quit program")
		wx.EVT_MENU(self, 7, self.quit)

		self.menubar.Append(fileitem, "&File")

		######################### For 'View' menu item #############################################		
		self.view = wx.Menu()

		#view genbank file
		self.view.Append(201, "View GenBank", "View GenBank")
		wx.EVT_MENU(self, 201, self.view_genbank)

		#view output panel
		self.view.Append(202, "View Output", "View Output")
		wx.EVT_MENU(self, 202, self.view_output)
		

		self.menubar.Append(self.view, "&View")

		######################### For 'Edit DNA' menu item #############################################
		self.edit = wx.Menu()
		#undo
		self.edit.Append(9, "Undo\tCtrl+Z", "Undo")
		wx.EVT_MENU(self, 9, self.Undo)

		#redo
		self.edit.Append(10, "Redo\tShift+Ctrl+Z", "Redo")
		wx.EVT_MENU(self, 10, self.Redo)
		self.edit.AppendSeparator() #________________________devider

		#cut
		self.edit.Append(11, "cut\tCtrl+X", "cut selected DNA")
		wx.EVT_MENU(self,11, self.cut)

		#copy
		self.edit.Append(12, "copy\tCtrl+C", "copy selected DNA")
		wx.EVT_MENU(self, 12, self.copy)

		#paste
		self.edit.Append(13, "paste\tCtrl+V", "paste DNA")
		wx.EVT_MENU(self, 13, self.paste)
		
		#cut reverse complement
		self.edit.Append(111, "cut Rev-Comp\tCtrl+Shift+X", "cut reverse-complement of selected DNA")
		wx.EVT_MENU(self,111, self.cut_reverse_complement)

		#copy reverse complement
		self.edit.Append(121, "copy Rev-Comp\tCtrl+Shift+C", "copy reverse-complement of selected DNA")
		wx.EVT_MENU(self, 121, self.copy_reverse_complement)		
		
		#paste reverse complement
		self.edit.Append(131, "paste Rev-Comp\tCtrl+Shift+V", "paste reverse-complement of DNA")
		wx.EVT_MENU(self, 131, self.paste_reverse_complement)

		#reverse-complement selection
		self.edit.Append(141, "Rev-Comp selection\tCtrl+R", "Reverse-complement seleected DNA")
		wx.EVT_MENU(self,141, self.reverse_complement_selection)
		self.edit.AppendSeparator() #________________________devider
		
		#select all
		self.edit.Append(14, "Select all\tCtrl+A", "Select all text")
		wx.EVT_MENU(self, 14, self.select_all)
		self.edit.AppendSeparator() #________________________devider

		#uppercase
		self.edit.Append(34, "uppercase\tCtrl+W", "Convert selected text to uppercase")
		wx.EVT_MENU(self, 34, self.uppercase)

		#lowercase
		self.edit.Append(35, "lowercase\tCtrl+E", "Convert selected text to lowercase")
		wx.EVT_MENU(self, 35, self.lowercase)
		self.edit.AppendSeparator() #________________________devider

		#Find mixed base codon
		self.edit.Append(135, "Find mixed base codon", "Find mixed base codon")
		wx.EVT_MENU(self, 135, self.mixed_base_codon_dlg)
		



		self.menubar.Append(self.edit, "Edit DNA")
		
			
	
		######## Features menu item ########
		self.features = wx.Menu()
		
		self.features.Append(40, "List Features", "List all features in file")
		wx.EVT_MENU(self,40, self.list_features)
		
#		self.features.Append(41, "Edit Features", "Edit features in file")
#		wx.EVT_MENU(self,41, self.edit_features)		
		
		self.menubar.Append(self.features, "Features")
		


		######## Protein menu item #########
		self.protein = wx.Menu()
		
		#translate
		self.protein.Append(30, "Translate\tCtrl+T", "Translate DNA to protein")
		wx.EVT_MENU(self,30, self.translate_selection)

		#translate reverse complement
		self.protein.Append(31, "Translate Rev-Comp\tCtrl+Shift+T", "Translate DNA to protein")
		wx.EVT_MENU(self,31, self.translate_selection_reverse_complement)

		#translate feature
		self.protein.Append(32, "Translate feature", "Translate DNA feature to protein")
		wx.EVT_MENU(self,32, self.translate_feature)

		self.menubar.Append(self.protein, "Protein")



		########## For 'Help' menu item #############
		self.help = wx.Menu()
		#about program
		self.help.Append(21, "About", "About this program")
		wx.EVT_MENU(self, 21, self.info)

		#print IUPAC codes for dna and AA
		self.help.Append(22, "IUPAC codes", "IUPAC codes of DNA and amino acids")
		wx.EVT_MENU(self, 22, self.IUPAC_codes)


		self.menubar.Append(self.help, "Help")

		self.SetMenuBar(self.menubar)










##### main loop
class MyApp(wx.App):
	def OnInit(self):
		frame = MyFrame(None, -1, "DNApy")
		frame.Show(True)
		self.SetTopWindow(frame)
		return True



if __name__ == '__main__': #if script is run by itself and not loaded	
	app = MyApp(0)
	app.MainLoop()
