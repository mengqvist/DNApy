#!/usr/bin/python


#DNApy is a DNA editor written purely in python. 
#The program is intended to be an intuitive, fully featured, 
#extendable, editor for molecular and synthetic biology.  
#Enjoy!
#
#copyright (C) 2014  Martin Engqvist 
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

from base_class import DNApyBaseClass
import wx
import re
import string
import genbank
import math
import collections



class EnzymeSelector(DNApyBaseClass):
	"""
	Class to select restriction enzymes.
	"""
	def __init__(self, parent, enzymeClass, id):
		self.parent = parent
		wx.Panel.__init__(self, parent)
		
		# make enzymes accesible
		self.enzymeClass = enzymeClass
		
		# enzyme selector GUI
		# we have three columns, each is in a flexGrid
		font             = wx.Font(9, wx.DEFAULT, wx.NORMAL, wx.BOLD)
		# veryfirst column - came after nameing the first first
		veryfirstColumn      = wx.BoxSizer(wx.VERTICAL)
		self.txt1     = wx.StaticText(self, -1, 'choose enzymeset to display')
		self.txt1.SetFont(font)
		
		self.radio1      = wx.RadioButton(self,id=-1,label="all")
		self.radio2      = wx.RadioButton(self,id=-1,label="singlecutters")
		self.radio3      = wx.RadioButton(self,id=-1,label="doublecutters")
		self.radio4      = wx.RadioButton(self,id=-1,label="commercial selection")
		
		self.radio1.Bind(wx.EVT_RADIOBUTTON,lambda event: self.showOnly(event, "all"))
		self.radio2.Bind(wx.EVT_RADIOBUTTON,lambda event:self.showOnly(event, 1))
		self.radio3.Bind(wx.EVT_RADIOBUTTON,lambda event: self.showOnly(event, 2))	
		self.radio4.Bind(wx.EVT_RADIOBUTTON,lambda event: self.showOnly(event, "commercial"))	
	
		veryfirstColumn.Add(item=self.txt1)
		veryfirstColumn.Add(item=self.radio1)
		veryfirstColumn.Add(item=self.radio4)
		veryfirstColumn.Add(item=self.radio2)
		veryfirstColumn.Add(item=self.radio3)

		

		
		# first column:
		firstColumn      = wx.BoxSizer(wx.VERTICAL)
		self.labellb     = wx.StaticText(self, -1, '')
		self.labellb.SetFont(font)
		self.lb          = wx.ListBox(self,
		                              size=(180, 300))
		self.lb.Bind(wx.EVT_LISTBOX, self.onSelect)
		
		# textfield for automatet text input
		self.txt = wx.TextCtrl(self, -1, size=(180,-1),style=wx.TE_PROCESS_ENTER)
		self.txt.Bind(wx.EVT_TEXT_ENTER, self.selectEnter)
		self.txt.Bind(wx.EVT_TEXT, self.keystroke)
		self.txt.SetFocus()
		
		# add items to sizer "firstColumn"
		firstColumn.Add(item=self.labellb)
		firstColumn.Add(item=self.lb)
		firstColumn.Add(item=self.txt)
		
		
		# second column 
		secondColumn     = wx.BoxSizer(wx.VERTICAL)
		self.oneadd      = wx.Button(self,-1, "add")		# add enzyme
		self.oneadd.Bind(wx.EVT_BUTTON, self.addOne)
		self.lb.Bind(wx.EVT_LISTBOX_DCLICK, self.addOne)
		
		self.remove      = wx.Button(self,-1, "remove")		# remove enzyme
		self.remove.Bind(wx.EVT_BUTTON, self.removeOne)

		
		self.removeAllB   = wx.Button(self,-1, "remove all")	# remove all
		self.removeAllB.Bind(wx.EVT_BUTTON,self.removeAll)
		
		# add buttons to sizer "secondColumn"
		secondColumn.Add(item=self.oneadd)
		secondColumn.Add(item=self.remove)
		secondColumn.Add(item=self.removeAllB)
		

		# thirt column
		thirdColumn      = wx.BoxSizer(wx.VERTICAL)
		self.labellb2     = wx.StaticText(self, -1, 'selected enzymes')
		self.labellb2.SetFont(font)
		self.lb2 = wx.ListBox(self,
		        	 size=(180, 300))
		self.lb2.Bind(wx.EVT_LISTBOX_DCLICK, self.removeOne)


		thirdColumn.Add(item=self.labellb2)
		thirdColumn.Add(item=self.lb2)

		
		# button to cancel or finish editing:
		self.cancel   = wx.Button(self,wx.ID_CANCEL)
		self.ok       = wx.Button(self,wx.ID_OK)
		

		# every column now has to be added to the grid
		#Use sizers add content in the correct arrangement
		sizerClose  = wx.BoxSizer(wx.HORIZONTAL)
		sizerClose.Add(item=self.cancel)
		sizerClose.Add(item=self.ok)	
					
		# grid settings
		hbox      = wx.BoxSizer(wx.HORIZONTAL)
		gridsizer = wx.FlexGridSizer(rows=2, cols=4, vgap=3, hgap=10)
		gridsizer.AddGrowableCol(0)					# make cols growable
		gridsizer.AddGrowableCol(1)					# make cols growable
		gridsizer.AddGrowableCol(2)					# make cols growable
		gridsizer.AddGrowableCol(3)					# make cols growable
		hbox.Add(gridsizer, 1, wx.EXPAND|wx.ALL, 15)
		
		# add the elements to the grid for flexible display
		gridsizer.Add(veryfirstColumn)      				# row 1, col 1
		gridsizer.Add(firstColumn)      				# row 1, col 1
		gridsizer.Add(secondColumn, 0, wx.ALIGN_CENTER_VERTICAL) 	# row 1, col 2
		gridsizer.Add(thirdColumn)     					# row 1, col 3
		gridsizer.AddSpacer(1) 			 			# row 2, col 1
		gridsizer.AddSpacer(1) 						# row 2, col 2
		gridsizer.AddSpacer(1) 						# row 2, col 2
		gridsizer.Add(item=sizerClose)   				# row 2, col 3
		
		#set sizer
		self.SetSizer(hbox)

		# load the enzymes from emboss
		self.showOnly(None, "all")
		
	
	
	def onButton(self, event):
		"""
		This method is fired when its corresponding button is pressed
		"""
		btn = event.GetEventObject()
		btn.Select()
		label = btn.GetLabel()
		message = "You just selected %s" % label
		dlg = wx.MessageDialog(None, message, 'Message', 
		                       wx.OK|wx.ICON_EXCLAMATION)
		dlg.ShowModal()
		dlg.Destroy()
		
	# Function to change the displayed enzyme in the self.lb first box
	def showOnly(self, event, n):
		
		# ordered dic to store objects
		enzymes2show = collections.OrderedDict()
		
		if n == "all":
			# show all
			enzymes2show = self.enzymeClass.enzymeObj
		elif n == 1:
			# show singlecutters
			enzymes2show = self.findCutters(1)
		elif n == 2:
			# show doublecutters
			enzymes2show = self.findCutters(2)
		elif n == "commercial":
			# show commercial
			with open ("resources/commercialEnzymes.lst", "r") as myfile:
				data=myfile.read().splitlines()
				# check if we know this enzyme
				for e in data:
					if e in self.enzymeClass.enzymeObj:
						enzymes2show[e] = self.enzymeClass.enzymeObj[e]
		else:
			print "what selection?"
	
		# clear list first
		self.lb.Clear()	
	
		# add the enzymes:
		for enzyme in enzymes2show:
			self.AddRestrictionEnzyme(enzymes2show[enzyme])	
			
		
		# reset focus on the text field 
		self.txt.SetFocus()
		return True
	
	####################################################
	# function to find singlecutters or similar
	def findCutters(self, n):
		cutter = collections.OrderedDict()
		dnaseq      = genbank.gb.gbfile["dna"]
		if dnaseq != None:
			# circular dna
			# we just inspect the region +-100 
			# to find enzyme, which cut near 0 in a ciclic plasmid
			if genbank.gb.gbfile['locus']['topology'] == 'circular':
				circularDnaHelper = dnaseq[:100]  # helper of 200bp from 0 to 100
			else:
				circularDnaHelper = ''
					
			wholeDNA2Inspect = '%s%s' % (dnaseq, circularDnaHelper)		     # dna width circular helper added

		
			# loop through every enzyme and look for match positions
			for enzyme in self.enzymeClass.enzymeObj:
				r         = self.enzymeClass.enzymeObj[enzyme].regex

				iterator  = r.finditer(wholeDNA2Inspect)      	# find in dnaseq and circularDnaHelper
				i         = 0 					# counter for occurence
				positions = []				     	# we remember the cut positions for every enzyme, so we do not get false doubles
				modulo    = len(dnaseq)			     	# all positions are modulo len(dnaseq)
				for match in iterator:
					start = (match.start()+1) % modulo
					if start not in positions:
						positions.append(start)
						i = i + 1 			# raise counter
			
				if i == n:					# compare length of hits to the given number
					cutter[enzyme] = self.enzymeClass.enzymeObj[enzyme]
	
	
		return cutter
	
	
	def keystroke(self, e):
		c = "^%s" % self.txt.GetValue()
		#n = chr(e.GetKeyCode())
		#self.txt.SetValue(("%s%s" %(c,n)))
		# serach for hits using regexp:
		regex = re.compile(c, re.IGNORECASE)
		allItems = self.lb.GetItems()
		for x in allItems:
			if re.search(regex, x):
				break
		if allItems.index(x) + 1 < len(allItems):
			self.lb.SetSelection(allItems.index(x))
	
	def selectEnter(self, e):
		self.addOne(e)
		self.txt.SetValue('')
	
	# on select focus back on txt -> not working, why?
	def onSelect(self, event):
		#self.txt.SetFocus()
        	return True


	# adds one item on buttonclick
	def addOne(self, event):
		item  		= self.lb.GetStringSelection()
		obj 		= self.lb.GetClientData(self.lb.GetSelection())
		allItems 	= self.lb2.GetItems()

		# just add the item, if its not already in there
		if item not in allItems:
			self.lb2.Append(item, obj)
			self.resort2list()
	
	# removes all items from 2
	def removeAll(self, event):
		self.lb2.Clear()
		return True
		
		

	# remove selected item from list
	def removeOne(self, event):
		n = self.lb2.GetSelection()
		if n >= 0: # to prevent error when none is selected
			self.lb2.Delete(n)
		

	# adds the given restriction enzyme to the first list
	def AddRestrictionEnzyme(self, item, i=1):
		if i == 1:
			self.lb.Append(item.name, item)
		if i == 2:
			self.lb2.Append(item.name, item)
			

	# can reorder the second list
	# is called by addOne(self)
	def resort2list(self):
		allitems = self.lb2.GetItems()	# get items
		self.lb2.Clear()				# clear lb2
		allSorted = sorted(allitems)	# sort items
		
		for enzyme in allSorted:		# readd items and objects
			self.lb2.Append(enzyme, self.enzymeClass.enzymeObj[enzyme])
	

	# gets called to retrive the sleected infos
	def getSelection(self):
		allItems 	= collections.OrderedDict()
		items 		=  self.lb2.GetItems()
		n = 0
		for i in items:
			allItems[i] = self.lb2.GetClientData(n)
			n = n +1

		return allItems


	def findRestrictionSites(self, selectedEnzymes):

		restrictionsitesList = [] # variable to return the restriction sites
					  # syntax:
					  # [[name, start, stop, cut 1, cut 2, sequence],[...]]
		dnaseq      = genbank.gb.gbfile["dna"]
		
		# circular dna
		# we just inspect the region +-100 
		# to find enzyme, which cut near 0 in a ciclic plasmid
		if genbank.gb.gbfile['locus']['topology'] == 'circular':
			circularDnaHelper = dnaseq[:100]  # helper of 200bp from 0 to 100
		else:
			circularDnaHelper = ''
					
		wholeDNA2Inspect = '%s%s' % (dnaseq, circularDnaHelper)		     # dna width circular helper added

		# loop all the selected enzymes
		for enzyme in selectedEnzymes:
			
			# load the regexp of the enzyme
			r           = self.enzymes[enzyme]["regexp"]
			# get the cut position
			offset1     = self.enzymes[enzyme]["cut5_1"]
			offset2     = self.enzymes[enzyme]["cut5_2"] 	# offset 2 not yet implimented

			# handle the cuts
			iterator    = r.finditer(wholeDNA2Inspect)      # find in dnaseq and circularDnaHelper
			for match in iterator:
				newEnz = []
				modulo = len(dnaseq)		# all positions are modulo len(dnaseq)
				# add offset to match position
				newEnz.append(str(enzyme))
				newEnz.append((match.start()+1) % modulo)
				newEnz.append(match.end() % modulo)
				newEnz.append((match.start() + offset1) % modulo)	# cuts after this base

				if offset2 == 0: # if it just cuts once, we'll say None to the second cut
					newEnz.append(None)
				else:
					newEnz.append((match.start() + offset2) % modulo)

				newEnz.append(wholeDNA2Inspect[match.start():match.end()])
				
				# save the new Enzyme if its new
				if newEnz not in restrictionsitesList:
					restrictionsitesList.append(newEnz)
			
			
		return restrictionsitesList
	

	

	def handleCut(self, enzyme, n):
		return True

	def update_globalUI(self):
		'''
		Method should be modified as to update other panels in response to changes in own panel.
		'''
		# here we should add the code to draw the restricion sites
		pass


	def update_ownUI(self):
		'''
		Updates to own panel can be made here.
		'''
		pass
		
		
		
		
class EnzymeSelectorDialog(wx.Dialog):
	'''A class that puts the Enzyme Selector capabilities in a dialog.'''
	def __init__(self, parent, title, oldSelection, enzymeClass):
		super(EnzymeSelectorDialog, self).__init__(parent=parent,id=wx.ID_ANY, title=title, size=(750, 430)) 		

		#add the panel (containing all the buttons/lists/interactive elements
		self.content = EnzymeSelector(self, enzymeClass,id=wx.ID_ANY,)	#get the feature edit panel
		

		# add the old selection:
		for item in oldSelection:

			self.content.AddRestrictionEnzyme(enzymeClass.enzymeObj[item],2)

		
		#add sizer
		sizer = wx.BoxSizer(wx.VERTICAL)
		sizer.Add(item=self.content, proportion=0, flag=wx.EXPAND|wx.ALL)

		#set sizer
		self.SetSizer(sizer)
		
		# pass the old selection to the window
		
		# show the window (enables the ok and cancel buttons)
		#result = EnzymeSelector.Show(self.content)	


		

	# this function is called to retrive the selected enzymes
	def GetSelection(self):
		'''
		Get the enzyme selection.
		Used to actually extract info from the dialog.
		'''
		return self.content.getSelection()

	# this function is also called by the main.py
	# it returns the cutting position of the enzyme
	def drawRestriction(self, enzymes):
		return self.content.findRestrictionSites(enzymes)
		




##############################################################################################################
#
# this class can show a dialog for simulation of an agarose gel
# it is mainly for calculating fragments and showing a gel image
# 
# usefull for students and researchers to visualize gel patterns
#
# can only be called with loaded file
# 
class EnzymeDigestionDialog(wx.Dialog):
	'''A class that puts the Enzyme digestion capabilities in a dialog.'''
	def __init__(self, parent, title, Enzymes,enzymeClass):
		super(EnzymeDigestionDialog, self).__init__(parent=parent,id=wx.ID_ANY, title=title, size=(600, 430)) 		
		
		# get width and height to pass to child class
		width, height = self.GetClientSize()	
		#add the panel, containing all the buttons/lists/interactive elements
		self.content = EnzymeDigestion(self,Enzymes, enzymeClass, width, height, id=wx.ID_ANY)	
		
		#add sizer
		sizer = wx.BoxSizer(wx.VERTICAL)
		sizer.Add(item=self.content, proportion=0, flag=wx.EXPAND|wx.ALL)

		#set sizer
		self.SetSizer(sizer)
		
    


class EnzymeDigestion(DNApyBaseClass):
	"""
	Class to show enzyme digestion
	"""
	def __init__(self, parent,Enzymes,enzymeClass,  width, height, id):
		self.parent = parent
		wx.Panel.__init__(self, parent,size=(width, height))
		
		# get the restriction enzymes
		self.enzymeClass = enzymeClass

		# box to keep all the GUI
		hbox      = wx.BoxSizer(wx.HORIZONTAL)
		# style for font
		font      = wx.Font(9, wx.DEFAULT, wx.NORMAL, wx.BOLD)			
		
		# find the cut positions
		self.cutpositions = self.findCutPositons(Enzymes, genbank.gb.gbfile["dna"])
		# we already know the cut positions! They are now always there
		
		
		
		# calculate the fragments, new and old style
		self.fragments, self.fragmentsObj     = self.calculateFragments(self.cutpositions, genbank.gb.gbfile['locus']['topology'])
		
		
		# box for drawing the gel
		self.drawingSpace = wx.ScrolledWindow(self, -1,size=(280,height-50))
		drawingbox        = wx.BoxSizer(wx.VERTICAL)
		
		self.txtAg        = wx.StaticText(self, -1, 'agarose gel preview')
		self.txtAg.SetFont(font)
		drawingbox.Add(self.txtAg)
		drawingbox.Add(self.drawingSpace)
		
		
		# add controls to control space
		controlbox       = wx.BoxSizer(wx.VERTICAL)

		self.txt1        = wx.StaticText(self, -1, 'Choose ladder to display')
		self.txt1.SetFont(font)
		
		# load the avaliable ladders from source
		self.ladders = self.loadLadders()
		# each ladder is selectable, so add it to the selector:
		selectableLadders = []
		for l in self.ladders:
			# only add the name
			selectableLadders.append(l[0])
		# initialise ladder selector
        	self.LadderSelect = wx.ComboBox(self,
				size=wx.DefaultSize,
				choices=selectableLadders)
		self.LadderSelect.Bind(wx.EVT_COMBOBOX, self.onSelect)
		# select the first ladder:
		self.LadderSelect.SetStringSelection(selectableLadders[0])
		
		# button to add a new custom ladder
		self.customButton   = wx.Button(self,-1, "Add custom ladder")	# add custom ladder
		self.customButton.Bind(wx.EVT_BUTTON,self.addcustomLadder)
		
		# add the controls
		controlbox.Add(self.txt1)		# Choose ladder to display
		controlbox.Add(self.LadderSelect)	# selector combobox
		controlbox.AddSpacer(20)		# spacer
		controlbox.Add(self.customButton)	# button "Add custom ladder"
		

		
		# grid settings
		hbox       = wx.BoxSizer(wx.HORIZONTAL)
		gridsizer  = wx.FlexGridSizer(rows=2, cols=2, vgap=3, hgap=10)
		gridsizer.AddGrowableCol(0)					# make cols growable
		hbox.Add(gridsizer, 1, wx.EXPAND|wx.ALL, 15)
		
		# add the elements to the grid for flexible display
		gridsizer.Add(controlbox)      				# row 1, col 1
		gridsizer.Add(drawingbox)      				# row 1, col 1
		
		#set sizer
		self.SetSizer(hbox)
		
		# draw the gel
		self.drawingSpace.Bind(wx.EVT_PAINT, lambda event: self.drawLines(event, self.LadderSelect.GetStringSelection()))
		
		# catch any mouse movement and click in the gel
		self.ClickHighlight = None 				# at first nothing is highlighted
		self.highlight      = None				# not even shortly based on movement 
		self.drawingSpace.Bind(wx.EVT_MOTION, self.OnMotion)
		self.drawingSpace.Bind(wx.EVT_LEFT_UP, self.OnClick)
		
		# catch scrolling
		self.drawingSpace.Bind(wx.EVT_MOUSEWHEEL, self.OnScroll)
		self.drawingSpace.Bind(wx.EVT_SCROLLWIN, self.OnDragScroll)
		self.zoom = 0
		self.drawingSpace.init_pos_x = 0
		self.drawingSpace.init_pos_y = 0
		self.scroll_y                = 0

	

	# respond to a change in ladder selection
	def onSelect(self, event):
		self.drawLines(event, self.LadderSelect.GetStringSelection())
	
	# paints the agarose gel
	def drawLines(self, event, wantedLadder=None, MoveHighlight=None):
		


		width, height = self.drawingSpace.GetClientSize()		# height and width of the parent
		lines = sorted(self.fragments)					# sort the fragments if not ordered yet
		
		# zoom the image
		zoom = 1 + (0.1 * self.zoom)
		self.drawingSpace.SetVirtualSize((width,height*zoom))
		self.drawingSpace.SetScrollRate(0,1)
		self.drawingSpace.Scroll(100, self.scroll_y)
		#print "scroll: ",self.scroll_y
		
		#print "zoom: ", zoom
		#print "old: ", height
		width, height = self.drawingSpace.GetVirtualSize()
		#print "new: ", height
		
		self.dc = wx.PaintDC(self.drawingSpace)				# initialise drawing space
		self.dc.Clear()							# clear
		self.gcdc = wx.GCDC(self.dc) 					# gcdc for nicer images
		wx.EmptyBitmap(self.drawingSpace.ClientSize[0], self.drawingSpace.ClientSize[1])

				
		laneSpacer = 12				# space between lanes
		
		#####################################
		# draw the ladder
		ladder = self.ladder(wantedLadder)
		ladder.pop(0) 				# remove name from ladder
		
		if len(lines) and len(lines) > 1:
			# we leave 10 pixel up and below the scale for purposes of aestetic
			positionSmall = height - 50
			positionLarge = 30
			
			# get the largest and smalles dna fragment
			largestDNA  = len(lines[0][0])
			smallestDNA = len(lines[0][0])
			for fragement in lines:
				l = len(fragement[0])
				if l < smallestDNA:
					smallestDNA = l
				if l > largestDNA:
					largestDNA = l
			# decide if ladder or dna is smaller or bigger, to show really every part of the gel
			if ladder[0] < smallestDNA:
				small = ladder[0]
			else:
				small = smallestDNA

			if ladder[(len(ladder)-1)] > largestDNA:
				large = ladder[(len(ladder)-1)]	
			else:
				large = largestDNA

		
			# we need m and n of the linear function y = m*log(x)+n
			m = (positionSmall - positionLarge)/(math.log10(small)-math.log10(large))
			n = positionLarge - m * math.log10(large)
			
			
			self.gcdc.SetPen(wx.Pen('#595959',3))		# pen for drawing the ladder
			for l in ladder:
				text = "%d bp" % l
				position = m * math.log10(l) + n							# position of the line
				unimport, position = self.drawingSpace.CalcScrolledPosition(0,position)
				self.gcdc.DrawLine(laneSpacer/2 + 55, position, (width-laneSpacer)/2 , position)	# actual line
				self.gcdc.DrawText(text, 10,position-7) 



			# draw the digested plasmid using the objects
			self.gcdc.SetPen(wx.Pen('#000000',3))
			for i, fragment in enumerate(self.fragmentsObj):
				length   = fragment.length
				
		
				text     = "%d bp" % length
				position = m * math.log10(length) + n
				unimport, position = self.drawingSpace.CalcScrolledPosition(0,position) # scrolling only
				
				# check if something needs highlighting
				if (MoveHighlight != None and MoveHighlight.length == length) or (self.ClickHighlight != None and self.ClickHighlight.length == length):
					self.gcdc.SetPen(wx.Pen('#007CE8',5))
					fragment.setPosition((width+laneSpacer)/2 + 55,width-laneSpacer/2,position+2,position-2) 
				else:
					self.gcdc.SetPen(wx.Pen('#000000',3))
					# save position to object +-1 in height because linewidth =3
					fragment.setPosition((width+laneSpacer)/2 + 55,width-laneSpacer/2,position+1,position-1) 
				# save new item:
				self.fragmentsObj[i] = fragment
				
				# draw the line:
				self.gcdc.DrawLine((width+laneSpacer)/2 + 55, position, width-laneSpacer/2, position)
				self.gcdc.DrawText(text, (width+laneSpacer)/2 + laneSpacer/2  ,position-7) 
				

			
				
		return True
	
	
	# two functions to enable vertical scrolling
	def OnScroll(self, event):
		
		# if wheel is spinned
		if event.GetWheelRotation() > 0:
			self.zoom = self.zoom + 1
		elif event.GetWheelRotation() < 0 and self.zoom >= 1:
			self.zoom = self.zoom - 1
		
		self.drawingSpace.prev_pos_y = self.drawingSpace.init_pos_y
	
		x, y = self.drawingSpace.ScreenToClient(wx.GetMousePosition())
	
		self.drawingSpace.init_pos_x, self.drawingSpace.init_pos_y = self.drawingSpace.GetViewStart()
	
		# now we make some calculation:
		width, height = self.drawingSpace.GetClientSize()
		xunit, yunit  = self.drawingSpace.GetScrollPixelsPerUnit()
		
		self.LastScroll_y = self.scroll_y
		addScroll_y       = (float(y)/float(height))*0.1*height  # add pixels per zoom
		#print x,y, height,addScroll_y
		if event.GetWheelRotation() > 0:
			self.scroll_y     = self.LastScroll_y + addScroll_y 
		elif self.LastScroll_y >= addScroll_y:
			self.scroll_y     = self.LastScroll_y - addScroll_y
		else:
			self.scroll_y     = 0

	
		self.drawLines(event, self.LadderSelect.GetStringSelection())
		return True 
	
	def OnDragScroll(self, event):
		self.scroll_y = event.GetPosition()
		self.drawLines(event, self.LadderSelect.GetStringSelection())
		return True
	
	
	# uses the EnzymeSelector class to return cut positions
	def findCutPositons(self, Enzymes, dna):
		
		#cuts = EnzymeSelector(self, id=wx.ID_ANY).findRestrictionSites(Enzymes)	# find cut positions:		
		cuts = []
		for enzyme in Enzymes:
			for i in self.enzymeClass.enzymeObj[enzyme].restrictionSites:
				cuts.append(i)
		
		
		
		self.cutpositions = []							# get the cuts positions

		for c in cuts:
			if c[3] not in self.cutpositions:
				self.cutpositions.append(c[3])					# get the first 5' cut
			if c[4] != None and c[4] not in self.cutpositions:
				self.cutpositions.append(c[4])				# only if there is, add the second cut

		self.cutpositions = sorted(self.cutpositions)				# sort the positions
		return self.cutpositions
	
	
	

	#
	# function calculates fragments of dna from the selected restriction enzymes
	# using the calculated restriction set
	#
	def calculateFragments(self, positions=[], topology="linear"):
		fragments = []				# old list			
		fragementList = []			# list to store objects
		positions  = sorted(positions)		# sort them in cased they are not sorted jet
		
		# fragments as following or as objects (new):
		#	fragement = [[dnastring, start, stop], [...]]
		
		dna = genbank.gb.gbfile["dna"]		# the dna
		length    = len(dna)			# length of the dna
		
		# if there are no positions the gel is just one fragment
		if  positions == []:
			fragments = [[dna, 1, len(dna)+1]]
		else:
		# else we cut the string:
			if topology == "circular":
				# first fragment is las cut first cut:
				f = []
				p1 = dna[positions[len(positions)-1]:len(dna)]
				p2 = dna[:positions[0]]
				part = "%s%s" % (p1,p2)
				f.append(part)				# the dna
				f.append(positions[len(positions)-1])	# start
				f.append(positions[0])			# stop
				newFragement = dnaFragment(part,positions[len(positions)-1],positions[0], len(part)) # using the new class
				
				if len(p1) > 0:				# only fragments bigger than 0
					fragments.append(f)		# add the dna segment
					fragementList.append(newFragement)
			else:
				# add two fragemtns, the last and the first one:
				f1 = []
				p1 = dna[positions[len(positions)-1]:len(dna)]
				f1.append(p1)
				f1.append(0)
				f1.append(positions[0])

				if len(p1) > 0:				# only fragments bigger than 0
					fragments.append(f1)		# add the dna segment
					newFragement = dnaFragment(p1,0,positions[0], len(p1)) # new class
					fragementList.append(newFragement)
				
				f2 = []
				p2 = dna[positions[len(positions)-1]:len(dna)]
				f2.append(p2)
				f2.append(positions[len(positions)-1])
				f2.append(len(dna))
				if len(p2) > 0:				# only fragments bigger than 0
					fragments.append(f2)		# add the dna segment
				
					newFragement = dnaFragment(p2,positions[len(positions)-1],len(dna), len(p2)) # new class
					fragementList.append(newFragement)
			
			# now we have to go from first cut to second, rom second to thirt etc.
			# our stop is the last one
			i = 1
			while i < len(positions):
				f = []
				p = dna[positions[i-1]:positions[i]]	# get the dna fragment
				f.append(p)				# append the dna
				f.append(positions[i-1])		# start
				f.append(positions[i])			# stop
				
				if len(p) > 0:				# only fragments bigger than 0
					fragments.append(f)		# add the dna segment
				
					newFragement = dnaFragment(p,positions[i-1],positions[i], len(p)) # new class
					fragementList.append(newFragement)
				
				i = i + 1 				# counter + 1
				
		return (fragments, fragementList)
	
	
	def ladder(self, wanted):
		# get the selected ladder
		wantedLadder = None
		
		# search avialiable ladders for the wanted
		allLadders = []
		for l in self.loadLadders():
			if l[0] == wanted:
				wantedLadder = l
			
		return wantedLadder
	
	def loadLadders(self):
		Ladders  = []

		# open the dna ladder file
		with open ("resources/dnaLadders.lst", "r") as myfile:
				data=myfile.read().splitlines()
		
		# generate a list of ladders by splitting the line
		for line in data:
			item = line.split(",")
			
			newline = []
			# loop to make integers out of strings if possible
			for i in item:
				try:
					i = int(i)
				except ValueError:
					i = i
				newline.append(i)
			Ladders.append(newline)
		

		
		return Ladders
	
	def reloadLadders(self):
		self.ladders = 	self.loadLadders()
		selectableLadders = []
		for l in self.ladders:
			selectableLadders.append(l[0])
		
		self.LadderSelect.Clear()
		for l in selectableLadders:
			self.LadderSelect.Append(l)
		
	
	def addcustomLadder(self, event):
		customLadderDi =  AddLadderDialog(None,'add custom ladder')
		customLadderDi.Center()
		res = customLadderDi.ShowModal()
		if res == wx.ID_OK:
			ladder = customLadderDi.returnLadder()
			# now save the ladder to the file and reload the ladders:
			strLadder = ''						# string to save to file
			for i in ladder:
				# first item does not need a komma
				if len(strLadder) < 1:
					strLadder = "\n%s" % (i)
				else:
					strLadder = "%s,%s" % (strLadder, i)
			# save to file
			with open("resources/dnaLadders.lst", "a") as myfile:
				myfile.write(strLadder)
			
			# reload ladders for selection
			self.reloadLadders()
			
		customLadderDi.Destroy()	
	
	
	# method to determine if a movement or a click hits an fragment
	def HitTest(self):
		'''Tests whether the mouse is over any fragment'''
		dc    = wx.ClientDC(self.drawingSpace) #get the client dc
		x, y  = self.drawingSpace.ScreenToClient(wx.GetMousePosition()) #get coordinate of mouse event
		pixel = dc.GetPixel(x,y)
		hit   = None
		# if not white
		if pixel != (255, 255, 255, 255):
			# check if its in a line:
			for f in self.fragmentsObj:
				left   = f.position[0]
				right  = f.position[1]
				top    = f.position[2]
				bottom = f.position[3]
				# check if x,y is in square build by left, right, top, bottom
				if x >= left and x <= right and y <= top and y >= bottom:
					hit = f
					break

				
					
		return hit

	
	def OnMotion(self, event):
		x, y = self.drawingSpace.ScreenToClient(wx.GetMousePosition())
		
		oldhighlight = self.highlight
		
		if x >= 0 and y >= 0:
			self.highlight = self.HitTest()
		
		# just redraw if something changed
		if oldhighlight != self.highlight:
			self.drawLines(event, self.LadderSelect.GetStringSelection(), self.highlight)
	
	
	def OnClick(self, event):
		x, y = self.drawingSpace.ScreenToClient(wx.GetMousePosition())
		highlight = self.HitTest()
		
		self.ClickHighlight = highlight
					
		self.drawLines(event, self.LadderSelect.GetStringSelection(), None)
		
		# select dna on click
		if highlight != None:
			self.set_dna_selection((highlight.start, highlight.stop))
		else:
			self.set_dna_selection((-1,0))							    
			
	
	def set_dna_selection(self, selection):
		'''Receives requests for DNA selection and then sends it.'''
		assert type(selection) == tuple, 'Error, dna selection must be a tuple'
		selection = (int(selection[0]+1), int(selection[1]))
		genbank.dna_selection = selection
		
		# update UI
		self.update_globalUI()


	def update_globalUI(self):
		'''
		Method should be modified as to update other panels in response to changes in own panel.
		'''
		
		self.GetParent().GetParent().tab_list[self.GetParent().GetParent().current_tab].update_ownUI()

		pass
		
		


	def update_ownUI(self):
		'''
		Updates to own panel can be made here.
		'''
		pass
		
	



# class to make dna fragments possible
# should be used troughout this file to store and handle dna fragments
# to be implemented in the next days
class dnaFragment():
	def __init__(self, dna, start, stop, length):
		self.dna    = dna
		self.start  = int(start)
		self.stop   = int(stop)
		self.length = int(length)
	
	def setPosition(self, left, right, top, bottom):
		self.position = [left, right, top, bottom]






##############################################################################################################
# Add custom ladder 
#
# simple dialog to add a ladder with name and diffrent fragments
#


class AddLadderDialog(wx.Dialog):
	"""
	Class to make custom ladder
	"""
	def __init__(self, parent, title):
		super(AddLadderDialog, self).__init__(parent=parent,id=wx.ID_ANY, title=title, size=(460, 300)) 		

		#add the panel (containing all the buttons/lists/interactive elements
		self.content = AddLadder(self, id=wx.ID_ANY)	#get the feature edit panel
		
		#add sizer
		sizer = wx.BoxSizer(wx.VERTICAL)
		sizer.Add(item=self.content, proportion=0, flag=wx.EXPAND|wx.ALL)

		#set sizer
		self.SetSizer(sizer)
	
	# this function is called to retrive the ladder
	def returnLadder(self):
		'''
		Get the enzyme selection.
		Used to actually extract info from the dialog.
		'''
		return self.content.returnLadder()	
						
		
		

class AddLadder(DNApyBaseClass):
	def __init__(self, parent, id):
		self.parent = parent
		wx.Panel.__init__(self, parent)
		
		
		
		font         = wx.Font(9, wx.DEFAULT, wx.NORMAL, wx.BOLD)	# style for font
		
		
		# different sizer
		column1     = wx.BoxSizer(wx.VERTICAL)
		column2     = wx.BoxSizer(wx.VERTICAL)
		okcancelColumn     = wx.BoxSizer(wx.HORIZONTAL)
		
		# button to cancel or finish editing:
		self.cancel   = wx.Button(self,wx.ID_CANCEL)
		self.ok       = wx.Button(self,wx.ID_OK)
		
		okcancelColumn.Add(self.cancel)
		okcancelColumn.Add(self.ok)

		
		# elements for column1
		self.txt1     = wx.StaticText(self, -1, 'laddername:')
		self.txt1.SetFont(font)
		
		self.nameField = wx.TextCtrl(self, -1, size=(180,-1),style=wx.TE_PROCESS_ENTER)

		
		self.txt2        = wx.StaticText(self, -1, 'add fragment:')
		self.txt2.SetFont(font)
		self.sizeInput   = wx.TextCtrl(self, -1, size=(180,-1),style=wx.TE_PROCESS_ENTER)
		self.sizeInput.Bind(wx.EVT_TEXT_ENTER, self.selectEnter)
		self.oneadd      = wx.Button(self,-1, "add")		# add enzyme
		self.oneadd.Bind(wx.EVT_BUTTON, self.addOne)

		column1.Add(self.txt1)
		column1.Add(self.nameField)
		column1.AddSpacer(20)
		column1.Add(self.txt2)
		column1.Add(self.sizeInput) 
		column1.Add(self.oneadd) 
		
		
		# elements for column2	
		self.fragmentList	= wx.ListBox(self, size=(200, 240))		
		self.fragmentList.Bind(wx.EVT_LISTBOX_DCLICK, self.remove)

		column2.Add(self.fragmentList)
		

		
		# grid settings
		hbox      = wx.BoxSizer(wx.HORIZONTAL)
		gridsizer = wx.FlexGridSizer(rows=2, cols=2, vgap=3, hgap=10)
		gridsizer.AddGrowableCol(0)					# make cols growable
		gridsizer.AddGrowableCol(1)					# make cols growable
		gridsizer.AddGrowableCol(2)					# make cols growable
		hbox.Add(gridsizer, 1, wx.EXPAND|wx.ALL, 15)
		
		# add the elements to the grid for flexible display
		gridsizer.Add(column1)      				# row 1, col 1
		gridsizer.Add(column2)      				# row 1, col 1
		gridsizer.AddSpacer(1)
		gridsizer.Add(okcancelColumn) 

		
		#set sizer
		self.SetSizer(hbox)
		
		# focus name field:
		self.nameField.SetFocus()


	def selectEnter(self, e):
		self.addOne(e)
		self.sizeInput.SetValue('')


	def addOne(self, event):
		# just integers can be added
		try:
			item     = int(self.sizeInput.GetValue())
			item     = str(item)
			allItems = self.fragmentList.GetItems()
			# just add the item, if its not already in there
			if item not in allItems:
				self.fragmentList.Append(item)
				self.resort2list()
			self.sizeInput.SetValue('')
		except ValueError:
			print "no integer"
			self.sizeInput.SetValue('')
			
		
	# order fragments by size
	def resort2list(self):
		allitems = self.fragmentList.GetItems()
		allIntItems = []
		for i in allitems:
			allIntItems.append(int(i))
		allSorted = sorted(allIntItems)	
		
		# remoformatt all to str again, cause SetItem requires us:
		allStrItems = []
		for i in allSorted:
			allStrItems.append(str(i))
		allSorted = allStrItems
			
		self.fragmentList.SetItems(allSorted)


	# remove selected item from list
	def remove(self, event):
		n = self.fragmentList.GetSelection()
		if n >= 0: # to prevent error when none is selected
			self.fragmentList.Delete(n)
	
	
	# returns the ladder if required
	def returnLadder(self):
		ladder   = []
		name     = self.nameField.GetValue()
		allItems = self.fragmentList.GetItems()
		
		ladder.append(name)
		for i in allItems:
			ladder.append(i)
		
		return ladder
