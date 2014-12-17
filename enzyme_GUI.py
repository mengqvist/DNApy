

from base_class import DNApyBaseClass
import wx
import re
import string

class EnzymeSelector(DNApyBaseClass):
	"""
	Class to select restriction enzymes.
	"""
	def __init__(self, parent, id):
		self.parent = parent
		wx.Panel.__init__(self, parent)
		
						
		# enzyme selector GUI
		self.lb       = wx.ListBox(self,
		                size=(180, 300))
		self.oneadd   = wx.Button(self,-1, "add")
		#self.multiadd = wx.Button(self, -1,"add all",pos=(10, 345))
		self.remove   = wx.Button(self,-1, "remove")


		# restriction enzymes of first list added to second:
		self.oneadd.Bind(wx.EVT_BUTTON, self.addOne)
		self.lb.Bind(wx.EVT_LISTBOX_DCLICK, self.addOne)

		# second list
		self.lb2 = wx.ListBox(self,
		        	 size=(180, 300))

		# remove enzymes from second by doubleclick
		self.remove.Bind(wx.EVT_BUTTON, self.removeOne)
		self.lb2.Bind(wx.EVT_LISTBOX_DCLICK, self.removeOne)
	 
		self.lb.Bind(wx.EVT_LISTBOX, self.onSelect)

		
		# button to cancel or finish editing:
		self.cancel   = wx.Button(self,wx.ID_CANCEL	)
		self.ok       = wx.Button(self,wx.ID_OK)
		



		#Use sizers add content in the correct arrangement
		sizer1      = wx.BoxSizer(wx.VERTICAL)
		sizerClose  = wx.BoxSizer(wx.VERTICAL)
		
		# add and remove
		sizer1.Add(item=self.oneadd)
		sizer1.Add(item=self.remove)

		sizerClose.Add(item=self.ok)		
		sizerClose.Add(item=self.cancel)
		

		hbox      = wx.BoxSizer(wx.HORIZONTAL)
		gridsizer = wx.FlexGridSizer(rows=2, cols=3, vgap=3)
		hbox.Add(gridsizer, proportion=1, flag=wx.ALL|wx.EXPAND, border=15)
		#gridsizer.SetFlexibleDirection( wx.BOTH )
		
		
		# add the elements to the grid for flexible display
		gridsizer.Add(item=self.lb)      # row 1, col 1
		gridsizer.Add(item=sizer1)       # row 1, col 2
		gridsizer.Add(item=self.lb2)     # row 1, col 3
		gridsizer.AddSpacer(1) 			 # row 2, col 1
		gridsizer.AddSpacer(1) 			 # row 2, col 2
		gridsizer.Add(item=sizerClose)   # row 2, col 3
		
		#set sizer
		self.SetSizer(hbox)

		# load the enzymes from emboss
		self.loadEnzymes()
	
	
	# on select nothing happens
	def onSelect(self, event):
        	return False


	# adds one item on buttonclick
	def addOne(self, event):
		#print self.lb.GetSelection()
		item = self.lb.GetStringSelection()
		allItems=self.lb2.GetItems()
		if item not in allItems:
			self.lb2.Append(item)
			self.resort2list()

	# remove selected item from list
	def removeOne(self, event):
		n = self.lb2.GetSelection()
		if n >= 0: # to prevent error when none is selected
			self.lb2.Delete(n)
		

	# adds the given restriction enzyme to the first list
	def AddRestrictionEnzyme(self, item, i):
		if i == 1:
			self.lb.Append(item)
		if i == 2:
			self.lb2.Append(item)

	# can reorder the second list
	# is called by addOne(self)
	def resort2list(self):
		allitems = self.lb2.GetItems() 
		allSorted = sorted(allitems)		
		self.lb2.SetItems(allSorted)
	
	####################################################
	# this function is going to transform the file 
	# "emboss_e.txt into an object of restriction
	# enzmyes
	def loadEnzymes(self):
		# first load the enzyme db

		# File strukture for line
		# name = name of enzyme
		# pattern = recognition site
		# len = length of pattern
		# ncuts = number of cuts made by enzyme
		#         Zero represents unknown
		# blunt = true if blunt end cut, false if sticky
		# c1 = First 5' cut
		# c2 = First 3' cut
		# c3 = Second 5' cut
		self.enzymes             = {}
		self.enzymes["names"]    = []
		self.enzymes["regexp"]   = []
		self.enzymes["ncuts"]    = []
		self.enzymes["cut5_1"]   = []
		self.enzymes["cut5_2"]   = []


		# for this we have enzymes in folder /resources
		with open('resources/emboss_e.txt') as f:
		    for line in f:	
		    	# each line is one enzyme, except the header
				if ( line[:1] != "#" ):	
			
					# split the line
					lineparts = re.split(r'\t+', line)	
				
					# list of all possible enzymes for selector
					self.enzymes["names"].append(lineparts[0])

					# new object for each enzmye, makes searching possible
					self.enzymes[lineparts[0]]        = {}
					self.enzymes[lineparts[0]]["enz"] = lineparts[0]
				

					# form the regexp:
					#Code	Meaning			Etymology	Complement	Opposite
					#A	A			Adenosine	T	B
					#T/U	T			Thymidine/Uridine	A	V
					#G	G			Guanine	C	H
					#C	C			Cytidine	G	D
					#K	G or T			Keto	M	M
					#M	A or C			Amino	K	K
					#R	A or G			Purine	Y	Y
					#Y	C or T			Pyrimidine	R	R
					#S	C or G			Strong	S	W
					#W	A or T			Weak	W	S
					#B	C or G or T		not A (B comes after A)	V	A
					#V	A or C or G		not T/U (V comes after U)	B	T/U
					#H	A or C or T		not G (H comes after G)	D	G
					#D	A or G or T		not C (D comes after C)	H	C
					#X/N	G or A or T or C	any	N	.
					#.	not G or A or T or C	.	N
					#-	gap of indeterminate length	
					regpattern = lineparts[1]
					regpattern = string.replace(regpattern,"U","T")
					regpattern = string.replace(regpattern,"K","(G|T)")
					regpattern = string.replace(regpattern,"M","(A|C)")
					regpattern = string.replace(regpattern,"R","(A|G)")
					regpattern = string.replace(regpattern,"Y","(C|T)")
					regpattern = string.replace(regpattern,"S","(C|G)")
					regpattern = string.replace(regpattern,"W","(A|T)")
					regpattern = string.replace(regpattern,"B","[CGT]")
					regpattern = string.replace(regpattern,"V","[ACG]")
					regpattern = string.replace(regpattern,"H","[ACT]")
					regpattern = string.replace(regpattern,"D","[AGT]")
					regpattern = string.replace(regpattern,"N","[AGTC]")
					regpattern = string.replace(regpattern,"X","[AGTC]")
	
					# save the info to the object
					self.enzymes[lineparts[0]]["regexp"] = re.compile(regpattern, re.IGNORECASE) # regexped, compiled and ready to use
					self.enzymes[lineparts[0]]["cut5_1"] = int(lineparts[5])
					self.enzymes[lineparts[0]]["cut5_2"] = int(lineparts[7])				# add the enzyme to the list of possible enzymes:
					self.AddRestrictionEnzyme(lineparts[0],1)
		

	# gets called to retrive the sleected infos
	def getSelection(self):
		return self.lb2.GetItems() 

	
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
	def __init__(self, parent, title, oldSelection):
		super(EnzymeSelectorDialog, self).__init__(parent=parent,id=wx.ID_ANY, title=title, size=(600, 400)) 		

		#add the panel (containing all the buttons/lists/interactive elements
		self.content = EnzymeSelector(self, id=wx.ID_ANY)	#get the feature edit panel
		
		# add the old selection:
		for item in oldSelection:
			self.content.AddRestrictionEnzyme(item,2)
		
		#add sizer
		sizer = wx.BoxSizer(wx.VERTICAL)
		sizer.Add(item=self.content, proportion=0, flag=wx.EXPAND|wx.ALL)

		#set sizer
		self.SetSizer(sizer)
		
		# pass the old selection to the window
		
		# show the window (enables the ok and cancel buttons)
		#result = EnzymeSelector.Show(self.content)	


		

		
	def GetSelection(self):
		'''
		Get the enzyme selection.
		Used to actually extract info from the dialog.
		'''
		return self.content.getSelection()
		

