

from base_class import DNApyBaseClass
import wx

class EnzymeSelector(DNApyBaseClass):
	"""
	Class to select restriction enzymes.
	"""
	def __init__(self, parent, id):
		self.parent = parent
		wx.Panel.__init__(self, parent)
		
		
		#set variable to hold the enzyme selection, this should be modified by user interaction
		self.selection = ['BamHI', 'EcoRI']
		
		#make content
		temp_text = wx.StaticText(self, id=wx.ID_ANY)
		textfont = wx.Font(14, wx.DECORATIVE, wx.ITALIC, wx.NORMAL)
		temp_text.SetFont(textfont)
		temp_text.SetLabel('I am an enzyme selection dialog. All my yummy stuff goes in here.')
				
		self.button = wx.Button(self, 1, 'Important button')

		#content 3
		#content 4
		#content n


		#Use sizers add content in the correct arrangement
		sizer1 = wx.BoxSizer(wx.VERTICAL)
		sizer1.Add(item=temp_text, proportion=-1, flag=wx.EXPAND)		
		sizer1.Add(item=self.button, proportion=0)
		#sizer.Add(item=content 3)
		#sizer.Add(item=content 4)
		#sizer.Add(item=content n)

		#set sizer
		self.SetSizer(sizer1)

		
		
		
	#def function 1...
	
	#def function 2...
	
	
	def update_globalUI(self):
		'''
		Method should be modified as to update other panels in response to changes in own panel.
		'''
		pass

	def update_ownUI(self):
		'''
		Updates to own panel can be made here.
		'''
		pass
		
		
		
		
class EnzymeSelectorDialog(wx.Dialog):
	'''A class that puts the Enzyme Selector capabilities in a dialog.'''
	def __init__(self, parent, title):
		super(EnzymeSelectorDialog, self).__init__(parent=parent,id=wx.ID_ANY, title=title, size=(700, 300)) 		

		#add the panel (containing all the buttons/lists/interactive elements
		self.content = EnzymeSelector(self, id=wx.ID_ANY)	#get the feature edit panel
		
		#add sizer
		sizer = wx.BoxSizer(wx.VERTICAL)
		sizer.Add(item=self.content, proportion=-1, flag=wx.EXPAND)

		#set sizer
		self.SetSizer(sizer)	
		

		
	def GetSelection(self):
		'''
		Get the enzyme selection.
		Used to actually extract info from the dialog.
		'''
		return self.content.selection
		

