#!/usr/bin/env python
import ast
import wx
import genbank
import sys, os
from string import *

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



class FeatureView(wx.Panel):
	def __init__(self, parent, id):
		wx.Panel.__init__(self, parent)
		self.listview = wx.ListCtrl(self, id=3001, style=wx.LC_REPORT|wx.LC_SINGLE_SEL)
		self.listview.InsertColumn(0, "Feature", format=wx.LIST_FORMAT_LEFT, width=200)
		self.listview.InsertColumn(1, "Type", format=wx.LIST_FORMAT_LEFT, width=70)
		self.listview.InsertColumn(2, "Location on DNA", format=wx.LIST_FORMAT_LEFT, width=200)
		self.listview.InsertColumn(3, "Strand", format=wx.LIST_FORMAT_LEFT, width=100)
		self.listview.InsertColumn(4, "Qualifiers", format=wx.LIST_FORMAT_LEFT, width=0)
		
		newfeature = wx.Button(self, 1, 'New Feature')
		deletefeature = wx.Button(self, 2, 'Delete Feature')
		moveup = wx.Button(self, 4, 'Move Up')
		movedown = wx.Button(self, 5, 'Move Down')
#		copytranslation = wx.Button(self, 5, 'Copy Translation')
		
		
		sizer = wx.BoxSizer(wx.HORIZONTAL)
		sizer.Add(newfeature)
		sizer.Add(deletefeature)
		sizer.Add(moveup)
		sizer.Add(movedown)
#		sizer.Add(copytranslation)

		
		sizer2 = wx.BoxSizer(wx.VERTICAL)
		sizer2.Add(self.listview, 3, wx.EXPAND)
		sizer2.Add(sizer, 0, wx.EXPAND)

		self.SetSizer(sizer2)


class EditFeatureView(wx.Panel):
	def __init__(self, parent, id):
		wx.Panel.__init__(self, parent)



class QualifierView(wx.Panel):
	def __init__(self, parent, id):
		wx.Panel.__init__(self, parent)
		self.listview = wx.ListCtrl(self, -1, style=wx.LC_REPORT)
		self.listview.InsertColumn(0, "Qualifier", format=wx.LIST_FORMAT_LEFT, width=70)
		self.listview.InsertColumn(1, "Tag", format=wx.LIST_FORMAT_LEFT, width=wx.LIST_AUTOSIZE)
		
		addqual = wx.Button(self, 7, 'Add Qualifier')
		deletequal = wx.Button(self, 8, 'Remove Qualifier')
		qualup = wx.Button(self, 9, 'Move Up')
		qualdown = wx.Button(self, 10, 'Move Down')		

		sizer = wx.BoxSizer(wx.HORIZONTAL)
		sizer.Add(addqual)
		sizer.Add(deletequal)
		sizer.Add(qualup)
		sizer.Add(qualdown)
		
		sizer2 = wx.BoxSizer(wx.VERTICAL)
		sizer2.Add(self.listview, 3, wx.EXPAND)
		sizer2.Add(sizer, 0, wx.EXPAND)

		self.SetSizer(sizer2)

		#bind qualifier buttions
		self.Bind(wx.EVT_BUTTON, self.OnAddQualifier, id=7)
		self.Bind(wx.EVT_BUTTON, self.OnRemoveQualifier, id=8)
		self.Bind(wx.EVT_BUTTON, self.OnMoveQualifierUp, id=9)
		self.Bind(wx.EVT_BUTTON, self.OnMoveQualifierDown, id=10)


	def OnAddQualifier(self, event):
		index = self.get_selection()
		feature = genbank.gb.get_feature(index)
		featurenumber = self.feature_list.listview.GetFirstSelected()
		qualifier = '/label=testing' #change this 
		genbank.gb.add_qualifier(feature, qualifier)
		self.updateUI()

		number = self.qualifier_list.listview.GetItemCount()-1
		self.update_qualifier_selection(featurenumber, number)
	

	def OnRemoveQualifier(self, event):
		index = self.get_selection()
		feature = genbank.gb.get_feature(index)
		featurenumber = self.feature_list.listview.GetFirstSelected()
		number = self.qualifier_list.listview.GetFirstSelected()
		if self.qualifier_list.listview.GetItemCount() != 1: #don't delete last qualifier
			genbank.gb.remove_qualifier(feature, number)
			self.updateUI()

			#set highlight, focus and selection
			if number == self.qualifier_list.listview.GetItemCount():
				number = number-1
			self.update_qualifier_selection(featurenumber, number)


	def OnMoveQualifierUp(self, event):
		index = self.get_selection()
		feature = genbank.gb.get_feature(index)
		featurenumber = self.feature_list.listview.GetFirstSelected()
		number = self.qualifier_list.listview.GetFirstSelected()
		genbank.gb.move_qualifier(feature, number, 'u')
		self.updateUI()
		if number != 0:
			number = number-1
		self.update_qualifier_selection(featurenumber, number)	


	def OnMoveQualifierDown(self, event):
		index = self.get_selection()
		feature = genbank.gb.get_feature(index)
		number = self.qualifier_list.listview.GetFirstSelected()
		featurenumber = self.feature_list.listview.GetFirstSelected()
		genbank.gb.move_qualifier(feature, number, 'd')
		self.updateUI()
		if number != self.qualifier_list.listview.GetItemCount()-1:
			number = number+1
		self.update_qualifier_selection(featurenumber, number)


	def update_qualifier_selection(self, featurenumber, number):
		'''Updates which feature is selected'''
		self.update_feature_selection(featurenumber) #make sure the right feature is selected
		self.qualifier_list.listview.SetItemState(number, wx.LIST_STATE_SELECTED,wx.LIST_STATE_SELECTED) #for the highlight
		self.qualifier_list.listview.Select(number, True) #to select it
		self.qualifier_list.listview.Focus(number) #to focus it



		

class MyPanel(wx.Panel):
	def __init__(self, parent):
		wx.Panel.__init__(self, parent)
		
		self.current_highlight_color = '#FFFFFF'
		
		#bind feature buttons
		self.Bind(wx.EVT_BUTTON, self.OnNew, id=1)
		self.Bind(wx.EVT_BUTTON, self.OnDelete, id=2)
		self.Bind(wx.EVT_BUTTON, self.OnMoveFeatureUp, id=4)
		self.Bind(wx.EVT_BUTTON, self.OnMoveFeatureDown, id=5)
#		self.Bind(wx.EVT_BUTTON, self.OnCopyTranslation, id=5)


		
		splitter1 = wx.SplitterWindow(self, -1, style=wx.SP_3D)
		splitter2 = wx.SplitterWindow(splitter1, -1, style=wx.SP_3D)


##################### basic setup

		#first panel, for showing feature overview
		self.feature_list = FeatureView(splitter1, id=wx.ID_ANY)
		font = wx.Font(pointSize=10, family=wx.FONTFAMILY_DEFAULT, style=wx.FONTSTYLE_NORMAL, weight=wx.FONTWEIGHT_NORMAL, underline=False, faceName='Monospace', encoding=wx.FONTENCODING_DEFAULT)
		self.feature_list.SetFont(font)
		self.feature_list.Bind(wx.EVT_LIST_ITEM_FOCUSED, self.ListOnSelect)


		#second panel, for editing features
		self.feature_dlg = wx.Panel(splitter2, id=wx.ID_ANY)


		#third panel, for showing qualifiers
		self.qualifier_list = QualifierView(splitter2, -1)
		font = wx.Font(pointSize=10, family=wx.FONTFAMILY_DEFAULT, style=wx.FONTSTYLE_NORMAL, weight=wx.FONTWEIGHT_NORMAL, underline=False, faceName='Monospace', encoding=wx.FONTENCODING_DEFAULT)
		self.qualifier_list.SetFont(font)

		splitter1.SplitVertically(self.feature_list, splitter2)
		splitter2.SplitHorizontally(self.feature_dlg, self.qualifier_list)


################    Feature edit panel #######################

		
		#feature label control
#		allfeatures = []
#		for entry in genbank.gb.gbfile['features']:
#			allfeatures.append(entry['qualifiers'][0].split('=')[1])
		self.featuretext = wx.StaticText(self.feature_dlg, id=wx.ID_ANY, label='Feature')
		textfont = wx.Font(18, wx.DECORATIVE, wx.ITALIC, wx.NORMAL)
		self.featuretext.SetFont(textfont)
#		self.features_combobox = wx.ComboBox(self.feature_dlg, id=2001, size=(300, -1), choices=allfeatures, style=wx.CB_READONLY)
		featuretext1 = wx.StaticText(self.feature_dlg, id=wx.ID_ANY, label='')
		featuretext2 = wx.StaticText(self.feature_dlg, id=wx.ID_ANY, label='')

		
		#feature type control
		typetext = wx.StaticText(self.feature_dlg, id=wx.ID_ANY, label='Type:')
		typetext.SetFont(wx.Font(11, wx.DECORATIVE, wx.ITALIC, wx.NORMAL))
		featuretypes = ["modified_base", "variation", "enhancer", "promoter", "-35_signal", "-10_signal", "CAAT_signal", "TATA_signal", "RBS", "5'UTR", "CDS", "gene", "exon", "intron", "3'UTR", "terminator", "polyA_site", "rep_origin", "primer_bind", "protein_bind", "misc_binding", "mRNA", "prim_transcript", "precursor_RNA", "5'clip", "3'clip", "polyA_signal", "GC_signal", "attenuator", "misc_signal", "sig_peptide", "transit_peptide", "mat_peptide", "STS", "unsure", "conflict", "misc_difference", "old_sequence", "LTR", "repeat_region", "repeat_unit", "satellite", "mRNA", "rRNA", "tRNA", "scRNA", "snRNA", "snoRNA", "misc_RNA", "source", "misc_feature", "misc_binding", "misc_recomb", "misc_structure", "iDNA", "stem_loop", "D-loop", "C_region", "D_segment", "J_segment", "N_region", "S_region", "V_region", "V_segment"]

		#can I create submenus with these???		
#		Genes: "promoter", "CDS", "exon", "intron", "gene", "5'UTR", "3'UTR", "polyA_site", "mRNA", "prim_transcript", "precursor_RNA", "5'clip", "3'clip"] 
#		Signals: "rep_origin", "promoter", "enhancer", "polyA_site", "polyA_signal", "terminator", "CAAT_signal", "TATA_signal", "-35_signal", "-10_signal", "GC_signal", "RBS", "attenuator", "misc_signal", "sig_peptide", "transit_peptide", "mat_peptide",
#		Binding: "primer_bind", "protein_bind", "misc_binding",
#		Variation: "variation", "STS", "unsure", "conflict", "modified_base", "misc_difference", "old_sequence",
#		Repeats: "LTR", "repeat_region", "repeat_unit", "satellite", 
#		RNA: "mRNA", "rRNA", "tRNA", "scRNA", "snRNA", "snoRNA", "misc_RNA",
#		Misc.: "source", "misc_feature", "misc_binding", "misc_recomb", "misc_structure", "iDNA", "stem_loop", "D-loop",
#		Ig: "C_region", "D_segment", "J_segment", "N_region", "S_region", "V_region", "V_segment"		

		self.type_combobox = wx.ComboBox(self.feature_dlg, id=1002, size=(150, -1), choices=featuretypes, style=wx.CB_READONLY)
		self.type_combobox.Bind(wx.EVT_COMBOBOX, self.TypeComboboxOnSelect)


		typetext2 = wx.StaticText(self.feature_dlg, id=wx.ID_ANY, label='')

		#location
		locationtext = wx.StaticText(self.feature_dlg, id=wx.ID_ANY, label='Location:')
		locationtext.SetFont(wx.Font(11, wx.DECORATIVE, wx.ITALIC, wx.NORMAL))
		self.location = wx.TextCtrl(self.feature_dlg, id=1003, size=(300,-1))
		self.location.Bind(wx.EVT_TEXT, self.LocationFieldOnText)
		locationtext2 = wx.StaticText(self.feature_dlg, id=wx.ID_ANY, label='')
		
		#complement or not
		self.complementbox = wx.CheckBox(self.feature_dlg, label="Complement?")
		self.complementbox.Bind(wx.EVT_CHECKBOX, self.ComplementCheckboxOnSelect)	

		#Join or not
		self.joinbox = wx.CheckBox(self.feature_dlg, label="Join?")
		self.joinbox.Bind(wx.EVT_CHECKBOX, self.JoinCheckboxOnSelect)

		#Order or not
		self.orderbox = wx.CheckBox(self.feature_dlg, label="Order?")
		self.orderbox.Bind(wx.EVT_CHECKBOX, self.OrderCheckboxOnSelect)
		#need to add logic that only makes join and order available if there are more than one location for a feature...

		hbox = wx.BoxSizer()
		sizer = wx.GridSizer(rows=4, cols=3, )

		sizer.AddMany([self.featuretext, featuretext1, featuretext2])
		sizer.AddMany([typetext, self.type_combobox, typetext2])
		sizer.AddMany([locationtext, self.location, locationtext2])
		sizer.AddMany([self.complementbox, self.joinbox, self.orderbox])


		hbox.Add(sizer, wx.EXPAND, wx.EXPAND)
		self.feature_dlg.SetSizer(hbox)
		
		
		globsizer = hbox = wx.BoxSizer()
		globsizer.Add(splitter1, -1, wx.EXPAND)
		self.SetSizer(globsizer)


		self.Centre()


################ Feature list ctrl methods
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
		self.feature_list.listview.DeleteAllItems()
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
			col4 = entry['qualifiers']
			
			self.feature_list.listview.Append([col0, col1, col2, col3, col4])	
		
			#coloring
			self.get_feature_color(entry)
			color = self.current_highlight_color
			item = n
			self.feature_list.listview.SetItemBackgroundColour(item, color)	
			n += 1
		#self.autosize()	

	def get_selection(self): #change this to actually accessing the file
		#get table content


#		listctrlcontent = []
#		for row in range(self.feature_list.listview.GetItemCount()):
#			featurename, featuretype, location, strand, qualifiers = [self.feature_list.listview.GetItem(row, col).GetText() for col in range(self.feature_list.listview.GetColumnCount())]
#			location = '[' + location + ']' #add brackets back
#			location =  ast.literal_eval(location) #convert to list
#			qualifiers = ast.literal_eval(qualifiers) #convert to list
#
#			if strand == 'leading': strand = False
#			elif strand == 'complement': strand = True
#			listctrlcontent.append(dict(feature=featurename, key=featuretype, location=location, complement=strand, qualifiers=qualifiers))

		#get the selected feature
		featureposition = self.feature_list.listview.GetFocusedItem()
#		feature = genbank.gb.get_feature(featureposition)
		return featureposition


	def ListOnSelect(self, event):
		'''Updates all fields depending on which feature is chosen'''
	
		#get selected feature
		index = self.get_selection()
		
		#update fields
		self.featuretext.SetLabel(genbank.gb.get_feature_label(index))
		self.type_combobox.SetStringSelection(genbank.gb.get_feature_type(index)) #update type
		self.location.ChangeValue(str(genbank.gb.get_feature_location(index))) #update location
		self.complementbox.SetValue(genbank.gb.get_feature_complement(index)) #update complement
		self.joinbox.SetValue(genbank.gb.get_feature_join(index)) #update join
		self.orderbox.SetValue(genbank.gb.get_feature_order(index)) #update order
		
		#update qualifier field
		self.qualifier_list.listview.DeleteAllItems()

#		for entry in genbank.gb.gbfile['features']:
#			if entry['qualifiers'][0].split('=')[1] == feature['feature']:

		index = self.get_selection()
		feature = genbank.gb.get_feature(index)
		for qualifier in feature['qualifiers']:
			#set content
			col0, col1 = qualifier.split('=')
			col0 = col0[1:]
			self.qualifier_list.listview.Append([col0, col1])
			self.qualifier_list.listview.SetColumnWidth(col=0, width=wx.LIST_AUTOSIZE)
			self.qualifier_list.listview.SetColumnWidth(col=1, width=wx.LIST_AUTOSIZE)					


	def ComplementCheckboxOnSelect(self, event):
		'''Toggle whether the feature is on the complement strand or not'''
		newcomplement = self.complementbox.GetValue()
		index = self.get_selection()
		feature = genbank.gb.get_feature(index)
		genbank.gb.set_feature_complement(feature, newcomplement)
		self.updateUI()

	def JoinCheckboxOnSelect(self, event):
		'''Toggle whether a feature with multiple locations should be joined or not'''
		newjoin = self.joinbox.GetValue()
		index = self.get_selection()
		feature = genbank.gb.get_feature(index)
		genbank.gb.set_feature_join(feature, newjoin)
		self.updateUI()

	def OrderCheckboxOnSelect(self, event):
		'''Toggle whether a feature with ultiple locations should be indicated as being in a certain order or not'''
		neworder = self.orderbox.GetValue()
		index = self.get_selection()
		feature = genbank.gb.get_feature(index)
		genbank.gb.set_feature_order(feature, neworder)
		self.updateUI()
		
	def TypeComboboxOnSelect(self, event):
		newkey = self.type_combobox.GetValue()
		index = self.get_selection()
		feature = genbank.gb.get_feature(index)
		genbank.gb.set_feature_type(feature, newkey)
		self.updateUI()

	def LocationFieldOnText(self, event): # fix  this! maybe use a different event to call it...
		newlocation = self.location.GetLineText(0) #get location
		print(type(newlocation))
		print(newlocation)
		index = self.get_selection()
		feature = genbank.gb.get_feature(index)
		genbank.gb.set_feature_location(feature, newlocation)
		self.updateUI()

	def OnNew(self, event):
		'''Make new feature'''
		#make feature and update interface
		genbank.gb.add_feature() #add arguments here!!!!!!!!!
		self.updateUI()

		number = self.feature_list.listview.GetItemCount()-1
		self.update_feature_selection(number)
		

	def OnDelete(self, event):
		'''Delete selected feature'''
		#identify feature, remove it and update interface
		index = self.get_selection()
		feature = genbank.gb.get_feature(index)
		number = self.feature_list.listview.GetFirstSelected()
		genbank.gb.remove_feature(feature)
		self.updateUI()
	
		#set highlight, focus and selection
		if number == self.feature_list.listview.GetItemCount():
			number = number-1
		self.update_feature_selection(number)


	def OnMoveFeatureUp(self, event):
		'''Move feature up one step'''
		index = self.get_selection()
		feature = genbank.gb.get_feature(index)
		number = self.feature_list.listview.GetFirstSelected()
		genbank.gb.move_feature(feature, 'u')	
		self.updateUI()
		if number != 0:
			number = number-1
		self.update_feature_selection(number)
		

	def OnMoveFeatureDown(self, event):
		'''Move feature up down step'''
		index = self.get_selection()
		feature = genbank.gb.get_feature(index)
		number = self.feature_list.listview.GetFirstSelected()
		genbank.gb.move_feature(feature, 'd')
		self.updateUI()
	
		if number != self.feature_list.listview.GetItemCount()-1:
			number = number+1
		self.update_feature_selection(number)


	def update_feature_selection(self, number):
		'''Updates which feature is selected'''
		self.feature_list.listview.SetItemState(number, wx.LIST_STATE_SELECTED,wx.LIST_STATE_SELECTED) #for the highlight
		self.feature_list.listview.Select(number, True) #to select it
		self.feature_list.listview.Focus(number) #to focus it


	def OnCopyFASTA(self, event):
		pass

	def OnCopyDNA(self, event):
		pass
		
	def OnCopyTranslation(self, event):
		pass

			




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



	
