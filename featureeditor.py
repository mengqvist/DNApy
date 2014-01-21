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
#add search function, narrow list based on keyword
#add menubar with load and save options
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
		self.listview = wx.ListCtrl(self, id=3001, style=wx.LC_REPORT)
		self.listview.InsertColumn(0, "Feature", format=wx.LIST_FORMAT_LEFT, width=200)
		self.listview.InsertColumn(1, "Type", format=wx.LIST_FORMAT_LEFT, width=70)
		self.listview.InsertColumn(2, "Location on DNA", format=wx.LIST_FORMAT_LEFT, width=200)
		self.listview.InsertColumn(3, "Strand", format=wx.LIST_FORMAT_LEFT, width=100)
		self.listview.InsertColumn(4, "Qualifiers", format=wx.LIST_FORMAT_LEFT, width=0)
		
		newfeature = wx.Button(self, 1, 'New Feature')
		deletefeature = wx.Button(self, 2, 'Delete Feature')
#		copyfasta = wx.Button(self, 5, 'Copy FASTA')
#		copydna = wx.Button(self, 4, 'Copy DNA')
#		copytranslation = wx.Button(self, 5, 'Copy Translation')
		
		
		sizer = wx.BoxSizer(wx.HORIZONTAL)
		sizer.Add(newfeature)
		sizer.Add(deletefeature)
#		sizer.Add(copyfasta)
#		sizer.Add(copydna)
#		sizer.Add(copytranslation)

		
		sizer2 = wx.BoxSizer(wx.VERTICAL)
		sizer2.Add(self.listview, 3, wx.EXPAND)
		sizer2.Add(sizer, 0, wx.EXPAND)

		self.SetSizer(sizer2)




class QualifierView(wx.Panel):
	def __init__(self, parent, id):
		wx.Panel.__init__(self, parent)
		self.listview = wx.ListCtrl(self, -1, style=wx.LC_REPORT)
		self.listview.InsertColumn(0, "Qualifier", format=wx.LIST_FORMAT_LEFT, width=70)
		self.listview.InsertColumn(1, "Tag", format=wx.LIST_FORMAT_LEFT, width=wx.LIST_AUTOSIZE)
		
		addqual = wx.Button(self, 7, 'Add Qualifier')
		deletequal = wx.Button(self, 8, 'Remove Qualifier')
		
		sizer = wx.BoxSizer(wx.HORIZONTAL)
		sizer.Add(addqual)
		sizer.Add(deletequal)
		
		sizer2 = wx.BoxSizer(wx.VERTICAL)
		sizer2.Add(self.listview, 3, wx.EXPAND)
		sizer2.Add(sizer, 0, wx.EXPAND)

		self.SetSizer(sizer2)
		

class MyPanel(wx.Panel):
	def __init__(self, parent):
		wx.Panel.__init__(self, parent)
		
		self.current_highlight_color = '#FFFFFF'
		
		#bind feature buttons
		self.Bind(wx.EVT_BUTTON, self.OnNew, id=1)
		self.Bind(wx.EVT_BUTTON, self.OnDelete, id=2)
#		self.Bind(wx.EVT_BUTTON, self.OnCopyFASTA, id=3)
#		self.Bind(wx.EVT_BUTTON, self.OnCopyDNA, id=4)
#		self.Bind(wx.EVT_BUTTON, self.OnCopyTranslation, id=5)

		#bind qualifier buttions
		self.Bind(wx.EVT_BUTTON, self.OnAddQualifier, id=7)
		self.Bind(wx.EVT_BUTTON, self.OnRemoveQualifier, id=8)


		
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


################    Feature edit panel

		
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
		complementtext = wx.StaticText(self.feature_dlg, id=1004, label='')
		self.complementbox = wx.CheckBox(self.feature_dlg, label="Complement?")
		self.complementbox.Bind(wx.EVT_CHECKBOX, self.ComplementCheckboxOnSelect)
		complementtext2 = wx.StaticText(self.feature_dlg, id=wx.ID_ANY, label='')
		
		


		hbox = wx.BoxSizer()
		sizer = wx.GridSizer(rows=4, cols=3, )

		sizer.AddMany([self.featuretext, featuretext1, featuretext2])
		sizer.AddMany([typetext, self.type_combobox, typetext2])
		sizer.AddMany([locationtext, self.location, locationtext2])
		sizer.AddMany([complementtext, self.complementbox, complementtext2])


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

	def get_selection(self):
		#get table content
		listctrlcontent = []
		for row in range(self.feature_list.listview.GetItemCount()):
			featurename, featuretype, location, strand, qualifiers = [self.feature_list.listview.GetItem(row, col).GetText() for col in range(self.feature_list.listview.GetColumnCount())]
			location = '[' + location + ']' #add brackets back
			location =  ast.literal_eval(location) #convert to list
			qualifiers = ast.literal_eval(qualifiers) #convert to list

			if strand == 'leading': strand = False
			elif strand == 'complement': strand = True
			listctrlcontent.append(dict(feature=featurename, key=featuretype, location=location, complement=strand, qualifiers=qualifiers))

		#get the selected feature
		featureposition = self.feature_list.listview.GetFocusedItem()
		feature = listctrlcontent[featureposition]

		return feature

	def ListOnSelect(self, event):
		'''Updates all fields depending on which feature is chosen'''
	
		#get selected feature
		feature = self.get_selection()
		
		#update fields
		self.featuretext.SetLabel(feature['feature'])
		self.type_combobox.SetStringSelection(feature['key']) #update type
		self.location.ChangeValue(str(feature['location'])) #update location
		self.complementbox.SetValue(feature['complement']) #update complement
		
		#update qualifier field
		self.qualifier_list.listview.DeleteAllItems()

		for entry in genbank.gb.gbfile['features']:
			if entry['qualifiers'][0].split('=')[1] == feature['feature']:
				for qualifier in entry['qualifiers']:
					#set content
					col0, col1 = qualifier.split('=')
					col0 = col0[1:]
					self.qualifier_list.listview.Append([col0, col1])
					self.qualifier_list.listview.SetColumnWidth(col=0, width=wx.LIST_AUTOSIZE)
					self.qualifier_list.listview.SetColumnWidth(col=1, width=wx.LIST_AUTOSIZE)					


	def ComplementCheckboxOnSelect(self, event):
		newcomplement = self.complementbox.GetValue()
		feature = self.get_selection()
		genbank.gb.change_feature_complement(feature, newcomplement)
		self.updateUI()
		
	def TypeComboboxOnSelect(self, event):
		newkey = self.type_combobox.GetValue()
		feature = self.get_selection()
		genbank.gb.change_feature_type(feature, newkey)
		self.updateUI()

	def LocationFieldOnText(self, event): # fix  this! maybe use a different event to call it...
		newlocation = self.location.GetLineText(0) #get location
		print(type(newlocation))
		print(newlocation)
		feature = self.get_selection()
		genbank.gb.set_location(feature, newlocation)
		self.updateUI()

	def OnNew(self, event):
		genbank.gb.add_empty_feature()
		self.updateUI()
		
	def OnDelete(self, event):
		'''Delete selected feature'''
		feature = self.get_selection()
		genbank.gb.remove_feature(feature)
		self.updateUI()
	

	def OnCopyFASTA(self, event):
		pass

	def OnCopyDNA(self, event):
		pass
		
	def OnCopyTranslation(self, event):
		pass


	## qualifiers ##
	
	def OnAddQualifier(self, event):
		feature = self.get_selection()
		qualifier = '/label=testing' #change this 
		genbank.gb.add_qualifier(feature, qualifier)
		self.updateUI()
	
	def OnRemoveQualifier(self, event):
		feature = self.get_selection()
		number = self.qualifier_list.listview.GetFocusedItem()
		genbank.gb.remove_qualifier(feature, number)
		self.updateUI()


################# #################### ###############

	### Buttons ###
		

			
			




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



	
