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
#change how qualifiers are edited


import ast
import wx
from wx.lib.agw import ultimatelistctrl as ULC
#from wx.lib.pubsub import pub

import genbank


import sys, os
import string

import copy

#for asserts and tests
import types
from base_class import DNApyBaseClass


files={}   #list with all configuration files
files['default_dir'] = os.path.abspath(os.path.dirname(sys.argv[0]))+"/"
files['default_dir']=string.replace(files['default_dir'], "\\", "/")
files['default_dir']=string.replace(files['default_dir'], "library.zip", "")
settings=files['default_dir']+"settings"   ##path to the file of the global settings
execfile(settings) #gets all the pre-assigned settings



class QualifierEdit(wx.Dialog):
	'''This class is for making a panel with two fields for choosing a qualifier type and to edit the tag associated with it.'''
	def __init__(self, qualifier, tag):
		wx.Dialog.__init__(self, None, id=wx.ID_ANY, title="Edit Qualifier", size=(420,250))
		self.qualifier = qualifier
		self.tag = tag


		#'qualifier'
		self.qualifier_header = wx.StaticText(self, id=wx.ID_ANY)
		textfont = wx.Font(18, wx.DECORATIVE, wx.ITALIC, wx.NORMAL)
		self.qualifier_header.SetFont(textfont)
		self.qualifier_header.SetLabel('Qualifier')

		choicelist = ['allele', 'altitude', 'anticodon', 'artificial_location', 'bio_material', 'bound_moiety', 'cell_line', 'cell_type', 'chromosome', 'citation', 'clone', 'clone_lib', 'codon_start', 'collected_by', 'collection_date', 'compare', 'country', 'cultivar', 'culture_collection', 'db_xref', 'dev_stage', 'direction', 'EC_number', 'ecotype', 'environmental_sample', 'estimated_length', 'exception', 'experiment', 'focus', 'frequency', 'function', 'gap_type', 'gene', 'gene_synonym', 'germline', 'haplogroup', 'haplotype', 'host', 'identified_by', 'inference', 'isolate', 'isolation_source', 'label', 'lab_host', 'lat_lon', 'linkage_evidence', 'locus_tag', 'macronuclear', 'map', 'mating_type', 'mobile_element_type', 'mod_base', 'mol_type', 'ncRNA_class', 'note', 'number', 'old_locus_tag', 'operon', 'organelle', 'organism', 'partial', 'PCR_conditions', 'PCR_primers', 'phenotype', 'plasmid', 'pop_variant', 'product', 'protein_id', 'proviral', 'pseudo', 'pseudogene', 'rearranged', 'replace', 'ribosomal_slippage', 'rpt_family', 'rpt_type', 'rpt_unit_range', 'rpt_unit_seq', 'satellite', 'segment', 'serotype', 'serovar', 'sex', 'specimen_voucher', 'standard_name', 'strain', 'sub_clone', 'sub_species', 'sub_strain', 'tag_peptide', 'tissue_lib', 'tissue_type', 'transgenic', 'translation', 'transl_except', 'transl_table', 'trans_splicing', 'type_material', 'variety']
		self.qualifier_combobox = wx.ComboBox(self, id=1002, size=(150, -1), choices=choicelist, style=wx.CB_READONLY)

		assert self.qualifier in choicelist, "Not a valid qualifier: %s" % str(self.qualifier)
		self.qualifier_combobox.SetStringSelection(self.qualifier) #set selection as current qualifier

		#spacer
		self.spacer1 = wx.StaticText(self, id=wx.ID_ANY)
		textfont = wx.Font(18, wx.DECORATIVE, wx.ITALIC, wx.NORMAL)
		self.spacer1.SetFont(textfont)
		self.spacer1.SetLabel('')


		#'tag'
		self.tag_header = wx.StaticText(self, id=wx.ID_ANY)
		textfont = wx.Font(18, wx.DECORATIVE, wx.ITALIC, wx.NORMAL)
		self.tag_header.SetFont(textfont)
		self.tag_header.SetLabel('Tag')

		#make text field that displays the qualifier tag
		textfont = wx.Font(12, wx.DECORATIVE, wx.NORMAL, wx.NORMAL)
		self.tag_field = wx.stc.StyledTextCtrl(self, size=(300,100))
		self.tag_field.StyleSetBackground(style=wx.stc.STC_STYLE_DEFAULT, back='#FFFFFF') #set background color of everything that is not text
		self.tag_field.StyleSetBackground(style=0, back='#FFFFFF') #set background color of text
		self.tag_field.StyleSetBackground(style=wx.stc.STC_STYLE_LINENUMBER, back='#FFFFFF') #sets color of left margin
		face = textfont.GetFaceName()
		size = textfont.GetPointSize()
		self.tag_field.StyleSetSpec(wx.stc.STC_STYLE_DEFAULT,"face:%s,size:%d" % (face, size))
		self.tag_field.SetUseHorizontalScrollBar(True)
		self.tag_field.SetUseVerticalScrollBar(True) 

		self.tag_field.SetText(self.tag) #put tag in the text field

		#ok button
		self.button_ok = wx.Button(self, 1, 'OK')
		self.button_ok.Bind(wx.EVT_BUTTON, self.OnOK, id=1)

		#cancel button
		self.button_cancel = wx.Button(self, 2, 'Cancel')
		self.button_cancel.Bind(wx.EVT_BUTTON, self.OnCancel, id=2)	

		#sizers
		buttonsizer = wx.BoxSizer(wx.HORIZONTAL)
		buttonsizer.Add(item=self.button_ok)
		buttonsizer.Add(item=self.button_cancel)

		globsizer = wx.BoxSizer(wx.VERTICAL)
		globsizer.Add(item=self.qualifier_header, flag=wx.LEFT, border=10)
		globsizer.Add(item=self.qualifier_combobox, flag=wx.LEFT, border=10)
		globsizer.Add(item=self.spacer1, flag=wx.LEFT, border=10)
		globsizer.Add(item=self.tag_header, flag=wx.LEFT, border=10)
		globsizer.Add(item=self.tag_field, flag=wx.LEFT, border=10)
		globsizer.Add(item=buttonsizer, flag=wx.LEFT, border=10)

		self.SetSizer(globsizer)
		self.Center()
		
	def OnOK(self, evt):
		'Return wx.ID_OK if ok button was pressed'
		self.EndModal(wx.ID_OK)

	def OnCancel(self, evt):
		'Return wx.ID_CANCEL if cancel button was pressed'
		self.EndModal(wx.ID_CANCEL)
		







class FeatureEdit(DNApyBaseClass):
	'''This class makes a panel with a field for editing feature properties and a list for displaying and editing qualifiers'''
	def __init__(self, parent, id):
		self.parent = parent
		wx.Panel.__init__(self, parent)

		##
		# first panel, for editing feature
		##
		self.feature_dlg = wx.Panel(self, id=wx.ID_ANY)

		#feature label control
		featuretext = wx.StaticText(self.feature_dlg, id=wx.ID_ANY, label='Feature:')
		featuretext.SetFont(wx.Font(11, wx.DECORATIVE, wx.ITALIC, wx.NORMAL))
		self.featuretext = wx.TextCtrl(self.feature_dlg, id=wx.ID_ANY, size=(200,-1))
		self.featuretext.Bind(wx.EVT_TEXT, self.FeatureFieldOnText)

		
		#feature type control
		typetext = wx.StaticText(self.feature_dlg, id=wx.ID_ANY, label='Type:')
		typetext.SetFont(wx.Font(11, wx.DECORATIVE, wx.ITALIC, wx.NORMAL))
		featuretypes = ["modified_base", "variation", "enhancer", "promoter", "-35_signal", "-10_signal", "CAAT_signal", "TATA_signal", "RBS", "5'UTR", "CDS", "gene", "exon", "intron", "3'UTR", "terminator", "polyA_site", "rep_origin", "primer_bind", "protein_bind", "misc_binding", "prim_transcript", "precursor_RNA", "5'clip", "3'clip", "polyA_signal", "GC_signal", "attenuator", "misc_signal", "sig_peptide", "transit_peptide", "mat_peptide", "STS", "unsure", "conflict", "misc_difference", "old_sequence", "LTR", "repeat_region", "repeat_unit", "satellite", "mRNA", "rRNA", "tRNA", "ncRNA","scRNA", "snRNA", "snoRNA", "misc_RNA", "source", "misc_feature", "misc_binding", "misc_recomb", "misc_structure", "iDNA", "stem_loop", "D-loop", "C_region", "D_segment", "J_segment", "N_region", "S_region", "V_region", "V_segment"]
		self.type_combobox = wx.ComboBox(self.feature_dlg, id=1002, size=(200, -1), choices=featuretypes, style=wx.CB_READONLY)
		self.type_combobox.Bind(wx.EVT_COMBOBOX, self.TypeComboboxOnSelect)

		#can I create submenus with these???		
#		Genes: "promoter", "CDS", "exon", "intron", "gene", "5'UTR", "3'UTR", "polyA_site", "mRNA", "prim_transcript", "precursor_RNA", "5'clip", "3'clip"] 
#		Signals: "rep_origin", "promoter", "enhancer", "polyA_site", "polyA_signal", "terminator", "CAAT_signal", "TATA_signal", "-35_signal", "-10_signal", "GC_signal", "RBS", "attenuator", "misc_signal", "sig_peptide", "transit_peptide", "mat_peptide",
#		Binding: "primer_bind", "protein_bind", "misc_binding",
#		Variation: "variation", "STS", "unsure", "conflict", "modified_base", "misc_difference", "old_sequence",
#		Repeats: "LTR", "repeat_region", "repeat_unit", "satellite", 
#		RNA: "mRNA", "rRNA", "tRNA","ncRNA", "scRNA", "snRNA", "snoRNA", "misc_RNA",
#		Misc.: "source", "misc_feature", "misc_binding", "misc_recomb", "misc_structure", "iDNA", "stem_loop", "D-loop",
#		Ig: "C_region", "D_segment", "J_segment", "N_region", "S_region", "V_region", "V_segment"		


		#location control
		locationtext = wx.StaticText(self.feature_dlg, id=wx.ID_ANY, label='Location:')
		locationtext.SetFont(wx.Font(11, wx.DECORATIVE, wx.ITALIC, wx.NORMAL))
		self.location = wx.TextCtrl(self.feature_dlg, id=1003, size=(200,-1), style=wx.TE_RICH)
		self.location.Bind(wx.EVT_TEXT, self.LocationFieldOnText)
		
		#complement or not
		self.complementbox = wx.CheckBox(self.feature_dlg, id=3001, label="Complement (DNA on complement strand)")
		self.complementbox.Bind(wx.EVT_CHECKBOX, self.ComplementCheckboxOnSelect)	

		#Join or not
		self.joinbox = wx.CheckBox(self.feature_dlg, id=3002, label="Join (location ranges should be joined)")
		self.joinbox.Bind(wx.EVT_CHECKBOX, self.JoinCheckboxOnSelect)

		#Order or not
		self.orderbox = wx.CheckBox(self.feature_dlg, id=3003, label="Order (location ranges indicate order)")
		self.orderbox.Bind(wx.EVT_CHECKBOX, self.OrderCheckboxOnSelect)
		#need to add logic that only makes join and order available if there are more than one location for a feature...


		#sizers
		label_sizer = wx.BoxSizer(wx.HORIZONTAL)
		label_sizer.Add(item=featuretext, flag=wx.ALL, border=5)
		label_sizer.Add(item=self.featuretext, flag=wx.LEFT, border=80-featuretext.GetSize()[0])

		type_sizer = wx.BoxSizer(wx.HORIZONTAL)
		type_sizer.Add(item=typetext, flag=wx.ALL, border=5)
		type_sizer.Add(item=self.type_combobox, flag=wx.LEFT, border=80-typetext.GetSize()[0])

		location_sizer = wx.BoxSizer(wx.HORIZONTAL)		
		location_sizer.Add(item=locationtext, flag=wx.ALL, border=5)
		location_sizer.Add(item=self.location, flag=wx.LEFT, border=80-locationtext.GetSize()[0])

		complement_sizer = wx.BoxSizer(wx.VERTICAL)
		complement_sizer.Add(item=self.complementbox)
		complement_sizer.Add(item=self.joinbox)
		complement_sizer.Add(item=self.orderbox)

		main_sizer = wx.BoxSizer(wx.VERTICAL)
		main_sizer.Add(item=label_sizer)
		main_sizer.Add(item=type_sizer)
		main_sizer.Add(item=location_sizer)
		main_sizer.Add(item=(0,40))
		main_sizer.Add(item=complement_sizer)
		self.feature_dlg.SetSizer(main_sizer)

		##
		#second panel, for qualifiers
		##
		self.qualifier_list = ULC.UltimateListCtrl(self, -1, agwStyle=ULC.ULC_REPORT|ULC.ULC_HAS_VARIABLE_ROW_HEIGHT|ULC.ULC_SINGLE_SEL)
		self.qualifier_list.InsertColumn(0, "Qualifier", format=wx.LIST_FORMAT_LEFT, width=120)
		self.qualifier_list.InsertColumn(1, "Tag", format=wx.LIST_FORMAT_LEFT, width=250)
#		self.qualifier_list.Bind(wx.EVT_LISTBOX_DCLICK, self.OnEditQualifier(None))

		imageFile = files['default_dir']+"/icon/new_small.png"
		image1 = wx.Image(imageFile, wx.BITMAP_TYPE_ANY).ConvertToBitmap()
		addqual = wx.BitmapButton(self, id=7, bitmap=image1, size = (image1.GetWidth()+15, image1.GetHeight()+15), name = "new")
		self.Bind(wx.EVT_BUTTON, self.OnAddQualifier, id=7)

		imageFile = files['default_dir']+"/icon/remove_small.png"
		image1 = wx.Image(imageFile, wx.BITMAP_TYPE_ANY).ConvertToBitmap()
		deletequal = wx.BitmapButton(self, id=8, bitmap=image1, size = (image1.GetWidth()+15, image1.GetHeight()+15), name = "remove")
		self.Bind(wx.EVT_BUTTON, self.OnRemoveQualifier, id=8)

		imageFile = files['default_dir']+"/icon/move_up.png"
		image1 = wx.Image(imageFile, wx.BITMAP_TYPE_ANY).ConvertToBitmap()
		qualup = wx.BitmapButton(self, id=9, bitmap=image1, size = (image1.GetWidth()+15, image1.GetHeight()+15), name = "move up")
		self.Bind(wx.EVT_BUTTON, self.OnMoveQualifierUp, id=9)

		imageFile = files['default_dir']+"/icon/move_down.png"
		image1 = wx.Image(imageFile, wx.BITMAP_TYPE_ANY).ConvertToBitmap()
		qualdown = wx.BitmapButton(self, id=10, bitmap=image1, size = (image1.GetWidth()+15, image1.GetHeight()+15), name = "move down")
		self.Bind(wx.EVT_BUTTON, self.OnMoveQualifierDown, id=10)

		imageFile = files['default_dir']+"/icon/edit.png"
		image1 = wx.Image(imageFile, wx.BITMAP_TYPE_ANY).ConvertToBitmap()
		qualedit = wx.BitmapButton(self, id=11, bitmap=image1, size = (image1.GetWidth()+15, image1.GetHeight()+15), name = "edit")
		self.Bind(wx.EVT_BUTTON, self.OnEditQualifier, id=11)	


		#sizers
		sizer = wx.BoxSizer(wx.VERTICAL)
		sizer.Add(item=addqual)
		sizer.Add(item=deletequal)
		sizer.Add(item=qualup)
		sizer.Add(item=qualdown)
		sizer.Add(item=qualedit)
		
		sizer2 = wx.BoxSizer(wx.HORIZONTAL)
		sizer2.Add(item=self.qualifier_list, proportion=3, flag=wx.EXPAND)
		sizer2.Add(item=sizer, proportion=0, flag=wx.EXPAND)

		main_sizer = wx.BoxSizer(wx.HORIZONTAL)
		main_sizer.Add(item=self.feature_dlg, proportion=-1, flag=wx.EXPAND)

		main_sizer.Add(item=sizer2, proportion=-1, flag=wx.EXPAND)
		self.SetSizer(main_sizer)
		self.Centre()

		self.update_ownUI()







####### Modify methods from base calss to fit current needs #########

	def update_globalUI(self):
		'''
		Method should be modified as to update other panels in response to changes in own panel.
		'''
		pass

	def update_ownUI(self):
		'''Updates all fields depending on which feature is chosen'''
		if len(genbank.gb.get_all_features()) == 0:
			pass
		else:
			index = int(copy.copy(genbank.feature_selection))
	
			#update fields
			self.featuretext.ChangeValue(genbank.gb.get_feature_label(index))
			self.type_combobox.SetStringSelection(genbank.gb.get_feature_type(index)) #update type

			locationstring = ''
			for location in genbank.gb.get_feature_location(index):
				if locationstring != '':
					locationstring += ', '
				locationstring += str(location)
			self.location.ChangeValue(locationstring) #update location

			self.complementbox.SetValue(genbank.gb.get_feature_complement(index)) #update complement
			self.joinbox.SetValue(genbank.gb.get_feature_join(index)) #update join
			self.orderbox.SetValue(genbank.gb.get_feature_order(index)) #update order

			#update qualifier field
			self.update_qualifiers()

		order = self.orderbox.GetValue()
		join = self.joinbox.GetValue()
		assert (join == True and order == True) == False, 'Error, order and join cannot both be True'
		if order == True:
			self.joinbox.Enable(False)
		elif join == True:
			self.orderbox.Enable(False)
		elif join == False and order == False:
			self.joinbox.Enable(True)
			self.orderbox.Enable(True)

#####################################################################




	def OnAddQualifier(self, event):
		'''For adding new qualifier'''
		index = int(copy.copy(genbank.feature_selection))
		feature = genbank.gb.get_feature(index)
		qualifier = 'label'
		tag = 'None'

		#make popup window with qualifier editing capabilities
		dlg = QualifierEdit(qualifier, tag) # creation of a panel in the frame
		output = dlg.ShowModal()
		if output == wx.ID_OK:
			qualifier = str(dlg.qualifier_combobox.GetValue())
			tag = str(dlg.tag_field.GetText())
		
			genbank.gb.add_qualifier(feature, '/%s=%s' % (qualifier, tag))
			self.update_ownUI()

			number = self.qualifier_list.GetItemCount()-1
			self.update_qualifier_selection(index, number)

		dlg.Destroy()
	

	def OnRemoveQualifier(self, event):
		'''For removing existing qualifier'''
		index = int(copy.copy(genbank.feature_selection)) #feature selected
		feature = genbank.gb.get_feature(index)
		number = self.qualifier_list.GetFirstSelected() #qualifier selected
		items = self.qualifier_list.GetItemCount() #number of qualifiers

		if items != 1: #don't delete last qualifier
			genbank.gb.remove_qualifier(feature, number)
			self.update_ownUI()
			self.update_globalUI()

			#set highlight, focus and selection
			if number == items-1: #if last item was deleted, selection is set as -1
				number = number-1
			self.update_qualifier_selection(index, number)


	def OnMoveQualifierUp(self, event):
		'''Move qualifier up one step in the list'''
		index = int(copy.copy(genbank.feature_selection))
		feature = genbank.gb.get_feature(index)
		number = self.qualifier_list.GetFirstSelected()
		genbank.gb.move_qualifier(feature, number, 'u')
		self.update_ownUI()
		self.update_globalUI()

		if number != 0:
			number = number-1
		self.update_qualifier_selection(index, number)	


	def OnMoveQualifierDown(self, event):
		'''Move qualifier down one step in the list'''
		index = int(copy.copy(genbank.feature_selection))
		feature = genbank.gb.get_feature(index)
		number = self.qualifier_list.GetFirstSelected()
		genbank.gb.move_qualifier(feature, number, 'd')
		self.update_ownUI()
		self.update_globalUI()

		if number != self.qualifier_list.GetItemCount():
			number = number+1
		self.update_qualifier_selection(index, number)


	def OnEditQualifier(self, event):
		'''Make popup window where qualifier and tag can be edited'''
		index = int(copy.copy(genbank.feature_selection))
		feature = genbank.gb.get_feature(index)
		number = self.qualifier_list.GetFirstSelected()
		qualifier, tag = genbank.gb.get_qualifier(index, number) #get actual info for that qualifier

		#make popup window with qualifier editing capabilities
		dlg = QualifierEdit(qualifier, tag) # creation of a panel in the frame
		output = dlg.ShowModal()
		if output == wx.ID_OK:
			qualifier = str(dlg.qualifier_combobox.GetValue())
			tag = str(dlg.tag_field.GetText())
		dlg.Destroy()

		#update file
		genbank.gb.set_qualifier(index, number, qualifier, tag)
		self.update_ownUI()
		self.update_globalUI()
		self.update_qualifier_selection(index, number)


	def update_qualifier_selection(self, index, number):
		'''Updates which feature and qualifier is selected'''
		self.qualifier_list.SetItemState(number, wx.LIST_STATE_SELECTED,wx.LIST_STATE_SELECTED) #for the highlight
		self.qualifier_list.Select(number, True) #to select it
		self.qualifier_list.Focus(number) #to focus it

	
	def ComplementCheckboxOnSelect(self, event):
		'''Toggle whether the feature is on the complement strand or not'''
		newcomplement = self.complementbox.GetValue()
		index = int(copy.copy(genbank.feature_selection))
		feature = genbank.gb.get_feature(index)
		genbank.gb.set_feature_complement(feature, newcomplement)
		self.update_ownUI()
		self.update_globalUI()

		

	def JoinCheckboxOnSelect(self, event):
		'''Toggle whether a feature with multiple locations should be joined or not'''
		newjoin = self.joinbox.GetValue()
		index = int(copy.copy(genbank.feature_selection))
		feature = genbank.gb.get_feature(index)
		genbank.gb.set_feature_join(feature, newjoin)
		self.update_ownUI()
		self.update_globalUI()


	def OrderCheckboxOnSelect(self, event):
		'''Toggle whether a feature with ultiple locations should be indicated as being in a certain order or not'''
		neworder = self.orderbox.GetValue()
		index = int(copy.copy(genbank.feature_selection))
		feature = genbank.gb.get_feature(index)
		genbank.gb.set_feature_order(feature, neworder)
		self.update_ownUI()
		self.update_globalUI()

		
	def TypeComboboxOnSelect(self, event):
		newkey = self.type_combobox.GetValue()
		index = int(copy.copy(genbank.feature_selection))
		feature = genbank.gb.get_feature(index)
		genbank.gb.set_feature_type(feature, newkey)
		self.update_ownUI()
		self.update_globalUI()	
		

	def LocationFieldOnText(self, event):
		'''When the location field is changed. Check if the entry is valid, if yes, then update features'''
		newlocation = str(self.location.GetLineText(0)) #get location as string
		newlocation.replace(' ', '') #get rid of spaces

		#convert to list
		if ',' in newlocation:
			locationlist = newlocation.split(',')
		else:
			locationlist = [newlocation]

		index = int(copy.copy(genbank.feature_selection))
		feature = genbank.gb.get_feature(index)
		
		is_valid = genbank.gb.IsValidLocation(locationlist) #test if location entry is valid
		if is_valid == True:
			self.location.SetForegroundColour(wx.BLACK)
			genbank.gb.set_feature_location(feature, locationlist)
			self.update_globalUI()

		elif is_valid == False:
			self.location.SetForegroundColour(wx.RED)
		

	def FeatureFieldOnText(self, event):
		newfeaturetext = str(self.featuretext.GetLineText(0)) #get location as string
	
		index = int(copy.copy(genbank.feature_selection))
		number = 0
		qualifier = genbank.gb.get_qualifier(index, number)[0]
		tag = newfeaturetext
		
		if ('=' in newfeaturetext) == False: # '=' is reserved for parsing, so it can't be in the string'
			self.featuretext.SetForegroundColour(wx.BLACK)
			genbank.gb.set_qualifier(index, number, qualifier, tag)
			self.update_globalUI()
			self.update_qualifiers()

		elif ('=' in newfeaturetext) == True:
			self.featuretext.SetForegroundColour(wx.RED)	

	def update_qualifiers(self):
		index = int(copy.copy(genbank.feature_selection))
		self.qualifier_list.DeleteAllItems()
		qualifiers = genbank.gb.get_qualifiers(index)
		for qualifier in qualifiers:
			col0, col1 = qualifier.split('=')
			col0 = col0[1:]
			self.qualifier_list.Append([col0, col1])	


######################################

class FeatureEditDialog(wx.Dialog):
	'''A class that puts the feature edit capabilities in a dialog'''
	def __init__(self, parent, title):
		super(FeatureEditDialog, self).__init__(parent=parent,id=wx.ID_ANY, title=title, size=(700, 300)) 		


		self.feature_edit = FeatureEdit(self, id=wx.ID_ANY)	#get the feature edit panel
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
		self.feature_edit.update_globalUI()
		self.Destroy()

## how to fix the cancel option????? ##
#	def OnCancel(self, event):
#		'''Reject new feeature from the "new feature" popup'''
#		genbank.gb.remove_feature(genbank.gb.get_feature(index=-1))
#		genbank.feature_selection = 0
#		self.Destroy()


######################################
######################################


##### main loop
class MyApp(wx.App):
	def OnInit(self):
		frame = wx.Frame(None, -1, title="Feature Edit", size=(700,300))
		panel =	FeatureEdit(frame, -1)
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
	assert len(sys.argv) == 3, 'Error, this script requires a path to a genbank file and a feature index (integer) as an argument.'
	print('Opening %s' % str(sys.argv[1]))

	genbank.gb = genbank.gbobject(str(sys.argv[1])) #make a genbank object and read file
	genbank.feature_selection = sys.argv[2]

	app = MyApp(0)
	app.MainLoop()


