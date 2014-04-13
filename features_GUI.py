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
#match features with qualifiers (mandatory and optional)
#add a function that checks that everything is ok


#add edit buttons to location and the first qualifier entry (which is used for name of feature)

#add a way of actually adding a new feature...


#make changes to how genbank handles /qualifier=xyz, the '=' is not always there...

import ast
import wx
import genbank
import sys, os
from string import *

#for asserts and tests
import types
import unittest

files={}   #list with all configuration files

files['default_dir'] = os.path.abspath(os.path.dirname(sys.argv[0]))+"/"
files['default_dir']=replace(files['default_dir'], "\\", "/")
files['default_dir']=replace(files['default_dir'], "library.zip", "")

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
		buttonsizer.Add(self.button_ok)
		buttonsizer.Add(self.button_cancel)

		globsizer = wx.BoxSizer(wx.VERTICAL)
		globsizer.Add(self.qualifier_header, flag=wx.LEFT, border=10)
		globsizer.Add(self.qualifier_combobox, flag=wx.LEFT, border=10)
		globsizer.Add(self.spacer1, flag=wx.LEFT, border=10)
		globsizer.Add(self.tag_header, flag=wx.LEFT, border=10)
		globsizer.Add(self.tag_field, flag=wx.LEFT, border=10)
		globsizer.Add(buttonsizer, flag=wx.LEFT, border=10)

		self.SetSizer(globsizer)
		self.Center()
		
	def OnOK(self, evt):
		'Return wx.ID_OK if ok button was pressed'
		self.EndModal(wx.ID_OK)

	def OnCancel(self, evt):
		'Return wx.ID_CANCEL if cancel button was pressed'
		self.EndModal(wx.ID_CANCEL)


class FeatureEdit(wx.Panel):
	def __init__(self, parent, id):
		wx.Panel.__init__(self, parent)

#		splitter2 = wx.SplitterWindow(self, -1, style=wx.SP_3D)

		self.feature_dlg = wx.Panel(self, id=wx.ID_ANY)

		#feature label control
		featuretext = wx.StaticText(self.feature_dlg, id=wx.ID_ANY, label='Feature:')
		featuretext.SetFont(wx.Font(11, wx.DECORATIVE, wx.ITALIC, wx.NORMAL))


		#make box that displays the mixed base codon for editing
#		textfont = wx.Font(18, wx.DECORATIVE, wx.ITALIC, wx.NORMAL)
#		self.featuretext = wx.stc.StyledTextCtrl(self, size=(200,30), style = wx.NO_FULL_REPAINT_ON_RESIZE)
#		self.featuretext.StyleSetBackground(style=wx.stc.STC_STYLE_DEFAULT, back='#FFFFFF') #set background color of everything that is not text
#		self.featuretext.StyleSetBackground(style=0, back='#FFFFFF') #set background color of text
##		self.featuretext.StyleSetBackground(style=wx.stc.STC_STYLE_LINENUMBER, back='#FFFFFF') #sets color of left margin
#		self.featuretext.SetUseHorizontalScrollBar(show=False) #hide scrollbar
#		face = textfont.GetFaceName()
#		size = textfont.GetPointSize()
#		self.featuretext.StyleSetSpec(wx.stc.STC_STYLE_DEFAULT,"face:%s,size:%d" % (face, size))

		self.featuretext = wx.TextCtrl(self.feature_dlg, id=wx.ID_ANY, size=(200,-1))
		self.featuretext.Bind(wx.EVT_TEXT, self.FeatureFieldOnText)

		
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

		self.type_combobox = wx.ComboBox(self.feature_dlg, id=1002, size=(200, -1), choices=featuretypes, style=wx.CB_READONLY)
		self.type_combobox.Bind(wx.EVT_COMBOBOX, self.TypeComboboxOnSelect)




		#location
		locationtext = wx.StaticText(self.feature_dlg, id=wx.ID_ANY, label='Location:')
		locationtext.SetFont(wx.Font(11, wx.DECORATIVE, wx.ITALIC, wx.NORMAL))
		self.location = wx.TextCtrl(self.feature_dlg, id=1003, size=(200,-1), style=wx.TE_RICH)
		self.location.Bind(wx.EVT_TEXT, self.LocationFieldOnText)
		
		#complement or not
		self.complementbox = wx.CheckBox(self.feature_dlg, id=3001, label="Complement")
		self.complementbox.Bind(wx.EVT_CHECKBOX, self.ComplementCheckboxOnSelect)	

		#Join or not
		self.joinbox = wx.CheckBox(self.feature_dlg, id=3002, label="Join (location ranges should be joined)")
		self.joinbox.Bind(wx.EVT_CHECKBOX, self.JoinCheckboxOnSelect)

		#Order or not
		self.orderbox = wx.CheckBox(self.feature_dlg, id=3003, label="Order (location ranges indicate order)")
		self.orderbox.Bind(wx.EVT_CHECKBOX, self.OrderCheckboxOnSelect)
		#need to add logic that only makes join and order available if there are more than one location for a feature...


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
		#second panel
		self.qualifier_list = wx.ListCtrl(self, -1, style=wx.LC_REPORT)
		self.qualifier_list.InsertColumn(0, "Qualifier", format=wx.LIST_FORMAT_LEFT, width=120)
		self.qualifier_list.InsertColumn(1, "Tag", format=wx.LIST_FORMAT_LEFT, width=250)
#		self.qualifier_list.Bind(wx.EVT_LISTBOX_DCLICK, self.OnEditQualifier(None))

		imageFile = files['default_dir']+"/icon/new_small.png"
		image1 = wx.Image(imageFile, wx.BITMAP_TYPE_ANY).ConvertToBitmap()
		addqual = wx.BitmapButton(self, id=7, bitmap=image1, size = (image1.GetWidth()+15, image1.GetHeight()+15), name = "new")

		imageFile = files['default_dir']+"/icon/remove_small.png"
		image1 = wx.Image(imageFile, wx.BITMAP_TYPE_ANY).ConvertToBitmap()
		deletequal = wx.BitmapButton(self, id=8, bitmap=image1, size = (image1.GetWidth()+15, image1.GetHeight()+15), name = "remove")

		imageFile = files['default_dir']+"/icon/move_up.png"
		image1 = wx.Image(imageFile, wx.BITMAP_TYPE_ANY).ConvertToBitmap()
		qualup = wx.BitmapButton(self, id=9, bitmap=image1, size = (image1.GetWidth()+15, image1.GetHeight()+15), name = "move up")

		imageFile = files['default_dir']+"/icon/move_down.png"
		image1 = wx.Image(imageFile, wx.BITMAP_TYPE_ANY).ConvertToBitmap()
		qualdown = wx.BitmapButton(self, id=10, bitmap=image1, size = (image1.GetWidth()+15, image1.GetHeight()+15), name = "move down")

		imageFile = files['default_dir']+"/icon/edit.png"
		image1 = wx.Image(imageFile, wx.BITMAP_TYPE_ANY).ConvertToBitmap()
		qualedit = wx.BitmapButton(self, id=11, bitmap=image1, size = (image1.GetWidth()+15, image1.GetHeight()+15), name = "edit")
		
		

		sizer = wx.BoxSizer(wx.VERTICAL)
		sizer.Add(addqual)
		sizer.Add(deletequal)
		sizer.Add(qualup)
		sizer.Add(qualdown)
		sizer.Add(qualedit)
		
		sizer2 = wx.BoxSizer(wx.HORIZONTAL)
		sizer2.Add(self.qualifier_list, 3, wx.EXPAND)
		sizer2.Add(sizer, 0, wx.EXPAND)

#		self.SetSizer(sizer2)

		#bind qualifier buttions
		self.Bind(wx.EVT_BUTTON, self.OnAddQualifier, id=7)
		self.Bind(wx.EVT_BUTTON, self.OnRemoveQualifier, id=8)
		self.Bind(wx.EVT_BUTTON, self.OnMoveQualifierUp, id=9)
		self.Bind(wx.EVT_BUTTON, self.OnMoveQualifierDown, id=10)
		self.Bind(wx.EVT_BUTTON, self.OnEditQualifier, id=11)

		##
		#third panel, for showing qualifiers
#		self.qualifier_list = FeatureEdit(splitter2, -1)
#		font = wx.Font(pointSize=10, family=wx.FONTFAMILY_DEFAULT, style=wx.FONTSTYLE_NORMAL, weight=wx.FONTWEIGHT_NORMAL, underline=False, faceName='Monospace', encoding=wx.FONTENCODING_DEFAULT)
#		self.qualifier_list.SetFont(font)

#		splitter2.SplitVertically(self.feature_dlg, self.qualifier_list)

		main_sizer = wx.BoxSizer(wx.HORIZONTAL)
		main_sizer.Add(self.feature_dlg, -1, wx.EXPAND)

		main_sizer.Add(sizer2, -1, wx.EXPAND)
		self.SetSizer(main_sizer)
		self.Centre()

		self.updateUI()

	def OnAddQualifier(self, event):
		index = genbank.gb.get_feature_selection()
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
			self.updateUI()

			number = self.qualifier_list.GetItemCount()-1
			self.update_qualifier_selection(index, number)

		dlg.Destroy()


	

	def OnRemoveQualifier(self, event):
		index = genbank.gb.get_feature_selection() #feature selected
		feature = genbank.gb.get_feature(index)
		number = self.qualifier_list.GetFirstSelected() #qualifier selected
		items = self.qualifier_list.GetItemCount() #number of qualifiers

		if items != 1: #don't delete last qualifier
			genbank.gb.remove_qualifier(feature, number)
			self.updateUI()
			try:
				self.GetParent().GetParent().feature_list.updateUI() #update feature viewer
			except:
				pass

			#set highlight, focus and selection
			if number == items-1: #if last item was deleted, selection is set as -1
				number = number-1
			self.update_qualifier_selection(index, number)


	def OnMoveQualifierUp(self, event):
		index = genbank.gb.get_feature_selection()
		feature = genbank.gb.get_feature(index)
		number = self.qualifier_list.GetFirstSelected()
		genbank.gb.move_qualifier(feature, number, 'u')
		self.updateUI()
		try:
			self.GetParent().GetParent().feature_list.updateUI() #update feature viewer
		except:
			pass

		if number != 0:
			number = number-1
		self.update_qualifier_selection(index, number)	


	def OnMoveQualifierDown(self, event):
		index = genbank.gb.get_feature_selection()
		feature = genbank.gb.get_feature(index)
		number = self.qualifier_list.GetFirstSelected()
		genbank.gb.move_qualifier(feature, number, 'd')
		self.updateUI()
		try:
			self.GetParent().GetParent().feature_list.updateUI() #update feature viewer
		except:
			pass

		if number != self.qualifier_list.GetItemCount():
			number = number+1
		self.update_qualifier_selection(index, number)

	def OnEditQualifier(self, event):
		'''Make popup window where qualifier and tag can be edited'''
		index = genbank.gb.get_feature_selection()
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
		self.updateUI()
		try:
			self.GetParent().GetParent().feature_list.updateUI() #update feature viewer
		except:
			raise IOError('Error updating feature viewer UI') 
		self.update_qualifier_selection(index, number)

	def update_qualifier_selection(self, index, number):
		'''Updates which feature and qualifier is selected'''
		try:
			self.GetParent().GetParent().feature_list.update_feature_selection(index) #make sure the right feature is selected
		except:
			pass

		###change this! it is not right. #edit# I cannot see what is wrong...

		self.qualifier_list.SetItemState(number, wx.LIST_STATE_SELECTED,wx.LIST_STATE_SELECTED) #for the highlight
		self.qualifier_list.Select(number, True) #to select it
		self.qualifier_list.Focus(number) #to focus it
	
	def ComplementCheckboxOnSelect(self, event):
		'''Toggle whether the feature is on the complement strand or not'''
		newcomplement = self.complementbox.GetValue()
		index = genbank.gb.get_feature_selection()
		feature = genbank.gb.get_feature(index)
		genbank.gb.set_feature_complement(feature, newcomplement)
		self.updateUI()
		try:
			self.GetParent().GetParent().feature_list.updateUI() #update feature viewer
			self.GetParent().GetParent().feature_list.update_feature_selection(index) #re-select the feature
		except:
			pass
		

	def JoinCheckboxOnSelect(self, event):
		'''Toggle whether a feature with multiple locations should be joined or not'''
		newjoin = self.joinbox.GetValue()
		index = genbank.gb.get_feature_selection()
		feature = genbank.gb.get_feature(index)
		genbank.gb.set_feature_join(feature, newjoin)
		self.updateUI()
		try:
			self.GetParent().GetParent().feature_list.updateUI() #update feature viewer
			self.GetParent().GetParent().feature_list.update_feature_selection(index) #re-select the feature
		except:
			pass

	def OrderCheckboxOnSelect(self, event):
		'''Toggle whether a feature with ultiple locations should be indicated as being in a certain order or not'''
		neworder = self.orderbox.GetValue()
		index = genbank.gb.get_feature_selection()
		feature = genbank.gb.get_feature(index)
		genbank.gb.set_feature_order(feature, neworder)
		self.updateUI()
		try:
			self.GetParent().GetParent().feature_list.updateUI() #update feature viewer
			self.GetParent().GetParent().feature_list.update_feature_selection(index) #re-select the feature
		except:
			pass
		
	def TypeComboboxOnSelect(self, event):
		newkey = self.type_combobox.GetValue()
		index = genbank.gb.get_feature_selection()
		feature = genbank.gb.get_feature(index)
		genbank.gb.set_feature_type(feature, newkey)
		self.updateUI()
		try:
			self.GetParent().GetParent().feature_list.updateUI() #update feature viewer
			self.GetParent().GetParent().feature_list.update_feature_selection(index) #re-select the feature
		except:
			pass
		
	def TestValidLocation(self, locationlist):
		'''Takes a location list and tests whether it is valid'''
		result = True		
		try:
			assert type(locationlist) == types.ListType
			for location in locationlist:
				start, finish = location.split('..')
				start = int(start)
				finish = int(finish)
				assert start < finish
				assert finish <= len(genbank.gb.get_dna())
		except:
			result = False
		return result



	def LocationFieldOnText(self, event):
		'''When the location field is changed. Check if the entry is valid, if yes, then update features'''
		newlocation = str(self.location.GetLineText(0)) #get location as string
		newlocation.replace(' ', '') #get rid of spaces

		#convert to list
		if ',' in newlocation:
			locationlist = newlocation.split(',')
		else:
			locationlist = [newlocation]

		index = genbank.gb.get_feature_selection()
		feature = genbank.gb.get_feature(index)
		
		
		is_valid = self.TestValidLocation(locationlist) #test if location entry is valid
		if is_valid == True:
			print(locationlist)
			print(type(locationlist))


			self.location.SetForegroundColour(wx.BLACK)
			genbank.gb.set_feature_location(feature, locationlist)

######
#Needs fixing! Changing location does not currently work the correct way. Almost there though
#####

#			self.updateUI()
			try:
				self.GetParent().GetParent().feature_list.updateUI() #update feature viewer
#				self.GetParent().GetParent().feature_list.update_feature_selection(index) #re-select the feature
			except:
				pass

		elif is_valid == False:
			self.location.SetForegroundColour(wx.RED)


	def FeatureFieldOnText(self, event):
		print('feature')

	def updateUI(self):
		'''Updates all fields depending on which feature is chosen'''
		try:
			#get selected feature
			index = genbank.gb.get_feature_selection()
		
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



#			if self.joinbox:
#				self.orderbox.EnableTool(512, 0)

			#update qualifier field
			self.qualifier_list.DeleteAllItems()
			qualifiers = genbank.gb.get_qualifiers(index)
			for qualifier in qualifiers:
				print(qualifier)
				col0, col1 = qualifier.split('=')
				col0 = col0[1:]
				self.qualifier_list.Append([col0, col1])
#				self.qualifier_list.SetColumnWidth(col=0, width=wx.LIST_AUTOSIZE)
#				self.qualifier_list.SetColumnWidth(col=1, width=wx.LIST_AUTOSIZE)			

		except:
			pass




class FeatureView(wx.Panel):
	"""Class for viewing features as a list"""
	def __init__(self, parent, id):
		wx.Panel.__init__(self, parent)
		self.feature_list = wx.ListCtrl(self, id=3001, style=wx.LC_REPORT|wx.LC_SINGLE_SEL)

		self.feature_list.InsertColumn(0, "Feature", format=wx.LIST_FORMAT_LEFT, width=200)
		self.feature_list.InsertColumn(1, "Type", format=wx.LIST_FORMAT_LEFT, width=100)
		self.feature_list.InsertColumn(2, "Location on DNA", format=wx.LIST_FORMAT_LEFT, width=200)
		self.feature_list.InsertColumn(3, "Strand", format=wx.LIST_FORMAT_LEFT, width=120)

#		font = wx.Font(pointSize=10, family=wx.FONTFAMILY_DEFAULT, style=wx.FONTSTYLE_NORMAL, weight=wx.FONTWEIGHT_NORMAL, underline=False, faceName='Monospace', encoding=wx.FONTENCODING_DEFAULT)
#		self.feature_list.SetFont(font)
		self.feature_list.Bind(wx.EVT_LIST_ITEM_FOCUSED, self.ListOnSelect)
		self.feature_list.Bind(wx.EVT_LIST_ITEM_ACTIVATED, self.ListOnActivate)

		
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

		#arrange buttons vertically		
		sizer = wx.BoxSizer(wx.VERTICAL)
		sizer.Add(newfeature)
		sizer.Add(deletefeature)
		sizer.Add(moveup)
		sizer.Add(movedown)

		#add feature list and buttons horisontally	
		sizer2 = wx.BoxSizer(wx.HORIZONTAL)
		sizer2.Add(self.feature_list, 3, wx.EXPAND)
		sizer2.Add(sizer, 0, wx.EXPAND)

		self.SetSizer(sizer2)
		self.updateUI()

	def get_selection(self):
		"""Get currently selected feature"""
		return self.feature_list.GetFocusedItem()

	def ListOnSelect(self, event):	
		'''Updates selection depending on which feature is chosen'''
		index = self.get_selection()
		genbank.gb.set_feature_selection(index)
		try:
			self.GetParent().GetParent().feature_edit.updateUI() #update feature editor
		except:
			pass

	def ListOnActivate(self, event):
		'''Updates feature AND DNA selection based on which feature is chosen'''
		index = self.get_selection()
		genbank.gb.set_feature_selection(index)

		locations = genbank.gb.get_feature_location(index)
		start = genbank.gb.get_location(locations[0])[0]
		finish = genbank.gb.get_location(locations[-1])[1]
		genbank.gb.set_dna_selection((start-1, finish))  #how to I propagate this to the DNA view???
		self.GetTopLevelParent().dnaview.gbviewer.SetSelection(start-1, finish) #update DNA selection
		self.GetTopLevelParent().dnaview.gbviewer.ShowPosition(start) #show the selection
######
######
######
######
	def OnNew(self, event):
		'''Make new feature'''
		#make feature and update interface
		self.GetTopLevelParent().match_selection()

		self.NewFeatureFrame = wx.Frame(None, title="New Feature", size=(700, 200)) # creation of a Frame with a title
		self.feature_edit = FeatureEdit(self.NewFeatureFrame, id=wx.ID_ANY)		

#		self.OK = wx.Button(self.NewFeatureFrame, 7, 'OK')
#		self.sizer = wx.BoxSizer(wx.VERTICAL)
#		self.sizer.Add(self.feature_edit)
#		self.sizer.Add(self.OK)

		self.NewFeatureFrame.Show()
		self.NewFeatureFrame.Center()
		

	def OnDelete(self, event):
		'''Delete selected feature'''
		#identify feature, remove it and update interface
		index = self.get_selection()
		feature = genbank.gb.get_feature(index)
		index = self.feature_list.GetFirstSelected()
		genbank.gb.remove_feature(feature)
		self.updateUI()
	
		#set highlight, focus and selection
		if index == self.feature_list.GetItemCount():
			index = index-1
		self.update_feature_selection(index)


	def OnMoveFeatureUp(self, event):
		'''Move feature up one step'''
		index = self.get_selection()
		feature = genbank.gb.get_feature(index)
		index = self.feature_list.GetFirstSelected()
		genbank.gb.move_feature(feature, 'u')	
		self.updateUI()
		if index != 0:
			index = index-1
		self.update_feature_selection(index)
		

	def OnMoveFeatureDown(self, event):
		'''Move feature up down step'''
		index = self.get_selection()
		feature = genbank.gb.get_feature(index)
		index = self.feature_list.GetFirstSelected()
		genbank.gb.move_feature(feature, 'd')
		self.updateUI()
	
		if index != self.feature_list.GetItemCount()-1:
			index = index+1
		self.update_feature_selection(index)


	def update_feature_selection(self, index):
		'''Updates which feature is selected'''
		genbank.gb.set_feature_selection(index)
		self.feature_list.SetItemState(index, wx.LIST_STATE_SELECTED,wx.LIST_STATE_SELECTED) #for the highlight
		self.feature_list.Select(index, True) #to select it
		self.feature_list.Focus(index) #to focus it


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
		try:
			self.feature_list.DeleteAllItems()
			n = 0 #for feautrecolor
			for entry in genbank.gb.gbfile['features']:
				col0 = entry['qualifiers'][0].split('=')[1]
		#		col0 = 'T7\terminator'
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
				color = self.current_highlight_color
				item = n
				self.feature_list.SetItemBackgroundColour(item, color)	
				n += 1
			#self.autosize()
		except:
			pass

class FeatureCreate(wx.Panel):
	def __init__(self, parent, id, editor):
		wx.Panel.__init__(self, parent)
		

		if editor == True: #if feature editor should be included
			splitter1 = wx.SplitterWindow(self, -1, style=wx.SP_3D)

			##
			#first panel, for showing feature overview
			self.feature_list = FeatureView(splitter1, id=wx.ID_ANY)

			##
			#second panel, for editing
			self.feature_edit = FeatureEdit(splitter1, id=wx.ID_ANY)

			splitter1.SplitHorizontally(self.feature_list, self.feature_edit)

			#global sizer		
			globsizer = wx.BoxSizer(wx.HORIZONTAL)
			globsizer.Add(splitter1, -1, wx.EXPAND)

		elif editor == False: #if feature editor should not be included
			##
			#first panel, for showing feature overview
			self.feature_list = FeatureView(self, id=wx.ID_ANY)			

			#global sizer		
			globsizer = wx.BoxSizer(wx.HORIZONTAL)
			globsizer.Add(self.feature_list, -1, wx.EXPAND)

		self.SetSizer(globsizer)
		self.Centre()

	def updateUI(self):
		"""Update feature list content"""
		self.feature_list.updateUI()
		try:
			self.feature_edit.updateUI()
		except:
			pass


