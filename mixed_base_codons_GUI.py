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

import wx
import mixed_base_codons as mbc

class MixedBaseCodon(wx.Panel):
	def __init__(self, parent):
		wx.Panel.__init__(self, parent)
#		self.dlg = wx.Panel(self, id=wx.ID_ANY)

		self.Ala = wx.ToggleButton(self, 1, 'Ala (A)')
		self.Ala.Bind(wx.EVT_TOGGLEBUTTON, self.OnToggle, id=1)
		self.Val = wx.ToggleButton(self, 2, 'Val (V)')
		self.Val.Bind(wx.EVT_TOGGLEBUTTON, self.OnToggle, id=2)
		self.Ile = wx.ToggleButton(self, 3, 'Ile (I)')
		self.Ile.Bind(wx.EVT_TOGGLEBUTTON, self.OnToggle, id=3)
		self.Leu = wx.ToggleButton(self, 4, 'Leu (L)')
		self.Leu.Bind(wx.EVT_TOGGLEBUTTON, self.OnToggle, id=4)
		self.Met = wx.ToggleButton(self, 5, 'Met (M)')
		self.Met.Bind(wx.EVT_TOGGLEBUTTON, self.OnToggle, id=5)
		sizer1 = wx.BoxSizer(wx.HORIZONTAL)
		sizer1.Add(self.Ala)
		sizer1.Add(self.Val)
		sizer1.Add(self.Ile)
		sizer1.Add(self.Leu)
		sizer1.Add(self.Met)

		self.Phe = wx.ToggleButton(self, 6, 'Phe (F)')
		self.Phe.Bind(wx.EVT_TOGGLEBUTTON, self.OnToggle, id=6)
		self.Tyr = wx.ToggleButton(self, 7, 'Tyr (Y)')
		self.Tyr.Bind(wx.EVT_TOGGLEBUTTON, self.OnToggle, id=7)
		self.Trp = wx.ToggleButton(self, 8, 'Trp (W)')
		self.Trp.Bind(wx.EVT_TOGGLEBUTTON, self.OnToggle, id=8)
		self.Asp = wx.ToggleButton(self, 9, 'Asp (D)')
		self.Asp.Bind(wx.EVT_TOGGLEBUTTON, self.OnToggle, id=9)
		self.Glu = wx.ToggleButton(self, 10, 'Glu (E)')
		self.Glu.Bind(wx.EVT_TOGGLEBUTTON, self.OnToggle, id=10)
		sizer2 = wx.BoxSizer(wx.HORIZONTAL)
		sizer2.Add(self.Phe)
		sizer2.Add(self.Tyr)
		sizer2.Add(self.Trp)
		sizer2.Add(self.Asp)
		sizer2.Add(self.Glu)

		self.Asn = wx.ToggleButton(self, 11, 'Asn (N)')
		self.Asn.Bind(wx.EVT_TOGGLEBUTTON, self.OnToggle, id=11)
		self.Gln = wx.ToggleButton(self, 12, 'Gln (Q)')
		self.Gln.Bind(wx.EVT_TOGGLEBUTTON, self.OnToggle, id=12)
		self.Ser = wx.ToggleButton(self, 13, 'Ser (S)')
		self.Ser.Bind(wx.EVT_TOGGLEBUTTON, self.OnToggle, id=13)
		self.Thr = wx.ToggleButton(self, 14, 'Thr (T)')
		self.Thr.Bind(wx.EVT_TOGGLEBUTTON, self.OnToggle, id=14)
		self.Arg = wx.ToggleButton(self, 15, 'Arg (R)')
		self.Arg.Bind(wx.EVT_TOGGLEBUTTON, self.OnToggle, id=15)
		sizer3 = wx.BoxSizer(wx.HORIZONTAL)
		sizer3.Add(self.Asn)
		sizer3.Add(self.Gln)
		sizer3.Add(self.Ser)
		sizer3.Add(self.Thr)
		sizer3.Add(self.Arg)

		self.His = wx.ToggleButton(self, 16, 'His (H)')
		self.His.Bind(wx.EVT_TOGGLEBUTTON, self.OnToggle, id=16)
		self.Lys = wx.ToggleButton(self, 17, 'Lys (K)')
		self.Lys.Bind(wx.EVT_TOGGLEBUTTON, self.OnToggle, id=17)
		self.Cys = wx.ToggleButton(self, 18, 'Cys (C)')
		self.Cys.Bind(wx.EVT_TOGGLEBUTTON, self.OnToggle, id=18)
		self.Gly = wx.ToggleButton(self, 19, 'Gly (G)')
		self.Gly.Bind(wx.EVT_TOGGLEBUTTON, self.OnToggle, id=19)
		self.Pro = wx.ToggleButton(self, 20, 'Pro (P)')
		self.Pro.Bind(wx.EVT_TOGGLEBUTTON, self.OnToggle, id=20)
		sizer4 = wx.BoxSizer(wx.HORIZONTAL)		
		sizer4.Add(self.His)
		sizer4.Add(self.Lys)
		sizer4.Add(self.Cys)
		sizer4.Add(self.Gly)
		sizer4.Add(self.Pro)

		self.Reset = wx.Button(self, 21, 'Reset')
		self.Reset.Bind(wx.EVT_BUTTON, self.OnReset, id=21)
		sizer5 = wx.BoxSizer(wx.HORIZONTAL)		
		sizer5.Add(self.Reset)

		spacer0 = wx.StaticText(self, id=wx.ID_ANY, label='codon')
		textfont = wx.Font(18, wx.DECORATIVE, wx.ITALIC, wx.NORMAL)
		spacer0.SetFont(textfont)
		spacer0.SetLabel('')

		self.codontext = wx.StaticText(self, id=wx.ID_ANY, label='codon')
		textfont = wx.Font(18, wx.DECORATIVE, wx.ITALIC, wx.NORMAL)
		self.codontext.SetFont(textfont)
		self.codontext.SetLabel('None')

#		self.codontext = wx.TextCtrl(self, id=wx.ID_ANY)
#		textfont = wx.Font(18, wx.DECORATIVE, wx.ITALIC, wx.NORMAL)
#		self.codontext.SetFont(textfont)
#		self.codontext.SetValue('None')
		
		spacer1 = wx.StaticText(self, id=wx.ID_ANY, label='')
		textfont = wx.Font(11, wx.DECORATIVE, wx.ITALIC, wx.NORMAL)
		spacer1.SetFont(textfont)
		spacer1.SetLabel('')

		self.target_hits = wx.StaticText(self, id=wx.ID_ANY, label='')
		textfont = wx.Font(11, wx.DECORATIVE, wx.ITALIC, wx.NORMAL)
		self.target_hits.SetFont(textfont)
		self.target_hits.SetLabel('Target amino acids: ')

		spacer2 = wx.StaticText(self, id=wx.ID_ANY, label='')
		textfont = wx.Font(11, wx.DECORATIVE, wx.ITALIC, wx.NORMAL)
		spacer2.SetFont(textfont)
		spacer2.SetLabel('')

		self.offtarget_hits = wx.StaticText(self, id=wx.ID_ANY, label='')
		textfont = wx.Font(11, wx.DECORATIVE, wx.ITALIC, wx.NORMAL)
		self.offtarget_hits.SetFont(textfont)
		self.offtarget_hits.SetLabel('Off-target amino acids: ')

		globsizer = wx.BoxSizer(wx.VERTICAL)
		globsizer.Add(sizer1)
		globsizer.Add(sizer2)
		globsizer.Add(sizer3)
		globsizer.Add(sizer4)
		globsizer.Add(sizer5)
		globsizer.Add(spacer0)
		globsizer.Add(self.codontext)
		globsizer.Add(spacer1)
		globsizer.Add(self.target_hits)
		globsizer.Add(spacer2)
		globsizer.Add(self.offtarget_hits)
		self.SetSizer(globsizer)

	def OnToggle(self, evt):
		'''When any togglebutton is pressed'''
		AA = []
		clicked_color = '#A2CD5A'
		unclicked_color = '#EEEEEE'
		if self.Ala.GetValue(): 
			AA.append('A')
			self.Ala.SetBackgroundColour(clicked_color)
		else:
			self.Ala.SetBackgroundColour(unclicked_color)

		if self.Val.GetValue(): 
			AA.append('V')
			self.Val.SetBackgroundColour(clicked_color)
		else:
			self.Val.SetBackgroundColour(unclicked_color)

		if self.Ile.GetValue(): 
			AA.append('I')
			self.Ile.SetBackgroundColour(clicked_color)
		else:
			self.Ile.SetBackgroundColour(unclicked_color)

		if self.Leu.GetValue(): 
			AA.append('L')
			self.Leu.SetBackgroundColour(clicked_color)
		else:
			self.Leu.SetBackgroundColour(unclicked_color)

		if self.Met.GetValue(): 
			AA.append('M')
			self.Met.SetBackgroundColour(clicked_color)
		else:
			self.Met.SetBackgroundColour(unclicked_color)

		if self.Phe.GetValue(): 
			AA.append('F')
			self.Phe.SetBackgroundColour(clicked_color)
		else:
			self.Phe.SetBackgroundColour(unclicked_color)

		if self.Tyr.GetValue(): 
			AA.append('Y')
			self.Tyr.SetBackgroundColour(clicked_color)
		else:
			self.Tyr.SetBackgroundColour(unclicked_color)

		if self.Trp.GetValue(): 
			AA.append('W')
			self.Trp.SetBackgroundColour(clicked_color)
		else:
			self.Trp.SetBackgroundColour(unclicked_color)

		if self.Asp.GetValue(): 
			AA.append('D')
			self.Asp.SetBackgroundColour(clicked_color)
		else:
			self.Asp.SetBackgroundColour(unclicked_color)

		if self.Glu.GetValue(): 
			AA.append('E')
			self.Glu.SetBackgroundColour(clicked_color)
		else:
			self.Glu.SetBackgroundColour(unclicked_color)

		if self.Asn.GetValue(): 
			AA.append('N')
			self.Asn.SetBackgroundColour(clicked_color)
		else:
			self.Asn.SetBackgroundColour(unclicked_color)

		if self.Gln.GetValue(): 
			AA.append('Q')
			self.Gln.SetBackgroundColour(clicked_color)
		else:
			self.Gln.SetBackgroundColour(unclicked_color)

		if self.Ser.GetValue(): 
			AA.append('S')
			self.Ser.SetBackgroundColour(clicked_color)
		else:
			self.Ser.SetBackgroundColour(unclicked_color)

		if self.Thr.GetValue(): 
			AA.append('T')
			self.Thr.SetBackgroundColour(clicked_color)
		else:
			self.Thr.SetBackgroundColour(unclicked_color)

		if self.Arg.GetValue(): 
			AA.append('R')
			self.Arg.SetBackgroundColour(clicked_color)
		else:
			self.Arg.SetBackgroundColour(unclicked_color)

		if self.His.GetValue(): 
			AA.append('H')
			self.His.SetBackgroundColour(clicked_color)
		else:
			self.His.SetBackgroundColour(unclicked_color)

		if self.Lys.GetValue(): 
			AA.append('K')
			self.Lys.SetBackgroundColour(clicked_color)
		else:
			self.Lys.SetBackgroundColour(unclicked_color)

		if self.Cys.GetValue(): 
			AA.append('C')
			self.Cys.SetBackgroundColour(clicked_color)
		else:
			self.Cys.SetBackgroundColour(unclicked_color)

		if self.Gly.GetValue(): 
			AA.append('G')
			self.Gly.SetBackgroundColour(clicked_color)
		else:
			self.Gly.SetBackgroundColour(unclicked_color)

		if self.Pro.GetValue(): 
			AA.append('P')
			self.Pro.SetBackgroundColour(clicked_color)
		else:
			self.Pro.SetBackgroundColour(unclicked_color)

		if len(AA) > 0:
			codon, target, offtarget, possibleAA = mbc.run(AA)
			self.codontext.SetLabel(codon)
			targetstring = ''
			for entry in target: #make string with the target AA
				targetstring += entry+' '
			offtargetstring = ''
			for entry in offtarget: #make string with the offtarget AA
				offtargetstring += entry+' '
			self.target_hits.SetLabel('Target amino acids: '+targetstring)
			self.offtarget_hits.SetLabel('Off-target amino acids: '+offtargetstring)	
	
	
			#set the yellow buttons to indicate which AA are possible without further off-target hits	
			possible_color = '#FFFF99'	
			if 'A' in possibleAA and self.Ala.GetValue() == False: 
				self.Ala.SetBackgroundColour(possible_color)

			if  'V' in possibleAA and self.Val.GetValue() == False: 
				self.Val.SetBackgroundColour(possible_color)

			if  'I' in possibleAA and self.Ile.GetValue() == False: 
				self.Ile.SetBackgroundColour(possible_color)

			if  'L' in possibleAA and self.Leu.GetValue() == False: 
				self.Leu.SetBackgroundColour(possible_color)

			if  'M' in possibleAA and self.Met.GetValue() == False: 
				self.Met.SetBackgroundColour(possible_color)

			if  'F' in possibleAA and self.Phe.GetValue() == False: 
				self.Phe.SetBackgroundColour(possible_color)

			if  'Y' in possibleAA and self.Tyr.GetValue() == False: 
				self.Tyr.SetBackgroundColour(possible_color)

			if  'W' in possibleAA and self.Trp.GetValue() == False: 
				self.Trp.SetBackgroundColour(possible_color)

			if  'D' in possibleAA and self.Asp.GetValue() == False: 
				self.Asp.SetBackgroundColour(possible_color)

			if  'E' in possibleAA and self.Glu.GetValue() == False: 
				self.Glu.SetBackgroundColour(possible_color)

			if  'N' in possibleAA and self.Asn.GetValue() == False: 
				self.Asn.SetBackgroundColour(possible_color)

			if  'Q' in possibleAA and self.Gln.GetValue() == False: 
				self.Gln.SetBackgroundColour(possible_color)

			if  'S' in possibleAA and self.Ser.GetValue() == False: 
				self.Ser.SetBackgroundColour(possible_color)

			if  'T' in possibleAA and self.Thr.GetValue() == False: 
				self.Thr.SetBackgroundColour(possible_color)

			if  'R' in possibleAA and self.Arg.GetValue() == False: 
				self.Arg.SetBackgroundColour(possible_color)

			if  'H' in possibleAA and self.His.GetValue() == False: 
				self.His.SetBackgroundColour(possible_color)

			if  'K' in possibleAA and self.Lys.GetValue() == False: 
				self.Lys.SetBackgroundColour(possible_color)

			if  'C' in possibleAA and self.Cys.GetValue() == False: 
				self.Cys.SetBackgroundColour(possible_color)

			if  'G' in possibleAA and self.Gly.GetValue() == False: 
				self.Gly.SetBackgroundColour(possible_color)

			if  'P' in possibleAA and self.Pro.GetValue() == False: 
				self.Pro.SetBackgroundColour(possible_color)

		else:
			self.codontext.SetLabel('None')
			self.target_hits.SetLabel('Target amino acids: ')
			self.offtarget_hits.SetLabel('Off-target amino acids: ')

	
	def OnReset(self, evt):
		'''Reset all buttons so that they are unclicked'''
		self.Ala.SetValue(False)
		self.Val.SetValue(False)
		self.Ile.SetValue(False)
		self.Leu.SetValue(False)
		self.Met.SetValue(False)
		self.Phe.SetValue(False)
		self.Tyr.SetValue(False)
		self.Trp.SetValue(False)
		self.Asp.SetValue(False)
		self.Glu.SetValue(False)
		self.Asn.SetValue(False)
		self.Gln.SetValue(False)
		self.Ser.SetValue(False)
		self.Thr.SetValue(False)
		self.Arg.SetValue(False)
		self.His.SetValue(False)
		self.Lys.SetValue(False)
		self.Cys.SetValue(False)
		self.Gly.SetValue(False)
		self.Pro.SetValue(False)
		self.codontext.SetLabel('None')
		self.target_hits.SetLabel('Target amino acids: ')
		self.offtarget_hits.SetLabel('Off-target amino acids: ')	
		self.OnToggle("")

if __name__ == '__main__': #if script is run by itself and not loaded	
	app = wx.App() # creation of the wx.App object (initialisation of the wxpython toolkit)
	frame = wx.Frame(None, title="Mixed Base Codons") # creation of a Frame with a title
	frame.MBC = MixedBaseCodon(frame) # creation of a richtextctrl in the frame
	frame.Show() # frames are invisible by default so we use Show() to make them visible
	app.MainLoop() # here the app enters a loop waiting for user input
