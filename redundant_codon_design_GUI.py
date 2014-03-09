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
import redundant_codon_design as rcd

class RedundantCodon(wx.Panel):
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

		self.codontext = wx.StaticText(self, id=wx.ID_ANY, label='codon')
		textfont = wx.Font(18, wx.DECORATIVE, wx.ITALIC, wx.NORMAL)
		self.codontext.SetFont(textfont)
		self.codontext.SetLabel('None')
		
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
		globsizer.Add(self.codontext)
		globsizer.Add(spacer1)
		globsizer.Add(self.target_hits)
		globsizer.Add(spacer2)
		globsizer.Add(self.offtarget_hits)
		self.SetSizer(globsizer)

	def OnToggle(self, evt):
		'''When any togglebutton is pressed'''
		AA = []
		if self.Ala.GetValue(): AA.append('A')
		if self.Val.GetValue(): AA.append('V')
		if self.Ile.GetValue(): AA.append('I')
		if self.Leu.GetValue(): AA.append('L')
		if self.Met.GetValue(): AA.append('M')
		if self.Phe.GetValue(): AA.append('F')
		if self.Tyr.GetValue(): AA.append('Y')
		if self.Trp.GetValue(): AA.append('W')
		if self.Asp.GetValue(): AA.append('D')
		if self.Glu.GetValue(): AA.append('E')
		if self.Asn.GetValue(): AA.append('N')
		if self.Gln.GetValue(): AA.append('Q')
		if self.Ser.GetValue(): AA.append('S')
		if self.Thr.GetValue(): AA.append('T')
		if self.Arg.GetValue(): AA.append('R')
		if self.His.GetValue(): AA.append('H')
		if self.Lys.GetValue(): AA.append('K')
		if self.Cys.GetValue(): AA.append('C')
		if self.Gly.GetValue(): AA.append('G')
		if self.Pro.GetValue(): AA.append('P')
		if len(AA) > 0:
			codon, target, offtarget = rcd.run(AA)
			self.codontext.SetLabel(codon)
			targetstring = ''
			for entry in target:
				targetstring += entry+' '
			offtargetstring = ''
			for entry in offtarget:
				offtargetstring += entry+' '
			self.target_hits.SetLabel('Target amino acids: '+targetstring)
			self.offtarget_hits.SetLabel('Off-target amino acids: '+offtargetstring)	
		else:
			self.OnReset('')
	
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

if __name__ == '__main__': #if script is run by itself and not loaded	
	app = wx.App() # creation of the wx.App object (initialisation of the wxpython toolkit)
	frame = wx.Frame(None, title="Redundant Codon Design") # creation of a Frame with a title
	frame.RCD = RedundantCodon(frame) # creation of a richtextctrl in the frame
	frame.Show() # frames are invisible by default so we use Show() to make them visible
	app.MainLoop() # here the app enters a loop waiting for user input
