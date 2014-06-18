#!/usr/bin/python


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

from output import create
import genbank
from base_class import DNApyBaseClass
from wx.lib.pubsub import pub

class MyPanel(create):
	'''Class to add update behavior to the output panel'''
	def __init__(self, parent, style):
		super(MyPanel, self).__init__(parent, style) 
	


		#determing which listening group from which to recieve messages about UI updates
		self.listening_group = 'from_dna_edit' 		
		pub.Publisher.subscribe(self.listen_to_updateUI, self.listening_group)

		self.listening_group2 = 'from_feature_edit'	
		pub.Publisher.subscribe(self.listen_to_updateUI, self.listening_group2)

		self.listening_group3 = 'from_feature_list'		
		pub.Publisher.subscribe(self.listen_to_updateUI, self.listening_group3)




####### Modify methods from base calss to fit current needs #########

	def update_globalUI(self):
		'''Method should be modified as to update other panels in response to changes in own panel.
		Preferred use is through sending a message using the pub module.
		Example use is: pub.Publisher.sendMessage('feature_list_updateUI', '').
		The first string is the "listening group" and deterimines which listeners get the message. 
		The second string is the message and is unimportant for this implementation.
		The listening group assigned here (to identify recipients) must be different from the listening group assigned in __init__ (to subscribe to messages).'''
		pass


	def update_ownUI(self):
		'''Updates all fields depending on which feature is chosen'''
		genbankfile = genbank.gb.make_gbstring()
		self.clear()
		self.write(genbankfile, 'Text')

#####################################################################







			
