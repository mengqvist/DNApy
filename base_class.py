#!/usr/bin/env python


#This file is part of DNApy. DNApy is a DNA editor written purely in python. 
#The program is intended to be an intuitive, fully featured, 
#extendable, editor for molecular and synthetic biology.  
#Enjoy!
#
#Copyright (C) 2014  Martin K. M. Engqvist | 
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
import os, sys
from wx.lib.pubsub import Publisher as pub


class DNApyBaseClass(wx.Panel):
	'''This is a base class which should serve as the base for any panel used in DNApy.'''
	def __init__(self, parent, id):
		wx.Panel.__init__(self, parent)

		self.listening_group = NotImplementedError #needs to be assigned a listening group name (a string) to recieve requests for own panel UI updates.
		pub.subscribe(self.listen_to_updateUI, self.listening_group)


#### own UI updates ####
	def update_ownUI(self):
		'''Method should be modified as to update graphical elements of own panel only. 
			Preferred use is through a self.update_ownUI() call.'''
		raise NotImplementedError

	def listen_to_updateUI(self, msg):
		'''For recieving requests for UI updates.
			This is the listening method for recieving requests for own UI updates.
			Its purpuse is to update the UI of own panel upon recieveing an UI update message.
			'''
		self.update_ownUI()


#### other UI updates ####
	def update_globalUI(self):
		'''Method should be modified as to update other panels in response to changes in own panel.
			Preferred use is through sending a message using the pub module.
			Example use is: pub.sendMessage('feature_list_updateUI', '').
			The first string is the "listening group" and deterimines which listeners get the message. 
			The second string is the message and is unimportant for this implementation.
			The listening group assigned here (to identify recipients) must be different from the listening group assigned in __init__ (to subscribe to messages).'''
		raise NotImplementedError

		

	


