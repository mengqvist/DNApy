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


# file to store options!


class options():
		
		# plasmid view settings
		pRadius 			= 26					# radius of plasmid
		pArrowThick   		= 2.4					# thickness of arrow
		pArrowHeadLength 	= 3.2 					# length of arrow head
		pRadChange    		= pArrowThick/2			# initial space between arrow and plasmid
		pLevAdd             = pArrowThick + 0.8		# level distance
		
		labelRadius			= pRadius+pRadChange+3*pArrowThick+pLevAdd # radius for labels and stuff
		
		
		pFeatureBC      	 = [0, 0, 0 ]			# bordercolor of arrow 
		pBW      			 = 0.1					# thickness of border
		
		pFeatureBChigh  	 = [1,0.2, 0]			# bordercolor of arrow when highlighted
		pBWhigh  			 = 0.35					# thickness of border when highlighted
		
		
		
		#####################
		
		
