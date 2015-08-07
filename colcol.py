#!/usr/bin/env python


#DNApy is a DNA editor written purely in python.
#The program is intended to be an intuitive and fully featured 
#editor for molecular and synthetic biology.
#Enjoy!
#
#copyright (C) 2014-2015  Martin Engqvist |
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#LICENSE:
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
#Get source code at: https://github.com/mengqvist/DNApy
#


#This file mainly deals with color conversions, color transformations, and generating color scales.

import re
import wsvg
import colorsys

def is_rgb(in_col):
	'''
	Check whether input is a valid RGB color.
	Return True if it is, otherwise False.
	'''
	if len(in_col) == 3 and type(in_col) == tuple:
		if 0<=in_col[0]<=255 and 0<=in_col[1]<=255 and 0<=in_col[2]<=255:
			return True
		else:
			return False
	else:
		return False
	
def rgb_to_hex(rgb):
	'''
	Convert RGB colors to hex.
	Input should be a tuple of integers (R, G, B).
	Output is a string representing a hex numer. For instance '#FFFFFF'.
	''' 
	#make sure input is ok
	assert is_rgb(rgb) is True, 'Error, %s is not a valid RGB color.' % str(rgb)
	
	#make conversion
	return "#%02x%02x%02x" % rgb

	
	
def is_hex(in_col):
	'''
	Check whether an input string is a valid hex value.
	Return True if it is, otherwise False.
	'''
	if type(in_col) is not str:
		return False
		
	regular_expression = re.compile(r'''^							#match beginning of string
										[#]{1} 						#exactly one hash
										[0-9a-fA-F]{6}				#exactly six of the hex symbols  0 to 9, a to f (big or small)
										$							#match end of string
										''', re.VERBOSE)	
	
	if regular_expression.match(in_col) == None:
		return False
	else:
		return True

	
def hex_to_rgb(in_col):
	'''
	Convert a hex color to RGB.
	Input should be a string. For example '#FFFFFF'.
	Output is a tuple of integers (R, G, B).
	'''
	#make sure input is ok
	assert is_hex(in_col) is True, 'Error, %s is not a valid hex color.' % in_col
	
	#make the conversion
	in_col = in_col.lstrip('#')
	return tuple([int(in_col[s:s+2], 16) for s in range(0, len(in_col), 2)])
	

def scale(col1, col2, white_mid=False):
	'''
	Function makes a color scale from 0 to 100 using the supplied colors.
	The variable white_mid is boolean and determines whether the colors 
	should transition over white in the middle or just smoothly run into each other.
	The function returns a dictionary of integer keys and hex color values corresponding to the scores 0 to 100.
	'''
	assert is_hex(col1) or is_rgb(col1), 'Error, the first color is not valid.'
	assert is_hex(col2) or is_rgb(col2), 'Error, the second color is not valid.'
	assert type(white_mid) is bool, 'Error, white_mid must be a boolean.'
	
	if is_hex(col1):
		col1 = hex_to_rgb(col1)
	if is_hex(col2):
		col2 = hex_to_rgb(col2)
	
	color_dict={}
	if white_mid is False:
		#determine how many points R, G and B should change with each score
		R_diff = (col2[0]-col1[0])/100.0
		G_diff = (col2[1]-col1[1])/100.0
		B_diff = (col2[2]-col1[2])/100.0

		#set starting values
		R, G, B = col1 
		for i in range(101):
			color_dict[i] = rgb_to_hex((R+int(R_diff*i), G+int(G_diff*i), B+int(B_diff*i)))
			
	elif white_mid is True:
		first_half = scale((col1), (255,255,255))
		for key in range(0, 100, 2):
			color_dict[key/2] = first_half[key]
			
		second_half = scale((255,255,255), col2)
		for key in range(0, 102, 2):
			color_dict[50+key/2] = second_half[key]
		
	return color_dict
	
	
def mix_colors(col1, col2):
	'''
	Mix two colors and return the result.
	The input colors can be either RGB or hex.
	It is also possible to use one RGB value and one hex value.
	'''
	color_dict = scale(col1, col2, white_mid=False)
	return color_dict[50]
	

def complementary(in_col):
	'''
	Returns input color and its complementary color as a list of hex or rgb values, depending on what was submitted.

		O
	   x x
	  x   x
	 x     x
	  x   x
	   x x
		O

	'''
	assert (is_rgb(in_col) or is_hex(in_col)), 'Error, the input must be a hex color string or an RGB tuple.'

	hex_color = is_hex(in_col) #check whether input is hex
	if hex_color: #if it is, make it RGB
		in_col = hex_to_rgb(in_col)

	#assign to  variables
	r, g, b = in_col

	# Convert to [0, 1]
	r, g, b = map(lambda x: x/255., [r, g, b]) 

	# RGB -> HLS
	h, l, s = colorsys.rgb_to_hls(r, g, b)     

	# Rotation by 180 degrees
	h = (h+0.5)
	color = map(lambda x: int(round(x*255)), colorsys.hls_to_rgb(h, l, s)) # H'LS -> new RGB

	if hex_color:
		colors = [rgb_to_hex(in_col), rgb_to_hex(tuple(color))]
	else:
		colors = [in_col, tuple(color)]

	return colors



def split_complementary(in_col):
	'''
	Returns input color and its split complementary colors (those adjecent to the complement) 
	as a list of of hex or rgb values, depending on what was submitted.

		x
	   O O
	  x   x
	 x     x
	  x   x
	   x x
		O


	'''
	assert (is_rgb(in_col) or is_hex(in_col)), 'Error, the input must be a hex color string or an RGB tuple.'

	hex_color = is_hex(in_col) #check whether input is hex
	if hex_color: #if it is, make it RGB
		in_col = hex_to_rgb(in_col)

	#assign to  variables
	r, g, b = in_col

	# Convert to [0, 1]
	r, g, b = map(lambda x: x/255., [r, g, b]) 

	# RGB -> HLS
	h, l, s = colorsys.rgb_to_hls(r, g, b)     

	# Rotation by 90 degrees
	angle = 150/360.0
	h_list = [(h+ang) % 1 for ang in (-angle, angle)]
	print(h)
	analagous = [map(lambda x: int(round(x*255)), colorsys.hls_to_rgb(h, l, s)) for h in h_list] # H'LS -> new RGB

	#add all the colors together
	colors = [tuple(analagous[0]), in_col, tuple(analagous[1])]

	#if the input was hex, convert it back
	if hex_color:
		colors = [rgb_to_hex(s) for s in colors]

	return colors




def triadic(in_color):
	'''
	Returns input color as well as the two triadic colors as a list of hex or rgb values, depending on what was submitted.

		x
	   x x
	  O   O
	 x     x
	  x   x
	   x x
		O

	'''
	rgb = False

	#make sure it's in hex
	if is_rgb(in_color):
		color = rgb_to_hex(in_color)
		rgb = True
	else:
		color = in_color	


	#add first triadic color
	color_list = ['#' + color.strip('#')[-2:] + color.strip('#')[:4]]

	#add the input color
	color_list.append(color)

	#add second tridadic color
	color_list.append('#' + color_list[-1].strip('#')[-2:] + color_list[-1].strip('#')[:4])

	if rgb is True:
		color_list = [hex_to_rgb(i) for i in color_list]	

	return color_list





def square(in_col):
	'''
		O
	   x x
	  x   x
	 O     O
	  x   x
	   x x
		O
	'''
	assert (is_rgb(in_col) or is_hex(in_col)), 'Error, the input must be a hex color string or an RGB tuple.'

	hex_color = is_hex(in_col) #check whether input is hex
	if hex_color: #if it is, make it RGB
		in_col = hex_to_rgb(in_col)

	#assign to  variables
	r, g, b = in_col

	# Convert to [0, 1]
	r, g, b = map(lambda x: x/255., [r, g, b]) 

	# RGB -> HLS
	h, l, s = colorsys.rgb_to_hls(r, g, b)     

	# Rotation by 90 degrees
	angle = 90/360.0
	h_list = [(h+ang) % 1 for ang in (-angle, angle, angle*2)]
	print(h)
	analagous = [map(lambda x: int(round(x*255)), colorsys.hls_to_rgb(h, l, s)) for h in h_list] # H'LS -> new RGB

	#add all the colors together
	colors = [tuple(analagous[0]), in_col, tuple(analagous[1]), tuple(analagous[2])]

	#if the input was hex, convert it back
	if hex_color:
		colors = [rgb_to_hex(s) for s in colors]

	return colors


def tetradic(in_col):
	'''


		O
	   x x
	  x   O
	 x     x
	  O   x
	   x x
		O
	'''
	assert (is_rgb(in_col) or is_hex(in_col)), 'Error, the input must be a hex color string or an RGB tuple.'

	hex_color = is_hex(in_col) #check whether input is hex
	if hex_color: #if it is, make it RGB
		in_col = hex_to_rgb(in_col)

	#assign to  variables
	r, g, b = in_col

	# Convert to [0, 1]
	r, g, b = map(lambda x: x/255., [r, g, b]) 

	# RGB -> HLS
	h, l, s = colorsys.rgb_to_hls(r, g, b)     

	# Rotation by 90 degrees
	angle = 30/360.0
	h_list = [(h+ang) % 1 for ang in (-angle*2, angle*4, angle*6)]
	print(h)
	analagous = [map(lambda x: int(round(x*255)), colorsys.hls_to_rgb(h, l, s)) for h in h_list] # H'LS -> new RGB

	#add all the colors together
	colors = [tuple(analagous[0]), in_col, tuple(analagous[1]), tuple(analagous[2])]

	#if the input was hex, convert it back
	if hex_color:
		colors = [rgb_to_hex(s) for s in colors]

	return colors


def analagous(in_col):
	'''
	Returns the input color as well as its analagous colors.

		x
	   x x
	  x   x
	 x     x
	  x   x
	   O O
		O

	'''
	assert (is_rgb(in_col) or is_hex(in_col)), 'Error, the input must be a hex color string or an RGB tuple.'

	hex_color = is_hex(in_col) #check whether input is hex
	if hex_color: #if it is, make it RGB
		in_col = hex_to_rgb(in_col)

	#assign to  variables
	r, g, b = in_col

	# Convert to [0, 1]
	r, g, b = map(lambda x: x/255., [r, g, b]) 

	# RGB -> HLS
	h, l, s = colorsys.rgb_to_hls(r, g, b)     

	# Rotation by 30 degrees
	degree = 30/360.0
	h = [(h+angle) % 1 for angle in (-degree, degree)]
	analagous = [map(lambda x: int(round(x*255)), colorsys.hls_to_rgb(hi, l, s)) for hi in h] # H'LS -> new RGB

	#add all the colors together
	colors = [tuple(analagous[0]), in_col, tuple(analagous[1])]

	#if the input was hex, convert it back
	if hex_color:
		colors = [rgb_to_hex(tuple(s)) for s in colors]

	return colors



def similar():
	'''
	Returns the input color as well as similar colors.
	'''
	pass



def monochromatic():
	'''
	Returns the input color as well as ....
	'''
	pass

def tints():
	'''
	Returns input color as well as 9 tints.
	'''
	pass


def shades():
	'''
	Returns input color as well as 9 shades.
	'''
	pass


def saturate():
	'''
	Returns the input color as well as 9 more saturated versions of that color.
	'''
	pass


def desaturate():
	'''
	Returns the input color as well as 9 less saturated versions of that color.
	'''
	pass





def visualize():
	'''
	Takes a list of colors and visualizes them.
	'''
	raise NotImplementedError


def NextRGB(color = (0,0,0)):
	'''
	Function for generating unique RGB colors. 
	The input is a tuple of an RGB color, for example (124,1,34), and the method returns the "next" color.
	When R reaches 255 one is added to G and R is reset.
	When R and G both reach 255 one is added to B and R and G are reset.
	This should generate over 1.6 million colors (255*255*255)
	'''
	assert is_rgb(color), 'Error, the input must be a tuple of three integers between 0 and 255'
	R, G, B = color

	if R == 255 and G == 255 and B == 255:
		raise ValueError, 'R, G and B all have the value 255, no further colors are available.'
	elif  R == 255 and G == 255:
		R = 0
		G = 0
		B += 1
	elif R == 255:
		R = 0
		G += 1
	else:
		R += 1
	return (R, G, B)
	
	
	
def test_scale():
	'''
	For testing the color scale, make an svg file.
	'''
	
	name = 'testing'
	scene = wsvg.Scene(name=name, size=(400, 400))

	def paintit(x, y, col1, col2, white_mid):
		col_dict = scale(col1, col2, white_mid)
		for n in range(0, 101):
			y += 0.5
			scene.add(wsvg.Rectangle(origin=(x, y), height=0.5, width=10, fill_color=col_dict[n], line_color=col_dict[n], line_width=0.5))	

	#purple to tiel
	col1 = '#800080'
	col2 = '#008080'
	x = 10
	y = 10
	paintit(x, y, col1, col2, white_mid=True)
	
	#orange to blue
	col1 = '#ffc500'
	col2 = '#056efa'
	x = 30
	y = 10
	paintit(x, y, col1, col2, white_mid=True)
	
	#blue to red
	col1 = '#817cbb'
	col2 = '#c12133'
	x = 50
	y = 10
	paintit(x, y, col1, col2, white_mid=True)
	scene.write_svg()
	
	#another blue to red
	col1 = '#30acdf'
	col2 = '#ef292b'
	x = 70
	y = 10
	paintit(x, y, col1, col2, white_mid=True)
	scene.write_svg()	

	#orange to green
	col1 = '#ff6600'
	col2 = '#2c9082'
	x = 90
	y = 10
	paintit(x, y, col1, col2, white_mid=True)
	scene.write_svg()		
	
	#orange to dark blue
	col1 = '#ff6600'
	col2 = '#18567d'
	x = 110
	y = 10
	paintit(x, y, col1, col2, white_mid=True)
	scene.write_svg()	
	
	#grey to orange
	col1 = '#666666'
	col2 = '#ff6600'
	x = 130
	y = 10
	paintit(x, y, col1, col2, white_mid=True)
	scene.write_svg()	

	#dark yellow to blue
	col1 = '#cba916'
	col2 = '#8eb6d5'
	x = 130
	y = 10
	paintit(x, y, col1, col2, white_mid=True)
	scene.write_svg()	

	
if __name__ == "__main__": 
#	test_scale()
	print(triadic('#ffa2b5'))
	print(triadic((255, 162, 181)))
	 	
	
