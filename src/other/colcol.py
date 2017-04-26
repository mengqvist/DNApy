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


#This file deals with color conversions, color transformations, and generating color scales.

import re
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


def is_hsl(in_col):
	'''
	Check whether an input is a valid HSL color.
	Return True if it is, otherwise False.
	'''
	if len(in_col) == 3 and type(in_col) == tuple:
		if 0<=in_col[0]<=360 and 0<=in_col[1]<=100 and 0<=in_col[2]<=100:
			return True
		else:
			return False
	else:
		return False


	
def rgb_to_hex(rgb):
	'''
	Convert RGB colors to hex.
	Input should be a tuple of integers (R, G, B) where each is between 0 and 255.
	Output is a string representing a hex numer. For instance '#FFFFFF'.
	''' 
	#make sure input is ok
	assert is_rgb(rgb) is True, 'Error, %s is not a valid RGB color.' % str(rgb)
	
	#make conversion
	return "#%02x%02x%02x" % rgb

	
	
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
	

def rgb_to_hsl(in_col):
	'''
	Convert RGB colors to HSL.
	Input should be a tuple of integers (R, G, B) where each is between 0 and 255.
	Output is a tuple of integers where h is between 0 and 360 and s and l are between 0 and 100. For instance: (162, 47, 39).
	'''
	assert is_rgb(in_col), 'Error, %s is not a valid RGB color.' % str(in_col)

	# Convert each RGB integer to a float between 0 and 1
	r, g, b = [x/255. for x in in_col] 

	# RGB -> HLS
	h, l, s = colorsys.rgb_to_hls(r, g, b)

	#convert it back to the appropriate integers
	h = int(round(360*h))
	l = int(round(100*l))
	s = int(round(100*s))

	return (h, s, l) 


def hsl_to_rgb(in_col):
	'''
	Convert HSL colors to RGB.
	Input should be a tuple of integers where h is between 0 and 360 and s and l are between 0 and 100. For instance: (162, 47, 39).
	Output is a tuple of integers (R, G, B) where each is between 0 and 255.
	'''
	assert is_hsl(in_col), 'Error, %s is not a valid HSL color.' % str(in_col)

	#assign to variables
	h, s, l = in_col

	# Convert each HSL integer to a float between 0 and 1
	h = h/360.0
	s = s/100.0
	l = l/100.0

	# RGB -> HLS
	r, g, b = colorsys.hls_to_rgb(h, l, s)

	#convert it back to the appropriate integers
	r = int(round(255*r))
	g = int(round(255*g))
	b = int(round(255*b))

	return (r, g, b) 	



def hex_to_hsl():
	'''
	Convert a hex color to RGB.
	Input should be a string. For example '#FFFFFF'.
	Output is a tuple of ...... . For instance: (h, s, l).
	'''
	pass


def hsl_to_hex():
	'''
	'''
	pass

def scale(col1, col2, number=101, mode='rgb', white_mid=False):
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
	r, g, b = [x/255. for x in [r, g, b]] 

	# RGB -> HLS
	h, l, s = colorsys.rgb_to_hls(r, g, b)     

	# Rotation by 180 degrees
	h = (h+0.5)
	color = [int(round(x*255)) for x in colorsys.hls_to_rgb(h, l, s)] # H'LS -> new RGB

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
	r, g, b = [x/255. for x in [r, g, b]] 

	# RGB -> HLS
	h, l, s = colorsys.rgb_to_hls(r, g, b)     

	# Rotation by 90 degrees
	angle = 150/360.0
	h_list = [(h+ang) % 1 for ang in (-angle, angle)]
	print(h)
	analagous = [[int(round(x*255)) for x in colorsys.hls_to_rgb(h, l, s)] for h in h_list] # H'LS -> new RGB

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

	#first color is wrong!

	'''
	rgb = False

	#make sure it's in hex
	if is_rgb(in_color):
		color = rgb_to_hex(in_color)
		rgb = True
	else:
		color = in_color	


	#add first triadic color
	#this one is wrong!!
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
	r, g, b = [x/255. for x in [r, g, b]] 

	# RGB -> HLS
	h, l, s = colorsys.rgb_to_hls(r, g, b)     

	# Rotation by 90 degrees
	angle = 90/360.0
	h_list = [(h+ang) % 1 for ang in (-angle, angle, angle*2)]
	print(h)
	analagous = [[int(round(x*255)) for x in colorsys.hls_to_rgb(h, l, s)] for h in h_list] # H'LS -> new RGB

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
	r, g, b = [x/255. for x in [r, g, b]] 

	# RGB -> HLS
	h, l, s = colorsys.rgb_to_hls(r, g, b)     

	# Rotation by 90 degrees
	angle = 30/360.0
	h_list = [(h+ang) % 1 for ang in (-angle*2, angle*4, angle*6)]
	print(h)
	analagous = [[int(round(x*255)) for x in colorsys.hls_to_rgb(h, l, s)] for h in h_list] # H'LS -> new RGB

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
	r, g, b = [x/255. for x in [r, g, b]] 

	# RGB -> HLS
	h, l, s = colorsys.rgb_to_hls(r, g, b)     

	# Rotation by 30 degrees
	degree = 30/360.0
	h = [(h+angle) % 1 for angle in (-degree, degree)]
	analagous = [[int(round(x*255)) for x in colorsys.hls_to_rgb(hi, l, s)] for hi in h] # H'LS -> new RGB

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


def tints(in_col, number=10):
	'''
	Returns input color as well as tints of that color (lighter colors).
	number specifies how many new ones to return.
	'''
	assert (is_rgb(in_col) or is_hex(in_col)), 'Error, the input must be a hex color string or an RGB tuple.'
	assert type(number) is int, 'Error, the input number must be an integer.'
	assert (2 <= number and number <= 1000), 'Error, the input number must be between 2 and 1000'

	#check whether input is hex, if it is, make it RGB
	hex_color = is_hex(in_col) 
	if hex_color:
		in_col = hex_to_rgb(in_col)

	#assign to  variables
	r, g, b = in_col

	# Convert to [0, 1]
	r, g, b = [x/255. for x in [r, g, b]] 

	# RGB -> HLS
	hue, lightness, saturation = colorsys.rgb_to_hls(r, g, b) 

	#what is the difference of 100% lightness and the current value	
	diff = 1.0-lightness 
	
	#devide the difference on a step size
	step = diff/float(number)

	#use that step size to generate the 10 increasing lightness values	
	lightness_list = [lightness + step*s for s in range(1, number+1)]

	#add the input color to a list, then build the 10 new HSL colors, convert to RGB and save in the same list
	colors = [in_col]
	colors.extend([[int(round(x*255)) for x in colorsys.hls_to_rgb(hue, l, saturation)] for l in lightness_list])

	#if the input was hex, convert it back
	if hex_color:
		colors = [rgb_to_hex(tuple(s)) for s in colors]

	return colors




def shades(in_col, number=10):
	'''
	Returns input color as well as shades of that color (darker colors).
	number specifies how many new ones to return.
	'''
	assert (is_rgb(in_col) or is_hex(in_col)), 'Error, the input must be a hex color string or an RGB tuple.'
	assert type(number) is int, 'Error, the input number must be an integer.'
	assert (2 <= number and number <= 1000), 'Error, the input number must be between 2 and 1000'

	#check whether input is hex, if it is, make it RGB
	hex_color = is_hex(in_col) 
	if hex_color:
		in_col = hex_to_rgb(in_col)

	#assign to  variables
	r, g, b = in_col

	# Convert to [0, 1]
	r, g, b = [x/255. for x in [r, g, b]] 

	# RGB -> HLS
	hue, lightness, saturation = colorsys.rgb_to_hls(r, g, b) 
	
	#devide the difference on a step size
	step = lightness/float(number)

	#use that step size to generate the 10 increasing lightness values	
	lightness_list = [lightness - step*s for s in range(1, number+1)]

	#add the input color to a list, then build the 10 new HSL colors, convert to RGB and save in the same list
	colors = [in_col]
	colors.extend([[int(round(x*255)) for x in colorsys.hls_to_rgb(hue, l, saturation)] for l in lightness_list])

	#if the input was hex, convert it back
	if hex_color:
		colors = [rgb_to_hex(tuple(s)) for s in colors]

	return colors


def saturate(in_col, number=10):
	'''
	Returns the input color as well as more saturated versions of that color.
	number specifies how many new ones to return.
	'''
	assert (is_rgb(in_col) or is_hex(in_col)), 'Error, the input must be a hex color string or an RGB tuple.'
	assert type(number) is int, 'Error, the input number must be an integer.'
	assert (2 <= number and number <= 1000), 'Error, the input number must be between 2 and 1000'

	#check whether input is hex, if it is, make it RGB
	hex_color = is_hex(in_col) 
	if hex_color:
		in_col = hex_to_rgb(in_col)

	#assign to  variables
	r, g, b = in_col

	# Convert to [0, 1]
	r, g, b = [x/255. for x in [r, g, b]] 

	# RGB -> HLS
	hue, lightness, saturation = colorsys.rgb_to_hls(r, g, b) 

	#what is the difference of 100% saturation and the current value	
	diff = 1.0-saturation 
	
	#devide the difference on a step size
	step = diff/float(number)

	#use that step size to generate the 10 increasing saturation values	
	saturation_list = [saturation + step*s for s in range(1, number+1)]

	#add the input color to a list, then build the 10 new HSL colors, convert to RGB and save in the same list
	colors = [in_col]
	colors.extend([[int(round(x*255)) for x in colorsys.hls_to_rgb(hue, lightness, s)] for s in saturation_list])

	#if the input was hex, convert it back
	if hex_color:
		colors = [rgb_to_hex(tuple(s)) for s in colors]

	return colors


def desaturate(in_col, number=10):
	'''
	Returns the input color as well as less saturated versions of that color.
	number specifies how many new ones to return.
	'''
	assert (is_rgb(in_col) or is_hex(in_col)), 'Error, the input must be a hex color string or an RGB tuple.'
	assert type(number) is int, 'Error, the input number must be an integer.'
	assert (2 <= number and number <= 1000), 'Error, the input number must be between 2 and 1000'

	#check whether input is hex, if it is, make it RGB
	hex_color = is_hex(in_col) 
	if hex_color:
		in_col = hex_to_rgb(in_col)

	#assign to  variables
	r, g, b = in_col

	# Convert to [0, 1]
	r, g, b = [x/255. for x in [r, g, b]] 

	# RGB -> HLS
	hue, lightness, saturation = colorsys.rgb_to_hls(r, g, b) 
	
	#devide the difference on a step size
	step = saturation/float(number)

	#use that step size to generate the 10 increasing saturation values	
	saturation_list = [saturation - step*s for s in range(1, number+1)]

	#add the input color to a list, then build the 10 new HSL colors, convert to RGB and save in the same list
	colors = [in_col]
	colors.extend([[int(round(x*255)) for x in colorsys.hls_to_rgb(hue, lightness, s)] for s in saturation_list])

	#if the input was hex, convert it back
	if hex_color:
		colors = [rgb_to_hex(tuple(s)) for s in colors]

	return colors


def continuum(in_col, number=10):
	'''
	Returns the entire color wheel starting with the input color.
	number specifies how many new ones to return.
	'''
	assert (is_rgb(in_col) or is_hex(in_col)), 'Error, the input must be a hex color string or an RGB tuple.'
	assert type(number) is int, 'Error, the input number must be an integer.'
	assert (2 <= number and number <= 1000), 'Error, the input number must be between 2 and 1000'

	#check whether input is hex, if it is, make it RGB
	hex_color = is_hex(in_col) 
	if hex_color:
		in_col = hex_to_rgb(in_col)

	#assign to  variables
	r, g, b = in_col

	# Convert to [0, 1]
	r, g, b = [x/255. for x in [r, g, b]] 

	# RGB -> HLS
	hue, lightness, saturation = colorsys.rgb_to_hls(r, g, b) 
	
	# Rotation by defined number of degrees
	angle = (360/float(number))/360.0
	h_list = [(hue+ang) % 1 for ang in [angle*s for s in range(number)]]

	colors = [[int(round(x*255)) for x in colorsys.hls_to_rgb(h, lightness, saturation)] for h in h_list] # H'LS -> new RGB
	colors = [tuple(s) for s in colors]

	#if the input was hex, convert it back
	if hex_color:
		colors = [rgb_to_hex(s) for s in colors]

	return colors


def visualize(color_list=['#acc123','#ffffff','#000000', '#1ccf9c']):
	'''
	Visualizes a list of colors.
	Useful to see what one gets out of the different functions.
	'''

	#need to adapt box size depending on how many colors are in list

	#asserts.... here....


	from tkinter import Tk, Canvas, Frame, BOTH

	#check whether input is hex, if it is, make it RGB
	rgb_color = is_rgb(color_list[0]) 
	if rgb_color:
		color_list = [rgb_to_hex(s) for s in color_list]
	 

	class Example(Frame):
	  
		def __init__(self, parent, cl):
			Frame.__init__(self, parent)   

			self.parent = parent   
			self.color_list = cl     
			self.initUI()
		    
		def initUI(self):

			self.parent.title("Colors")        
			self.pack(fill=BOTH, expand=1)

			canvas = Canvas(self)

			#modify rectangle size based on how many colors there are
			rect_size = 700/float(len(color_list))

			for i in range(len(self.color_list)):
				canvas.create_rectangle(10+rect_size*i, 10, 10+rect_size*(i+1), 110, outline=self.color_list[i], fill=self.color_list[i])
			canvas.pack(fill=BOTH, expand=1)


	def main():

		root = Tk()
		ex = Example(root, color_list)
		root.geometry("720x120+250+300")
		root.mainloop()  

	main() 


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
		raise ValueError('R, G and B all have the value 255, no further colors are available.')
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
	
	
	
#good scales

#purple to tiel
#col1 = '#800080'
#col2 = '#008080'

#orange to blue
#col1 = '#ffc500'
#col2 = '#056efa'

#blue to red
#col1 = '#817cbb'
#col2 = '#c12133'

#another blue to red
#col1 = '#30acdf'
#col2 = '#ef292b'

#orange to green
#col1 = '#ff6600'
#col2 = '#2c9082'

#orange to dark blue
#col1 = '#ff6600'
#col2 = '#18567d'

#grey to orange
#col1 = '#666666'
#col2 = '#ff6600'

#dark yellow to blue
#col1 = '#cba916'
#col2 = '#8eb6d5'

	

	 	
	
