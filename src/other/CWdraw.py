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


import wsvg

def Draw(table=1):
	'''Draw an SVG codon wheel'''
	scene = wsvg.Scene(name="Wheel", size=(400,400))
	xc = scene.size[0]/2 #centre of codon circle in x
	yc = scene.size[1]/2 #centre of codon circle in y
	Radius = min(xc, yc)/1.2


	nucleotide_color = "#28220b"
	line_color = "#8B835F"
	first_segment_color = "#ffe7ab"
	second_segment_color = "#ffd976"
	third_segment_color = "#ffc700"
	aa_segment_color = "#fff2d1"

	first_nucleotide_thickness = Radius/3
	second_nucleotide_thickness = 2*(Radius/3)/3
	third_nucleotide_thickness = 1*(Radius/3)/3
	amino_acid_thickness = Radius/2.25

	
	#draw first nucleotide
	radius = 0
	thickness = first_nucleotide_thickness
	nucleotides = ['U', 'C', 'A', 'G']
	for i in range(len(nucleotides)):
		start_angle = 0 + 90*i
		finish_angle = 90+90*i

		#draw the segment (background) for the current base
		scene.add(wsvg.DoubleArc(origin=(xc, yc), radius=radius, width=thickness, start_ang=start_angle, end_ang=finish_angle, fill_color=first_segment_color, line_color=first_segment_color, line_width=0.5))

		#draw the base
		x1, y1 = wsvg.PolarToCartesian(xc, yc, (radius+thickness)/2, finish_angle-(finish_angle-start_angle)/2)
		fontsize = int(thickness/1.5)
		scene.add(wsvg.Text(text=nucleotides[i], origin=(x1,y1+fontsize/3), angle=0, size=fontsize, weight="bold", color=nucleotide_color, anchor="middle"))


	#draw second nucleotide
	radius = first_nucleotide_thickness
	thickness = second_nucleotide_thickness
	nucleotides = ['UU', 'UC', 'UA', 'UG','CU', 'CC', 'CA', 'CG','AU', 'AC', 'AA', 'AG', 'GU', 'GC', 'GA', 'GG']
	for i in range(len(nucleotides)):
		start_angle = 0 + 22.5*i
		finish_angle = 22.5+22.5*i
		nucleotide = nucleotides[i][1]

		#draw the segment (background) for the current nucleotide
		scene.add(wsvg.DoubleArc(origin=(xc, yc), radius=radius, width=thickness, start_ang=start_angle, end_ang=finish_angle, fill_color=second_segment_color, line_color=second_segment_color, line_width=0.5))

		#draw the nucleotide
		x1, y1 = wsvg.PolarToCartesian(xc, yc, radius+thickness/2, finish_angle-(finish_angle-start_angle)/2)
		fontsize = int(thickness/1.5)
		scene.add(wsvg.Text(text=nucleotide, origin=(x1,y1+fontsize/3), angle=0, size=fontsize, weight="bold", color=nucleotide_color, anchor="middle"))



	#draw third nucleotide
	radius = first_nucleotide_thickness+second_nucleotide_thickness
	thickness = third_nucleotide_thickness
	nucleotides = ['UUU', 'UUC', 'UUA', 'UUG','UCU', 'UCC', 'UCA', 'UCG','UAU', 'UAC', 'UAA', 'UAG', 'UGU', 'UGC', 'UGA', 'UGG',\
				'CUU', 'CUC', 'CUA', 'CUG','CCU', 'CCC', 'CCA', 'CCG','CAU', 'CAC', 'CAA', 'CAG', 'CGU', 'CGC', 'CGA', 'CGG',\
				'AUU', 'AUC', 'AUA', 'AUG','ACU', 'ACC', 'ACA', 'ACG','AAU', 'AAC', 'AAA', 'AAG', 'AGU', 'AGC', 'AGA', 'AGG',\
				'GUU', 'GUC', 'GUA', 'GUG','GCU', 'GCC', 'GCA', 'GCG','GAU', 'GAC', 'GAA', 'GAG', 'GGU', 'GGC', 'GGA', 'GGG']
	for i in range(len(nucleotides)):
		start_angle = 0 + 5.625*i
		finish_angle = 5.625+5.625*i
		nucleotide = nucleotides[i][2]

		#draw the segment (background) for the current nucleotide
		scene.add(wsvg.DoubleArc(origin=(xc, yc), radius=radius, width=thickness, start_ang=start_angle, end_ang=finish_angle, fill_color=third_segment_color, line_color=third_segment_color, line_width=0.5))

		#draw the nucleotide
		x1, y1 = wsvg.PolarToCartesian(xc, yc, radius+thickness/2, finish_angle-(finish_angle-start_angle)/2)
		fontsize = int(thickness/1.5)
		scene.add(wsvg.Text(text=nucleotide, origin=(x1,y1+fontsize/3), angle=0, size=fontsize*0.9, weight="bold", color=nucleotide_color, anchor="middle"))



	#draw amino acids
	radius = first_nucleotide_thickness+second_nucleotide_thickness+third_nucleotide_thickness
	thickness = amino_acid_thickness

	AA = ['F', 'L', 'S', 'Y', 'stop', 'C', 'stop2', 'W', 'L2', 'P', 'H', 'Q', 'R', 'I', 'M', 'T', 'N', 'K', 'S2', 'R2', 'V', 'A', 'D', 'E', 'G']
	AA_width = {'F':2, 'L':2, 'S':4, 'Y':2, 'stop':2, 'C':2, 'stop2':1, 'W':1, 'L2':4, 'P':4, 'H':2, 'Q':2, 'R':4, 'I':3, 'M':1, 'T':4, 'N':2, 'K':2, 'S2':2, 'R2':2, 'V':4, 'A':4, 'D':2, 'E':2, 'G':4}
	AA_full = {'F':'Phenylalanine', 'L':'Leucine', 'S':'Serine', 'Y':'Tyrosine', 'stop':'Stop', 'C':'Cysteine', 'stop2':'Stop', 'W':'Tryptophan', 'L2':'Leucine', 'P':'Proline', 'H':'Histidine', 'Q':'Glutamine', 'R':'Arginine', 'I':'Isoleucine', 'M':'Methionine', 'T':'Threonine', 'N':'Asparagine', 'K':'Lysine', 'S2':'Serine', 'R2':'Arginine', 'V':'Valine', 'A':'Alanine', 'D':'Aspartic acid', 'E':'Glutamic acid', 'G':'Glycine'}
	finish_angle = 0
	for i in range(len(AA)):
		start_angle = finish_angle
		finish_angle = start_angle+5.625*AA_width[AA[i]]

		#draw the segment (background) for the current amino acid
		scene.add(wsvg.DoubleArc(origin=(xc, yc), radius=radius, width=thickness, start_ang=start_angle, end_ang=finish_angle, fill_color=aa_segment_color, line_color="#000000", line_width=0))
		
		#draw the amino acid
		fontsize = third_nucleotide_thickness/1.5
		text_angle = finish_angle-(finish_angle-start_angle)/2
		text_extent = fontsize*0.66
		text_radius = (first_nucleotide_thickness+second_nucleotide_thickness+third_nucleotide_thickness)*1.05
		if finish_angle <= 180:
			tx, ty = wsvg.PolarToCartesian(xc, yc, text_radius, finish_angle-(finish_angle-start_angle)/2)
			scene.add(wsvg.Text(text=AA_full[AA[i]], origin=(tx,ty+fontsize/3), angle=text_angle-90, size=fontsize*0.9, color="#000000", anchor="left"))
		else:
			tx, ty = wsvg.PolarToCartesian(xc, yc, text_radius, finish_angle-(finish_angle-start_angle)/2)
			scene.add(wsvg.Text(text=AA_full[AA[i]], origin=(tx,ty+fontsize/3), angle=text_angle+90, size=fontsize*0.9, color="#000000", anchor="end"))

		#draw lines
		angle = start_angle
		if angle in [0,90,180,270]:
			line_radius = 0
		elif angle % 22.5 == 0:
			line_radius = first_nucleotide_thickness
		elif angle % 5.625 ==0:
			line_radius = first_nucleotide_thickness+second_nucleotide_thickness

		x1, y1 = wsvg.PolarToCartesian(xc, yc, line_radius, angle)
		line_radius = first_nucleotide_thickness+second_nucleotide_thickness+third_nucleotide_thickness+amino_acid_thickness
		x2, y2 = wsvg.PolarToCartesian(xc, yc, line_radius, angle)
		scene.add(wsvg.Line(start=(x1,y1), end=(x2,y2), color=line_color, width=1))


	scene.write_svg()


if __name__ == "__main__": 
	Draw()
