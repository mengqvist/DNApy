#!/usr/bin/env python
"""\
wsvg.py (Write SVG) - Construct/display SVG scenes.

The following code is a lightweight wrapper around SVG files. The metaphor
is to construct a scene, add objects to it, and then write it to a file
to display it.

This program uses ImageMagick to display the SVG files. ImageMagick also 
does a remarkable job of converting SVG files into other formats.

Isendrak Skatasmid (http://code.activestate.com/recipes/578123-draw-svg-images-in-python-python-recipe-enhanced-v/) created an enhanced Version of Rick Muller's Code from http://code.activestate.com/recipes/325823-draw-svg-images-in-python/
This was in turn enhanced by Martin Engqvist to contain code for generating arcs and double arcs.
"""

#import os
import math
#display_prog = "display"

class Scene:
    def __init__(self, name="svg", size=(400,400)):
        self.name = name
        self.items = []
        self.height = size[0]
        self.width = size[1]
        return

    def add(self,item): 
        self.items.append(item)

    def strarray(self):
        var = ["<?xml version=\"1.0\"?>\n",
               "<svg height=\"%r\" width=\"%r\" >\n" % (self.height,self.width),
               " <g>\n"]
        for item in self.items: var += item.strarray()            
        var += [" </g>\n</svg>\n"]
        return var

    def write_svg(self,filename=None):
        if filename:
            self.svgname = filename
        else:
            self.svgname = self.name + ".svg"
        file = open(self.svgname,'w')
        file.writelines(self.strarray())
        file.close()
        return

#    def display(self,prog=display_prog):
#        os.system("%s %s" % (prog,self.svgname))
#        return        

class Line:
    def __init__(self,start,end,color,width):
        self.start = start
        self.end = end
        self.color = color
        self.width = width
        return

    def strarray(self):
        return ["""  <line 
    x1=\"%r\" 
    y1=\"%r\" 
    x2=\"%r\" 
    y2=\"%r\" 
    style=\"stroke:%s;stroke-width:%r\"/>\n""" %\
                (self.start[0],self.start[1],self.end[0],self.end[1],colorstr(self.color),self.width)]

class Circle:
    def __init__(self,center,radius,fill_color,line_color,line_width):
        self.center = center
        self.radius = radius
        self.fill_color = fill_color
        self.line_color = line_color
        self.line_width = line_width
        return

    def strarray(self):
        return ["""  <circle 
    cx=\"%r\" 
    cy=\"%r\" 
    r=\"%r\"\n""" %\
                (self.center[0],self.center[1],self.radius),
                "    style=\"fill:%s;stroke:%s;stroke-width:%r\"/>\n" % (colorstr(self.fill_color),colorstr(self.line_color),self.line_width)]

class Ellipse:
    def __init__(self, center, radius_x, radius_y, fill_color, line_color, line_width):
        self.center = center
        self.radius_x = radius_x
        self.radius_y = radius_y
        self.fill_color = fill_color
        self.line_color = line_color
        self.line_width = line_width
    def strarray(self):
        return ["""  <ellipse 
    cx=\"%r\" 
    cy=\"%r\" 
    rx=\"%r\" 
    ry=\"%r\"\n""" %\
                (self.center[0],self.center[1],self.radius_x,self.radius_y),
                "    style=\"fill:%s;stroke:%s;stroke-width:%r\"/>\n" % (colorstr(self.fill_color),colorstr(self.line_color),self.line_width)]

class Polygon:
    def __init__(self, points, fill_color, line_color, line_width):
        self.points = points
        self.fill_color = fill_color
        self.line_color = line_color
        self.line_width = line_width
    def strarray(self):
        polygon=""
        for point in self.points:
            polygon+="%r,%r " % (point[0],point[1])
        return ["  <polygon\n",
    			"    points=\"%s\"\n" % str(polygon),
				"    style=\"fill:%s;stroke:%s;stroke-width:%r\"/>\n" % (colorstr(self.fill_color),colorstr(self.line_color),self.line_width)]

class Rectangle:
    def __init__(self,origin,height,width,fill_color,line_color,line_width):
        self.origin = origin
        self.height = height
        self.width = width
        self.fill_color = fill_color
        self.line_color = line_color
        self.line_width = line_width
        return

    def strarray(self):
        return ["""  <rect 
    x=\"%r\" 
    y=\"%r\" 
    height=\"%r\"\n""" %\
                (self.origin[0],self.origin[1],self.height),
                "    width=\"%r\" style=\"fill:%s;stroke:%s;stroke-width:%r\"/>\n" %\
                (self.width,colorstr(self.fill_color),colorstr(self.line_color),self.line_width)]

class Arc:
	def __init__(self, origin, radius, start_ang, end_ang, fill_color, line_color, line_width):
		self.origin = origin
		self.radius = radius
		self.start_ang = start_ang
		self.end_ang = end_ang
		self.fill_color = fill_color
		self.line_color = line_color
		self.line_width = line_width
		return

	def makeArc(self, radius):
		'''Draws an arc in the forward (clockwise) direction'''
		radius = radius
		angle = self.start_ang
		path = ''
		for i in range(90, 450, 90):
			if self.start_ang + i <= self.end_ang:
				ang_increment = 90
			elif self.start_ang + i > self.end_ang:
				ang_increment = self.end_ang - angle
			angle += ang_increment

			#empirically the bezier length is about 0.35 of the arc length (an arc is maximally 1/4 of a circle)
			#more accurately it seems to be (1.104568/pi)*a	    where a is the arc length
			bezier_length =  (1.104568/math.pi)*(radius*2*math.pi)*(ang_increment/360.0) 

			#I need to calculate where the beziers go. These are basically defined length tangents of a circle 
			#the radius of the circle is adjacent
			#the bezier length is opposite	
			#z is the radius of a new circle on which the bezier point lies
			z = math.sqrt(math.pow(radius,2) + math.pow(bezier_length,2)) #determine length of z (the hypotenuse)

			#now I need the angle of so that I can calculate the x,y coordinate of the point that the hypotenuse leads to
			tanangle = bezier_length/radius #tan(angle)
			radians = math.atan(tanangle) #negate the Tan part and get radians
			degrees = radians*(180/math.pi)	#convert radians to degrees

			#I'm calculating the internal angle of the right triangle used to find the bezier point, I needs to do the real (non-bezier) angle minus this to get the real angle
			x, y = PolarToCartesian(self.origin[0], self.origin[1], radius, angle) #next point
			x1, y1 = PolarToCartesian(self.origin[0], self.origin[1], z, angle-ang_increment+degrees) #bezier on previous point
			x2, y2 = PolarToCartesian(self.origin[0], self.origin[1], z, angle-degrees) #bezier on previous point
			path += "%r,%r %r,%r %r,%r " % (x1, y1, x2, y2, x, y) #c is curve in the form (x1 y1 x2 y2 x y). Draws line from current point to xy. Uses x1 y1 and x2 y2 as control points for beziers.
		return path	


	def makeRevArc(self, radius):
		'''Draws an arc in the reverse (anti-clockwise) direction'''
		radius = radius
		angle = self.end_ang
		path = ''
		while angle != self.start_ang:
			if (angle-self.start_ang) % 90 != 0: #if the span is not evenly devisible by 90
				ang_increment = -((self.end_ang-self.start_ang) % 90)
			else:
				ang_increment = -90
			angle += ang_increment

			#empirically the bezier length is about 0.35 of the arc length (an arc is maximally 1/4 of a circle)
			#more accurately it seems to be (1.104568/pi)*a	    where a is the arc length
			bezier_length =  (1.104568/math.pi)*(radius*2*math.pi)*(ang_increment/360.0) 

			#I need to calculate where the beziers go. These are basically defined length tangents of a circle 
			#the radius of the circle is adjacent
			#the bezier length is opposite	
			#z is the radius of a new circle on which the bezier point lies
			z = math.sqrt(math.pow(radius,2) + math.pow(bezier_length,2)) #determine length of z (the hypotenuse)

			#now I need the angle of so that I can calculate the x,y coordinate of the point that the hypotenuse leads to
			tanangle = bezier_length/radius #tan(angle)
			radians = math.atan(tanangle) #negate the Tan part and get radians
			degrees = radians*(180/math.pi)	#convert radians to degrees

			#I'm calculating the internal angle of the right triangle used to find the bezier point, I needs to do the real (non-bezier) angle minus this to get the real angle
			x, y = PolarToCartesian(self.origin[0], self.origin[1], radius, angle) #next point
			x1, y1 = PolarToCartesian(self.origin[0], self.origin[1], z, angle-ang_increment+degrees) #bezier on previous point
			x2, y2 = PolarToCartesian(self.origin[0], self.origin[1], z, angle-degrees) #bezier on previous point
			path += "%r,%r %r,%r %r,%r " % (x1, y1, x2, y2, x, y) #c is curve in the form (x1 y1 x2 y2 x y). Draws line from current point to xy. Uses x1 y1 and x2 y2 as control points for beziers.
		return path			

	def strarray(self):
		arc_path = '    d="' #for storing the points when drawing the arc
		
		angle = self.end_ang
		x, y = PolarToCartesian(self.origin[0], self.origin[1], self.radius, angle)
		arc_path += "M %r,%r " % (x, y) #m is moveto, basically lift pen and start as designated point
		arc_path += "C " #C for 
		arc_path += self.makeArc(self.radius) #make the actual points and beziers

		if self.start_ang == self.end_ang or self.start_ang == self.end_ang-360: #complete circle
			arc_path += ' z "\n'
		else:
			arc_path += '"\n'
		

		return ["""  <path 
    cx=\"%r\" 
    cy=\"%r\"\n""" %\
			(self.origin[0], self.origin[1]),
			arc_path,
			"    style=\"fill:%s;stroke:%s;stroke-width:%r\"/>\n" % (colorstr(self.fill_color), colorstr(self.line_color), self.line_width)]



class DoubleArc(Arc):
	def __init__(self, origin, radius, width, start_ang, end_ang, fill_color, line_color, line_width):
		self.origin = origin
		self.radius = radius
		self.width = width
		self.start_ang = start_ang
		self.end_ang = end_ang
		self.fill_color = fill_color
		self.line_color = line_color
		self.line_width = line_width
		return



	def strarray(self):
		arc_path = '    d="' #for storing the points when drawing the arc
		angle = self.start_ang

		if self.start_ang == self.end_ang or self.start_ang == self.end_ang-360: #complete circle
			x, y = PolarToCartesian(self.origin[0], self.origin[1], self.radius, angle)
			arc_path += "M %r,%r " % (x, y) #move pen to start of inner arc
			arc_path += "C " #C for curveto
			inner_arc = self.makeArc(self.radius) #make the actual points and beziers
			arc_path += inner_arc #add the inner arc
			arc_path += 'z ' #close it

			x, y = PolarToCartesian(self.origin[0], self.origin[1], self.radius+self.width, angle)
			arc_path += "M %r,%r " % (x, y) #move pen to start of outer arc
			arc_path += "C " #C for curveto
			outer_arc = self.makeRevArc(self.radius+self.width) #make the actual points and beziers
			arc_path += outer_arc  #add outer arc
			arc_path += 'z"\n' #close it

		else: #partial circle
			x, y = PolarToCartesian(self.origin[0], self.origin[1], self.radius+self.width, angle)
			arc_path += "M %r,%r " % (x, y) #move pen to start of outer arc

			if self.radius == 0:
				x, y = self.origin
				arc_path += "L %r,%r " % (x, y) #draw line to start of inner arc
			elif self.radius != 0:
				x, y = PolarToCartesian(self.origin[0], self.origin[1], self.radius, angle)
				arc_path += "L %r,%r " % (x, y) #draw line to start of inner arc
				inner_arc = "C " + self.makeArc(self.radius) #make the actual points and beziers
				arc_path += inner_arc #add the inner arc

			x, y = PolarToCartesian(self.origin[0], self.origin[1], self.radius+self.width, self.end_ang)
			arc_path += "L %r,%r " % (x, y) #draw line to end of outer arc
			
			outer_arc = "C " + self.makeRevArc(self.radius+self.width) #draw the points and beziers of the outer arc in reverse
			arc_path += outer_arc  #add outer arc
			arc_path += 'z"\n' #close it
			
		return ["""  <path 
    cx=\"%r\" 
    cy=\"%r\"\n""" %\
			(self.origin[0], self.origin[1]),
			arc_path,
			"    style=\"fill:%s;stroke:%s;stroke-width:%r\"/>\n" % (colorstr(self.fill_color), colorstr(self.line_color), self.line_width)]

		


class Text:
    def __init__(self,origin,text,size,color):
        self.origin = origin
        self.text = text
        self.size = size
        self.color = color
        return

    def strarray(self):
        return ["""  <text 
    x=\"%r\" 
    y=\"%r\" 
    font-size=\"%r\" 
    fill=\"%s\">\n""" %\
                (self.origin[0],self.origin[1],self.size,colorstr(self.color)),
                "   %s\n" % self.text,
                "  </text>\n"]

def colorstr(rgb):
	'''Convert RGB colors to hex''' 
	R, G, B = rgb
	assert 0<=R<=255
	assert 0<=G<=255
	assert 0<=B<=255
	return "#%02x%02x%02x" % (R,G,B)

def PolarToCartesian(centre_x, centre_y, radius, angle):	
	'''Takes the centre of a circle, an angle (in degrees) and a radius and calculates the correspoinding XY coordinate on the circle'''
	assert type(centre_x) is int and type(centre_y) is int, 'Error, the input coordinates need to be integers.'
	assert (type(angle) is int or type(angle) is float) and (type(radius) is int or type(radius) is float), 'Error, the input angle and radius need to be integers or floats.'
	assert 0<=angle<=360, 'Error, the input angle is %r, but needs to be between 0 and 360.' % angle
	x = centre_x + radius * math.cos((angle-90)*(math.pi/180))
	y = centre_y + radius * math.sin((angle-90)*(math.pi/180))
	return (x, y)


def test():
	'''Example usage'''
	scene = Scene("test", size=(400,400))
	scene.add(Rectangle((100,100),200,200,(0,255,255),(0,0,0),1))
	scene.add(Polygon([(5,5),(5,10),(15,10),(20,20)], (0,255,255),(0,0,0),1))
	scene.add(Line((200,200),(200,300),(0,0,0),1))
	scene.add(Circle((200,200),30,(0,0,255),(0,0,0),1))
	scene.add(Text((50,50),"Testing SVG",24,(0,0,0)))
	scene.add(Ellipse((300,300), 40, 25, (255,0,255),(0,255,0),1))
	scene.add(Arc(origin=(10,10), radius=50, start_ang=0, end_ang=90, fill_color=(0,0,0), line_color=(255,0,0), line_width=3))
	scene.add(Arc(origin=(100,100), radius=15, start_ang=1, end_ang=89, fill_color=(0,0,0), line_color=(255,0,0), line_width=3))
	scene.add(DoubleArc(origin=(100,200), radius=0, width=14, start_ang=0, end_ang=90, fill_color=(0,0,0), line_color=(255,0,0), line_width=0))
	scene.add(DoubleArc(origin=(100,300), radius=15, width=14, start_ang=190, end_ang=200, fill_color=(0,0,0), line_color=(255,0,0), line_width=0))
	scene.add(DoubleArc(origin=(200,100), radius=15, width=5, start_ang=0, end_ang=360, fill_color=(0,0,0), line_color=(255,0,0), line_width=0))
	scene.write_svg()
#	scene.display()
	return

if __name__ == "__main__": 
	test()
