from buffer_test import BufferedWindow


import time

import wx
import wx.lib.graphics

import math

class Test(BufferedWindow):
	def __init__(self, *args, **kwargs):	
		self.x = 0
		self.y = 0

		self.selection_down = 36
		self.selection_up = 49

		self.centre_x = 0
		self.centre_y = 0
		BufferedWindow.__init__(self, *args, **kwargs)
	


#		wx.EVT_PAINT(self, self.OnPaint)
#		wx.EVT_SIZE(self, self.OnSize)
		wx.EVT_LEFT_DOWN(self, self.OnLeftDown)
		wx.EVT_LEFT_UP(self, self.OnLeftUp)
		wx.EVT_RIGHT_UP(self, self.OnRightUp)
		wx.EVT_MOTION(self, self.OnMotion)

#		self.Centre()
#		self.Show(True)



	def Draw(self, dc):
		# make a path that contains a circle and some lines
		self.centre_x = self.size[0]/2 #centre of window in x
		self.centre_y = self.size[1]/2 #centro of window in y

		dc.SetBackground(wx.Brush("White") )
		dc.Clear() # make sure you clear the bitmap!
		Pen = wx.Pen(colour='#444444', width=5)
		dc.SetPen(Pen)
		dc.SetBrush(wx.Brush("White"))
#		path = dc.CreatePath()
#		path.AddCircle(50.0, 50.0, 50.0)
#		path.MoveToPoint(0.0, 50.0)
#		path.AddLineToPoint(100.0, 50.0)
#		path.MoveToPoint(50.0, 0.0)
#		path.AddLineToPoint(50.0, 100.0)
#		path.CloseSubpath()
#		path.AddRectangle(25.0, 25.0, 50.0, 50.0)
#		dc.StrokePath(path)

#		dc.SetPen(wx.BLACK_PEN)
#		path = dc.CreatePath()
#		path.AddCircle(self.x, self.y, 50.0)
#		dc.StrokePath(path)

		#DNA circle
		if self.centre_x > self.centre_y:
			self.Radius = self.centre_y/2
		else:
			self.Radius = self.centre_x/2
		dc.DrawCircle(x=self.centre_x, y=self.centre_y, radius=self.Radius) #DNA circle
		#features

		#selection
		dc.SetBrush(wx.Brush("Orange"))
		if self.selection_down != None:
			xc=self.centre_x
			yc=self.centre_y

			x1 = xc + self.Radius * math.cos((self.selection_up-90)*(math.pi/180)) #the latter needs to be first as the arc draws backwards
			y1 = yc + self.Radius * math.sin((self.selection_up-90)*(math.pi/180))
			x2 = xc + self.Radius * math.cos((self.selection_down-90)*(math.pi/180))
			y2 = yc + self.Radius * math.sin((self.selection_down-90)*(math.pi/180))


			dc.DrawArc(x1, y1, x2, y2, xc, yc)
			
		#enzymes
		#names

	def calc_location(self):
		'''Calculates the DNA location corresponding to mouse clicks (as a fraction float)'''
		z = math.sqrt(math.pow((self.x-self.centre_x),2) + math.pow((self.y-self.centre_y),2)) #determine z
#		cosangle = (self.x-centre_x)/z
#		cos_radians = math.acos(cosangle)
#		degrees = cos_radians*(180/math.pi)
#		print('cos', degrees)

		sinangle = (self.y-self.centre_y)/z 	
		sin_radians = math.asin(sinangle)
		degrees = sin_radians*(180/math.pi)	
		#now I have the angle of the triangle. Based on where it is placed I have to calculate the 'real' angle

		#for these calculations it is important to remember that y increases as you go down...
		x_difference = self.x-self.centre_x
		y_difference = self.y-self.centre_y
		if x_difference > 0 and y_difference > 0: #triangle is in bottom right of circle
			angle = 90+degrees
#			print('bottom right')
		elif  x_difference < 0 and y_difference < 0: #triangle is in top left of circle
			angle = 180+90-degrees
#			print('top left')
		elif x_difference > 0 and y_difference < 0: #triangle is in in top right of circle
			angle = 90+degrees
#			print('top right')
		elif x_difference < 0 and y_difference > 0: #triangle is in bottom left of circle
			angle = 180+90-degrees
#			print('bottom left')
		print('angle', angle)

		#now calculate what percentage of the circle that is
		#length = (diameter*pi*angle)/360
#		length = (self.Radius*2*math.pi*angle)/360
#		percent_length = length/(self.Radius*2*math.pi)
#		print('percent', percent_length)
		return angle


	def OnLeftDown(self, event):
		self.centre_x = self.size[0]/2 #centre of window in x
		self.centre_y = self.size[1]/2 #centro of window in y
		x, y = self.ScreenToClient(wx.GetMousePosition())	
		self.x = x
		self.y = y

		angle = self.calc_location()
		self.selection_down = angle
		self.selection_up = angle +1
		self.UpdateDrawing()

	def OnLeftUp(self, event):
		self.centre_x = self.size[0]/2 #centre of window in x
		self.centre_y = self.size[1]/2 #centro of window in y
		x, y = self.ScreenToClient(wx.GetMousePosition())	
		self.x = x
		self.y = y

		angle = self.calc_location()
		if (self.selection_down-1)<=angle<=(self.selection_down+1):
			self.selection_up = self.selection_down +1
		else:
			self.selection_up = angle
		self.UpdateDrawing()




	def OnMotion(self, event):
	#	time.sleep(0.01)
		if event.Dragging() and event.LeftIsDown():
			x, y = self.ScreenToClient(wx.GetMousePosition())	
			self.x = x
			self.y = y

			angle = self.calc_location()
			self.selection_up = angle

			self.UpdateDrawing()

	def OnRightUp(self, event):
		self.Close()



class TestFrame(wx.Frame):
	def __init__(self, parent=None):
		wx.Frame.__init__(self, parent, size = (500,500), title="Double Buffered Test", style=wx.DEFAULT_FRAME_STYLE)

		self.Window = Test(self)
		self.Show()



class DemoApp(wx.App):
    def OnInit(self):
        frame = TestFrame()
        self.SetTopWindow(frame)

        return True

if __name__ == "__main__":
    app = DemoApp(0)
    app.MainLoop()
