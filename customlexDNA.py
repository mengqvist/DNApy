# Chapter 10: Creating Components and Extending Functionality
# Recipe 6: StyledTextCtrl Custom Highlighting
#
import wx
import wx.stc

class BaseLexer(object):
    """Defines simple interface for custom lexer objects"""
    def __init__(self):
        super(BaseLexer, self).__init__()

    def StyleText(self, event):
        raise NotImplementedError

class DNALexer(BaseLexer):
    """Simple lexer to highlight vowels"""
    # Define some style IDs
    STC_STYLE_DEFAULT = 0
    STC_STYLE_ADENINE = 1
    STC_STYLE_THYMINE = 2
    STC_STYLE_GUANINE = 3
    STC_STYLE_CYTOSINE = 4

    def __init__(self):
        super(DNALexer, self).__init__()
        
        # Attributes
        self.adenine = [ord(char) for char in "aA"]
        self.thymine = [ord(char) for char in "tT"]
        self.guanine = [ord(char) for char in "gG"]
        self.cytosine = [ord(char) for char in "cC"]

    def StyleText(self, event):
        """Handle the EVT_STC_STYLENEEDED event"""
        stc = event.GetEventObject()
        # Last correctly styled character
        last_styled_pos = stc.GetEndStyled()
        # Get styling range for this call
        line = stc.LineFromPosition(last_styled_pos)
        start_pos = stc.PositionFromLine(line)
        end_pos = event.GetPosition()
        # Walk the line and find all the vowels to style
        # Note: little inefficient doing one char at a time
        #       but just to illustrate the process.
        while start_pos < end_pos:
            stc.StartStyling(start_pos, 0x1f)
            char = stc.GetCharAt(start_pos)
            if char in self.adenine:
                style = DNALexer.STC_STYLE_ADENINE
            elif char in self.thymine:
                style = DNALexer.STC_STYLE_THYMINE
            elif char in self.guanine:
                style = DNALexer.STC_STYLE_GUANINE
            elif char in self.cytosine:
                style = DNALexer.STC_STYLE_CYTOSINE
            else:
                style = DNALexer.STC_STYLE_DEFAULT	
            # Set the styling byte information for 1 char from
            # current styling position (start_pos) with the
            # given style.
            length = 1
            stc.SetStyling(length, style)
            start_pos += 1

class CustomSTC(wx.stc.StyledTextCtrl):
    def __init__(self, *args, **kwargs):
        super(CustomSTC, self).__init__(*args, **kwargs)

        # Attributes
        self.custlex = None

        # Event Handlers
        self.Bind(wx.stc.EVT_STC_STYLENEEDED, self.OnStyle)

 
    def OnStyle(self, event):
        # Delegate to custom lexer object if one exists
        if self.custlex:
            self.custlex.StyleText(event)
        else:
            event.Skip()

    def SetLexer(self, lexerid, lexer=None):
        """Overrides StyledTextCtrl.SetLexer
        Adds optional param to pass in custom container
        lexer object.
        """
        self.custlex = lexer
        super(CustomSTC, self).SetLexer(lexerid)



#---- End Recipe Code ----#

class StyledTextApp(wx.App):
    def OnInit(self):
        self.frame = StcFrame(None, title="Custom Lexer")
        self.frame.Show()
        return True




class StcFrame(wx.Frame):
    """Main application window"""
    def __init__(self, parent, *args, **kwargs):
        super(StcFrame, self).__init__(parent,
                                       *args,
                                       **kwargs)

        # Attributes
        self.stc = CustomSTC(self)
        self.stc.SetLayoutCache(wx.stc.STC_CACHE_DOCUMENT) #cache laout calculations and only draw when something changes
        self.stc.SetWrapMode(wx.stc.STC_WRAP_WORD) #enable word wrap
        self.stc.SetText('''aaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttcaaaaaaaaacacacaacccccccccccccccccccccccccccccccccccccccgggggccttttc''')

        #set whether letters or background should be colored
        self.DNAmode = 'face' 

        # Setup STC for DNA highlighting
        if self.DNAmode == 'face':
			#colored letters
		    style = DNALexer.STC_STYLE_DEFAULT
		    self.stc.StyleSetSpec(style, "fore:#FF00FF,face:Mono")
		    style = DNALexer.STC_STYLE_ADENINE
		    self.stc.StyleSetSpec(style, "fore:#00CC33,face:Mono") #green
		    style = DNALexer.STC_STYLE_THYMINE
		    self.stc.StyleSetSpec(style, "fore:#FF0000,face:Mono") #red
		    style = DNALexer.STC_STYLE_CYTOSINE
		    self.stc.StyleSetSpec(style, "fore:#0033FF,face:Mono") #blue
		    style = DNALexer.STC_STYLE_GUANINE
		    self.stc.StyleSetSpec(style, "fore:#000000,face:Mono") #black
        elif self.DNAmode == 'back':
			#colored back
		    style = DNALexer.STC_STYLE_DEFAULT
		    self.stc.StyleSetSpec(style, "fore:#000000,back:#00CC33,face:Mono")
		    style = DNALexer.STC_STYLE_ADENINE
		    self.stc.StyleSetSpec(style, "fore:#000000,back:#00CC33,face:Mono")
		    style = DNALexer.STC_STYLE_THYMINE
		    self.stc.StyleSetSpec(style, "fore:#000000,back:#FF6666,face:Mono")
		    style = DNALexer.STC_STYLE_CYTOSINE
		    self.stc.StyleSetSpec(style, "fore:#000000,back:#3366FF,face:Mono")
		    style = DNALexer.STC_STYLE_GUANINE
		    self.stc.StyleSetSpec(style, "fore:#000000,back:#666666,face:Mono")


        self.stc.SetLexer(wx.stc.STC_LEX_CONTAINER, DNALexer())

        # Layout
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.stc, 1, wx.EXPAND)
        self.SetSizer(sizer)
        self.SetInitialSize((300, 300))

if __name__ == '__main__':
    app = StyledTextApp(False)
    app.MainLoop()
