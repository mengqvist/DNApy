
	def Cut(self, start, finish, complement=False):
		'''Cuts selection and place it in clipboard together with any features present on that DNA'''
		assert (type(start) == int and type(finish) == int), 'Function requires two integers.'
		assert start <= finish, 'Startingpoint must be before finish'
		self.RichCopy(start, finish, complement) # changed to RichCopy
		deletedsequence = self.GetDNA(start, finish)
		self.changegbsequence(start, finish, 'd', deletedsequence)
		self.add_file_version()

	def CutRC(self, start, finish):
		'''Cuts the reverese-complement of a selection and place it in clipboard together with any features present on that DNA'''
		assert (type(start) == int and type(finish) == int), 'Function requires two integers.'
		assert start <= finish, 'Startingpoint must be before finish'
		self.Cut(start, finish, True)


	def Paste(self, ip, DNA=None):
		'''Paste DNA present in clipboard and any features present on that DNA
			If a string is passed to DNA then this will over-ride anything that is present in the clipboard'''
		assert type(ip) == int, 'The insertion point must be an integer.'

		if DNA == None:
			system_clip = re.sub(r'\s+', '', pyperclip.paste()) #remove whitespace from item in system clipboard
			if self.clipboard['dna'] != system_clip: #if internal and system clipboard is not same, then system clipboard takes presidence
				self.clipboard['dna'] = system_clip
				self.clipboard['features'] = []
			DNA = copy.copy(self.clipboard['dna'])
			self.changegbsequence(ip, ip, 'i', DNA) #change dna sequence
			for i in range(len(self.clipboard['features'])): #add features from clipboard
				self.paste_feature(self.clipboard['features'][i], ip-1)
		else:
			self.changegbsequence(ip, ip, 'i', DNA) #change dna sequence
		self.add_file_version()


	def PasteRC(self, ip):
		'''Paste reverse complement of DNA in clipboard'''
		assert type(ip) == int, 'The insertion point must be an integer.'
		self.RichPaste(ip,True)

	def CopyRC(self, start, finish):
		'''Copy the reverse complement of DNA and all the features for a certain selection'''
		assert (type(start) == int and type(finish) == int), 'Function requires two integers.'
		assert start <= finish, 'Startingpoint must be before finish'
		self.RichCopy(start, finish, True)


	def RichCopy(self, start, finish, complement=False):
		'''a method to copy not only dna to clipboard but exchange features and dna
		   beetween two windows
		   this depreciates the Copy method from this file
		'''
		assert (type(start) == int and type(finish) == int), 'Function requires two integers.'

		dna = self.GetDNA(start, finish)

		dnapyClipboard             = {}
		dnapyClipboard['dna']      = self.GetDNA(start, finish) #copy to internal clipboard
		dnapyClipboard['features'] = []
		self.allgbfeatures_templist = copy.deepcopy(self.gbfile['features'])
		for i in range(len(self.gbfile['features'])): #checks to match dna change
			if len(self.gbfile['features'][i]['location']) == 1:
				featurestart, featurefinish = self.get_location(self.gbfile['features'][i]['location'][0])
			else:
				n = 0
				featurestart = self.get_location(self.gbfile['features'][i]['location'][n])[0]
				n = len(self.gbfile['features'][i]['location'])-1
				featurefinish = self.get_location(self.gbfile['features'][i]['location'][n])[1]

			if (start<=featurestart and featurefinish<=finish) == True: #if change encompasses whole feature
				dnapyClipboard['features'].append(self.allgbfeatures_templist[i])
				for n in range(len(self.gbfile['features'][i]['location'])):
					newlocation = self.add_or_subtract_to_locations(self.gbfile['features'][i]['location'][n], -start+1, 'b')
					dnapyClipboard['features'][-1]['location'][n] = newlocation

		# if complement copy is hoped for, reverse it here:
		if complement == True:
			dnapyClipboard = self.reverse_complement_clipboard(dnapyClipboard)

		# save simple dna as txt
		self.richClipboard = wx.DataObjectComposite()
		text               = wx.TextDataObject(dnapyClipboard['dna'])

		# save clipboard to json string
		JSONdnapyClipboard     = json.dumps(dnapyClipboard)



		dnapy              = wx.CustomDataObject("application/DNApy")
		dnapy.SetData(JSONdnapyClipboard)

		self.richClipboard.Add(text, True)	# add clear dna text as preffered object
		self.richClipboard.Add(dnapy)		# add rich dna info as json



		# save to the clipboard
		if wx.TheClipboard.Open():
			wx.TheClipboard.SetData(self.richClipboard)
			wx.TheClipboard.Close()


		return True



	def RichPaste(self, ip, complement=False):
		'''Paste DNA present in clipboard and any features present on that DNA
			If a string is passed to DNA then this will over-ride anything that is present in the clipboard
			this depreciates the Paste method from this file
			'''
		assert type(ip) == int, 'The insertion point must be an integer.'



		self.textClipboard      = wx.TextDataObject()
		self.dnaClipboard       = wx.CustomDataObject("application/DNApy")

		# retrive the clipboard content
		if wx.TheClipboard.Open():
			containsDNApy  = wx.TheClipboard.GetData(self.dnaClipboard)
			containsText   = wx.TheClipboard.GetData(self.textClipboard)
			wx.TheClipboard.Close()

			if containsDNApy:
				pasteRich = True
				# get the json object
				pasteDNA  = json.loads(self.dnaClipboard.GetData())

			elif containsText:
				pasteRich = False
				# paste the pure text
				paste     = self.textClipboard.GetDataHere()
			else:
				# else there is nothing to paste
				pasteRich = False
				paste     = ''

		# if complement paste is wished we need to reverse somethings
		if complement == True and pasteRich == True:
			pasteDNA = self.reverse_complement_clipboard(pasteDNA) 	# reverse the clipboard content
		elif complement == True and pasteRich == False:
			paste = dna.reverse_complement(paste) 									# reverse the dna txt

		# here we paste
		if pasteRich:
			DNA = pasteDNA['dna']
			self.changegbsequence(ip, ip, 'i', DNA) #change dna sequence
			for i in range(len(pasteDNA['features'])): #add features from clipboard
				self.paste_feature(pasteDNA['features'][i], ip-1)
		else:
			DNA = paste
			self.changegbsequence(ip, ip, 'i', DNA) #change dna sequence
		self.add_file_version()


		return True



	def reverse_complement_clipboard(self, clipboard):
		'''Reverse-complements the DNA and all features in clipboard'''
		clipboard['dna'] = dna.reverse_complement(clipboard['dna']) #change dna sequence

		for i in range(len(clipboard['features'])): #checks self.allgbfeatures to match dna change
			if clipboard['features'] [i]['complement'] == True:
				clipboard['features'] [i]['complement'] = False
			elif clipboard['features'] [i]['complement'] == False:
				clipboard['features'] [i]['complement'] = True

			for n in range(len(clipboard['features'][i]['location'])):
				start, finish = self.get_location(clipboard['features'][i]['location'][n])
				featurelength = finish - start    #how much sequence in feature?
				trail = len(clipboard['dna']) - finish     # how much sequence after feature?

				clipboard['features'][i]['location'][n] = self.add_or_subtract_to_locations(clipboard['features'][i]['location'][n], -finish+(trail+featurelength+1), 'f')
				clipboard['features'][i]['location'][n] = self.add_or_subtract_to_locations(clipboard['features'][i]['location'][n], -start+trail+1, 's')
			clipboard['features'][i]['location'].reverse() #reverse order of list elements

		return clipboard


	def RCselection(self, start, finish):
		'''Reverse-complements current DNA selection'''
		assert (type(start) == int and type(finish) == int), 'Function requires two integers.'
		assert start <= finish, 'Startingpoint must be before finish'
		self.RichCopy(start, finish, True) # RichCopy
		#self.reverse_complement_clipboard()
		self.Delete(start, finish, visible=False)
		self.RichPaste(start)


	def Upper(self, start=1, finish=-1):
		'''Change DNA selection to uppercase characters.
			Start and finish are optional arguments (integers). If left out change applies to the entire DNA molecule.
			If specified, start and finish determines the range for the change.'''
		if finish == -1:
			finish = len(self.GetDNA())
		assert (type(start) == int and type(finish) == int), 'Function requires two integers.'
		assert start <= finish, 'Startingpoint must be before finish'
		string = self.GetDNA(start, finish)
		self.changegbsequence(start, finish, 'r', string.upper())
		self.add_file_version()

	def Lower(self, start=1, finish=-1):
		'''Change DNA selection to lowercase characters.
			Start and finish are optional arguments (integers). If left out change applies to the entire DNA molecule.
			If specified, start and finish determines the range for the change.'''
		if finish == -1:
			finish = len(self.GetDNA())
		assert (type(start) == int and type(finish) == int), 'Function requires two integers.'
		assert start <= finish, 'Startingpoint must be before finish'
		string = self.GetDNA(start, finish)
		self.changegbsequence(start, finish, 'r', string.lower())
		self.add_file_version()

	def delete(self, start, finish, visible=True):
		'''Deletes current DNA selection.
			Start and finish should be integers.
			The optional variable 'hidden' can be set to True or False.
			If set to True, no file versions are added to the undo/redo record.
			If set to False, it does add file versions to the undo/redo record.'''
		assert (type(start) == int and type(finish) == int), 'Function requires two integers.'
		assert start <= finish, 'Startingpoint must be before finish'
		deletedsequence = self.GetDNA(start, finish)
		self.changegbsequence(start, finish, 'd', deletedsequence)
		if visible == True:
			self.add_file_version()
