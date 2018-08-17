#!/usr/bin/env python3


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


from dnapy.resources import bioseq, fasta
from os.path import exists

# Data format expeced as below (without double ## and tab-separated)
# number is total number of that codon-optimized
# fraction is the percent the codons for each amino acid are used
# /1000 is the average usage of that codon per 1000 codons

## #organism=escherichia_coli
## #translation_table=1
## triplet  amino_acid  number  fraction  /1000
## UUU      F           78285   12.77     0.13
## UUC      F           145663  23.76     0.26
## ...      ...         ...     ...       ...


class FreqTable(object):
	"""
	Class for assembling a codon frequency table from fasta file.
	"""

	def __init__(self):
		self.translation_table = None
		self.table = None
		self.organism = None
		self.decimals = 2


	def calculate_ratio(self, codons):
		'''
		Calculate the ratio with wich each codon is used for each amino acid.
		'''
		codon_info = self.codon_table.get_codons()

		data = {}
		for AA in codon_info.keys():
			total = sum([codons[codon] for codon in codon_info[AA]])
			for codon in codon_info[AA]:
				data[codon] = codons[codon] / total
		return data


	def make_codon_freq_table(self, organism, translation_table, filepath, silent=False):
		'''
		Input is a file path.
		Counts the usage of each codon in a FASTA file of DNA sequences.
		Then converts that as codon usage per 1000 codons.
		Good for generating codon tables.
		Output is a dictionary of codon frequencies per 1000 codons and the total number in brackets.
		'''
		self.organism = organism.lower().replace(' ', '_')
		self.translation_table = int(translation_table)
		self.codon_table = bioseq.CodonTable(self.translation_table)
		print(organism)

		skipped_nondiv = 0
		skipped_noncomp = 0
		used = 0
		codons = {}

		records = fasta.parse_file(filepath, seqtype='dna')
		for record in records:
			header, sequence = record

			# Skip all that aren't divisible by 3
			if len(sequence) % 3 != 0:
				skipped_nondiv +=1
				continue

			# Count codons while skipping all that has other characters than ATCG
			try:
				cds = bioseq.DNA(sequence)
				gene_codons = cds.count_codons()

			except ValueError as e:
				skipped_noncomp +=1
				continue

			# Use these codons
			used += 1
			if codons == {}:
				codons = gene_codons
			else:
				for k in codons.keys():
					codons[k] += gene_codons[k]

		#sum codons
		codon_sum = sum(codons.values())

		#divide each by the sum and multiply by 1000
		ratio_data = self.calculate_ratio(codons)
		self.table = {}
		for triplet in list(codons.keys()):
			per_1000 = 1000*(codons[triplet]/codon_sum)
			number = codons[triplet]
			fraction = ratio_data[triplet]
			amino_acid = bioseq.RNA(triplet).translate(table=self.translation_table)

			self.table[triplet] = {'amino_acid':amino_acid,
					'number':int(number),
					'fraction':round(float(fraction), self.decimals),
					'per_1000':round(float(per_1000), self.decimals)}

		if not silent:
			print('Skipped %s sequences that were not divisible by 3' % str(skipped_nondiv))
			print('Skipped %s sequences that had non-compliant characters' % str(skipped_nondiv))
			print('Used %s sequences that were ok' % str(used))


	def read_table(self, filepath):
		'''
		Read a codon frequency table.
		'''
		assert exists(filepath)

		self.table = {}
		with open(filepath, 'r') as f:
			line = f.readline()
			assert line.startswith('#organism=')
			self.organism = int(line.strip().split('=')[1])

			line = f.readline()
			assert line.startswith('#translation_table=')
			self.translation_table = int(line.strip().split('=')[1])

			for line in f:
				line = line.strip()
				if line == '':
					continue
				triplet, amino_acid, number, fraction, per_1000 = line.split()
				self.table[triplet] = {'amino_acid':amino_acid,
									'number':int(number),
									'fraction':round(float(fraction), self.decimals),
									'per_1000':round(float(per_1000), self.decimals)}


	def write_table(self, filepath):
		'''
		Write codon table to file.
		'''
		with open(filepath, 'w') as f:
			f.write('#organism=%s\n' % str(self.organism))
			f.write('#translation_table=%s\n' % str(self.translation_table))
			f.write('triplet\tamino_acid\tnumber\tfraction\t/1000\n')
			for triplet in sorted(self.table.keys()):
				amino_acid = self.table[triplet]['amino_acid']
				number = str(self.table[triplet]['number'])
				fraction = '%.2f' % round(self.table[triplet]['fraction'], self.decimals)
				per_1000 = '%.2f' % round(self.table[triplet]['per_1000'], self.decimals)
				f.write('%s\t%s\t%s\t%s\t%s\n' % (triplet, amino_acid, number, fraction, per_1000))


	def frequency_table(self):
		'''
		Return the codon frequency table to user.
		'''
		return self.table


	def translation_table(self):
		'''
		Return the translation table to user.
		'''
		return self.translation_table


	def organism(self):
		'''
		Return the organism to user.
		'''
		return self.organism


def create(fasta_file, outfile, organism, translation_table):
	'''
	Create a codon frequency table.
	'''
	table_obj = FreqTable()
	table_obj.make_codon_freq_table(organism, translation_table, fasta_file)
	table_obj.write_table(outfile)


def load(filepath):
	'''
	Read codon table from file.
	'''
	table_obj = FreqTable()
	table_obj.read_table(filepath)
	return table_obj.frequency_table()

# How do I send on the organism and translation table data?
