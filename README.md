DNApy
=====


## This repo is currently broken. Efforts are underway to restore it and move to a development and deployment branch system.
## Re-vamping the underlying scripts first, then the GUI... Please be patient....

A free and open source GUI toolkit for DNA editing - written in python

This project aims to provide a powerful codebase for viewing, editing and creating DNA in the GenBank format. The code is free to use, modify and re-distribute under a GPL license. Contributions in the form of improvements and new functions are welcome and encouraged!



## Installation
Download repository and unzip. cd to the project base folder and execute the command below:

```
pip3 install -e .
```

## Sequence manipulations (bioseq module)
For all commands below, first do this:

```
# Import the library
From dnapy.resources import bioseq
```


```
# Create DNA sequence object
x = bioseq.DNA('ATGGGATGGTAA')

# Reverse sequence
x.reverse()

# Get complement of sequence
x.complement()

# Get reverse-complement of sequence
x.reverse_complement()

# Transcribe to RNA
x.transcribe()

# Transcribe the reverse complement to RNA
x.transcribe_reverse_complement()

# Translate to protein
x.translate()

# Translate reverse-complement to protein
x.translate_reverse_complement()

# Randomize sequence
x.randomize()

# Get molecular weight of sequence
x.mass()

# Count number of each of the nucleotides in sequence
x.count_bases()

# Count codons in sequence
x.count_codons()
```


```
# Create RNA sequence object
x = bioseq.RNA('AUGGGAUGGUAA')

# Reverse sequence
x.reverse()

# Get complement of sequence
x.complement()

# Get reverse-complement of sequence
x.reverse_complement()

# Reverse-transcribe from RNA to DNA
x.reverse_transcribe()

# Reverse-transcribe reverse complement from RNA to DNA
x.reverse_transcribe_reverse_complement()

# Translate to protein
x.translate()

# Translate reverse-complement to protein
x.translate_reverse_complement()

# Randomize sequence
x.randomize()

# Get molecular weight of sequence
x.mass()

# Count number of each of the nucleotides in sequence
x.count_bases()

# Count codons in sequence
x.count_codons()
```


```
# Create Protein sequence object
x = bioseq.Protein('MGW*')

#
x.code_one()

#
x.code_three()

#
x.code_full()

# Count how many of each of the amino acids there are in the sequence
x.count_aa()

# Reverse translate to RNA using randomly chosen codon
x.reverse_translate()

#
x.reverse_translate_ambiguous()

# Get molecular weight of sidechains
x.sidechain_mass()

# Get molecular weight of entire protein
x.mass()

# Get hydropathy for entire protein
x.hydropathy()

# Get hydrophobicity for entire protein
x.hydrophobicity()

```

A description regarding what else one can do with these is forthcoming.


## Genbank file manipulations
# here


## Parse fasta files (fasta module)
# here
