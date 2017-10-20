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

## Sequence manipulations
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

# Transcribe the reverse-complement to RNA
x.reverse_transcribe()

# Translate to protein
x.translate()

# Translate reverse-complement to protein
x.reverse_translate()

# Randomize sequence
x.randomize()

# Get molecular weight of sequence
x.mass()

# Count nucleotides in sequence
x.count_bases()

# Count codons in sequence
x.count_codons()
```


```
# Create RNA sequence object
x = bioseq.RNA('AUGGGAUGGUAA')
```


```
# Create Protein sequence object
x = bioseq.Protein('MGW*')
```

A description regarding what else one can do with these is forthcoming.
