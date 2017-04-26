DNApy
=====


## This repo is currently broken. Efforts are underway to restore it and move to a development and deployment branch system.

A free and open source GUI toolkit for DNA editing - written in python

This project aims to provide a powerful codebase for viewing, editing and creating DNA in the GenBank format. The code is free to use, modify and re-distribute under a GPL license. Contributions in the form of improvements and new functions are welcome and encouraged!

The software is being developed on a Linux machine and works well in that environment. As the software matures testing will start on Windows and Mac to make sure it is cross-platform. The only external library that is needed is wxPython. Launch main.py to give the software a go (it is still under heavy development though).

## Testing
To start testing the software you have to install python, wxpython and pycairo:

```
sudo apt-get install python2.7 python-wxgtk2.8 python-cairo
```

Then you can download the software and run it:
```
cd ~
git clone https://github.com/mengqvist/DNApy.git
cd DNApy
python main.py
```

### Known problems
On Ubuntu 15.04 the following message might occure:
```
/usr/bin/xsel
Segmentation fault (core dumped)
```
This can be resolved by additionaly installing the package xclip:
```
sudo apt-get install xclip
```

![DNApy GUI](/Screenshot.png?raw=true "DNApy")

Implemented Software features
=====

* Visualization of DNA sequence with sequence features

* Plasmid view visualization

* DNA editing, copy, paste, reverse complement

* Unlimited undo/redo

* Easy search for nucleotide or amino acid positions or sequence

* Easy mutation by nucleotide or amino acid position (this one is pretty awesome!)

* Design of mixed base codons for libraries

* DNA translation to protein

* Addition/removal/modification of genbank features

* Addition/removal/modification of genbank qualifiers



Not Yet Implemented Software features (and priority list)
=====

* Analysis of sequence reads in the .ab1 format

* Restriction enzyme finder
 - [done] located restrictionsites in dna
 - [partly] display restrictionsites in plasmid (missing zoom)
 - [todo] improve dna editor to visualise cut location

* Addition/removal/modification of genbank header entries
 - [todo] improve genbank parser to allow parsing of corrupted genbank files from ApE, Serial Cloner, SnapGen Viewer

* DNA codon optimization

* Fetch genes/plasmids from NCBI

* Primer design

* Calculation of ribosome binding strength

* NCBI blast for homologous genes

* Simulate PCR

* (Multiple) Sequence alignment
