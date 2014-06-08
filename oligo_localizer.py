#!/usr/bin/python2

##################################################################################
#   An oligo or motif finder written by Andrea Cabibbo           #      ######   #
#   Feel free to use/modify/redistribute the code               # #     #        #
#   If you use significant parts of this code please preserve  #   #    #        #
#   this header.                                              #######   #        # 
#   If you find bugs or have suggestions, please contact     #       #  #        #
#   the author at andrea.cabibbo@uniroma2.it                #         # ######   #
#   The tool is available online at                                              #
#   http://www.cellbiol.com/scripts/oligo/oligo_motif_sequence_finder.php        #
##################################################################################

#import cgitb; cgitb.enable()
import sys
import os
sys.path.insert(0, os.getcwd())
#import cgi
import re
import string
from copy import deepcopy

def cleaner(seq):
    '''Function for cleaning DNA of non-DNA characters'''
    i=0
    so=[] #seq_out
    for char in seq:
        if char in string.digits or char==('\\' or '/'):
            pass
            #print char,'character number %s is digit!! deleted'%i
        #elif char==' ' or char=='  ' or (char in string.whitespace):
            #print char,'is whitespace!! deleted'
        if char not in 'GATCRYWSMKHBVDN':   #IUPAC ambiguous_dna_values
            pass
            #so.append(dna_seq[i])
            # print 'WARNING: INVALID DNA CHARACTER <B>%s</B> DETECTED IN POSITION <B>%s</B>. Character NOT deleted'%(char,str(i+1))
        else:
            so.append(seq[i])
            #print char, 'normal character, not a number, not deleted'  
        i+=1
    so=string.join(so,'')
    return so



def re_const(oligo):
    '''This is the actual search function'''
    oligo=oligo.upper()
    wo=list(oligo)
    re_out=[]
    for char in oligo:
        if char=='G':
            re_out.append('[GRSKVDBXN]')
        elif char=='A':
            re_out.append('[AMRWVHDXN]')
        elif char=='T':
            re_out.append('[TWYKHDBXN]')
        elif char=='C':
            re_out.append('[CMSVHBXN]')
        elif char=='M':
            re_out.append('[MACVHXN]')
        elif char=='R':
            re_out.append('[RAGVDXN]')
        elif char=='W':
            re_out.append('[WATHDXN]')
        elif char=='S':
            re_out.append('[SCGVBXN]')
        elif char=='Y':
            re_out.append('[YCTHBXN]')
        elif char=='K':
            re_out.append('[KGTDBXN]')
        elif char=='V':
            re_out.append('[VACGXNMRS]')
        elif char=='H':
            re_out.append('[HACTXNMWY]')
        elif char=='D':
            re_out.append('[DAGTXNRWK]')
        elif char=='B':
            re_out.append('[BCGTXNSYK]')
        elif char=='X':
            re_out.append('[XNMRWSYKVHDBGATC]')
        elif char=='N':
            re_out.append('[NXMRWSYKVHDBGATC]')
    re_out=string.join(re_out,'')
    re_out="(?=(%s))" % re_out
    re_out_comp=re.compile(re_out,re.IGNORECASE)
    return re_out_comp


def match_oligo(seq,oligo,mismatches=0):
	'''Function for searching for a certain oligo'''
	re_oligo=re_const(oligo)
	L_out=[]
	for match in re_oligo.finditer(seq):
		L_out.append([match.end(),match.start(),match.group(1)])
	return L_out

     
if __name__=='__main__':
	dna_seq='cacccaaaaccaaaasdddc'
	dna_seq=dna_seq.upper()
	oligo='nnnnnnn'
	oligo=oligo.upper()
	dna_seq=list(dna_seq)
	myseq=cleaner(dna_seq)

	if oligo=='':
		print 'The oligo is missing, please <a href="oligo_finder.html">try again</A>'
	else:
		if match_oligo(myseq,oligo)!=[]:
			print('The following matches were found: ')
			for match in match_oligo(myseq,oligo):
				lm=len(match[2])
				print('from %s to %s %s' % (match[0]+1,match[1]+lm,match[2]))
		else:
			print('Sorry, no matches were found')



