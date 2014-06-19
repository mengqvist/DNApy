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


#script from above source modified by Martin Engqvist 


import re
import string

def re_const(oligo):
    '''This is the actual search function'''
    oligo=oligo.upper()
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
		L_out.append([match.start()+1, match.end()+len(match.group(1))])
	return L_out

