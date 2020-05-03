#!/usr/bin/env python

# ##################################################################################################
#  Disclaimer                                                                                      #
#  This file is a python3 translation of AutoDockTools (v.1.5.7)                                   #
#  Modifications made by Valdes-Tresanco MS (https://github.com/Valdes-Tresanco-MS)                #
#  Tested by Valdes-Tresanco-MS and Valdes-Tresanco ME                                             #
#  There is no guarantee that it works like the original distribution,                             #
#  but feel free to tell us if you get any difference to correct the code.                         #
#                                                                                                  #
#  Please use this cite the original reference.                                                    #
#  If you think my work helps you, just keep this note intact on your program.                     #
#                                                                                                  #
#  Modification date: 2/5/20 19:51                                                                 #
#                                                                                                  #
# ##################################################################################################

#
# 
#
# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/Utilities24/prepare_ligand_vif.py,v 1.2 2012/01/31 17:57:37 rhuey Exp $ 
#
import os 

from MolKit import Read
from string import split, strip




if __name__ == '__main__':
    import sys
    import getopt


    def usage():
        "Print helpful, accurate usage statement to stdout."
        print("Usage: prepare_ligand_vif.py -l filename")
        print()
        print("    Description of command...")
        print("         -l     ligand_filename (.pdbqt format)")
        print("    Optional parameters:")
        print("        [-v]    verbose output")
        print("        [-o pdbqt_filename] (default output filename is ligand_filename_stem + '_L.pdbqt')")
        print("        [-P]    list of indicies of residues to write with LP dummy atoms '22,26,30,40,41,42' ")
        print("        [-S]    list of indicies of residues to write with LS dummy atoms '15,17' ")


    # process command arguments
    try:
        opt_list, args = getopt.getopt(sys.argv[1:], 'l:vo:P:S:h')
    except getopt.GetoptError as msg:
        print('prepare_ligand_vif.py: %s' %msg)
        usage()
        sys.exit(2)

    # initialize required parameters
    #-l: ligand
    ligand_filename =  None
    # optional parameters
    verbose = None

    #-o outputfilename
    outputfilename = None
    #-P LP_atom_residues
    LP_atom_residues = None
    #-S LS_atom_residues
    LS_atom_residues = None

    #'l:vo:P:S:'
    for o, a in opt_list:
        #print "o=", o, " a=", a
        if o in ('-l', '--l'):
            ligand_filename = a
            if verbose: print('set ligand_filename to ', a)
            outputfilename = ligand_filename.split('.')[0] + '_LLL.pdbqt'
            if verbose: print('set outputfilename from ligand_filename to ', outputfilename)
        if o in ('-v', '--v'):
            verbose = True
            if verbose: print('set verbose to ', True)
        if o in ('-o', '--o'):
            outputfilename = a
            if verbose: print('set outputfilename to ', outputfilename)
        if o in ('-P', '--P'):
            LP_atom_residues = a
            if LP_atom_residues.find(',')>-1:
                LP_atom_residues = list(map(int, split(a, ',')))
            if verbose: print('set LP_atom_residues to ', LP_atom_residues)
        if o in ('-S', '--S'):
            LS_atom_residues = a
            if LS_atom_residues.find(',')>-1:
                LS_atom_residues = list(map(int, split(a, ',')))
            if verbose: print('set LS_atom_residues to ', LS_atom_residues)
        if o in ('-h', '--'):
            usage()
            sys.exit()


    # check input 
    if not ligand_filename:
        print('prepare_ligand_vif: ligand filename must be specified.')
        usage()
        sys.exit()
    
    if not LP_atom_residues:
        print('prepare_ligand_vif: LP_atom_residues must be specified.')
        usage()
        sys.exit()

    if not LS_atom_residues:
        print('prepare_ligand_vif: LS_atom_residues must be specified.')
        usage()
        sys.exit()
    
    if verbose: 
        print('LP_atom_residues=', LP_atom_residues)
        print('LS_atom_residues=', LS_atom_residues)
        print("reading ", ligand_filename)

    # process ligand
    #??check that ligand in pdbqt format??
    ext = os.path.splitext(ligand_filename)[1]
    assert ext=='.pdbqt'
    fptr = open(ligand_filename)
    liglines = fptr.readlines()
    fptr.close()
    if verbose: print('read ', len(liglines), ' lines from ', ligand_filename) 
    
    optr = open(outputfilename, 'w')
    if verbose: print("writing ", outputfilename)

    # check whether already has ROOT/ENDROOT/TORSDOF
    i = 0
    has_root = 0
    if liglines[i]=='ROOT\n':
        optr.write('ROOT\n')
        has_root = 1
        i+=1
    else:
        optr.write('ROOT\n')

    for j in liglines[i:]:
        ll = split(j)
        if j[0:4]in ['ATOM', 'HETA']:
            atname = strip(j[12:16])
            resnum = strip(j[22:26])
            if strip(atname)=='CA':
                optr.write(j)
                #@@ guard against duplicates?
                #if j[-4:].find("RP")<0:
                optr.write(j[:-4]+' RP\n') 
                if verbose: print("wrote RP")
                resnum = strip(resnum)
                if int(resnum) in LP_atom_residues:
                    optr.write(j[:-4]+' LP\n')
                    if verbose: print(" LP ", resnum)
                if int(resnum) in LS_atom_residues:
                    optr.write(j[:-4]+' LS\n')
                    if verbose: print(" LS ", resnum)
            else:
                optr.write(j)
        else:
            optr.write(j)
    if not has_root:
        optr.write("ENDROOT\n")
        optr.write("TORSDOF 0 \n")

    optr.close()

# To execute this command type:
#prepare_ligand_vif -l 1vzf.pdbqt -P 22,26,30,40,41,42 -S 15,17 -o 1vzf_L.pdbqt
