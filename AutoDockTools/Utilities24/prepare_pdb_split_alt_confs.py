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
# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/Utilities24/prepare_pdb_split_alt_confs.py,v 1.4.2.1 2015/08/26 22:45:31 sanner Exp $
#

import os
from MolKit import Read
from MolKit.pdbParser import PdbParser

def hasAt(name):
    if '@' in name: 
        return 1
    else: 
        return 0

def getAT_types(m):
    AT_SET = m.allAtoms.get(lambda x: hasAt(x.name))
    AT_SET_SET = set(AT_SET.name)
    alt_items = {}
    for ent in AT_SET_SET:
        alt_items[ent.split("@")[1]] = 1
    print("@@returning ", list(alt_items.keys()), " @@")
    return list(alt_items.keys())
    

if __name__ == '__main__':
    import sys
    import getopt


    def usage():
        "Print helpful, accurate usage statement to stdout."
        print("Usage: prepare_pdb_split_alt_confs.py -r filename.pdb")
        print()
        print("    Description of command...")
        print("        -r   filename.pdb ")
        print("        create separate pdb files for file with alt_loc coords")
        print("    Optional parameters:")
        print("        [-n]  remove '   new' and '   flip' at end of lines (apparent issue from MolProbity)")
        print("        [-f]  filename to contain lines after new and flip removed (MolProbity)")
        print("        [-v]  verbose output (default is minimal output)")
        print("        [-o pdb_stem]  (default creates 'filename_A.pdb', 'filename_B.pdb' etc)")


    # process command arguments
    try:
        opt_list, args = getopt.getopt(sys.argv[1:], 'r:nf:vo:h')

    except getopt.GetoptError as msg:
        print('prepare_pdb_split_alt_confs.py: %s' %msg)
        usage()
        sys.exit(2)

    # initialize required parameters
    #-r: filename 
    filename =  None

    # optional parameters
    #-v verbose
    verbose = None
    #-o pdb_stem
    pdb_stem = None
    #-n remove_new
    remove_new = False
    #-f filename_nonews
    filename_nonews = None

    #'r:vo:nf:h'
    for o, a in opt_list:
        if o in ('-r', '--r'):
            filename = a
            if verbose: print('set filename to ', filename)
        if o in ('-n', '--n'):
            remove_new = True
            if verbose: print('set remove_new to ', remove_new)
        if o in ('-f', '--f'):
            filename_nonews = a
            if verbose: print('set filename for clean output ("new"+"flip" removed) to ', filename_nonews)
        if o in ('-v', '--v'):
            verbose = True
            if verbose: print('set verbose to ', verbose)
        if o in ('-o', '--o'):
            pdb_stem = a
            if verbose: print('set pdb_stem to ', pdb_stem)
        if o in ('-h', '--'):
            usage()
            sys.exit()


    if not filename:
        print('prepare_pdb_split_alt_confs: filename for cleanup must be specified.')
        usage()
        sys.exit()

    file_stem = filename.split('.')[0]
    if pdb_stem is not None:
        file_stem = pdb_stem

    if filename_nonews is None:
        filename_nonews = file_stem + "_nonews.pdb"
        if verbose: print("set filename_nonews to ", filename_nonews)

    if remove_new:
        optr = open(filename)
        lines = optr.readlines()
        if verbose: print("len(lines)==", len(lines))
        new_lines = []
        for ll in lines:
            if ll.find("   new")==78:
                ll = ll[:78] + "\n"
            if ll.find("   flip")==78:
                ll = ll[:78] + "\n"
            new_lines.append(ll)
        nptr = open(filename_nonews, 'w')
        if verbose: print("len(new_lines)==", len(new_lines))
        for ll in new_lines:
            nptr.write(ll)
        nptr.close()
        filename = filename_nonews
    mols = Read(filename)
    if not len(mols):
        msg = 'problem reading file: ' + filename + '! Unable to continue'
        print("ERROR reading " + filename)
        sys.exit()
    if verbose: print('read ', filename)
    mol = mols[0]
    if len(mols)>1:
        if verbose: print("more than one molecule in file using molecule with most atoms")
        #use the molecule with the most atoms
        ctr = 1
        for m in mols[1:]:
            ctr += 1
            if len(m.allAtoms)>len(mol.allAtoms):
                mol = m
                if verbose: print("mol set to ", ctr, "th molecule with", len(mol.allAtoms), "atoms")
    ats = mol.allAtoms.get("*@*")
    if not len(ats):
        print('Nothing to do:no alt loc atoms found in ', filename)
        sys.exit()
    list_to_write = getAT_types(mol)
    ATOMLINES = mol.parser.getAtomsLines(-2,0)
    
    for altT in list_to_write:
        fn = file_stem + '_' + altT + '.pdb'
        fptr = open(fn, 'w')
        ctr = 1
        for ll in ATOMLINES:
            if ll[16]==altT: #'B'
                newL = ll[:6] +"%5d" %(ctr) + ll[11:16]+" " + ll[17:] 
                ctr = ctr + 1
            elif ll[16]==' ':
                newL = ll[:6] +"%5d" %(ctr) + ll[11:16]+" " + ll[17:] 
                ctr = ctr + 1
            else:
                newL = ""
            if len(newL): fptr.write(newL)
        fptr.close()


# To execute this command type:
# prepare_pdb_split_alt_confs.py -r pdb_file -o pdb_stem 

