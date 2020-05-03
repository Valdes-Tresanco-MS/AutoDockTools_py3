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
# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/Utilities24/write_component_energies.py,v 1.4 2011/01/10 23:18:41 rhuey Exp $
#
import os 

from MolKit import Read
from PyAutoDock.InternalEnergy import InternalEnergy
from PyAutoDock.MolecularSystem import MolecularSystem
from PyAutoDock.AutoDockScorer import AutoDock4Scorer, AutoDock41Scorer




if __name__ == '__main__':
    import sys
    import getopt


# write_component_energies.py -r receptorfilename -l ligandfilename [-o outputfilename -a append_to_file] -v
    def usage():
        "Print helpful, accurate usage statement to stdout."
        print("Usage: write_component_energies.py -r receptorfilename -l ligandfilename")
        print()
        print("    Description of command...")
        print("         -r     receptorfilename")
        print("         -l     ligandfilename")
        print("    Optional parameters:")
        print("        [-o]    outputfilename (default is 'component_energies.txt')")
        print("        [-w]    overwrite_file (default mode is to append)")
        print("        [-n]    use AutoDock4 weights (default mode uses AutoDock41 weights)")
        print("        [-v]    verbose output")

    # initialize required parameters
    #-r: receptorfilename
    receptorfilename =  None
    #-l: ligandfilename
    ligandfilename =  None
    #-a append to outputfile
    append_to_outputfile = True
    # optional parameters
    verbose = None
    #-o outputfilename
    outputfilename = 'component_energies.txt'
    #-n use_autodock40_weights 
    use_autodock40_weights = False

    # process command arguments
    try:
        opt_list, args = getopt.getopt(sys.argv[1:], 'r:l:o:anvh')
    except getopt.GetoptError as msg:
        print('write_component_energies.py: %s' %msg)
        usage()
        sys.exit(2)

    #'r:l:o:wnvh'
    for o, a in opt_list:
        #print "o=", o, " a=", a
        if o in ('-r', '--r'):
            receptorfilename = a
            if verbose: print('set receptorfilename to ', a)
        if o in ('-l', '--l'):
            ligandfilename = a
            if verbose: print('set ligandfilename to ', a)
        if o in ('-o', '--o'):
            outputfilename = a
            if verbose: print('set outputfilename to ', a)
        if o in ('-w', '--w'):
            append_to_outputfile = False
            if verbose: print('set append_to_outputfile to ', append_to_outputfile)
        if o in ('-n', '--n'):
            use_autodock40_weights = True
            if verbose: print('set use_autodock41_weights to ', use_autodock41_weights)
        if o in ('-v', '--v'):
            verbose = True
            if verbose: print('set verbose to ', True)
        if o in ('-h', '--'):
            usage()
            sys.exit()


    if not receptorfilename:
        print('write_component_energies: receptorfilename must be specified.')
        usage()
        sys.exit()

    if not ligandfilename:
        print('write_component_energies: ligandfilename must be specified.')
        usage()
        sys.exit()

    receptor = Read(receptorfilename)[0]
    if verbose: print('read ', receptorfilename)
    receptor.buildBondsByDistance()

    ligand = Read(ligandfilename)[0]
    if verbose: print('read ', ligandfilename)
    ligand.buildBondsByDistance()
    mode = 'w'
    first = True
    if append_to_outputfile:
        mode = 'a'
        first = not os.path.exists(outputfilename)
    if verbose: print('first is ', first)
    optr = open(outputfilename, mode)
    if first:
        tstr = "RECEPTOR_LIGAND         #ESTAT       #HB      #VDW    #DSOLV  #TORSDOF #ATS #TORS\n"
        optr.write(tstr)
    #get the scores
    ms = MolecularSystem()
    ms.add_entities(receptor.allAtoms)
    ms.add_entities(ligand.allAtoms)
    if use_autodock40_weights:
        ad_scorer = AutoDock4Scorer()
    else:
        ad_scorer = AutoDock41Scorer()
    ad_scorer.set_molecular_system(ms)
    score_list = ad_scorer.get_score_per_term()
    ostr = "%8s_%6s," %(receptor.name,ligand.name)
    for score in score_list:
        ostr = ostr + " % 8.4f," %score
    ostr = ostr + " % 8d," %ligand.TORSDOF
    ostr = ostr + " % 8d," %len(ligand.allAtoms)
    ostr = ostr + " % 8d," %len(ligand.torTree.torsionMap)
    ostr = ostr + "\n"   
    optr.write(ostr)
    optr.close()
    if verbose:
        print('wrote %s' %outputfilename)


# To execute this command type:
# write_component_energies.py -r receptorfilename -l ligandfilename [-o outputfilename -a append_to_file -n use_autodock40_weights] -v




