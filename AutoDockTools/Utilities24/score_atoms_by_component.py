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

# SETUP input files:
#1. In your local MGLToolsPckgs directory:
#-create a pdbqt file containing the docked ligand conformation for the energy-breakdown.
#( I saved mine as "ind4.1.1best.pdbqt")
#-copy the receptor file here, too ('hsg1.pdbqt')
#3. In your local MGLToolsPckgs directory, start this script with mgltoolspckgs pythonsh
# to output to terminal:
#../bin/pythonsh ./script_score_atoms_by_component.py -l ind4.1.1best.pdbqt -r hsg1.pdbqt 
# to save output to a file:
#../bin/pythonsh ./script_score_atoms_by_component.py -l ind4.1.1best.pdbqt -r hsg1.pdbqt -o output.txt
# add '-d' to include an energy breakdown grouped as in a dlg:
#../bin/pythonsh ./script_score_atoms_by_component.py -l ind4.1.1best.pdbqt -r hsg1.pdbqt -o output.txt -d


#imports
import os
from MolKit import Read
from PyAutoDock.AutoDockScorer import AutoDock41Scorer, AutoDockTermWeights41
from PyAutoDock.MolecularSystem import MolecularSystem
import numpy


if __name__ == '__main__':
    import sys
    import getopt


    def usage():
        "Print helpful, accurate usage statement to stdout."
        print("Usage: score_atoms_by_component.py -r receptor_file -l ligand_file")
        print()
        print("    Description of command...")
        print("         -r     receptor_file")
        print("         -l     ligand_file")
        print("    Optional parameters:")
        print("        [-o]    output_file")
        print("                      (default is scores.txt)")
        print("        [-d]    also include grouping energies in dlg-format")
        print("        [-v]    verbose output")


    # process command arguments
    try:
        opt_list, args = getopt.getopt(sys.argv[1:], 'r:l:o:dvh')
    except getopt.GetoptError as msg:
        print('score_atoms_by_component.py: %s' %msg)
        usage()
        sys.exit(2)

    # initialize required parameters
    verbose = None
    #-r receptor_file
    receptor_file = None
    #-l ligand_file
    ligand_file = None
    #[-o output_file]
    output_file = None
    optr = None   
    #[-d group energies as in dlg 'Run']
    include_dlg_grouping = False
                                                                                                                                 
    #'r:l:o:dvh'
    for o, a in opt_list:
        #print "o=", o, " a=", a
        if o in ('-o', '--o'):
            output_file = a
            try:
                optr = open(output_file, 'w')
            except:
                 print("problem opening file to contain per-atom energy breakdown:", output_file)
                 sys.exit()
            if verbose: print('set output_file to ', a)
            optr = open(output_file, 'w')
        if o in ('-r', '--r'):
            receptor_file = a
            if not os.path.isfile(receptor_file):
                 print("receptor_file:", receptor_file, " does not exist!")
                 sys.exit()
            if verbose: print('set receptor_file to ', a)
        if o in ('-l', '--l'):
            ligand_file = a
            if not os.path.isfile(ligand_file):
                 print("ligand_file ", ligand_file, " does not exist!")
                 sys.exit()
            if verbose: print('set ligand_file for energy component output to ', a)
        if o in ('-d', '--d'):
            include_dlg_grouping = True
            if verbose: print('set include_dlg_grouping to ', True)
        if o in ('-v', '--v'):
            verbose = True
            if verbose: print('set verbose to ', True)
        if o in ('-h', '--'):
            usage()
            sys.exit()


    if not receptor_file:
        print('score_docked_atoms_by_component: receptor filename must be specified.')
        usage()
        sys.exit()

    if not ligand_file:
        print('score_docked_atoms_by_component: ligand file must be specified.')
        usage()
        sys.exit()


    #1: create the python molecules
    rec = Read(receptor_file)[0]
    recbnds = rec.buildBondsByDistance()
    lig = Read(ligand_file)[0]
    ligbnds = lig.buildBondsByDistance()

    #2: setup python scorers for your molecules:
    ms = MolecularSystem()
    ms.add_entities(rec.allAtoms)
    ms.add_entities(lig.allAtoms)

    #3: setup python scorers
    scorer = AutoDock41Scorer()
    scorer.set_molecular_system(ms)
    estat_scores = numpy.add.reduce(scorer.terms[0][0].get_score_array())
    dsolv_scores = numpy.add.reduce(scorer.terms[1][0].get_score_array())
    vdw_scores = numpy.add.reduce(scorer.terms[2][0].get_score_array())
    hbond_scores = numpy.add.reduce(scorer.terms[3][0].get_score_array())
    wts = AutoDockTermWeights41()

    #4: open outputfile for writing scores
    if output_file is not None:
        fptr = open(output_file, 'w')

    #5. score the docked pose outputting individual terms for each atom:
    all_estat = []
    all_VHD = []
    for i in range(len(lig.allAtoms)):
        estat = round(estat_scores[i]*wts.estat_weight,3)
        vdw = round(vdw_scores[i]*wts.vdw_weight,3)
        hbond = round(hbond_scores[i]*wts.hbond_weight,3)
        dsolv = round(dsolv_scores[i]*wts.dsolv_weight,3)
        all_estat.append(estat)
        all_VHD.append(vdw+hbond+dsolv)
        ostr= "%4d:%6.3f %6.3f %6.3f %6.3f = %6.3f\n"%(i+1,estat,vdw,hbond,dsolv,estat+vdw+hbond+dsolv)
        if optr is None:
            if i==0: print("\n PyAutoDock AD4 energies:\n at#  estat    vdw  hbond  dsolv    total")
            print(ostr[:-1])
        else:
            if i==0: optr.write( "\n PyAutoDock AD4 energies:\n at#  estat    vdw  hbond  dsolv    total\n")
            optr.write(ostr)

    #6. [optional: alternative grouping of output for comparison with numbers in dlg]
    #              output vdw+hb+ds ,estat for each atom:
    if include_dlg_grouping:
        for i in range(len(lig.allAtoms)):
            ostr ="%4d:%6.3f  %6.3f = %6.3f\n"%(i+1,all_estat[i], all_VHD[i], all_estat[i]+all_VHD[i])
            if optr is None:
                if i==0: print("\n energies grouped as in dlg:\n at#    vdW    Elec    total")
                print(ostr[:-1])
            else:
                if i==0: optr.write( "\n energies grouped as in dlg:\n at#  ~vdW~    Elec    total\n")
                optr.write(ostr)
    # !Done!
    if optr is not None:
        optr.close()


