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
# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/Utilities24/compute_AD42ScoreC.py,v 1.3.6.1 2015/08/26 22:44:52 sanner Exp $
#

import os 

from MolKit import Read
import MolKit.molecule
import MolKit.protein
#?@@?
from AutoDockFR import evaluateEnergy
from AutoDockFR.ADCscorer import AD42ScoreC
#rigidRecAtoms, ligAtoms, ligandTorTree, TORSDOF, flexRecAtoms=None,
#		cutoff=1.0, RR_L=True, FR_L=True, L_L=True, RR_RR=True, RR_FR=True, FR_FR=True,
#		RR_L_Fitness=True, FR_L_Fitness=True, L_L_Fitness=True, RR_RR_Fitness=True, 
#                RR_FR_Fitness=True, FR_FR_Fitness=True
#6 possible scorers:
#case1: rec-ligand
#CREATE SPECIFIED molecular systems:
#L_L (1)    T   vs dlgs w/flexres    w/out flexRes
#RR_RR (2)  T   @@ ==>False         False
#FR_FR (3)  T         True      ==> False
#FR_L (4)   T         True      ==> False
#RR_L (5)   T         True
#RR_FR (6)  T         True      ==> False
# RR_L, FR_L, RR_FR <=Intermolecular energies
#@@ *_Fitness internal energy 
#INCLUDE score in GA Fitness
#RR_L_Fitness  T 
#FR_L_Fitness  T
#L_L_Fitness   T @@
#RR_RR_Fitness T     False         False
#RR_FR_Fitness T     True          False
#FR_FR_Fitness T @@  True          True



#check these autodock comparisons:
#case1: rec-ligand
#case2: rec-ligand-flexres

if __name__ == '__main__':
    import sys
    import getopt


# compute_AD42ScoreC.py -r receptorfilename -l ligandfilename [-o outputfilename -v
# for comparison with docking output: these are set to False: RR_RR, 
    def usage():
        "Print helpful, accurate usage statement to stdout."
        print("Usage: compute_AD42ScoreC.py -r receptorfilename -l ligandfilename")
        print()
        print("    Description of command...")
        print("      Calculates the AD42ScoreC energy between two molecules and stores it in 'AD42ScoreC_scores.txt'")
        print("      By default....")
        print("         -r     receptorfilename")
        print("         -l     ligandfilename")
        print("    Optional parameters:")
        print("        [-o]    outputfilename (default is 'AD42ScoreC_scores.txt')")
        print("        [-x]    flexresfilename (default is none)")
        print("        [-w]    write_file_mode (default mode is 'a')")
        print("        [-p]    printAllScoreTerms (default is not to print them all)")
#defaults based on assumption using simplest case: receptor-ligand
#        print "        [-R]    include Rigid-Rec - Rigid-Rec (default is False)"
#        print "        [-I]    include Rigid-Rec - Flex-Rec (default is False)"
#        print "        [-f]    include FR_L_Fitness (default is False)"
#        print "        [-g]    include RR_RR_Fitness (default is True)"
#        print "        [-j]    include RR_FR_Fitness (default is False)"
#        print "        [-k]    include FR_FR_Fitness (default is False)"
#        print "        [-m]    flexRec Atoms  (default is None)"
#        print "        [-p]    printAllScoreTerms (default is False)"
        print("        [-v]    verbose output")

    FR_L = False  #no Flex-Rec - Ligand @@ True by default in AD42ScoreC
    L_L = True    # internal energy ligand
    RR_RR = False  # internal energy Rigid-Rec
    RR_L = True   # rigid-receptor ligand
    # if flexresfilename, turn these on: RR_FR, FR_FR, FR_L_Fitness, 
    RR_FR = False #no Rigid-Rec - Flex-Rec @@ True by default
    FR_FR = False  #no Flex-Rec - Flex-Rec @@ True by default
    RR_L_Fitness = True #Rigid-Rec - Ligand True by default
    L_L_Fitness = True  # ?
    FR_L_Fitness = False #no Flex-Rec - Ligand score @@ True by default
    RR_RR_Fitness = False #Rigid-Rec - Rigid-Rec Fitness
    RR_FR_Fitness = False #no Rigid-Rec - Flex-Rec Fitness @@ True by default
    FR_FR_Fitness = False #no Flex-Rec - Flex-Rec score @@ True by default
    cutoff = 1.00

    # process command arguments
    try:
        opt_list, args = getopt.getopt(sys.argv[1:], 'r:l:o:x:w:pvh')
        #opt_list, args = getopt.getopt(sys.argv[1:], 'r:l:o:x:wzRIfgjkmvh')
    except getopt.GetoptError as msg:
        print('compute_AD42ScoreC.py: %s' %msg)
        usage()
        sys.exit(2)

    # initialize required parameters
    #-r: receptorfilename
    receptorfilename =  None
    #-l: ligandfilename
    ligandfilename =  None
    # optional parameters
    #-o outputfilename
    outputfilename = 'AD42ScoreC_scores.txt'
    #-x flexresfilename 
    flexresfilename = None
    flexRecAtoms = None
    FR_L_Fitness = False
    #-w write new outputfile
    write_file_mode = 'a'
    #FR_L = False
    #RR_FR = False
    #-R include_RigidRec_RigidRec #@@ july3 sanity
    #RR_RR = False #@@True by default
    #-I include_RigidRec_FlexRec 
    #RR_FR = False
    #-f include_FR_L_Fitness
    #-g include_RR_RR_Fitness 
    #RR_RR_Fitness = False
    #-j include_RR_FR_Fitness 
    #RR_FR_Fitness = False
    #-k include_FR_FR_Fitness 
    #FR_FR_Fitness = False
    #-p printAllScoreTerms 
    printAllScoreTerms = False
    #-v verbose
    verbose = None
    # current args:
    # r:l:o:x:w:pvh
    # previously:
    #'r:l:o:x:w:RIfgjkmpvh'

    for o, a in opt_list: #r:l:o:x:w:pvh
        if o in ('-r', '--r'):
            receptorfilename = a
            if verbose: print('set receptorfilename to ', a)
        if o in ('-l', '--l'):
            ligandfilename = a
            if verbose: print('set ligandfilename to ', a)
        if o in ('-o', '--o'):
            outputfilename = a
            if verbose: print('set outputfilename to ', a)
        if o in ('-x', '--x'):
            flexresfilename = str(a)
            print("line 148: a=", a)
            if verbose: print('set flexresfilename to ', flexresfilename)
            flexres = Read(flexresfilename)[0]
            flexres.buildBondsByDistance()
            flexRecAtoms = flexres.allAtoms
            #print "now flexRecAtoms=", flexRecAtoms
            FR_L = True
            FR_L_Fitness= True
            #RR_FR = True awaiting fix of  BUG in AutoDockFR/ScoringFunction.py => recIE never initialized
            #FR_FR = True #@@ internal energy of flexres
            ##RR_FR_Fitness = True #1
            ##FR_FR_Fitness = True #2
        if o in ('-w', '--w'):
            write_file_mode = str(a)
            if verbose: print('set write_file_mode to ', a)
        if o in ('-p', '--p'):
            printAllScoreTerms = True
            if verbose: print("set printAllScoreTerms to ", printAllScoreTerms)
        if o in ('-v', '--v'):
            verbose = True
            if verbose: print('set verbose to ', True)
        if o in ('-h', '--'):
            usage()
            sys.exit()

    if not receptorfilename:
        print('compute_AD42ScoreC: receptorfilename must be specified.')
        usage()
        sys.exit()

    if not ligandfilename:
        print('compute_AD42ScoreC: ligandfilename must be specified.')
        usage()
        sys.exit()

    ##@@ L_L = True  #@@ ==  AD42ScoreC default: Ligand - Ligand
#    RR_L = True  #@@ == AD42ScoreC default: Rigid-Rec - Ligand
#    FR_L = False  #no Flex-Rec - Ligand @first pass: assume no flexres
#    RR_RR = False  #@@    Rigid-Rec - Rigid-Rec  internal energy of receptor :o 
#    RR_FR = False #no Rigid-Rec - Flex-Rec @first pass
#    FR_FR = False #no Flex-Rec - Flex-Rec  @first pass
#    FR_L_Fitness = False #no Flex-Rec - Ligand score
#    RR_L_Fitness = True #@@ 
#    L_L_Fitness = True #@@
#    RR_RR_Fitness = False #@@
#    RR_FR_Fitness = False #no Rigid-Rec - Flex-Rec score
#    FR_FR_Fitness = False #no Flex-Rec - Flex-Rec score
 

    receptor = Read(receptorfilename)[0]
    assert os.path.splitext(receptor.parser.filename)[-1]=='.pdbqt',"receptor file not in required '.pdbqt' format"
    if verbose: print('read ', receptorfilename)
    receptor.buildBondsByDistance()
    #check format of receptor ?@@?
    ligand = Read(ligandfilename)[0]
    assert os.path.splitext(ligand.parser.filename)[-1]=='.pdbqt',"ligand file not in required '.pdbqt' format"
    if verbose: print('read ', ligandfilename)
    ligand.buildBondsByDistance()    

    ad42sc = AD42ScoreC(receptor.allAtoms, ligand.allAtoms, ligand.torTree, ligand.TORSDOF,
                        flexRecAtoms=flexRecAtoms, RR_L=RR_L, FR_L=FR_L, RR_FR=RR_FR, RR_RR=RR_RR, FR_FR=FR_FR,
                        FR_L_Fitness=FR_L_Fitness, FR_FR_Fitness=FR_FR_Fitness, RR_RR_Fitness=RR_RR_Fitness,
                        RR_FR_Fitness=RR_FR_Fitness)

    # compute the score...
    score = ad42sc.score()
    print("score=", score)

    if printAllScoreTerms:
        #print "calling ad42sc.printAllScoreTerms"
        #print dir(ad42sc)
        ad42sc.printAllScoreTerms()
    mode = 'a'
    first = not os.path.exists(outputfilename)
    if write_file_mode:
        mode = 'w'
        first = True
    if verbose: print('first is ', first)
    optr = open(outputfilename, mode)
    if first:
        tstr = "    Receptor      Ligand   AD42Score\n" 
        optr.write(tstr)
    #setup the scorer:
    ostr = "%12s-%12s      %6.3f\n" %(receptor.name,ligand.name, score)
    optr.write(ostr)
    if verbose:
        print('wrote %s' %outputfilename)
    optr.close()

# To execute this command type:
# compute_AD42ScoreC.py -r receptorfilename -l ligandfilename [-o outputfilename  -w write_file_mode] -v




