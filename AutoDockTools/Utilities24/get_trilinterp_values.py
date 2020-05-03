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
# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/Utilities24/get_trilinterp_values.py,v 1.2 2011/12/21 21:37:42 rhuey Exp $
#
import os, glob
from MolKit import Read
from Volume.IO.AutoGridReader import ReadAutoGrid
from Volume.Grid3D import Grid3D
from Volume.Operators.trilinterp import trilinterp



if __name__ == '__main__':
    import sys
    import getopt

    def usage():
        print("Usage: get_trilinterp_energies.py -r receptor_filename -l ligand_filename -n atom_indicies")
        print("    Given filenames for receptor and ligand and a list of atom indicies")
        print("      this script: ")
        print("         (1) determines atom types + coordinates to interpolate")
        print("         (2) reads corresponding autogrid maps ")
        print("         (3) calculates interpolated energy from appropriate autogrid map ")
        print("         OPTIONS: ")
        print("         (4) write output to a file")
        print("         (5) output indices and energies of 8 gridpts used by trilinear interpolation for each atom  ")
        print("    ")
        print("    -r receptor_filename")
        print("    -l ligand_filename")
        print("    -i comma-separated list of indicies of atoms in ligand")
        print("    -s autogrid map stem")
        print("    ")
        print("    Optional parameters:")
        print("        [-o output filename] (default output to screen)")
        print("        [-t output total] (default output per atom index )")
        print("        [-c output 8pts] (default not to do this)")
        print("        [-v]    verbose output")

    # process command arguments
    try:
        opt_list, args = getopt.getopt(sys.argv[1:], 'r:l:i:s:o:tcvh')
    except getopt.GetoptError as msg:
        print('get_trilinterp_energies.py: %s' %msg)
        usage()
        sys.exit(2)

    # initialize required parameters
    #-r receptorfilename
    receptorfilename = None
    #-l ligandfilename
    ligandfilename = None
    #-i indicies
    indicies = ""
    #-s mapstem
    mapstem = None

    # optional parameters
    #-o outputfilename
    outputfilename = None
    #-t output_total
    output_total = 0
    #-c output_8pts
    output_8pts = 0
    #-v verbose
    verbose = False

    #'r:l:i:s:o:tvh'
    if verbose: print("opt_list=", opt_list)
    for o, a in opt_list:
        if verbose: print("o=", o, " a=", a)
        if o in ('-r', '--r'):
            receptorfilename = a
            if verbose: print('set receptorfilename to ', a)
        if o in ('-l', '--l'):
            ligandfilename = a
            if verbose: print('set ligandfilename to ', a)
        if o in ('-i', '--i'):
            indicies = a
            if verbose: print('set indicies to ', a)
        if o in ('-s', '--s'):
            mapstem = a
            if verbose: print('set mapstem to ', a)
        if o in ('-o', '--o'):
            outputfilename = a
            if verbose: print('write output to file:', a)
        if o in ('-t', '--t'):
            output_total = True
            if verbose: print('set output_total to ', True)
        if o in ('-c', '--c'):
            output_8pts = True
            if verbose: print('set output_8pts to ', True)
        if o in ('-v', '--v'):
            verbose = True
            if verbose: print('set verbose to ', True)
        if o in ('-h', '--'):
            usage()
            sys.exit()

    if verbose: print("receptorfilename=", receptorfilename)
    if receptorfilename is None:
        print("trilinterp_atoms.py: receptor filename must be specified.")
        usage()
        sys.exit()
    if ligandfilename is None:
        print("trilinterp_atoms.py: ligand filename must be specified.")
        usage()
        sys.exit()
    #OPTIONAL:
    if mapstem is None and verbose:
        print("trilinterp_atoms.py: map stem set from receptorfilename")
    if outputfilename is None: 
        if verbose: print("trilinterp_atoms.py: output to screen")
    else:
        print('writing output to %s' %outputfilename)

    #Start to work:
    rec = Read(receptorfilename)[0]
    rec.buildBondsByDistance() #??@@??
    #setup mapstem
    if mapstem is None: 
        mapstem = rec.name

    lig = Read(ligandfilename)[0]
    lig.buildBondsByDistance()
    #get atoms for trilinterp
    if indicies is None or indicies=='':
        indicies = list(range(len(lig.allAtoms)))
    else:
        indicies = list(map(int, indicies.split(',')))
    if verbose: print("now indicies=", indicies)

    total = 0
    pts = []
    maps = {}
    inv_spacing = -1
    if outputfilename: fptr = open(outputfilename, 'w')
    for IND in indicies:
        at = lig.allAtoms[IND]
        pts= [(at.coords)]
        mapname = "%s.%s.map"%(mapstem,at.autodock_element)
        if at.autodock_element not in list(maps.keys()):
            reader = ReadAutoGrid()
            maps[at.autodock_element] = reader.read(mapname,normalize=True) #@@ CHECK THIS!!
        grid3D = maps[at.autodock_element]
        if inv_spacing <0:
            inv_spacing = (1./grid3D.stepSize[0], 1./grid3D.stepSize[1], 1./grid3D.stepSize[2])
        if output_8pts:
            try:
                vals = trilinterp(pts, grid3D.data, inv_spacing, grid3D.getOriginReal(), output_8pts=1)
            except:
                print("IGNORING output_8pts: NOT supported in current Volume/Operators/trilinterp.py")
                vals = trilinterp(pts, grid3D.data, inv_spacing, grid3D.getOriginReal())
        else:
            vals = trilinterp(pts, grid3D.data, inv_spacing, grid3D.getOriginReal())
        for i in range(len(vals)):
            print(at.full_name(), ':', pts[i], '--', vals[i])
            if outputfilename:
                fptr.write("%s : [%6.3f, %6.3f, %6.3f] -- %14.11f\n" %(at.full_name(), pts[i][0],pts[i][1],pts[i][2], vals[i]))
            total += float(vals[i])
    if output_total or verbose:
        print("Total energy =", total)

# To execute this command type:
# get_trilinterp_energies.py  [-r receptorfilename -l ligandfilename -n indicies -o outputfilename] -v




