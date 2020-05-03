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
# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/Utilities24/write_lowest_energy_ligand.py,v 1.3 2012/01/31 18:00:39 rhuey Exp $
#
import os, glob

from MolKit import Read
from AutoDockTools.Docking import Docking
from mglutil.math.rmsd import RMSDCalculator




if __name__ == '__main__':
    import sys
    import getopt


    def usage():
        "Print helpful, accurate usage statement to stdout."
        print("Usage: write_lowest_energy_ligand.py -f dlgfilename")
        print()
        print("    Description of command...")
        print("         -f     dlgfilename")
        print("    Optional parameters:")
        print("        [-m]    multiple dlg files")
        print("        [-N]    include information about best conformation")
        print("        [-v]    verbose output")
        print("        [-o pdbqt_filename] (output filename)")

    # process command arguments
    try:
        opt_list, args = getopt.getopt(sys.argv[1:], 'f:mNvo:h')
    except getopt.GetoptError as msg:
        print('write_lowest_energy_ligand.py: %s' %msg)
        usage()
        sys.exit(2)

    # initialize required parameters
    #-f: dlgfilename
    dlgfilename =  None
    # optional parameters
    verbose = None
    #-o outputfilename
    outputfilename = None
    #-m multiple
    multiple = False
    #-N include_best_docking_info
    include_best_docking_info = False

    #'f:vmo:h'
    for o, a in opt_list:
        #print "o=", o, " a=", a
        if o in ('-f', '--f'):
            dlgfilename = a
            if verbose: print('set dlgfilename to ', a)
        if o in ('-v', '--v'):
            verbose = True
            if verbose: print('set verbose to ', True)
        if o in ('-o', '--o'):
            outputfilename = a
            if verbose: print('set outputfilename to ', a)
        if o in ('-m', '--m'):
            multiple = True
            if verbose: print('read all dlgs: set multiple to ', multiple)
        if o in ('-N', '--N'):
            include_best_docking_info = True
            if verbose: print('read all dlgs: set include_best_docking_info to ', include_best_docking_info)
        if o in ('-h', '--'):
            usage()
            sys.exit()


    if not  dlgfilename:
        print('write_lowest_energy_ligand: dlgfilename must be specified.')
        usage()
        sys.exit()

    #what about nucleic acids???

    rms_tolerance = 5.0
    d = Docking()
    if not multiple:
        d.readDlg(dlgfilename)
    else:    
        fl = glob.glob("*.dlg")
        for n in fl:
            d.readDlg(n)
    if verbose: print('read ', dlgfilename)
    if multiple:
        coords = d.ligMol.allAtoms.coords[:]
        d.clusterer.rmsTool = RMSDCalculator(coords)
        d.clusterer.make_clustering(rms_tolerance)
    key0 = list(d.clusterer.clustering_dict.keys())[0]
    conf0 = d.clusterer.clustering_dict[key0][0][0]
    #figure out its dlg filename
    for dlo in d.dlo_list:
        if conf0 in dlo.conformations:
            conf0_filename = dlo.parser.filename
    d.ch.set_conformation(conf0)
    #binding_energy, energy, intermol_energy: -427.09
    #electrostatic_energy +3.19
    #run 5
    #vdw_hb_desolv_energy -430.28
    #translation:(-30.497708000000003, -30.385196, -39.081098)
    #quaternion0: (-0.001463, -0.010084, -0.171606, -0.985113)
    #qtn_nx,qtn_ny, qtn_nz:(-0.008511,-0.058662,-0.998242)
    #qtn_ang_deg -19.797555
    #conf0.trn_x,  conf0.trn_y, conf0.trn_z: (-3.077708, 3.784804, -21.501098)
    #conf0.translation: (-30.497708000000003, -30.385196, -39.081098)
    parser = d.ligMol.parser
    lines = []
    #have to add newline character to lines read from dlg
    for l in parser.allLines:
        l+= '\n'
        lines.append(l)
    parser.allLines = lines
    coords = d.ligMol.allAtoms.coords
    if outputfilename is None:
        outputfilename = d.ligMol.name  + '_BE.pdbqt'
    parser.write_with_new_coords(coords, outputfilename) 
    if verbose:
        print('wrote %s' %outputfilename)
    if include_best_docking_info:
        fptr = open(outputfilename)
        lines = fptr.readlines()
        fptr.close()
        optr = open(outputfilename, 'w')
        ostr = "REMARK best energy run %d in file:%s\n" %(conf0.run, conf0_filename)
        optr.write(ostr)
        ostr = "REMARK binding_energy %10.5f\n" %(conf0.binding_energy)
        optr.write(ostr)
        ostr = "REMARK translation %10.4f %10.4f %10.4f\n" %(conf0.trn_x, conf0.trn_y, conf0.trn_z)
        optr.write(ostr)
        ostr = "REMARK quaternion0 %10.4f %10.4f %10.4f %10.4f\n" %(conf0.quaternion0)
        optr.write(ostr)
        for l in lines:
            optr.write(l)
        optr.close()
    print("wrote ", outputfilename)
# To execute this command type:
# write_lowest_energy_ligand.py -f dlgfilename [-o outputfilename] -v




