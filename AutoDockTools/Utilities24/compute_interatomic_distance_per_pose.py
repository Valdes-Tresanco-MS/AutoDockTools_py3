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
# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/Utilities24/compute_interatomic_distance_per_pose.py,v 1.4.4.2 2016/02/11 09:24:08 annao Exp $
#
import os, glob
import numpy
from math import sqrt

from MolKit import Read
from MolKit.molecule import MoleculeSet
from MolKit.stringSelector import CompoundStringSelector

from AutoDockTools.Docking import Docking



def dist(coords1, coords2):
    """return distance between two atoms, a and b.
    """
    npts = len(coords1)
    pt1 = numpy.add.reduce(coords1)/npts
    pt2 = numpy.add.reduce(coords2)/npts
    d = numpy.array(pt1) - numpy.array(pt2)
    return sqrt(numpy.sum(d*d))


if __name__ == '__main__':
    import sys
    import getopt


    def usage():
        "Print helpful, accurate usage statement to stdout."
        print("Usage: compute_interatomic_distance_per_pose.py -r hsg1:A:ASP29:N -l ind:I:IND201:N5 -d ind.dlg")
        print()
        print("  Description of command...")
        print("  Compute the distance between a specified receptor atom and a specified ligand atom ")
        print("      for each docked pose in specified dlg file:    ")
        print("      -r     full name of receptor atom ")
        print("      -l     full name of ligand atom ")
        print("      -d     dlg filename")
        print("  Optional parameters:")
        print("     [-o]    output filename")
        print("     [-v]    verbose output")


    # process command arguments
    try:
        opt_list, args = getopt.getopt(sys.argv[1:], 'r:l:d:o:vh')
    except getopt.GetoptError as msg:
        print('compute_interatomic_distance_per_pose.py: %s' %msg)
        usage()
        sys.exit(2)

    # initialize required parameters
    #-r: receptor atom fullname: eg 'hsg1:A:ASP29:N'
    ratom_name =  ""
    #-l: receptor atom fullname: eg 'ind:I:IND201:N5'
    latom_name =  ""
    #-d: dlg_filename
    dlg_filename =  ""
    #-o outputfilename
    outputfilename =  None #"rec_dockedLIG_"+dlg_filename +"_dists.txt"
    # optional parameters
    #-v detailed output 
    verbose = None

    #'r:l:d:o:vh'
    for o, a in opt_list:
        #print "o=", o, " a=", a
        if o in ('-r', '--r'):
            ratom_name = a
            if verbose: print('set receptor atom full_name to ', a)
        if o in ('-l', '--l'):
            latom_name = a
            if verbose: print('set ligand atom full_name to ', a)
        if o in ('-d', '--d'):
            dlg_filename = a
            if verbose: print('set dlg filename to ', a)
        if o in ('-o', '--o'):
            outputfilename = a
            if verbose: print('set outputfilename to ', a)
        if o in ('-v', '--v'):
            verbose = True
            if verbose: print('set verbose to ', True)
        if o in ('-h', '--'):
            usage()
            sys.exit()

    if not ratom_name:
        print('compute_interatomic_distance_per_pose: receptor atom full_name must be specified: eg: "hsg1:A:ASP29:N"')
        usage()
        sys.exit()

    if not latom_name:
        print('compute_interatomic_distance_per_pose: ligand atom full_name must be specified: eg: "ind:I:IND201:N5"')
        usage()
        sys.exit()

    if not dlg_filename:
        print('compute_interatomic_distance_per_pose: autodock dlg filename must be specified: eg: "ind.dlg"')
        usage()
        sys.exit()

    if outputfilename is None:
        outputfilename = "rec_dockedLIG_"+dlg_filename +"_dists.txt"

    #process docking in reference directory
    d = Docking()
    d.readDlg(dlg_filename)
    lig = d.ligMol
    #determine the receptor molecule to read in
    rec_file = ratom_name.split(":")[0] + ".pdbqt"
    msg = rec_file + " not found"
    assert os.path.exists(rec_file), msg 
    recs = Read(rec_file)
    rec = recs[0]
    # locate the receptor atom
    css = CompoundStringSelector()
    rec_ats = css.select(recs, ratom_name)[0]
    msg = ratom_name + " did not match exactly 1 atom in " + rec_file
    assert len(rec_ats)==1, msg
    rec_at = rec_ats[0]
    rec_at_coords = rec_at.coords
    if verbose: 
        print("found rec_at = ", rec_at.full_name(), end=' ') 
        print("with initial coords = ", rec_at.coords)
    # locate the ligand atom
    lig_ats = css.select(MoleculeSet([d.ligMol]), latom_name)[0]
    assert len(lig_ats)==1
    lig_at = lig_ats[0]
    print("found lig_at =>", lig_at.full_name())
    print("initial coords=", lig_at.coords)
    init_coords = lig_at.coords
    #open the outputfile
    if os.path.exists(outputfilename):
        fptr = open(outputfilename, 'a')
    else:
        fptr = open(outputfilename, 'w')
        ostr= "run     rec_at coords             lig_at coords        distance \n"
        fptr.write(ostr)
    print(" opened output file:", outputfilename)
    # set the pose to Run 1
    ctr = 1 
    print("about to step through ", len(d.ch.conformations), " docked conformations")
    for i in range(len(d.ch.conformations)):
        d.ch.set_conformation(d.ch.conformations[i])
        if verbose: print("now lig_at.coords=", lig_at.coords)
        new_dist = dist(rec_at.coords, lig_at.coords)
        if verbose: print(" new_dist =", new_dist)
        ostr = "%2d % 6.3f % 6.3f % 6.3f --% 8.3f % 8.3f % 8.3f = %8.4f \n" %(i+1, rec_at.coords[0], rec_at.coords[1], rec_at.coords[2], lig_at.coords[0], lig_at.coords[1], lig_at.coords[2], new_dist)
        fptr.write(ostr)
        ctr += 1
    if verbose: print("wrote ", ctr , " lines to ", outputfilename)
    fptr.close()
    

# To execute this command type:
# compute_interatomic_distance_per_pose.py -r receptor_atom_full_name -l ligand_atom_full_name -d dlg filename 
#           -o  output_filename  -v verbose
# NOTE: -r, -l and -d require arguments, as does the optional argument -o while -v sets a boolean

