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
# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/Utilities24/compute_interatomic_distance_per_vina_pose.py,v 1.3.4.2 2016/02/11 09:24:08 annao Exp $
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
        print("Usage: compute_interatomic_distance_per_vina_pose.py -r hsg1:A:ASP29:N -l ind:I:IND201:N5")
        print()
        print("  Description of command...")
        print("  Compute the distance between a specified receptor atom and a specified ligand atom ")
        print("      for each docked pose in specified vina result file:    ")
        print("      -r     full name of receptor atom, eg 'hsg1:A:ASP29:N'")
        print("      -l     full name of ligand atom in vina_result, eg 'ind_vina:I:IND201:C27'")
        print("  Optional parameters:")
        print("     [-o]    output filename default 'rec_dockedLIG_+vinafilename_dists.txt'")
        print("     [-v]    verbose output")


    # process command arguments
    try:
        opt_list, args = getopt.getopt(sys.argv[1:], 'r:l:o:vh')
    except getopt.GetoptError as msg:
        print('compute_interatomic_distance_per_vina_pose.py: %s' %msg)
        usage()
        sys.exit(2)

    # initialize required parameters
    #-r: receptor atom fullname: eg 'hsg1:A:ASP29:N'
    ratom_name =  "" 
    #-l: ligand atom fullname: eg 'ind:I:IND201:N5' @@includes vina_result filename
    latom_name =  ""
    # optional parameters
    #-o outputfilename
    outputfilename =  None  #"rec_dockedLIG_"+vina_result +"_dists.txt"
    #-v detailed output 
    verbose = None
    

    #'r:l:o:vh'
    for o, a in opt_list:
        #print "o=", o, " a=", a
        if o in ('-r', '--r'):
            ratom_name = a
            if verbose: print('set receptor atom full_name to ', a)
        if o in ('-l', '--l'):
            latom_name = a
            if verbose: print('set ligand atom full_name to ', a)
            try:
                vina_result = latom_name.split(':')[0] + ".pdbqt"
                if verbose: print('set vina_result filename to ', a)
            except:
                print("invalid ligand name ", a , "unable to find vina result_filename ", a + ".pdbqt") 
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
        print('compute_interatomic_distance_per_vina_pose: receptor atom full_name must be specified! eg: "hsg1:A:ASP29:N"')
        usage()
        sys.exit()

    if not latom_name:
        print('compute_interatomic_distance_per_vina_pose: ligand atom full_name must be specified! eg: "ind:I:IND201:N5"')
        usage()
        sys.exit()

    if outputfilename is None:
        outputfilename = "rec_dockedLIG_"+vina_result +"_dists.txt"

    #process vina_result in reference directory
    msg = 'invalid vina_result filename: ' + vina_result + " not found"
    assert os.path.exists(vina_result), msg 
    ligs = Read(vina_result, "conformations")
    msg = vina_result + " contains no valid molecules"
    assert len(ligs), msg
    lig = ligs[0]
    #determine the receptor molecule to read in
    rec_file = ratom_name.split(":")[0] + ".pdbqt"
    msg = rec_file + " file not found"
    assert os.path.exists(rec_file), msg 
    recs = Read(rec_file)
    msg = rec_file + " contains no valid molecules"
    assert len(recs), msg
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
    lig_ats = css.select(ligs, latom_name)[0]
    msg = latom_name + " did not match exactly 1 atom in " + ligs[0].name
    assert len(lig_ats)==1, msg
    lig_at = lig_ats[0]
    if verbose: 
        print("found lig_at =>", lig_at.full_name())
        print("with initial coords=", lig_at.coords)
    init_coords = lig_at.coords
    #open the outputfile
    printed_rec_coords = False
    if os.path.exists(outputfilename):
        fptr = open(outputfilename, 'a')
    else:
        fptr = open(outputfilename, 'w')
        ostr= "        rec_at coords             lig_at coords         distance\n"
        fptr.write(ostr)
        ostr= "run   %12s         %12s   (Angstroms)\n"%(ratom_name, latom_name)
        fptr.write(ostr)
    if verbose:
        print(" opened output file:", outputfilename)
    # set the pose to Run 1
    ctr = 0 
    if verbose:
        print("about to step through ", len(lig_at._coords), " docked conformations")
    for i in range(len(lig_at._coords)):
        lig.allAtoms.setConformation(i)
        if verbose: print("now lig_at.coords=", lig_at.coords)
        new_dist = dist(rec_at.coords, lig_at.coords)
        if verbose: print(" new_dist =", new_dist)
        if not printed_rec_coords:
            ostr = "%2d % 6.3f % 6.3f % 6.3f --% 8.3f % 8.3f % 8.3f = %8.4f \n" %(i+1, rec_at.coords[0], rec_at.coords[1], rec_at.coords[2], lig_at.coords[0], lig_at.coords[1], lig_at.coords[2], new_dist)
            printed_rec_coords = True
        else:
            ostr = "%2d                       --% 8.3f % 8.3f % 8.3f = %8.4f \n" %(i+1, lig_at.coords[0], lig_at.coords[1], lig_at.coords[2], new_dist)
        fptr.write(ostr)
        ctr += 1
    if verbose: print("wrote ", ctr , " distances to ", outputfilename)
    fptr.close()
    

# To execute this command type:
# compute_interatomic_distance_per_vina_pose.py -r receptor_atom_full_name -l ligand_atom_full_name 
#           [-o  output_filename  -v verbose] 
# NOTE: -r  and -l require arguments, as does the optional argument -o while -v sets a boolean

