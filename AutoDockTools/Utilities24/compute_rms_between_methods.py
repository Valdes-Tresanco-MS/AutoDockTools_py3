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
# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/Utilities24/compute_rms_between_methods.py,v 1.6.2.1 2016/02/11 09:24:08 annao Exp $
#
import os, glob
import numpy
from math import sqrt

from MolKit import Read
from AutoDockTools.Docking import Docking
from mglutil.math.rmsd import RMSDCalculator


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
        print("Usage: compute_rms_between_methods.py -f reference directory -d directory")
        print("  Compare two docking methods by computing rms of docked conformations:LE and LE_LC")
        print("   By default, output written to 'summary_rms_results'")
        print()
        print("    Description of command...")
        print("         -f     reference directory")
        print("         -d     directory")
        print("    Optional parameters:")
        print("        [-t]    rmsd tolerance (default is 1.5)")
        print("        [-o]    output filename")
        print("                (default is 'summary_rms_results')")
        print("        [-v]    verbose output")


    # process command arguments
    try:
        opt_list, args = getopt.getopt(sys.argv[1:], 'f:d:o:t:vh')
    except getopt.GetoptError as msg:
        print('compute_rms_between_methods.py: %s' %msg)
        usage()
        sys.exit(2)

    # initialize required parameters
    #-f: refdirectory
    refdirectory =  None
    #-d: directory
    directory =  None
    #-t: rms_tolerance
    rms_tolerance =  1.5
    # optional parameters
    verbose = None
    #-o outputfilename
    outputfilename = "summary_rms_results"

    #'f:d:o:t:vh'
    for o, a in opt_list:
        #print "o=", o, " a=", a
        if o in ('-f', '--f'):
            refdirectory = a
            if verbose: print('set refdirectory to ', a)
        if o in ('-d', '--d'):
            directory = a
            if verbose: print('set directory to ', a)
        if o in ('-o', '--o'):
            outputfilename = a
            if verbose: print('set outputfilename to ', a)
        if o in ('-t', '--t'):
            rms_tolerance = float(a)
            if verbose: print('set rms_tolerance to ', a)
        if o in ('-v', '--v'):
            verbose = True
            if verbose: print('set verbose to ', True)
        if o in ('-h', '--'):
            usage()
            sys.exit()

    msg = ""
    if not refdirectory:
        msg += 'reference directory must be specified!\n'
        if not directory: 
            msg = 'Both reference directory and directory must be specified!\n'
    elif not directory:
        msg += 'directory must be specified!'

    if not refdirectory or not directory:
        print(msg)
        usage()
        sys.exit()

    #process docking in reference directory
    #read all the docking logs in reference directory as one Docking
    ref_dlg_list = glob.glob(refdirectory + '/*.dlg')
    ref_d = Docking()
    for dlg in ref_dlg_list:
        ref_d.readDlg(dlg)
    #setup rmsd tool
    coords = ref_d.ligMol.allAtoms.coords[:]
    ref_d.clusterer.rmsTool = RMSDCalculator(coords)
    ref_d.clusterer.make_clustering(rms_tolerance) 
    clust0 = ref_d.clusterer.clustering_dict[rms_tolerance][0]
    c = clust0[0] #lowest energy overall in reference docking
    ref_d.ch.set_conformation(c)
    ref_LE_coords = ref_d.ligMol.allAtoms.coords[:]
    ref_largest = clust0
    for clust in ref_d.clusterer.clustering_dict[rms_tolerance]:
        if verbose: print("current largest cluster len= ", len(clust))
        if len(clust)>len(ref_largest): 
            if verbose: print("resetting largest clust: now len=", len(clust))
            ref_largest = clust
    ref_d.ch.set_conformation(ref_largest[0])
    ref_LC_coords = ref_d.ligMol.allAtoms.coords[:]


    #process docking in test directory
    dlg_list = glob.glob(directory + '/*.dlg')
    d = Docking()
    for dlg in dlg_list:
        d.readDlg(dlg)
    #setup rmsd tool
    #coords = d.ligMol.allAtoms.coords[:]
    d.clusterer.rmsTool = RMSDCalculator(coords)
    d.clusterer.make_clustering(rms_tolerance) 
    clust0 = d.clusterer.clustering_dict[rms_tolerance][0]
    c = clust0[0]
    #SET d.ligMol to LOWEST ENERGY conf
    d.ch.set_conformation(c)
    LE_coords = d.ligMol.allAtoms.coords[:]
    # lowest energy in test v. lowest energy in ref
    rms_LE_LE = c.getRMSD(ref_LE_coords)
    dist_LE_LE = dist(ref_LE_coords, LE_coords) 
    #repeat for ref LC coords
    # lowest energy in test v. largest cluster lowest energy in ref
    rms_LE_LC = c.getRMSD(ref_LC_coords)
    dist_LE_LC = dist(ref_LC_coords, LE_coords) 

    # repeat for largest cluster 
    largest = clust0
    for clust in d.clusterer.clustering_dict[rms_tolerance]:
        if verbose: print("current largest cluster len= ", len(clust))
        if len(clust)>len(largest): 
            if verbose: print("resetting largest clust: now len=", len(clust))
            largest = clust
    conf_LC = largest[0]  #lowest energy conformation in largest cluster
    #SET d.ligMol to LARGEST CLUSTER's LOWEST ENERGY conf
    d.ch.set_conformation(conf_LC)
    LC_coords = d.ligMol.allAtoms.coords[:]
    # largest cluster lowest energy in test v. lowest energy in ref
    rms_LC_LE = conf_LC.getRMSD(ref_LE_coords)
    dist_LC_LE = dist(LC_coords, ref_LE_coords)
    #repeat for ref LC coords
    # largest cluster lowest energy in test v. largest cluster lowest energy in ref
    rms_LC_LC = conf_LC.getRMSD(ref_LC_coords)
    dist_LC_LC = dist(LC_coords, ref_LC_coords)
    

    first = not os.path.exists(outputfilename)
    if first:
        fptr = open(outputfilename, 'w')
        ostr= "# (in human terms: test_le->ref_le     test_le ->ref_lc        test_lc->ref_lc    test_lc ->ref_le)\n"
        fptr.write(ostr)
        ostr = "# ligand      rms_LE_LE dist_LE_LE  rms_LE_LC dist_LE_LC  rms_LC_LC dist_LC_LC  rms_LC_LE dist_LC_LE\n"
        fptr.write(ostr)
    else:
        fptr = open(outputfilename, 'a')
   
    ostr = "% 10s,% 10.4f,% 10.4f,% 10.4f,% 10.4f,% 10.4f,% 10.4f,% 10.4f,% 10.4f,\n" %(d.ligMol.name,\
                rms_LE_LE, dist_LE_LE, rms_LE_LC, dist_LE_LC, \
                rms_LC_LC, dist_LC_LC, rms_LC_LE, dist_LC_LE)
    fptr.write(ostr)
    fptr.close()
    

# To execute this command type:
# compute_rms_between_methods.py -f refdirectory -d directory -t rmsd tolerance 
#                         -v verbose
# NOTE: -f, -d and -t require arguments whereas the other, -v, sets a boolean

