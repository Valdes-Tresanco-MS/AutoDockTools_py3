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
# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/Utilities24/write_each_cluster_LE_conf.py,v 1.3 2011/12/27 22:41:53 rhuey Exp $
#
import os, glob

from MolKit import Read
from AutoDockTools.Docking import Docking
from mglutil.math.rmsd import RMSDCalculator
from string import strip




if __name__ == '__main__':
    import sys
    import getopt


    def usage():
        print("Usage: write_largest_cluster_ligand.py ")
        print("    This script does the following: ")
        print("         (1) read all the files with extension '.dlg' into one Docking")
        print("         (2) compute a clustering at the specified rms tolerance ")
        print("         (3) write the ligand with the coordinates of the ")
        print("             lowest-energy conformation in each cluster to a separate file") 
        print("    ")
        print("    Optional parameters:")
        print("        [-t]    rms_tolerance (default 2.0)")
        print("        [-o pdbqt_filename] (default ligandstem_clust#.pdbqt)")
        print("        [-v]    verbose output")

    # process command arguments
    try:
        opt_list, args = getopt.getopt(sys.argv[1:], 't:o:vh')
    except getopt.GetoptError as msg:
        print('write_each_cluster_LE_conf.py: %s' %msg)
        usage()
        sys.exit(2)

    # initialize required parameters
    # optional parameters
    verbose = None
    #-o outputfilestem
    outputfilestem = None
    #-t rms_tolerance
    rms_tolerance = 2.0

    #'t:o:vh'
    for o, a in opt_list:
        #print "o=", o, " a=", a
        if o in ('-o', '--o'):
            outputfilestem = a
            if verbose: print('set stem for outputfile to ', a)
        if o in ('-t', '--t'):
            rms_tolerance = float(a)
            if verbose: print('set rms_tolerance to ', a)
        if o in ('-v', '--v'):
            verbose = True
            if verbose: print('set verbose to ', True)
        if o in ('-h', '--'):
            usage()
            sys.exit()

    dlg_list = glob.glob('*.dlg')
    d = Docking()
    for dlg in dlg_list: 
        d.readDlg(dlg)
    d.clusterer.rmsTool = RMSDCalculator(d.ligMol.allAtoms.coords[:])
    d.clusterer.make_clustering(rms_tolerance)
    clustering = d.clusterer.clustering_dict[rms_tolerance]
    for clust in clustering:
        clustStr = str(clustering.index(clust))
        if verbose: print("processing clust number ", clustStr)
        #update the coordinates to those of conf with lowest energy in current cluster
        d.ch.set_conformation(clust[0])
        parser = d.ligMol.parser
        lines = []
        #have to add newline character to lines read from dlg
        for l in parser.allLines:
            l = strip(l)
            l+= '\n'
            lines.append(l)
        parser.allLines = lines
        coords = d.ligMol.allAtoms.coords
        if outputfilestem is None:
            parser.write_with_new_coords(coords, d.ligMol.name  + "_LE_clust"+clustStr +'.pdbqt') 
        else:
            parser.write_with_new_coords(coords, outputfilestem  + clustStr +'.pdbqt')
        if verbose: print('wrote %s' %outputfilestem)


# To execute this command type:
# write_each_cluster_LE_conf.py  [-t rms_tolerance, -o outputfilestem] -v




