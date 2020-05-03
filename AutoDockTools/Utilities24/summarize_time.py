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
# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/Utilities24/summarize_time.py,v 1.4.2.1 2016/02/11 09:24:08 annao Exp $
#
import os, glob
from AutoDockTools.DlgParser import DlgParser


if __name__ == '__main__':
    import sys
    import getopt


    def usage():
        "Print helpful, accurate usage statement to stdout."
        print("Usage: summarize_time.py -d directory")
        print()
        print("    Description of command...")
        print("         -d     directory")
        print("    Optional parameters:")
        print("        [-o]    output filename")
        print("                      (default is 'summary_of_time')")
        print("        [-v]    verbose output")


    # process command arguments
    try:
        opt_list, args = getopt.getopt(sys.argv[1:], 'd:o:vh')
    except getopt.GetoptError as msg:
        print('summarize_time.py: %s' %msg)
        usage()
        sys.exit(2)

    # initialize required parameters
    #-d: directory
    directory =  None

    # optional parameters
    #-o outputfilename
    outputfilename = "summary_of_time.txt"
    #-v: verbose
    verbose = None
    #-h: help
    help = None


    #'d:o:vh'
    for o, a in opt_list:
        #print "o=", o, " a=", a
        if o in ('-d', '--d'):
            directory = a
            if verbose: print('set directory to ', a)
        if o in ('-o', '--o'):
            outputfilename = a
            if verbose: print('set outputfilename to ', a)
        if o in ('-v', '--v'):
            verbose = True
            if verbose: print('set verbose to ', True)
        if o in ('-h', '--'):
            usage()
            sys.exit()

    if not  directory:
        print('summarize_time: directory must be specified.')
        usage()
        sys.exit()

    p = DlgParser()
    dlg_list = glob.glob(directory + '/*.dlg')
    total_time = 0
    for dlg in dlg_list:
        p.parse(dlg)
        total_time = total_time + p.total_time
        if verbose: print("dlg time=", total_time)
    stem = dlg_list[0].split('.')[0]
    mins = total_time/60.
    if verbose: print("total time in minutes =", mins)
    ostr = '%s: % 12.4f m \n' %(directory, mins)
    fptr = open(outputfilename, 'a')
    fptr.write(ostr)
    fptr.close()

# To execute this command type:
# summarize_time.py -d directory -o outputfilename -v

