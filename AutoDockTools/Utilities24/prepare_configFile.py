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
#  Modification date: 1/2/22, 5:10 PM                                                              #
#                                                                                                  #
# ##################################################################################################

#
#
#
# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/Utilities24/prepare_configFile.py,v 1.1 2012/10/11 21:14:47 rhuey Exp $
# $Id: prepare_configFile.py,v 1.1 2012/10/11 21:14:47 rhuey Exp $
#

import os.path
from MolKit import Read #?detect center of molecule+bounding box???
from AutoDockTools.DockingParameters import ConfigFileMaker



def usage():
    print("Usage: prepare_configFile.py -l pdbqt_file -r pdbqt_file")
    print("Description of command...")
    print("Prepare a config file for AutoDock vina: config.txt")
    print("  containing input parameters")
    print("    -l ligand_filename")
    print("    -r receptor_filename")
    print()
    print("Optional parameters:")
    print("    [-a center_x]") #default [0,0,0]?
    print("    [-b center_y]")
    print("    [-c center_z]")
    print("    [-X size_x]")  #default is bounding box of receptor?
    print("    [-Y size_y]")
    print("    [-Z size_z]")
    print("    [-o output config filename (default is 'config.txt')]")
    print("    [-L optionally, write log file]")
    print("    [-i template config filename]")
    print("    [-x flex_filename]")
    print("    [-p parameter_name=new_value]")
    print("    [-S write score_only config file ]")
    print("    [-L write local_only config file ]")
    print("    [-R write randomize_only config file ]")
    print("    [-F filename for new config file ]")
    print("    [-v] verbose output")
    print()


if __name__ == '__main__':
    import getopt
    import sys

    try:
        opt_list, args = getopt.getopt(sys.argv[1:], 'l:r:a:b:c:X:Y:Z:o:L:i:x:p:F:SLRFvh')
    except getopt.GetoptError as msg:
        print('prepare_configFile.py: %s' % msg)
        usage()
        sys.exit(2)

    receptor_filename = ligand_filename = None
    out = None
    reference_config = None
    flex_filename = None
    kw = {}
    verbose = None
    config_filename = ""
    kw = {}
    kw['parameters'] = []
    for o, a in opt_list:
        if verbose: print("o=", o, ' a=', a)
        if o in ('-h', '--'):
            usage()
            sys.exit()
        if o in ('-v', '--v'):
            verbose = 1
            if verbose: print('verbose output')
        if o in ('-r', '--r'):   #receptor filename
            receptor_filename = a
            if verbose: print('receptor_filename =', receptor_filename)
        if o in ('-x', '--x'):   #flex_filename
            flex_filename = a
            kw['flex'] = a
            if verbose: print('flex =', flex_filename)
        if o in ('-l', '--l'):   #ligand filename
            ligand_filename = a
            if verbose: print('ligand_filename =', ligand_filename)
        if o in ('-a', '--a'):   #center_x
            center_x = a
            kw['center_x'] = a
            if verbose: print('center_x =', center_x)
        if o in ('-b', '--b'):   #center_y
            center_y = a
            kw['center_y'] = a
            if verbose: print('center_y =', center_y)
        if o in ('-c', '--c'):   #center_z
            center_z = a
            kw['center_z'] = a
            if verbose: print('center_z =', center_z)
        if o in ('-X', '--X'):   #size_x
            size_x = a
            kw['size_x'] = a
            if verbose: print('size_x =', size_x)
        if o in ('-Y', '--Y'):   #size_y
            size_y = a
            kw['size_y'] = a
            if verbose: print('size_y =', size_y)
        if o in ('-Z', '--Z'):   #size_z
            size_z = a
            kw['size_z'] = a
            if verbose: print('size_z =', size_z)
        #optional : 'l:r:a:b:c:X:Y:Z:o:i:x:p:SLRvh'
        if o in ('-o', '--o'):   #out filename (pdbqt)
            out = a
            kw['out'] = a
            if verbose: print('out =', out)
        if o in ('-L', '--L'):   #log filename
            log = a
            kw['log'] = a
            if verbose: print('log =', log)
        if o in ('-i', '--i'):   #reference_config filename
            reference_config = a
            if verbose: print('reference_config =', reference_config)
        if o in ('-p', '--p'):   #parameter
            kw['parameter'].append(a)
            if verbose: print('kw =', kw)
    #[-S write score_only config file ],[-L write local_only config file ],[-R write randomize_only config file ]
        if o in ('-S', '--S'):   #score_only
            score_only = 1
            #parameters.append((o,1))
            kw['score_only'] = True
            if verbose: print('1: kw =', kw)
        if o in ('-L', '--L'):   #parameter_list_to_write
            kw['local_only'] = True
            local_search = 1
            if verbose: print('2: kw =', kw)
        if o in ('-R', '--R'):   #randomize_only
            randomize_only = True
            kw['randomize_only'] = True
            if verbose: print('3: kw =', kw)
        if o in ('-F', '--F'):   # config filename
            config_filename = a
            kw['config_filename'] = a
            if verbose: print('3: kw =', kw)
        if o in ('-h', '--'):
            usage()
            sys.exit()

    if (not receptor_filename) or (not ligand_filename):
        print("prepare_configFile.py: ligand and receptor filenames")
        print("                    must be specified.")
        usage()
        sys.exit()

    if len(config_filename)==0:
        config_filename = "config.txt"
    cm = ConfigFileMaker(receptor=receptor_filename,ligand=ligand_filename, **kw)
    if verbose: print("writing new config file:", config_filename)
    cm.write(config_filename)
    if verbose: print("results will be written by vina in ", out)


#prepare_configFile.py -l indinavir.pdbqt -r 1hsg.pdbqt

