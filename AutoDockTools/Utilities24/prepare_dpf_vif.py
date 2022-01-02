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
# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/Utilities24/prepare_dpf_vif.py,v 1.4 2012/03/01 23:09:58 rhuey Exp $
#

import os.path
from math import floor
from MolKit import Read
from AutoDockTools.DockingParameters import DockingParameters, DockingParameter42FileMaker, genetic_algorithm_list, \
    genetic_algorithm_local_search_list4, local_search_list4, \
    simulated_annealing_list4


def usage():
    print("Usage: prepare_dpf_vif.py -l pdbqt_file  -i template_filename")
    print("    -l ligand_filename")
    print("    -i template_filename")
    print()
    print("Optional parameters:")
    print("    [-v] verbose output")
    print("    [-o] output dpf_filename")
    print()
    print("Prepare a docking parameter file (DPF) for docking constrained models with AutoDock4.")
    print()


if __name__ == '__main__':
    import getopt
    import sys

    # print "    [-C] use two atoms with autodock_element LS for setting about"
    #    opt_list, args = getopt.getopt(sys.argv[1:], 'sLShvl:r:i:o:x:p:k:eC')
    try:
        opt_list, args = getopt.getopt(sys.argv[1:], 'vl:i:o:h')
    except getopt.GetoptError as msg:
        print('prepare_dpf_vif.py: %s' % msg)
        usage()
        sys.exit(2)

    ligand_filename = None
    template_filename = None
    dpf_filename = None
    verbose = None
    for o, a in opt_list:
        if o in ('-v', '--v'):
            verbose = 1
            if verbose: print('verbose output')
        if o in ('-l', '--l'):  # ligand filename
            ligand_filename = a
            if verbose: print('ligand_filename =', ligand_filename)
        if o in ('-i', '--i'):  # input reference
            template_filename = a
            if verbose: print('template_filename =', template_filename)
        if o in ('-o', '--o'):  # output filename
            dpf_filename = a
            if verbose: print('output dpf_filename =', dpf_filename)
        if o in ('-h', '--'):
            usage()
            sys.exit()

    if (not ligand_filename) or (not template_filename):
        print("prepare_dpf_vif.py: ligand filename and template filename")
        print("                    must be specified.")
        usage()
        sys.exit()

    if template_filename is not None:  # setup values by reading dpf
        tptr = open(template_filename)
        TEMPLATE_DPFlines = tptr.readlines()
        tptr.close()
    # check for format, determine center from LS atoms
    ligand = Read(ligand_filename)
    if not len(ligand):
        print(" problem reading ligand ", ligand_filename)
        sys.exit()
    lig = ligand[0]
    lig_types = set(lig.allAtoms.autodock_element)
    if "LS" not in lig_types:  # verify format for constrained docking
        print("@@no LS atoms found in ligand ", ligand_filename, "@@")
        sys.exit()

    found_cen_res = False
    LS_ats = lig.allAtoms.get(lambda x: x.autodock_element == 'LS')
    if not len(LS_ats):
        print(ligand_filename, " does not have any 'LS' atoms")
        raise Exception('improperly formatted ligand file')  # @@ just exit??
    if len(LS_ats):
        cen_res = ""
        len_ls_ats = len(LS_ats)
        parent_indices = list(map(int, LS_ats.parent.number))
        if verbose: print("parent_indices=", parent_indices)
        first_parent_index = int(parent_indices[0])
        if verbose: print("first parent_index is ", first_parent_index)
        mid_index = int(len_ls_ats / 2.)
        if verbose: print("mid_index=", mid_index)
        if len_ls_ats % 2 == 1:
            mid_index = int(floor(len_ls_ats / 2.))
            if verbose: print("NOW mid_index=", mid_index)
        cen_res_number = first_parent_index + mid_index
        if verbose: print("cen_res_number = ", cen_res_number)
        cen_res = lig.chains.residues.get(lambda x: int(x.number) == cen_res_number)
        if len(cen_res):
            cen_at = cen_res[0].atoms.get(lambda x: x.autodock_element == 'RP')[0]
            if verbose: print('cen_at coords= ', cen_at.coords, ' cen_at.parent=', cen_at.parent.full_name())
            about = cen_at.coords
            if verbose: print("ABOUT= ", about)
            found_cen_res = True
            about_value = cen_at.coords
            if verbose: print("about is ", dm.dpo['about']['value'])
        else:
            raise ligand_filename + ':problem finding center residue!'

    # OUTPUT with new about value
    # dm.write_dpf(dpf_filename, parameter_list, pop_seed)
    if dpf_filename is None:
        dpf_filename = "3ir2mod_" + lig.name + ".dpf"
        print("set dpf_filename to ", dpf_filename)
    optr = open(dpf_filename, 'w')
    for l in TEMPLATE_DPFlines:
        if l.find('about') == 0:
            new_l = 'about  %5.2f %5.2f %5.2f             # small molecule center\n' % (about[0], about[1], about[2])
            print("wrote new about:" + new_l)
            optr.write(new_l)
        elif l.find('move') == 0:
            new_l = 'move  ' + lig.parser.filename + "               # small molecule\n"
            optr.write(new_l)
            print("wrote new move:" + new_l)
        else:
            optr.write(l)
    optr.close()

# To execute this command type:
# prepare_dpf_vif.py -l pdbqt_file  -o output_dpf_filename -i template dpf_filename
