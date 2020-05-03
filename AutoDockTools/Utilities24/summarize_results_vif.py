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
# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/Utilities24/summarize_results_vif.py,v 1.1.2.2 2016/03/11 00:52:51 annao Exp $
#
# $Id: summarize_results_vif.py,v 1.1.2.2 2016/03/11 00:52:51 annao Exp $ 
#
import os, glob

from MolKit import Read
from MolKit.hydrogenBondBuilder import HydrogenBondBuilder
from AutoDockTools.Docking import Docking
from mglutil.math.rmsd import RMSDCalculator
from PyAutoDock.InternalEnergy import InternalEnergy
from PyAutoDock.MolecularSystem import MolecularSystem
from PyAutoDock.AutoDockScorer import AutoDock4Scorer
import numpy


num_atoms = 0
num_hydrogens = 0
num_heavy_atoms = 0
be = 0.
be_lc = 0.
#vif
lp_ats = None
ls_ats = None
rp_ats = None
lp_energy = 0
ls_energy = 0
rp_energy = 0
report_energy_breakdown = False

def get_filename(d, conf):
    #find the filename of dlg containing this conformation
    for k in d.dlo_list:
        if conf in k.conformations:
            dlo_filename = os.path.splitext(k.parser.filename)[0]
            break
    return dlo_filename


def get_energy_breakdown(receptor_atoms, ligand_atoms, coords_to_use):
    ms = MolecularSystem()
    ms.add_entities(receptor_atoms)
    ms.add_entities(ligand_atoms)
    ms.set_coords(1, coords_to_use) 
    #ms.set_coords(1, c.getCoords()[:])
    ad_scorer = AutoDock4Scorer()
    ad_scorer.set_molecular_system(ms)
    score_list = ad_scorer.get_score_per_term()
    ostr = ""
    for score in score_list:
        ostr = ostr + "% 8.4f," %score
    return ostr[:-1]


def construct_hydrogen_bonds(receptor_atoms, ligand_atoms, verbose=False):
    atDict = HydrogenBondBuilder().build(receptor_atoms, ligand_atoms)
    num_hbonds = len(atDict)
    if verbose:
        print("found ", num_hbonds, " hydrogen bonds")
        for k,v in list(atDict.items()):
            print(k.full_name(), ' ',len(v), " hydrogen bond(s) to ", end=' ')
            for a in atDict[k]:
                print(a.full_name(), " ", end=' ')
            print()
    ostr = "%2d," % num_hbonds
    return ostr


def get_best_energy_info(d, rms_tolerance, build_hydrogen_bonds=False, \
                                report_energy_breakdown=False,
                                receptor_filename=None, refCoords=None,
                                subtract_internal_energy=False):
    global be
    ostr = "" 
    clust0 = d.clusterer.clustering_dict[rms_tolerance][0]
    c = clust0[0]
    d.ch.set_conformation(c)
    #find the filename of the best result
    dlo_filename = get_filename(d, c)
    #print "83: set dlo_filename to ", dlo_filename
    if verbose: print("set dlo_filename to ", dlo_filename)
    ostr = dlo_filename + ", "
    be = c.binding_energy
    #print "87: set be to ", be
    if subtract_internal_energy:
        be = c.binding_energy - c.total_internal
        #print "90: after subtracting_internal_energy be = ", be
    if refCoords:
        ostr += "%3d,%3d,%3d,% 8.4f,% 8.4f," %(len(d.clusterer.data), num_clusters, len(clust0),  be, c.getRMSD(refCoords=refCoords))
        #print "93: len(clust0)=", len(clust0), ", be=", be, " rmsd_LE=", c.getRMSD(refCoords=refCoords)
        #print "94: ostr=", ostr
        #ostr += "%3d,%3d,%3d,% 8.4f,% 8.4f," %(len(d.clusterer.data), num_clusters, len(clust0),  c.binding_energy, c.getRMSD(refCoords=refCoords))
    else:
        ostr += "%3d,%3d,%3d,% 8.4f,% 8.4f," %(len(d.clusterer.data), num_clusters, len(clust0),  be, c.getRMSD())
        print("99:  rmsd_LE=", c.getRMSD(), " @@@")
        #print "98: len(clust0)=", len(clust0), ", be=", be, " rmsd_LE=", c.getRMSD()
        #print "99: ostr=", ostr
        #ostr += "%3d,%3d,%3d,% 8.4f,% 8.4f," %(len(d.clusterer.data), num_clusters, len(clust0),  c.binding_energy, c.getRMSD())
    receptor = None
    if build_hydrogen_bonds:
        if not receptor_filename:
            print("receptor_filename must be specified in order to build_hydrogen_bonds")
            return 
        receptor = Read(receptor_filename)[0] 
        receptor.buildBondsByDistance()
        d.ligMol.buildBondsByDistance()
        ostr += construct_hydrogen_bonds(receptor.allAtoms, d.ligMol.allAtoms, verbose=verbose)
    if report_energy_breakdown:
        if not receptor_filename:
            print("receptor_filename must be specified in order to build_hydrogen_bonds")
        if not receptor:
            receptor = Read(receptor_filename)[0] 
            receptor.buildBondsByDistance()
        #build bonds if necessary 
        d.ligMol.buildBondsByDistance()
        #ostr += get_energy_breakdown(receptor.allAtoms, d.ligMol.allAtoms, c)
        ostr += get_energy_breakdown(receptor.allAtoms, d.ligMol.allAtoms, c.getCoords()[:])
    print("128: len(lp_ats)=", len(lp_ats))
    if len(lp_ats):
        try:
            ostr += " % 8.4f," % numpy.add.reduce(lp_ats.estat_energy+lp_ats.vdw_energy)
            print("131: lp_ats 1: ostr=", ostr, "\n")
        except:
            ostr += "         " 
    if len(ls_ats):
        ostr += " % 8.4f," % numpy.add.reduce(ls_ats.estat_energy+ls_ats.vdw_energy)
        print("135 ls_ats 2: ostr=", ostr, "\n")
    print("len(rp_ats)=", len(rp_ats))
    if len(rp_ats):
        ostr += " % 8.4f," % numpy.add.reduce(rp_ats.estat_energy+rp_ats.vdw_energy)
        print("139: ostr=", ostr, "\n")
    return ostr


def get_largest_cluster_info(d, rms_tolerance, build_hydrogen_bonds=False, \
                                    report_energy_breakdown=False,
                                    receptor_filename=None, refCoords=None,
                                    subtract_internal_energy=False):
    global be_lc
    largest = d.clusterer.clustering_dict[rms_tolerance][0]
    if verbose: print("set largest to ", len(largest))
    for clust in d.clusterer.clustering_dict[rms_tolerance]:
        if verbose: print("current largest cluster len= ", len(clust))
        if len(clust)>len(largest): 
            if verbose: print("resetting largest clust: now len=", len(clust))
            largest = clust
    c = largest[0]   #print info about the lowest energy member of this cluster
    d.ch.set_conformation(c)
    ostr = " "
    dlo_filename = get_filename(d, c)
    print("lc: dlo_filename= ", dlo_filename)
    ostr += dlo_filename + ","
    print("161 ostr=", ostr, "\n")
    be_lc = c.binding_energy
    print("163 set be_lc to ", be_lc, "\n")
    if subtract_internal_energy:
        be_lc = c.binding_energy - c.total_internal
    if refCoords:
        ostr += "%3d,%3d,%3d,% 8.4f,% 8.4f," %(len(d.clusterer.data), num_clusters, len(largest),  be_lc, c.getRMSD(refCoords=refCoords))
    else:
        ostr += "%3d,%3d,%3d,% 8.4f,% 8.4f," %(len(d.clusterer.data), num_clusters, len(largest),  be_lc, c.getRMSD())
    print("170 2c: ostr=", ostr, "\n")
    if verbose: print("set dlo_filename to ", dlo_filename)
    #update the coords only if you need to?
    receptor = None
    if build_hydrogen_bonds:
        if not receptor_filename:
            print("receptor_filename must be specified in order to build_hydrogen_bonds")
            return 
        receptor = Read(receptor_filename)[0] 
        receptor.buildBondsByDistance()
        ostr += construct_hydrogen_bonds(receptor.allAtoms, d.ligMol.allAtoms, verbose=verbose)
    if report_energy_breakdown:
        if not receptor_filename:
            print("receptor_filename must be specified in order to build_hydrogen_bonds")
            return 
        if not receptor:
            receptor = Read(receptor_filename)[0] 
            receptor.buildBondsByDistance()
        #ostr += get_energy_breakdown(receptor.allAtoms, d.ligMol.allAtoms, c)
        ostr += get_energy_breakdown(receptor.allAtoms, d.ligMol.allAtoms, c.getCoords()[:])
    if len(lp_ats):
        ostr += " % 8.4f," % numpy.add.reduce(lp_ats.estat_energy+lp_ats.vdw_energy)
        print("199  lp_ats:  ", ostr, "\n")
    if len(ls_ats):
        ostr += " % 8.4f," % numpy.add.reduce(ls_ats.estat_energy+ls_ats.vdw_energy)
        print("202 ls_ats:  ", ostr, "\n")
    if len(rp_ats):
        ostr += " % 8.4f," % numpy.add.reduce(rp_ats.estat_energy+rp_ats.vdw_energy)
        print("205 rp_ats:  ", ostr, "\n")
    return ostr


if __name__ == '__main__':
    import sys
    import getopt


    def usage():
        "Print helpful, accurate usage statement to stdout."
        print("Usage: summarize_results_vif.py -d directory")
        print()
        print("    Description of command...")
        print("         -d     directory")
        print("    Optional parameters:")
        print("        [-t]    rmsd tolerance (default is 1.0)")
        print("        [-f]    rmsd reference filename ")
        print("        (default is to use input ligand coordinates from docking log)")
        print("        [-b]    print best docking info only (default is print all)")
        print("        [-L]    print largest cluster info only (default is print all)")
        print("        [-B]    print best docking and largest cluster info only (default is print all)")
        print("        [-o]    output filename")
        print("                      (default is 'summary_of_results')")
        print("        [-a]    append to  output filename")
        print("                      (default is to open output filename 'w')")
        print("        [-k]    build hydrogen bonds and report number built")
        print("        [-e]    compute estat, vdw, hb + desolv energies and report breakdown")
        print("        [-r]    receptor filename")
        print("        [-i]    subtract internal energy")
        print("        [-p]    report depth of torsion tree")
        print("        [-v]    verbose output")


    # process command arguments
    try:
        opt_list, args = getopt.getopt(sys.argv[1:], 'd:o:t:f:r:bLBkaueipvh')
    except getopt.GetoptError as msg:
        print('summarize_results_vif.py: %s' %msg)
        usage()
        sys.exit(2)

    # initialize required parameters
    #-d: directory
    directory =  None
    #-t: rms_tolerance
    rms_tolerance =  1.0
    #-f: rms reference
    rms_reference = None
    # optional parameters
    verbose = None
    #-o outputfilename
    outputfilename = "summary_of_results"
    #-r receptor_filename
    receptor_filename = None
    #-b print best only
    print_best_only = False
    #-L print Largest cluster only
    print_Largest_only = False
    #-B print best and Largest cluster only
    print_best_and_Largest_only = False
    #-k build hydrogen bonds 
    build_hydrogen_bonds = False
    #-a append to outputfile
    append_to_outputfile = False
    #-e report_energy_breakdown report_energy_breakdown = False
    #-p report_torsion_tree_depth
    report_torsion_tree_depth = False
    #-i subtract_internal_energy
    subtract_internal_energy = False

    #'d:o:t:r:f:bLBkaueivh'
    for o, a in opt_list:
        #print "o=", o, " a=", a
        if o in ('-d', '--d'):
            directory = a
            if verbose: print('set directory to ', a)
        if o in ('-o', '--o'):
            outputfilename = a
            if verbose: print('set outputfilename to ', a)
        if o in ('-t', '--t'):
            rms_tolerance = float(a)
            if verbose: print('set rms_tolerance to ', a)
        if o in ('-f', '--f'):
            rms_reference = a
            if verbose: print('set rms_reference to ', a)
        if o in ('-r', '--r'):
            receptor_filename = a
            if verbose: print('set receptor_filename to ', a)
        if o in ('-b', '--b'):
            print_best_only = True
            if verbose: print('set print_best_only to ', True)
        if o in ('-L', '--L'):
            print_Largest_only = True
            if verbose: print('set print_Largest_only to ', True)
        if o in ('-B', '--B'):
            print_best_and_Largest_only = True
            print_best_only = True
            print_Largest_only = True
            if verbose: print('set print_best_and_Largest_only to ', True)
        if o in ('-k', '--k'):
            build_hydrogen_bonds = True
            if verbose: print('set build_hydrogen_bonds to ', True)
        if o in ('-a', '--a'):
            append_to_outputfile = True
            if verbose: print('set append_to_outputfile to ', True)
        #this option is apparently working currently:
        if o in ('-e', '--e'):
            report_energy_breakdown = True
            if verbose: print('set report_energy_breakdown to ', True)
        if o in ('-i', '--i'):
            subtract_internal_energy = True
            if verbose: print('set subtract_internal_energy to ', True)
        if o in ('-p', '--p'):
            report_torsion_tree_depth= True
            if verbose: print('set report_torsion_tree_depth to ', True)
        if o in ('-v', '--v'):
            verbose = True
            if verbose: print('set verbose to ', True)
        if o in ('-h', '--'):
            usage()
            sys.exit()

    if outputfilename=='summary_of_results':
        outputfilename = outputfilename + '_' + str(rms_tolerance)
    if build_hydrogen_bonds is True and receptor_filename is None:
        print('to build hydrogen bonds, receptor filename must be specified')
        sys.exit()
    if report_energy_breakdown is True and receptor_filename is None:
        print('to report energy breakdown, receptor filename must be specified')
        sys.exit()
    if not  directory:
        print('summarize_results_vif: directory must be specified.')
        usage()
        sys.exit()

    #read all the docking logs in as one Docking
    dlg_list = glob.glob(directory + '/*.dlg')
    d = Docking()
    for dlg in dlg_list:
        #code contributed by S.Forli
        try:  
            d.readDlg(dlg)
        except:
            print("Error reading the file \"%s\""  % dlg)
            #break
    if not d.dlo_list:
        print("Sorry, no valid DLG results found in the path : %s" % directory)
        sys.exit(1)

    #when >50 truncate the dockings at 50
    #d.dlo_list = d.dlo_list[:50]
    if receptor_filename is not None:
        receptor_filename = directory + '/' + receptor_filename
    #setup rmsd tool
    coords = d.ligMol.allAtoms.coords[:]
    lp_ats = d.ligMol.allAtoms.get(lambda x: x.autodock_element=='LP')
    ls_ats = d.ligMol.allAtoms.get(lambda x: x.autodock_element=='LS')
    rp_ats = d.ligMol.allAtoms.get(lambda x: x.autodock_element=='RP')
    refCoords = None
    if rms_reference is not None:
        file = directory + '/' + rms_reference
        ref = Read(file)
        if not len(ref): 
            print("unable to read ", rms_reference)
        else:
            ref = ref[0]
            coords = ref.allAtoms.coords
            refCoords = coords[:]
    d.clusterer.rmsTool = RMSDCalculator(coords)
    d.clusterer.make_clustering(rms_tolerance) 
    mode = 'w'
    first = True
    if append_to_outputfile:
        mode = 'a'
        first = not os.path.exists(outputfilename)
    fptr = open(outputfilename, mode)
    num_clusters = len(d.clusterer.clustering_dict[rms_tolerance])
    #output a header if file opened with mode 'w' or if this is first use of outputfile...
    if verbose: print("first is ", first)
    if first:
        tstr = ""
        if print_best_only:
            tstr += "lowestEnergy_dlgfn #runs #cl #LEC LE LP_E RP_E LS_E rmsd_LE "  #report energy for LP,RP,LS ats
            #if build_hydrogen_bonds:
            #    tstr+= "#hb "
            #if report_energy_breakdown:
            #    tstr+= "#ESTAT #HB #VDW #DSOLV "
        if print_Largest_only:
            if print_best_only:
                tstr += "   "
            tstr += "largestCl_dlgfn #runs #cl #LEC LE LP_E RP_E LS_E rmsd_LC "  #report energy for LP,RP,LS ats
            #if build_hydrogen_bonds:
            #    tstr+= "#hb "
            #if report_energy_breakdown:
            #    tstr+= "#ESTAT #HB #VDW #DSOLV "
        if not print_best_only and not print_Largest_only:
            tstr = "#dlgfn #incluster  #LE    LP_E   RP_E   LS_E #rmsd "
            #if build_hydrogen_bonds:
            #    tstr+= "#hb "
            #if report_energy_breakdown:
            #    tstr+= "#ESTAT #HB #VDW #DSOLV "
        tstr+= "#ats #hydrogens lig_eff" 
        #tstr+= "#ats #tors "
        #if report_torsion_tree_depth:
        #    tstr+= "#depthTT "
        #tstr+= "#h_ats #lig_eff "
        tstr+="\n"
        fptr.write(tstr)
        first = False
    #set up counts
    ligAts = d.ligMol.allAtoms
    total_num_atoms = len(ligAts)
    num_hydrogen_atoms = len(ligAts.get(lambda x: x.element=='H'))
    if verbose: print("num_hydrogen_atoms = ", num_hydrogen_atoms)
    num_heavy_atoms = total_num_atoms - num_hydrogen_atoms
    if verbose: print("num_heavy_atoms = ", num_heavy_atoms)
    ####lp_ats = ligAts.get(lambda x: x.autodock_element=='LP')
    ####ls_ats = ligAts.get(lambda x: x.autodock_element=='LS')
    ####rp_ats = ligAts.get(lambda x: x.autodock_element=='RP')
    #torTree_depth = d.ligMol.torTree.get_depth()
    ostr = ""
    if print_best_only or print_Largest_only:
        if print_best_only:
            if verbose: print('print_best_only is True')
            ostr += get_best_energy_info(d, rms_tolerance, build_hydrogen_bonds, report_energy_breakdown, receptor_filename, refCoords=refCoords, subtract_internal_energy=subtract_internal_energy)
        if print_Largest_only:
            if verbose: print('print_Largest_only is True')
            if len(ostr): ostr += ','
            ostr += get_largest_cluster_info(d, rms_tolerance, build_hydrogen_bonds, report_energy_breakdown, receptor_filename, refCoords=refCoords, subtract_internal_energy=subtract_internal_energy)
        #add number of atoms and number of torsions to output string and depth of torsion tree
        #ostr += "%4d,%2d" %(len(d.ligMol.allAtoms), d.ligMol.parser.keys.count('BRANCH'))
        ostr += " %4d," %(total_num_atoms)
        if  report_torsion_tree_depth:
            ostr += " %2d," %(torTree_depth)
        lig_eff = be/num_heavy_atoms
        ostr += " %2d,% 8.4f," %(num_hydrogen_atoms, lig_eff)
        ostr +='\n'
        try:
            fptr.write(ostr)
        except:
            print('except')
            fptr.write("\n")
    else:
        receptor = None
        if receptor_filename is not None:
            receptor = Read(receptor_filename)
            receptor = receptor[0]
            receptor.buildBondsByDistance()
        clust_ind = 0
        for clust in d.clusterer.clustering_dict[rms_tolerance]:
            #for each cluster in the clustering at this rms_tolerance:
            #output length of cluster and energies of the lowest energy result
            c = clust[0]
            d.ch.set_conformation(c)
            #write the dlg filename for the best in this cluster
            dlo_filename = get_filename(d, c)
            ostr = dlo_filename + "," 
            fptr.write(ostr)
            be = c.binding_energy
            if subtract_internal_energy:
                be = c.binding_energy - c.total_internal
            #print "477: clust_index=", clust_ind, "len(clust)=", len(clust), ", be=", be, " rmsd_LE=", c.getRMSD(refCoords=refCoords)
   
    
            #ostr = "%d,% 8.4f,% 8.4f" %(len(clust), be, c.getRMSD(refCoords=refCoords))
            ostr = "%d,% 8.4f," %(len(clust), be)
            if len(lp_ats):
                ostr += " % 8.4f," % numpy.add.reduce(lp_ats.estat_energy+lp_ats.vdw_energy)
            else:
                ostr += " %     0.0,"
            if len(rp_ats):
                ostr += "% 8.4f," % numpy.add.reduce(rp_ats.estat_energy+rp_ats.vdw_energy)
            else:
                ostr += " %     0.0,"
            #ostr = "%3d,% 8.4f,% 8.4f" %(len(clust), c.binding_energy, c.getRMSD(refCoords=refCoords))
            if len(ls_ats):
                ostr += " % 8.4f," % numpy.add.reduce(ls_ats.estat_energy+ls_ats.vdw_energy)
            else:
                ostr += " %     0.0,"
            #add rmsd
            ostr += " % 8.4f," % c.getRMSD()
            ostr += " %4d," %(total_num_atoms)
            lig_eff = be/num_heavy_atoms
            #ostr += ", %2d,% 8.4f" %(num_hydrogen_atoms, lig_eff)
            ostr += " %2d, % 8.4f" %(num_hydrogen_atoms, lig_eff)
            ostr += "\n"
            fptr.write(ostr)
            clust_ind+=1
    fptr.close()
    
    

# To execute this command type:
# summarize_results_vif.py -d directory -t rmsd tolerance 
# -b report info about lowest energy docked result 
# -f rms reference filename
# -L report largest cluster information 
# -B print best docking info  and largest cluster info only 
# -r receptor filename  -k build hydrogen bonds -o outputfilename 
# -a append to  outputfilename  -e report energy breakdown"
# -i subtract internal energy
# -p report depth of torsion tree 
# -v verbose output
