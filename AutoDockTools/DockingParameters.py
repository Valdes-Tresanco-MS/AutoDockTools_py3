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
#  Modification date: 1/2/22, 5:32 PM                                                              #
#                                                                                                  #
# ##################################################################################################

#############################################################################
#
# Author: William LINDSTROM
#
# Copyright: M. Sanner TSRI 2000
#
#############################################################################


import os.path
import sys
from collections import UserDict

from AutoDockTools.atomTypeTools import AutoDock4_AtomTyper
from MolKit import Read
from .energyConstants import Rij, epsij


class DockingParameters(UserDict):
    def __init__(self, receptor_filename='', ligand_filename='', flexres_filename=''):
        UserDict.__init__(self)

        # if the docking parameters have been read from or written
        # to a file otherthan stdout,
        # then the following instance variables will be set:
        self.dpf_filename = ''
        self.dpf_written_filename = ''
        self.file_params = []

        # begin dictionary
        self['about'] = {
            'keyword': 'about',
            'default': [0., 0., 0.],
            'comment': "small molecule center",
            'value': [0., 0., 0.]
        }
        self['accs'] = {
            'keyword': 'accs',
            'default': 25000,  # 5/2012 previously 100
            'comment': "maximum number of accepted steps per cycle",
            'value': 25000
        }
        self['analysis'] = {
            'keyword': 'analysis',
            'default': 1,
            'comment': 'perform a ranked cluster analysis',
            'value': 1,  # true or false
        }
        self['axisangle0'] = {
            'keyword': 'axisangle0',
            'default': 'random',
            'comment': "initial orientation",
            'value': 'random'
        }
        self['barrier'] = {
            'keyword': 'barrier',
            'default': '',
            'comment': 'Torsion barrier-energy maximum value',
            'value': '65535'
        }
        self['bin_energies_by_rmsd'] = {  # same as investigate
            'keyword': 'bin_energies_by_rmsd',
            'default': '1000 100000 10',  # OutputEveryNTests,maxTests,NumLocalTests
            'comment': 'bin energies by RMSD from ref',
            'value': '1000 100000 10',  # OutputEveryNTests,maxTests,NumLocalTests
        }
        self['charmap'] = {  # @@ OBSOLETE ??
            'keyword': 'charmap',
            'default': '',
            'comment': "charmap",
            'value': ''
        }
        self['cluster'] = {
            'keyword': 'cluster',
            'default': '',
            'comment': 'structure binning',
            'value': ''
        }
        self['coliny'] = {  # @@ OBSOLETE ??
            'keyword': 'coliny',
            'default': '',  # algname,&nruns
            'comment': "coliny",
            'value': ''  # algname,&nrums
        }
        self['compute_unbound_extended'] = {
            'keyword': 'compute_unbound_extended',
            'default': '',
            'comment': "compute extended ligand energy",
            'value': ''
        }
        self['compute_unbound_extended_flag'] = {
            'keyword': 'compute_unbound_extended_flag',
            'default': 1,
            'comment': "whether to compute unbound ligand energy",
            'value': 1  # True/False
        }
        self['custom_parameter_file'] = {
            'keyword': 'custom_parameter_file',
            'default': 0,
            'comment': "use custom parameter library",
            'value': 0,
        }
        self['cycles'] = {
            'keyword': 'cycles',
            'default': 50,
            'comment': "number of temperature reduction cycles",
            'value': 50
        }
        self['desolvmap'] = {
            'keyword': 'desolvmap',
            'default': '',
            'comment': "desolvation map",
            'value': ''
        }
        self['dihe0'] = {
            'keyword': 'dihe0',
            'default': 'random',
            'comment': "initial dihedrals (relative) or random",
            'value': 'random'
        }
        self['dihrf'] = {
            'keyword': 'dihrf',
            'default': 1.0,
            'comment': "per cycle reduction factor for dihedrals",
            'value': 1.0
        }
        self['do_cpso'] = {  # do_cpso, do_pso
            'keyword': 'do_cpso',
            'default': 50,
            'comment': "do this many cpso runs",
            'value': 50
        }
        self['do_pso'] = {  # do_pso
            'keyword': 'do_pso',
            'default': 50,
            'comment': "do this many pso runs",
            'value': 50
        }
        self['do_global_only'] = {  # ga_only_run
            'keyword': 'do_global_only',
            'default': 50,
            'comment': "do this many GA runs",
            'value': 50
        }
        self['do_local_only'] = {  # ls_run
            'keyword': 'do_local_only',
            'default': 50,
            'comment': "do this many LS runs",
            'value': 50
        }
        self['do_gals'] = {  # gals_run, do_gals
            'keyword': 'ga_run',
            'default': 10,
            'comment': "do this many hybrid GA-LS runs",
            'value': 10
        }
        self['dstep'] = {
            'keyword': 'dstep',
            'default': 5.0,  # 5/2012 previously 50.
            'comment': "torsion step/deg",
            'value': 5.0
        }
        self['elecmap'] = {
            'keyword': 'elecmap',
            'default': '',
            'comment': "electrostatics map",
            'value': ''
        }
        self['epdb'] = {
            'keyword': 'epdb',
            'default': "",
            'comment': "evaluate ligand specified with move command",
            'value': ""
        }
        self['epdb_flag'] = {
            'keyword': 'epdb_flag',
            'default': 0,
            'comment': "whether to include epdb keyword",
            'value': 0  # true/false
        }
        self['e0max'] = {
            'keyword': 'e0max',
            'default': [0.0, 10000],
            'comment': "max initial energy; max number of retries",
            'value': [0.0, 10000]
        }
        self['extnrg'] = {
            'keyword': 'extnrg',
            'default': 1000.0,
            'comment': "external grid energy",
            'value': 1000.0
        }
        self['fld'] = {
            'keyword': 'fld',
            'default': None,
            'comment': "grid_data_file",
            'value': None
        }
        self['flex'] = {
            'keyword': 'flex',
            'default': '',
            'comment': "flexible side-chains, cannot translate",
            'value': ''
        }
        self['flexres_flag'] = {
            'keyword': 'flexres_flag',
            'default': 0,
            'comment': "whether to include flexres file",
            'value': 0
        }
        self['flexres'] = {
            'keyword': 'flexres',
            'default': flexres_filename,
            'comment': "file containing flexible residues",
            'value': flexres_filename
        }
        self['flexible_residues'] = {
            'keyword': 'flexible_residues',
            'default': flexres_filename,
            'comment': "file containing flexible residues",
            'value': flexres_filename
        }
        self['fmap'] = {
            'keyword': 'fmap',
            'default': '',
            'comment': "floating map",
            'value': ''
        }
        self['ga_boltzman_selection'] = {
            'keyword': 'ga_boltzman_selection',
            'default': 0.0,
            'comment': "Boltzman_selection is not yet implemented",
            'value': 0.0
        }
        self['ga_cauchy_alpha'] = {
            'keyword': 'ga_cauchy_alpha',
            'default': 0.0,
            'comment': "Alpha parameter of Cauchy distribution",
            'value': 0.0
        }
        self['ga_cauchy_beta'] = {
            'keyword': 'ga_cauchy_beta',
            'default': 1.0,
            'comment': "Beta parameter Cauchy distribution",
            'value': 1.0
        }
        self['ga_crossover_rate'] = {
            'keyword': 'ga_crossover_rate',
            'default': 0.80,
            'comment': "rate of crossover",
            'value': 0.80
        }
        self['ga_crossover_mode_flag'] = {
            'keyword': 'ga_crossover_mode_flag',
            'default': 0,
            'comment': "mode of crossover",
            'value': 0  # False/True
        }
        self['ga_crossover_mode'] = {
            'keyword': 'ga_crossover_mode',
            'default': 'twopt',
            'comment': 'mode of crossover',
            'value': 'twopt'  # onept,twopt,uniform,arithmetic,branch
        }
        self['ga_elitism'] = {
            'keyword': 'ga_elitism',
            'default': 1,
            'comment': "number of top individuals to survive to next generation",
            'value': 1
        }
        self['ga_high'] = {
            'keyword': 'ga_high',
            'default': 1,
            'comment': " ga high is not yet implemented ",
            'value': 1
        }
        self['gausstorcon'] = {  # DPF_GAUSSTORCON, DPF_HARDTORCON
            'keyword': 'gausstorcon',
            'default': "0 0 0",
            'comment': " gausstorcon is not yet implemented ",
            'value': 1
        }
        self['ga_linear_ranking_selection'] = {  # _probability_ratio
            'keyword': 'ga_linear_ranking_selection',
            'default': 2.0,
            'comment': "linear_ranking_selection",
            'value': 2.0
        }
        self['ga_low'] = {  # sets low
            'keyword': 'ga_low',
            'default': 0,
            'comment': " ga low is not yet implemented ",
            'value': 0
        }
        self['ga_mutation_rate'] = {
            'keyword': 'ga_mutation_rate',
            'default': 0.02,
            'comment': "rate of gene mutation",
            'value': 0.02
        }
        self['ga_num_evals'] = {
            'keyword': 'ga_num_evals',
            'default': 2500000,
            'comment': "maximum number of energy evaluations",
            'value': 2500000
        }
        self['ga_num_generations'] = {
            'keyword': 'ga_num_generations',
            'default': 27000,
            'comment': "maximum number of generations",
            'value': 27000
        }
        self['ga_pop_size'] = {
            'keyword': 'ga_pop_size',
            'default': 150,
            'comment': "number of individuals in population",
            'value': 150
        }
        self['ga_proportional_selection'] = {  # gals_run, do_gals
            'keyword': 'ga_proportional_selection',
            'default': 0,  # enum Selection_Mode { Proportional=0, LinearRanking=1, Tournament=2, Boltzmann=3 };
            'comment': "selection mode used in GA and LGS searchs",
            'value': 0
        }
        self['ga_run'] = {  # gals_run, do_gals
            'keyword': 'ga_run',
            'default': 10,
            'comment': "do this many hybrid GA-LS runs",
            'value': 10
        }
        self['ga_termination'] = {  # not yet implemented
            'keyword': 'ga_termination',  # energy 0.1; evals 25000; or time 120 s
            'default': 0,
            'comment': "selection mode used in GA and LGS searchs",
            'value': 0
        }  # or ga_termination_criterion
        self['ga_tournament_selection'] = {  # @@not yet implemented@@
            'keyword': 'ga_tournament_selection',
            'default': 0,  # ?
            'comment': "selection mode used in GA and LGS searchs",
            'value': 0
        }
        self['ga_window_size'] = {
            'keyword': 'ga_window_size',
            'default': 10,
            'comment': '',
            'value': 10
        }
        self['gals_run'] = {  # gals_run, do_gals
            'keyword': 'ga_run',
            'default': 10,
            'comment': "do this many hybrid GA-LS runs",
            'value': 10
        }
        self['include_1_4_interactions'] = {
            'keyword': 'include_1_4_interactions',
            'default': 1.0,
            'comment': "include internal 1-4 interactions",
            'value': 1.0  # weight
        }
        self['include_1_4_interactions_flag'] = {
            'keyword': 'include_1_4_interactions_flag',
            'default': 0,
            'comment': "whether to include internal 1-4 interactions",
            'value': 0  # true/false
        }
        self['intelec'] = {
            'keyword': 'intelec',
            'default': 0,
            'comment': "calculate internal electrostatics",
            'value': 0  # true/false
        }
        self['intelec4'] = {
            'keyword': 'intelec4',
            'default': 0,
            'comment': "calculate internal electrostatics",
            'value': '0.1465'  # from latest force field
        }
        self['intnbp_r_eps'] = {
            'keyword': 'intnbp_r_eps',
            'default': None,
            'comment': '',  # to be set at write-time
            'value': None,
        }
        self['investigate'] = {
            'keyword': 'investigate',
            'default': '1000 100000 10',  # OutputEveryNTests,maxTests,NumLocalTests
            'comment': 'bin energies by RMSD from ref',
            'value': '1000 100000 10',  # OutputEveryNTests,maxTests,NumLocalTests
        }
        self['ligand_is_not_inhibitor'] = {
            'keyword': 'ligand_is_not_inhibitor',
            'default': 1,
            'comment': "Kd calculated if not inhibitor",  # main.cc ~line 3802
            'value': 1,
        }
        self['ligand_types'] = {
            'keyword': 'ligand_types',
            'default': 'C A HD OA NA',
            'comment': "atoms types in ligand",
            'value': 'C A HD OA NA',
        }
        self['linear_schedule'] = {
            'keyword': 'linear_schedule',
            'default': 1,
            'comment': "use linear, arithmetic temperature reduction",
            'value': 1  # true/false
        }
        self['linsched'] = {
            'keyword': 'linsched',
            'default': 1,
            'comment': "use linear, arithmetic temperature reduction",
            'value': 1  # true/false
        }
        self['ls_run'] = {  # do_local_only
            'keyword': 'ls_run',  # @@ change to do_local_only??
            'default': 50,
            'comment': "do this many LS runs",
            'value': 50
        }
        self['ls_search_freq'] = {
            'keyword': 'ls_search_freq',
            'default': 0.06,
            'comment': "probability of performing local search on individual",
            'value': 0.06
        }
        self['map'] = {
            'keyword': 'map',
            'default': '',
            'comment': "atom-specific affinity map",
            'value': ''
        }
        self['move'] = {  # 'ligand' is same @@
            'keyword': 'move',
            'default': ligand_filename,
            'comment': "small molecule",
            'value': ligand_filename
        }
        self['ndihe'] = {
            'keyword': 'ndihe',
            'default': 0,
            'comment': "number of active torsions",
            'value': 0
        }
        self['outlev'] = {
            'keyword': 'outlev',
            'default': 1,
            'comment': "diagnostic output level",
            'value': 1
        }
        self['output_pop_file'] = {
            'keyword': 'output_pop_file',
            'default': "gen_state.txt",
            'comment': "generation state file",
            'value': "gen_state.txt"
        }
        self['output_population_statistics'] = {  # must be after outlev to have effect
            'keyword': 'output_population_statistics',  # @@ ASK MP
            'default': "",  # c_mode_str, everyNgens, everyNevals
            'comment': "",
            'value': ""
        }
        self['output_resnum_as'] = {
            # pdbqt format for residues in dlgs
            # *  default is to keep the residue number string from input
            # *  possible values: 'resnum' and 'runnum'
            # *  intended to replace '-k' option in setflags.cc
            'keyword': 'output_resnum_as',  # @@ ASK MP
            'default': 'resnum',  # c_mode_str, everyNgens, everyNevals
            'comment': 'pdbqt format for residues in dlgs',
            'value': 'resnum'
        }
        self['parameter_file'] = {  # ??parameter_library??
            'keyword': 'parameter_file',
            'default': 'AD4.1_bound.dat',
            'comment': "parameter library filename ",
            'value': 'AD4.1_bound.dat',
        }
        self['pso_adaptive_velocity'] = {
            'keyword': 'pso_adaptive_velocity',
            'default': 0,  # False
            'comment': "pso_adaptive_velocity ",
            'value': 0,  # False
        }
        self['pso_c1'] = {
            'keyword': 'pso_c1',
            'default': 2.05,  # from autodock/pso.h
            'comment': "pso_c1 ",
            'value': 2.05,  # Falsefrom autodock/pso.h
        }
        self['pso_c2'] = {
            'keyword': 'pso_c2',
            'default': 2.05,  # from autodock/pso.h
            'comment': "pso_c2 ",
            'value': 2.05,  # Falsefrom autodock/pso.h
        }
        self['pso_interpolate_as_scalars'] = {
            'keyword': 'pso_interpolate_as_scalars',
            'default': 1,  # True in autodock/pso.h
            'comment': "pso_interpolate_as_scalars ",
            'value': 1,  # True in autodock/pso.h
        }
        self['pso_k'] = {
            'keyword': 'pso_k',  # s
            'default': 4,  # pso_K set to 4 in autodock/pso.h
            'comment': "pso_k ",
            'value': 4,  # pso_K set to 4 in autodock/pso.h
        }
        self['pso_neighbors'] = {
            'keyword': 'pso_neighbors',  #
            'default': 4,  # pso_K set to 4 in autodock/pso.h
            'comment': "pso_neighbors ",
            'value': 4,  # pso_K set to 4 in autodock/pso.h
        }
        self['pso_neighbors_dynamic'] = {
            'keyword': 'pso_neighbors_dynamic',  #
            'default': False,  # set in autodock/pso.h
            'comment': "pso_neighbors_dynamic ",
            'value': False,  # pso_K set to 4 in autodock/pso.h
        }
        self['pso_neighbors_symmetric'] = {  # MP not yet implemented
            'keyword': 'pso_neighbors_symmetric',  #
            'default': False,  # set in autodock/pso.h
            'comment': "pso_neighbors_symmetric ",
            'value': False,  # pso_K set to 4 in autodock/pso.h
        }
        self['pso_qvmax'] = {  #
            'keyword': 'pso_qvmax',  #
            'default': 1.0,  # set in main.cc
            'comment': "pso_qvmax ",
            'value': 1.0,  # set in main.cc
        }
        self['pso_random_by_dimension'] = {
            'keyword': 'pso_random_by_dimension',  #
            'default': True,  # set in pso.h
            'comment': "pso_random_by_dimension ",
            'value': True,  # set in pso.h
        }
        self['pso_regenerate_at_limit'] = {
            'keyword': 'pso_regenerate_at_limit',  #
            'default': True,  # set in pso.h
            'comment': "pso_regenerate_at_limit ",
            'value': True,  # set in pso.h
        }
        self['pso_rvmax'] = {
            'keyword': 'pso_rvmax',  #
            'default': 50.0,  # set in main.cc
            'comment': "pso_rvmax ",
            'value': 50.0,  # set in main.cc
        }
        self['pso_stage2constriction'] = {
            'keyword': 'pso_stage2constriction',  #
            'default': False,  # set in pso.h
            'comment': "pso_stage2constriction ",
            'value': False,  # set in pso.h
        }
        self['pso_tvmax'] = {  #
            'keyword': 'pso_tvmax',  #
            'default': 2.0,  # set in main.cc
            'comment': "pso_tvmax ",
            'value': 2.0,  # set in main.cc
        }
        self['pso_w_end'] = {  #
            'keyword': 'pso_w_end',  #
            'default': 0.4,  # set in pso.h
            'comment': "pso_w_end ",
            'value': 0.4,  # set in pso.h
        }
        self['pso_w_start'] = {  #
            'keyword': 'pso_w_start',  #
            'default': 0.9,  # set in pso.h
            'comment': "pso_w_start ",
            'value': 0.9,  # set in pso.h
        }
        self['psw_trans_scale'] = {
            'keyword': 'psw_trans_scale',
            'default': 1,
            'comment': "pseudo sw translation rho scale",
            'value': 1,
        }
        self['psw_rot_scale'] = {
            'keyword': 'psw_rot_scale',
            'default': 0.05,
            'comment': "pseudo sw rotation rho scale",
            'value': 0.05,
        }
        self['psw_tors_scale'] = {
            'keyword': 'psw_tors_scale',
            'default': 0.1,
            'comment': "pseudo sw torsion rho scale",
            'value': 0.1,
        }
        self['qstep'] = {
            'keyword': 'qstep',
            'default': 5.0,  # 5/2012 previously 50.
            'comment': "quaternion step/deg",
            'value': 5.0
        }
        self['quarf'] = {
            'keyword': 'quarf',
            'default': 1.0,
            'comment': "per cycle reduction factor for quaternions",
            'value': 1.0
        }
        self['quat0'] = {
            'keyword': 'quat0',
            'default': 'random',
            'comment': "initial quaternion",
            'value': 'random'
        }
        self['quaternion0'] = {
            'keyword': 'quaternion0',
            'default': 'random',
            'comment': "initial orientation",
            'value': 'random'
        }
        self['rejs'] = {
            'keyword': 'rejs',
            'default': 25000,  # 5/2012 previously 100
            'comment': "maximum number of rejected steps per cycle",
            'value': 25000
        }
        self['reorient_flag'] = {
            'keyword': 'reorient_flag',
            'default': 0,  # do not routinely include reorient
            'comment': "whether to reorient ligand at start of each run and how",
            'value': 0,
        }
        self['reorient'] = {
            'keyword': 'reorient',
            'default': 'random',
            'comment': "initial orientation of ligand",
            'value': 'random'
        }
        self['rmsatoms'] = {
            'keyword': 'rmsatoms',
            'default': 'ligand_only',  # other choice is 'all'
            'comment': "cluster reference pdbqt file",
            'value': 'ligand_only',
        }
        self['rmsatoms_flag'] = {
            'keyword': 'rmsatoms_flag',
            'default': 0,  # do not routinely include this keyword
            'comment': "what to use for cluster reference pdbqt file",
            'value': 0,
        }
        self['rmsref'] = {
            'keyword': 'rmsref',
            'default': os.path.basename(ligand_filename),
            'comment': "cluster reference file",
            'value': os.path.basename(ligand_filename)
        }
        self['rmsref_flag'] = {
            'keyword': 'rmsref_flag',
            'default': 0,
            'comment': "whether rmsref is present",
            'value': 0,
        }
        self['rmsmode'] = {  # atype|unique_pair|heavy_atoms_only
            'keyword': 'rmsmode',
            'default': 'atype',
            'comment': "all pairs|uniq_pair|heavy_atoms_only",
            'value': 'atype'  # all pairs
        }
        self['rmsnosym'] = {  # atype|unique_pair|heavy_atoms_only
            'keyword': 'rmsnosym',
            'default': 'atype',
            'comment': "B_symmetry_flag",  # @@
            'value': 'atype'  # all pairs
        }
        self['rmstol'] = {
            'keyword': 'rmstol',
            'default': 2.0,
            'comment': "cluster_tolerance/A",
            'value': 2.0
        }
        self['rt0'] = {
            'keyword': 'rt0',
            'default': 616.0,  # 5/2012 previously 1000.
            'comment': "initial annealing temperature (times gas constant)",
            'value': 616.0
        }
        self['rtrf'] = {
            'keyword': 'rtrf',
            'default': 0.95,
            'comment': "annealing temperature reduction factor",
            'value': 0.95
        }
        self['runs'] = {
            'keyword': 'runs',
            'default': 10,
            'comment': "",
            'value': 10
        }
        self['seed'] = {
            'keyword': 'seed',
            'default': ['pid', 'time'],
            'comment': "seeds for random generator",
            'value': ['pid', 'time']
        }
        self['select'] = {
            'keyword': 'select',
            'default': 'm',
            'comment': "state selection flag: (m)inimum or (l)ast state",
            'value': 'm'
        }
        self['set_ga'] = {
            'keyword': 'set_ga',
            'default': 1,
            'comment': "set the above parameters for GA or LGA",
            'value': 1  # true/false
        }
        self['set_sw1_flag'] = {
            'keyword': 'set_sw1_flag',
            'default': 0,
            'comment': "set ls to sw",
            'value': 0  # true/false
        }
        self['set_sw1'] = {
            'keyword': 'set_sw1',
            'default': 0,
            'comment': "set the above Solis & Wets parameters",
            'value': 0  # true/false
        }
        self['set_psw1_flag'] = {
            'keyword': 'set_psw1_flag',
            'default': 1,
            'comment': "set ls to psw",
            'value': 1  # true/false
        }
        self['set_psw1'] = {
            'keyword': 'set_psw1',
            'default': 1,
            'comment': "set the above pseudo-Solis & Wets parameters",
            'value': 1  # true/false
        }
        self['simanneal'] = {
            'keyword': 'simanneal',
            'default': 1,
            'comment': "do as many SA runs as set by runs keyword above",
            'value': 1  # true/false
        }
        self['sw_lb_rho'] = {
            'keyword': 'sw_lb_rho',
            'default': 0.01,
            'comment': "lower bound on rho",
            'value': 0.01
        }
        self['sw_max_fail'] = {
            'keyword': 'sw_max_fail',
            'default': 4,
            'comment': "consecutive failures before changing rho",
            'value': 4
        }
        self['sw_max_its'] = {
            'keyword': 'sw_max_its',
            'default': 300,
            'comment': "iterations of Solis & Wets local search",
            'value': 300
        }
        self['sw_max_succ'] = {
            'keyword': 'sw_max_succ',
            'default': 4,
            'comment': "consecutive successes before changing rho",
            'value': 4
        }
        self['sw_rho'] = {
            'keyword': 'sw_rho',
            'default': 1.0,
            'comment': "size of local search space to sample",
            'value': 1.0
        }
        self['torsdof'] = {
            'keyword': 'torsdof',
            'default': [0, 0.3113],
            'comment': "torsional degrees of freedom and coeffiecent",
            'value': [0, 0.3113]
        }
        self['torsdof4'] = {
            'keyword': 'torsdof4',
            'default': [0],
            'comment': "torsional degrees of freedom",
            'value': [0]
        }
        self['tran0'] = {
            'keyword': 'tran0',
            'default': 'random',
            'comment': "initial coordinates/A or random",
            'value': 'random'
        }
        self['trnrf'] = {
            'keyword': 'trnrf',
            'default': 1.0,
            'comment': "per cycle reduction factor for translation",
            'value': 1.0
        }
        self['tstep'] = {
            'keyword': 'tstep',
            'default': 0.2,  # 5/2012 previously 2.0
            'comment': "translation step/A",
            'value': 0.2
        }
        self['types'] = {
            'keyword': 'types',
            'default': 'ACONSH',
            'comment': "atom type names",
            'value': 'ACONSH'
        }
        self['unbound'] = {  # 4.2
            'keyword': 'unbound',
            'default': 0.0,
            'comment': "free energy of ligand's unbound state ",
            'value': 0.0
        }
        self['unbound_flag'] = {
            'keyword': 'unbound_flag',
            'default': 0,
            'comment': "whether to include a custom unbound value",
            'value': 0
        }
        self['unbound_energy'] = {  # 4.2
            'keyword': 'unbound_energy',
            'default': 0.0,
            'comment': "free energy of ligand's unbound state ",
            'value': 0.0
        }
        self['unbound_energy_flag'] = {
            'keyword': 'unbound_energy_flag',
            'default': 0,
            'comment': "whether to include a custom unbound_energy value",
            'value': 0
        }
        self['unbound_intnbp_coeffs'] = {
            'keyword': 'unbound_intnbp_coeffs',
            'default': "20. 0. 1 2",
            'comment': "intnbp ff coeffs",
            'value': "20. 0. 1 2"  # two floats and two integers:
        }
        self['unbound_intnbp_coeffs_flag'] = {
            'keyword': 'unbound_intnbp_coeffs_flag',
            'default': 0,
            'comment': "whether to include custom intnbp values",
            'value': 0
        }

        self['unbound_model'] = {
            'keyword': 'unbound_model',
            'default': "bound",  # 4.1 default
            'comment': "state of unbound ligand",
            'value': "bound",  # possible values: bound, extended (AD4.0), compact(n/a)
        }
        self['unbound_model_flag'] = {
            'keyword': 'unbound_model_flag',
            'default': 0,  #
            'comment': "whether to include unbound_model keyword",
            'value': 0,  #
        }

        self['autodock_parameter_version'] = {
            'keyword': 'autodock_parameter_version',
            'default': 4.2,
            'comment': "used by autodock to validate parameter set",
            'value': 4.2
        }

        self['write_all'] = {
            'keyword': 'write_all',
            'default': "",
            'comment': "write all conformations in a cluster",
            'value': ""  # true/false
        }
        self['write_all_flag'] = {
            'keyword': 'write_all_flag',
            'default': 0,
            'comment': "whether to include the write all keyword",
            'value': 0  # true/false
        }
        # end dictionary

        self.set_receptor(receptor_filename)  # also sets self.receptor_stem
        self.set_ligand(ligand_filename)
        self.boolean_param_list = [
            'analysis',
            'epdb_flag',
            'flexres_flag',
            'compute_unbound_extended_flag',
            'include_1_4_interactions_flag',
            'intelec',
            'linear_schedule',
            'reorient_flag',
            'rmsatoms_flag',
            'set_ga',
            'set_sw1',
            'set_psw1',
            'simanneal',
            'unbound_flag',
            'unbound_energy_flag',
            'unbound_model_flag',
            'unbound_intnbp_coeffs_flag',
            'write_all',
            'write_all_flag'
        ]
        # end __init__

    def set_version(self, version):
        if float(version) in [3.05, 4.0, 4.1, 4.2]:
            self['autodock_parameter_version']['value'] = float(version)
        else:
            print(version, " is not valid. Valid autodock versions are 3.05, 4.1, 4.2]")

    def set_ligand(self, ligand_filename):
        self.ligand_filename = os.path.basename(ligand_filename)
        basename = os.path.basename(ligand_filename)
        self['rmsref']['value'] = basename
        self['move']['value'] = basename  # @@RH

    def set_receptor(self, receptor_filename):
        self.receptor_filename = os.path.basename(receptor_filename)
        # self.receptor_filename = receptor_filename
        basename = os.path.basename(receptor_filename)
        self.receptor_stem = basename[:basename.rfind('.')]
        if receptor_filename != '':
            self['fld']['value'] = self.receptor_stem + '.maps.fld'

    def set_ligand_types_from_filename(self, ligand_filename):
        ligand = Read(ligand_filename)[0]
        if ligand is None:
            print("unable to read ", ligand_filename)
            return "ERROR"
        d = {a.autodock_element: 1 for a in ligand.allAtoms}
        keys = list(d.keys())
        type_str = keys[0]
        for t in keys[1:]:
            type_str = type_str + ' ' + t + ' '
        # eg "C A NA N OA SA HD"
        self['ligand_types']['value'] = type_str
        self['torsdof4']['value'][0] = ligand.TORSDOF

    def set_ligand_types_from_filename_v3(self, ligand_filename):
        ligand = Read(ligand_filename)
        if ligand is None:
            print("unable to read ", ligand_filename)
            return "ERROR"
        d = {a.autodock_element: 1 for a in ligand.allAtoms}
        keys = list(d.keys())
        keys.sort()
        type_str = keys[0]
        for t in keys[1:]:
            type_str = type_str + t
        # eg "CANOSH"
        self['types']['value'] = type_str
        self['torsdof']['value'][0] = ligand.TORSDOF

    def set_ligand_types3_from_ligand_types(self, ligand_types):
        d = {}
        for t in ligand_types:
            if len(t) == 1:
                d[t] = 1
            elif t[1] in ['A', 'D']:  # AD4 special cases: NA,SA,OA,HD
                d[t[0]] = 1
            elif t in ['Cl', 'CL', 'cl']:  # AD3 special case: chlorine
                d['c'] = 1
            elif t in ['Br', 'BR', 'br']:  # AD3 special case: bromine
                d['b'] = 1
            elif t in ['Fe', 'FE', 'fe']:  # AD3 special case: iron
                d['f'] = 1
            else:
                print("unrecognized ligand_atom_type:", t)
        all_types = list(d.keys())
        all_types.sort()
        type_str = all_types[0]
        for t in all_types[1:]:
            type_str = type_str + t
        self['types']['value'] = type_str

    #
    # read methods
    #
    def read(self, filename):
        """Read lines from the file and call _parse to set current state.
        """
        self.dpf_filename = os.path.basename(filename)
        self.dpf_written_filename = filename
        with open(filename) as dpf_ptr:
            lines = dpf_ptr.readlines()
        self._parse(lines)

    def _parse(self, lines):
        """set the current state according to lines.
        """
        self.file_params = []
        keys = list(self.keys())
        found_compute_unbound_extended = False
        found_include_1_4_interactions = False
        found_reorient = False
        for line in lines:
            p = ''
            words = line.replace('\t', ' ').split()
            if words != [] and words[0][0] != '#':
                p = words[0]
                if p not in keys:
                    print("WARNING: unrecognized parameter:\n", p)
                    continue
                # maintain a list of the parameters read from the file
                if self.file_params == [] or p != self.file_params[-1]:
                    self.file_params.append(p)
                # parse the line
                l = len(words)
                for i in range(l):
                    if words[i][0] == '#':
                        l = i
                values = words[1:l]
                if p == 'include_1_4_interactions':
                    self['include_1_4_interactions']['value'] = self._get_val(values[0])
                    self['include_1_4_interactions_flag']['value'] = True
                    found_include_1_4_interactions = True
                elif p == 'set_sw1':
                    self[p]['value'] = 1
                    self['set_sw1_flag']['value'] = 1
                    self['set_psw1']['value'] = 0
                    self['set_psw1_flag']['value'] = 0
                elif p == 'set_psw1':
                    self[p]['value'] = 1
                    self['set_psw1_flag']['value'] = 1
                    self['set_sw1_flag']['value'] = 0
                    self['set_sw1']['value'] = 0
                elif p == 'write_all':
                    self['write_all_flag']['value'] = True
                elif p in self.boolean_param_list:
                    self[p]['value'] = 1
                elif p == 'ligand_types':
                    self['ligand_types']['value'] = ''.join(words[1:l])
                elif p == 'types':
                    self['types']['value'] = ''.join(words[1:l])
                    self['autodock_parameter_version']['value'] = 3.05
                elif p == 'rmsref':
                    self['rmsref']['value'] = self._get_val(values[0])
                    self['rmsref_flag']['value'] = True
                elif p == 'reorient':
                    self['reorient']['value'] = self._get_val(values[0])
                    self['reorient_flag']['value'] = True
                    found_reorient = True
                elif p in ['epdb', 'epdb_flag']:
                    # if len(values)>0:
                    #    self['epdb']['value'] = self._get_val(values[0])
                    self['epdb_flag']['value'] = True
                elif p == 'compute_unbound_extended':
                    self['compute_unbound_extended_flag']['value'] = True
                    found_compute_unbound_extended = True
                elif p == 'ga_crossover_mode':
                    self['ga_crossover_mode']['value'] = self._get_val(values[0])
                    self['ga_crossover_mode_flag']['value'] = True
                elif p == 'unbound':
                    self['unbound_flag']['value'] = True
                    self['unbound']['value'] = self._get_val(values[0])
                elif p == 'unbound_intnbp_coeffs':
                    self['unbound_intnbp_coeffs_flag']['value'] = True
                    self['unbound_intnbp_coeffs']['value'] = ''.join(words[1:l])
                elif p == 'rmsatoms' and values[0] == 'all':
                    self['rmsatoms_flag']['value'] = True
                    self['rmsatoms']['value'] = values[0]  # this should be 'ligand_only' or'all'
                elif p == 'flexres':
                    self['flexres']['value'] = self._get_val(values[0])
                    self['flexres_flag']['value'] = True
                elif p == 'parameter_file':
                    self['custom_parameter_file']['value'] = 1
                    self['parameter_file']['value'] = values[0]
                elif ((len(values) == 1) and
                      (type(self[p]['default']) != list)):
                    self[p]['value'] = self._get_val(values[0])
                elif p in self:
                    self[p]['value'] = []
                    for v in values:
                        self[p]['value'].append(self._get_val(v))
                    if p == 'quat0' or p == 'quaternion0':  # @@RH
                        self['quaternion0']['value'] = self[p]['value']
                        # self['axisangle0']['value'] = self[p]['value']
                else:
                    print('WARNING: unknown keyword=', p)
            if p == 'fld':
                # so words[1] ends .maps.fld
                ind = words[1].index('maps') - 1
                self.receptor_stem = words[1][:ind]
            elif p == 'map':
                # problem in this case is that elements may be 1 or 2 char
                # so name could end .C.map or .CL.map
                name = words[1]
                if name[-6] == '.':
                    self.receptor_stem = name[:-6]
                elif name[-7] == '.':
                    self.receptor_stem = name[:-7]
        if not found_compute_unbound_extended:
            self['compute_unbound_extended_flag']['value'] = 0
        if not found_include_1_4_interactions:
            self['include_1_4_interactions_flag']['value'] = 0
        if not found_reorient:
            self['reorient_flag']['value'] = 0

    def _get_val(self, val_str):
        try:
            return int(val_str)
        except ValueError:
            pass
        try:
            return float(val_str)
        except ValueError:
            pass
        if type(val_str) == bytes:
            return val_str
        else:
            raise NotImplementedError("value: %s of unsupport type %s" % (val_str, type(val_str).__name__))

    #
    # write methods
    #
    def write(self, filename, param_list):
        """Write the current state to a file

        file is a pointer to a writeable file
        param_list is a list of parameter strings.
        For best results use the parameter_lists supplied by this class.
        """
        if filename == '':
            dpf_ptr = sys.stdout
            self.dpf_filename = ''
            self.dpf_written_filename = None
        else:
            dpf_ptr = open(filename, 'w')
            self.dpf_filename = os.path.basename(filename)
            self.dpf_written_filename = filename
            self.file_params = param_list

        for p in param_list:
            # maps are a special case
            if p == 'map':
                for a in self['types']['value']:
                    dpf_ptr.write(self.make_map_string(p, a))
                # write the electrostatics map; kluge comment
                tmp = self[p]['comment']
                self[p]['comment'] = "electrostatics map"
                dpf_ptr.write(self.make_map_string(p, 'e'))
                self[p]['comment'] = tmp
            elif p == 'fmap':
                dpf_ptr.write(self.make_map_string(p, 'f'))
            # intnbp_r_eps is a special case
            elif p == 'intnbp_r_eps':
                atoms = self['types']['value']
                # strings don't have index method in python1.5.2
                # for a1 in atoms:
                # for a2 in atoms[atoms.index(a1):]:
                # dpf_ptr.write( self.make_intnbp_r_eps_string(a1, a2))
                lenTypes = len(atoms)
                for i in range(lenTypes):
                    for j in range(i, lenTypes):
                        a1 = atoms[i]
                        a2 = atoms[j]
                        dpf_ptr.write(self.make_intnbp_r_eps_string(a1, a2))

            elif p == 'intelec' and self[p]['value']:
                dpf_ptr.write('intelec 0.1146                       # calculate internal electrostatics\n')
            elif p == 'rmsref_flag':
                if self['rmsref_flag']['value']:
                    dpf_ptr.write(self.make_param_string('rmsref'))
            elif p == 'rmsref':
                pass
            elif p == 'set_sw1':
                # print "in write with set_sw1:", self[p]['value']
                if int(self[p]['value']) == 0:
                    self['set_psw1']['value'] = 1
                    dpf_ptr.write(self.make_param_string('set_psw1'))
                else:
                    self[p]['value'] = 1
                    dpf_ptr.write(self.make_param_string(p))
            # all the other parameters handle themselves
            elif p == 'set_psw1':
                # print "in write with set_psw1:", self[p]['value']
                if int(self[p]['value']) == 0:
                    self['set_sw1']['value'] = 1
                    dpf_ptr.write(self.make_param_string('set_sw1'))
                else:
                    self[p]['value'] = 1
                    dpf_ptr.write(self.make_param_string(p))
            elif p == 'about':  # round to 3 decimal places
                about_v = self['about']['value']
                dpf_ptr.write(
                    "about %.3f %.3f %.3f            # small molecule center\n" % (about_v[0], about_v[1], about_v[2]))
            # all the other parameters handle themselves
            else:
                dpf_ptr.write(self.make_param_string(p))
        if dpf_ptr != sys.stdout:
            dpf_ptr.close()

    def make_param_string(self, param):
        """return the output string for the given param using the value
           and comment entries in its dictionary.
        """
        p = self[param]
        vt = type(p['value'])
        if param in self.boolean_param_list:
            if not p['value']:
                return "#\n"
            else:
                val_str = ""
        elif vt in [int, float, str]:
            val_str = str(p['value'])
        elif vt in [list, tuple]:
            val_str = ""
            for v in p['value']:
                val_str = val_str + str(v) + " "
        else:
            raise NotImplementedError("type (%s) of parameter %s unsupported" % (vt.__name__, param))
        return self._make_string(p, val_str)

    def make_intnbp_r_eps_string(self, atom1, atom2):
        p = self['intnbp_r_eps']
        index = "lj" + atom1 + atom2

        val_str = "%5.2f %9.7f 12 6" % (Rij[index], epsij[index])
        p['comment'] = "%s-%s lj" % (atom1, atom2)
        return self._make_string(p, val_str)

    def make_map_string(self, param, type):
        p = self[param]
        val_str = self.receptor_stem + ".%s.map" % (type)
        return self._make_string(p, val_str)

    def _make_string(self, p, val_str):
        # fix 1/2013 for bug report:
        # map bbbb_B99990001_mod_rigid.maps.fld# grid_data file
        # fix 4/2014 for bug report:
        #
        return "%s %s%s # %s\n" % (p['keyword'],
                                   val_str,
                                   " " * (35 - (len(p['keyword']) + len(val_str))),
                                   p['comment'])

    def write4(self, filename, param_list):
        """Write the current state to an AutoDock4 dpf file

        file is a pointer to a writeable file
        param_list is a list of parameter strings.
        For best results use the parameter_lists supplied by this class.
        """
        if filename == '':
            dpf_ptr = sys.stdout
            self.dpf_filename = ''
            self.dpf_written_filename = None
        else:
            dpf_ptr = open(filename, 'w')
            self.dpf_filename = os.path.basename(filename)
            self.dpf_written_filename = filename
            self.file_params = param_list

        # to write dpf4, set unbound_model to extended OR unbound_model extended float
        for p in param_list:
            if p == 'autodock_parameter_version':
                oldval = self['autodock_parameter_version']['value']
                self['autodock_parameter_version']['value'] = 4.2
                dpf_ptr.write(self.make_param_string('autodock_parameter_version'))
                # NEW 4/1/2009
                # also set unbound_model to extended
                self['unbound_model']['value'] = "extended"
                self['unbound_model_flag']['value'] = 1
                self['autodock_parameter_version']['value'] = oldval
            elif p == 'custom_parameter_file':
                if self['custom_parameter_file']['value']:
                    dpf_ptr.write(self.make_param_string('parameter_file'))
            elif p == 'map':
                # maps are a special case
                for a in self['ligand_types']['value'].split():
                    dpf_ptr.write(self.make_map_string(p, a))
                # write the electrostatics map
                dpf_ptr.write(self.make_map_string('elecmap', 'e'))
                # write the desolvation map
                dpf_ptr.write(self.make_map_string('desolvmap', 'd'))
            elif p == 'reorient_flag':
                if self['reorient_flag']['value']:
                    dpf_ptr.write(self.make_param_string('reorient'))
            elif p == 'reorient':
                pass
            elif p == 'set_psw1' or p == 'set_sw1':
                if self['set_psw1_flag']['value']:
                    dpf_ptr.write(self.make_param_string(p))
                elif self['set_sw1_flag']['value']:
                    dpf_ptr.write(self.make_param_string('set_sw1'))
                else:
                    pass
            elif p == 'fmap':
                self[p]['comment'] = "floating point map"
                dpf_ptr.write(self.make_map_string(p, 'f'))
            elif p == 'include_1_4_interactions_flag':
                if self['include_1_4_interactions_flag']['value']:
                    dpf_ptr.write(self.make_param_string('include_1_4_interactions'))
            elif p == 'include_1_4_interactions':
                pass
            # not support 4/1/2009->
            # elif p=='compute_unbound_extended_flag':
            #    if self['compute_unbound_extended_flag']['value']:
            #        dpf_ptr.write(self.make_param_string('compute_unbound_extended'))
            # elif p=='compute_unbound_extended':
            #    pass
            elif p == 'unbound_model_flag':
                self['unbound_model']['value'] = "extended"
                dpf_ptr.write(self.make_param_string('unbound_model'))
            elif p == 'unbound_model':
                pass
            # elif p=='unbound_model_flag':
            #    if self['unbound_model_flag']['value']:
            #        dpf_ptr.write(self.make_param_string('unbound model'))
            # elif p=='unbound_model':
            #    pass
            elif p == 'unbound_flag':
                if self['unbound_flag']['value']:
                    dpf_ptr.write(self.make_param_string('unbound'))
            elif p == 'unbound':
                pass
            elif p == 'unbound_intnbp_coeffs_flag':
                if self['unbound_intnbp_coeffs_flag']['value']:
                    dpf_ptr.write(self.make_param_string('unbound_intnbp_coeffs'))
            elif p == 'ga_crossover_mode_flag':
                if self['ga_crossover_mode_flag']['value']:
                    dpf_ptr.write(self.make_param_string('ga_crossover_mode'))
            elif p == 'ga_crossover_mode':
                pass
            elif p == 'unbound_intnbp_coeffs':
                pass
            elif p == 'rmsatoms_flag':
                if self['rmsatoms_flag']['value'] and self['rmsatoms']['value'] == 'all':
                    dpf_ptr.write(self.make_param_string('rmsatoms'))
            elif p == 'rmsatoms':
                pass
            elif p == 'flexres_flag':
                if self['flexres_flag']['value']:
                    dpf_ptr.write(
                        "flexres %s                  # file containing flexible residues\n" % self['flexres']['value'])
            elif p == 'flexres':
                pass
            elif p == 'unbound_model_flag':
                if self['unbound_model_flag']['value'] > 0:
                    dpf_ptr.write(self.make_param_string('unbound_model'))
            elif p == 'write_all_flag':
                if self['write_all_flag']['value']:
                    dpf_ptr.write("write_all                  # write all conformations in a cluster\n")
            elif p == 'write_all':
                pass
            elif p == 'epdb_flag':
                if self['epdb_flag']['value']:
                    dpf_ptr.write("epdb                     # evaluate ligand specified with move command\n")
            elif p == 'epdb':
                pass
            elif p == 'rmsref_flag':
                flag = self['rmsref_flag']['value']
                if type(flag) == type(''):
                    flag = eval(flag)
                    self['rmsref_flag']['value'] = flag
                if self['rmsref_flag']['value']:
                    if self['rmsref']['value'] == "":
                        dpf_ptr.write(
                            "rmsref %s                  # reference ligand conformation\n" % self.ligand_filename)
                    else:
                        dpf_ptr.write(
                            "rmsref %s                  # reference ligand conformation\n" % self['rmsref']['value'])
            elif p == 'rmsref':
                pass
            elif p == 'torsdof4':
                dpf_ptr.write('torsdof %d                            # torsional degrees of freedom\n' \
                              % (self['torsdof4']['value'][0]))
            elif p == 'intelec':  # always include internal electrostatics
                dpf_ptr.write('intelec                              # calculate internal electrostatics\n')
            elif p == 'about':  # round to 3 decimal places
                about_v = self['about']['value']
                dpf_ptr.write(
                    "about %.3f %.3f %.3f            # small molecule center\n" % (about_v[0], about_v[1], about_v[2]))
            # all the other parameters handle themselves
            else:
                dpf_ptr.write(self.make_param_string(p))
        if dpf_ptr != sys.stdout:
            dpf_ptr.close()

    def write41(self, filename, param_list):
        """Write the current state to an AutoDock4.1 dpf file
        file is a pointer to a writeable file
        param_list is a list of parameter strings.
        For best results use the parameter_lists supplied by this class.
        """
        if filename == '':
            dpf_ptr = sys.stdout
            self.dpf_filename = ''
            self.dpf_written_filename = None
        else:
            dpf_ptr = open(filename, 'w')
            self.dpf_filename = os.path.basename(filename)
            self.dpf_written_filename = filename
            self.file_params = param_list

        for p in param_list:
            if p == 'custom_parameter_file':
                self['custom_parameter_file']['value'] = 1
                # self['parameter_file']['value'] = 'AD4.1_bound.dat'
                if 'parameter_file' not in param_list:
                    oldval = self['parameter_file']['value']
                    self['parameter_file']['value'] = 'AD4.1_bound.dat'
                    dpf_ptr.write(self.make_param_string('parameter_file'))
                    self['parameter_file']['value'] = oldval
            elif p == 'map':
                # maps are a special case
                for a in self['ligand_types']['value'].split():
                    dpf_ptr.write(self.make_map_string(p, a))
                # write the electrostatics map
                dpf_ptr.write(self.make_map_string('elecmap', 'e'))
                # write the desolvation map
                dpf_ptr.write(self.make_map_string('desolvmap', 'd'))
            elif p == 'reorient_flag':
                if self['reorient_flag']['value']:
                    dpf_ptr.write(self.make_param_string('reorient'))
            elif p == 'reorient':
                pass
            elif p == 'set_psw1' or p == 'set_sw1':
                if self['set_psw1_flag']['value']:
                    dpf_ptr.write(self.make_param_string(p))
                elif self['set_sw1_flag']['value']:
                    dpf_ptr.write(self.make_param_string('set_sw1'))
                else:
                    pass
            elif p == 'fmap':
                self[p]['comment'] = "floating point map"
                dpf_ptr.write(self.make_map_string(p, 'f'))
            elif p == 'include_1_4_interactions_flag':
                if self['include_1_4_interactions_flag']['value']:
                    dpf_ptr.write(self.make_param_string('include_1_4_interactions'))
            elif p == 'include_1_4_interactions':
                pass
            elif p == 'unbound_energy_flag':
                # IF user specifies an unbound_energy, write it
                if self['unbound_energy_flag']['value']:
                    dpf_ptr.write(self.make_param_string('unbound_energy'))
            elif p == 'unbound_energy':
                pass
            elif p == 'unbound_model_flag':
                self['unbound_model']['value'] = "bound"
                dpf_ptr.write(self.make_param_string('unbound_model'))
            elif p == 'unbound_model':
                pass
            elif p == 'compute_unbound_extended_flag':
                self['compute_unbound_extended_flag']['value'] = 0
            elif p == 'compute_unbound_extended':
                pass
            elif p == 'unbound_flag':
                self['unbound_flag']['value'] = 0
            elif p == 'unbound':
                pass
            elif p == 'unbound_intnbp_coeffs_flag':
                self['unbound_intnbp_coeffs_flag']['value'] = 0
            elif p == 'ga_crossover_mode_flag':
                if self['ga_crossover_mode_flag']['value']:
                    dpf_ptr.write(self.make_param_string('ga_crossover_mode'))
            elif p == 'ga_crossover_mode':
                pass
            elif p == 'unbound_intnbp_coeffs':
                pass
            elif p == 'rmsatoms_flag':
                if self['rmsatoms_flag']['value'] and self['rmsatoms']['value'] == 'all':
                    dpf_ptr.write(self.make_param_string('rmsatoms'))
            elif p == 'rmsatoms':
                pass
            elif p == 'flexres_flag':
                if self['flexres_flag']['value']:
                    dpf_ptr.write(
                        "flexres %s                  # file containing flexible residues\n" % self['flexres']['value'])
            elif p == 'flexres':
                pass
            elif p == 'unbound_model_flag':
                if self['unbound_model_flag']['value'] > 0:
                    dpf_ptr.write(self.make_param_string('unbound_model'))
            elif p == 'unbound_model':
                pass
            elif p == 'write_all_flag':
                if self['write_all_flag']['value']:
                    dpf_ptr.write("write_all                  # write all conformations in a cluster\n")
            elif p == 'write_all':
                pass
            elif p == 'epdb_flag':
                if self['epdb_flag']['value']:
                    dpf_ptr.write("epdb                     # evaluate ligand specified with move command\n")
            elif p == 'epdb':
                pass
            elif p == 'rmsref_flag':
                flag = self['rmsref_flag']['value']
                if type(flag) == type(''):
                    flag = eval(flag)
                    self['rmsref_flag']['value'] = flag
                if self['rmsref_flag']['value']:
                    if self['rmsref']['value'] == "":
                        dpf_ptr.write(
                            "rmsref %s                  # reference ligand conformation\n" % self.ligand_filename)
                    else:
                        dpf_ptr.write(
                            "rmsref %s                  # reference ligand conformation\n" % self['rmsref']['value'])
            elif p == 'rmsref':
                pass
            elif p == 'torsdof4':
                dpf_ptr.write('torsdof %d                            # torsional degrees of freedom\n' % (
                    self['torsdof4']['value'][0]))
            elif p == 'intelec':  # always include internal electrostatics
                dpf_ptr.write('intelec                              # calculate internal electrostatics\n')
            # all the other parameters handle themselves
            else:
                dpf_ptr.write(self.make_param_string(p))
        if dpf_ptr != sys.stdout:
            dpf_ptr.close()

    def write42(self, filename, param_list):
        """Write the current state to an AutoDock4.2 dpf file
        file is a pointer to a writeable file
        param_list is a list of parameter strings.
        For best results use the parameter_lists supplied by this class.
        """
        local_only = 'do_local_only' in param_list
        if local_only:
            self['tran0']['value'] = self['about']['value']
            self['quaternion0']['value'] = "0 0 0 1"
            # self['dihe0']['value'] = "0 " *self.ligand.ndihe
        if filename == '':
            dpf_ptr = sys.stdout
            self.dpf_filename = ''
            self.dpf_written_filename = None
        else:
            dpf_ptr = open(filename, 'w')
            self.dpf_filename = os.path.basename(filename)
            self.dpf_written_filename = filename
            self.file_params = param_list

        dpf_ptr.write(self.make_param_string('autodock_parameter_version'))
        for p in param_list[1:]:
            # both custom_parameter_file + parameter_file must be set
            # if parameter_file is missing, defaults to 'AD4.1_bound.dat'
            if self['custom_parameter_file']['value'] == 1:
                # write parameter_file only once HERE
                self['custom_parameter_file']['value'] = 0
                dpf_ptr.write(self.make_param_string('parameter_file'))
            elif p == 'custom_parameter_file':  # handled above
                self['custom_parameter_file']['value'] = 1
            elif p == 'map':
                # maps are a special case
                for a in self['ligand_types']['value'].split():
                    dpf_ptr.write(self.make_map_string(p, a))
                # write the electrostatics map
                dpf_ptr.write(self.make_map_string('elecmap', 'e'))
                # write the desolvation map
                dpf_ptr.write(self.make_map_string('desolvmap', 'd'))
            elif p == 'reorient_flag':
                if self['reorient_flag']['value']:
                    dpf_ptr.write(self.make_param_string('reorient'))
            elif p == 'reorient':
                pass
            elif p == 'set_psw1' or p == 'set_sw1':
                if self['set_psw1_flag']['value']:
                    dpf_ptr.write(self.make_param_string(p))
                elif self['set_sw1_flag']['value']:
                    dpf_ptr.write(self.make_param_string('set_sw1'))
                else:
                    pass
            elif p == 'fmap':
                self[p]['comment'] = "floating point map"
                dpf_ptr.write(self.make_map_string(p, 'f'))
            elif p == 'include_1_4_interactions_flag':
                if self['include_1_4_interactions_flag']['value']:
                    dpf_ptr.write(self.make_param_string('include_1_4_interactions'))
            elif p == 'include_1_4_interactions':
                pass
            elif p == 'unbound_energy_flag':
                # IF user specifies an unbound_energy, write it
                if self['unbound_energy_flag']['value']:
                    dpf_ptr.write(self.make_param_string('unbound_energy'))
            elif p == 'unbound_energy':
                pass
            elif p == 'unbound_model_flag':
                self['unbound_model']['value'] = "bound"
                dpf_ptr.write(self.make_param_string('unbound_model'))
            elif p == 'unbound_model':
                pass
            elif p == 'compute_unbound_extended_flag':
                self['compute_unbound_extended_flag']['value'] = 0
            elif p == 'compute_unbound_extended':
                pass
            elif p == 'unbound_flag':
                self['unbound_flag']['value'] = 0
            elif p == 'unbound':
                pass
            elif p == 'unbound_intnbp_coeffs_flag':
                self['unbound_intnbp_coeffs_flag']['value'] = 0
            elif p == 'ga_crossover_mode_flag':
                if self['ga_crossover_mode_flag']['value']:
                    dpf_ptr.write(self.make_param_string('ga_crossover_mode'))
            elif p == 'ga_crossover_mode':
                pass
            elif p == 'unbound_intnbp_coeffs':
                pass
            elif p == 'rmsatoms_flag':
                if self['rmsatoms_flag']['value'] and self['rmsatoms']['value'] == 'all':
                    dpf_ptr.write(self.make_param_string('rmsatoms'))
            elif p == 'rmsatoms':
                pass
            elif p == 'flexres_flag':
                if self['flexres_flag']['value']:
                    dpf_ptr.write(
                        "flexres %s                  # file containing flexible residues\n" % self['flexres']['value'])
            elif p == 'flexres':
                pass
            elif p == 'unbound_model_flag':
                if self['unbound_model_flag']['value'] > 0:
                    dpf_ptr.write(self.make_param_string('unbound_model'))
            elif p == 'unbound_model':
                pass
            elif p == 'write_all_flag':
                if self['write_all_flag']['value']:
                    dpf_ptr.write("write_all                  # write all conformations in a cluster\n")
            elif p == 'write_all':
                pass
            elif p == 'epdb_flag':
                if self['epdb_flag']['value']:
                    dpf_ptr.write("epdb                     # evaluate ligand specified with move command\n")
            elif p == 'epdb':
                pass
            elif p == 'rmsref_flag':
                flag = self['rmsref_flag']['value']
                if type(flag) == type(''):
                    flag = eval(flag)
                    self['rmsref_flag']['value'] = flag
                if self['rmsref_flag']['value']:
                    if self['rmsref']['value'] == "":
                        dpf_ptr.write(
                            "rmsref %s                  # reference ligand conformation\n" % self.ligand_filename)
                    else:
                        dpf_ptr.write(
                            "rmsref %s                  # reference ligand conformation\n" % self['rmsref']['value'])
            elif p == 'rmsref':
                pass
            elif p == 'torsdof4':
                dpf_ptr.write('torsdof %d                            # torsional degrees of freedom\n' % (
                    self['torsdof4']['value'][0]))
            elif p == 'intelec':  # always include internal electrostatics
                dpf_ptr.write('intelec                              # calculate internal electrostatics\n')
            # all the other parameters handle themselves
            elif p == 'about':  # round to 3 decimal places
                about_v = self['about']['value']
                dpf_ptr.write(
                    "about %.3f %.3f %.3f            # small molecule center\n" % (about_v[0], about_v[1], about_v[2]))
            elif p == 'ligand_types':  # make sure there is at least one space before pound sign
                # new_val = 'ligand_types '
                ligTypes_val = self['ligand_types']['value']
                marker = ligTypes_val.find("#")
                if marker == -1:
                    dpf_ptr.write(self.make_param_string(p))
                else:
                    ss1 = ligTypes[:marker]
                    self['ligand_types']['value'] = ligTypes[:marker]
                    dpf_ptr.write(self.make_param_string(p))
            else:
                dpf_ptr.write(self.make_param_string(p))
        if dpf_ptr != sys.stdout:
            dpf_ptr.close()


# the following lists are class variables describing the keywords
# that must be output for a file of the specified type.

cluster_list = ['types', 'rmstol', 'cluster', 'analysis']  # 3.05

genetic_algorithm_list = [  # 3.05
    'outlev',
    'seed',
    'types',
    'fld',
    'map',
    'move',
    'about',
    'tran0',
    'quat0',
    'ndihe',
    'dihe0',
    'torsdof',
    'intnbp_r_eps',
    'intelec',
    'rmstol',
    'rmsref',
    'rmsref_flag',
    'extnrg',
    'e0max',
    'ga_pop_size',
    'ga_num_evals',
    'ga_num_generations',
    'ga_elitism',
    'ga_mutation_rate',
    'ga_crossover_rate',
    'ga_window_size',
    'ga_cauchy_alpha',
    'ga_cauchy_beta',
    'set_ga',
    'do_global_only',  # ga_only_run
    'analysis']

genetic_algorithm_local_search_list = [  # 3.05
    'outlev',
    'seed',
    'types',
    'fld',
    'map',
    'move',
    'about',
    'tran0',
    'quat0',
    'ndihe',
    'dihe0',
    'torsdof',
    'intnbp_r_eps',
    'intelec',
    'rmstol',
    'rmsref',
    'rmsref_flag',
    'extnrg',
    'e0max',
    'ga_pop_size',
    'ga_num_evals',
    'ga_num_generations',
    'ga_elitism',
    'ga_mutation_rate',
    'ga_crossover_rate',
    'ga_window_size',
    'ga_cauchy_alpha',
    'ga_cauchy_beta',
    'set_ga',
    'sw_max_its',
    'sw_max_succ',
    'sw_max_fail',
    'sw_rho',
    'sw_lb_rho',
    'ls_search_freq',
    'set_sw1',
    'ga_run',
    'analysis'
]

local_search_list = [  # 3.05
    'outlev',
    'seed',
    'types',
    'fld',
    'map',
    'move',
    'about',
    'tran0',
    'quat0',
    'ndihe',
    'dihe0',
    'torsdof',
    'intnbp_r_eps',
    'intelec',
    'rmstol',
    'rmsref',
    'rmsref_flag',
    'extnrg',
    'e0max',
    'sw_max_its',
    'sw_max_succ',
    'sw_max_fail',
    'sw_rho',
    'sw_lb_rho',
    'ls_search_freq',
    'set_sw1',
    'do_local_only',
    'analysis'
]

simulated_annealing_list = [  # 3.05
    'outlev',
    'seed',
    'types',
    'fld',
    'map',
    'move',
    'about',
    'tran0',
    'quat0',
    'ndihe',
    'dihe0',
    'tstep',
    'qstep',
    'dstep',
    'torsdof',
    'intnbp_r_eps',
    'intelec',
    'rmstol',
    'rmsref',
    'rmsref_flag',
    'extnrg',
    'e0max',
    'rt0',
    'trnrf',
    'quarf',
    'dihrf',
    'runs',
    'cycles',
    'accs',
    'rejs',
    'select',
    'linear_schedule',
    'simanneal',
    'analysis'
]

docking_parameter_list = [  # 3.05
    'about',
    'accs',
    'analysis',
    'cluster',
    'cycles',
    'dihe0',
    'dihrf',
    'do_global_only',
    'do_local_only',
    'dstep',
    'e0max',
    'extnrg',
    'fld',
    'fmap',
    'ga_cauchy_alpha',
    'ga_cauchy_beta',
    'ga_crossover_rate',
    'ga_elitism',
    'ga_mutation_rate',
    'ga_num_evals',
    'ga_num_generations',
    'ga_pop_size',
    'ga_run',
    'ga_window_size',
    'intnbp_r_eps',
    'intelec',
    'linear_schedule',
    'ls_search_freq',
    'map',
    'move',
    'ndihe',
    'outlev',
    'qstep',
    'quarf',
    'quat0',
    'rejs',
    'rmsref',
    'rmstol',
    'rmsref',
    'rmsref_flag',
    'rt0',
    'rtrf',
    'runs',
    'seed',
    'select',
    'set_ga',
    'set_sw1',
    'set_psw1',
    'simanneal',
    'sw_lb_rho',
    'sw_max_fail',
    'sw_max_its',
    'sw_max_succ',
    'sw_rho',
    'torsdof',
    'tran0',
    'trnrf',
    'tstep',
    'types',
    'write_all',
    'write_all_flag']

#####AUTODOCK4####

###???parameter_file???
cluster_list4 = ['autodock_parameter_version', 'custom_parameter_file', 'ligand_types', 'rmstol', 'cluster', 'analysis']

genetic_algorithm_list4 = [
    'autodock_parameter_version',
    'outlev',
    'custom_parameter_file',  # NEW
    'include_1_4_interactions',  # NEW
    'include_1_4_interactions_flag',  # NEW
    'intelec',
    'seed',
    'unbound_intnbp_coeffs',
    'unbound_intnbp_coeffs_flag',
    'ligand_types',  # NEW
    'fld',
    'map',
    'move',
    'flexres_flag',
    'flexres',
    'about',
    'reorient_flag',
    'tran0',
    # 'axisangle0',
    'quaternion0',
    'dihe0',
    'torsdof4',
    'rmstol',
    'rmsref',
    'rmsref_flag',
    'rmsatoms',
    'rmsatoms_flag',
    'extnrg',
    'e0max',
    'ga_pop_size',
    'ga_num_evals',
    'ga_num_generations',
    'ga_elitism',
    'ga_mutation_rate',
    'ga_crossover_rate',
    'ga_crossover_mode_flag',
    'ga_window_size',
    'ga_cauchy_alpha',
    'ga_cauchy_beta',
    'set_ga',
    'unbound',  # set to user value
    'unbound_flag',  # 0 by default
    'unbound_model',  # set to "extended"
    'unbound_model_flag',  # 1 by default
    'epdb',  # NEW ???
    'epdb_flag',  # NEW ???
    # 'compute_unbound_extended',
    # 'compute_unbound_extended_flag',
    'do_global_only',  # why isn't ga_run here???
    'write_all',
    'write_all_flag',
    'analysis']

genetic_algorithm_local_search_list4 = [
    'autodock_parameter_version',
    'outlev',
    'custom_parameter_file',  # NEW
    'include_1_4_interactions',  # NEW
    'include_1_4_interactions_flag',  # NEW
    'intelec',
    'seed',
    'unbound_intnbp_coeffs',
    'unbound_intnbp_coeffs_flag',
    'ligand_types',  # NEW
    'fld',
    'map',
    'move',
    'flexres_flag',
    'flexres',
    'about',
    'reorient_flag',
    'tran0',
    # 'axisangle0',
    'quaternion0',
    'dihe0',
    'torsdof4',
    'rmstol',
    'rmsref',
    'rmsref_flag',
    'rmsatoms',
    'rmsatoms_flag',
    'extnrg',
    'e0max',
    'ga_pop_size',
    'ga_num_evals',
    'ga_num_generations',
    'ga_elitism',
    'ga_mutation_rate',
    'ga_crossover_rate',
    'ga_crossover_mode_flag',
    'ga_window_size',
    'ga_cauchy_alpha',
    'ga_cauchy_beta',
    'set_ga',
    'sw_max_its',
    'sw_max_succ',
    'sw_max_fail',
    'sw_rho',
    'sw_lb_rho',
    'ls_search_freq',
    'set_psw1',
    'unbound',  # set to user value
    'unbound_flag',  # 0 by default
    'unbound_model',  # set to "extended"
    'unbound_model_flag',  # 1 by default
    'epdb',  # NEW ???
    'epdb_flag',  # NEW ???
    # 'compute_unbound_extended',
    # 'compute_unbound_extended_flag',
    'ga_run',
    'write_all',
    'write_all_flag',
    'analysis'
]

genetic_algorithm_local_search_list4_with_parameter_file = [
    'autodock_parameter_version',
    'outlev',
    'parameter_file',  # NEW
    'include_1_4_interactions',  # NEW
    'include_1_4_interactions_flag',  # NEW
    'intelec',
    'seed',
    'unbound_intnbp_coeffs',
    'unbound_intnbp_coeffs_flag',
    'ligand_types',  # NEW
    'fld',
    'map',
    'move',
    'flexres_flag',
    'flexres',
    'about',
    'reorient_flag',
    'tran0',
    # 'axisangle0',
    'quaternion0',
    'dihe0',
    'torsdof4',
    'rmstol',
    'rmsref',
    'rmsref_flag',
    'rmsatoms',
    'rmsatoms_flag',
    'extnrg',
    'e0max',
    'ga_pop_size',
    'ga_num_evals',
    'ga_num_generations',
    'ga_elitism',
    'ga_mutation_rate',
    'ga_crossover_rate',
    'ga_crossover_mode_flag',
    'ga_window_size',
    'ga_cauchy_alpha',
    'ga_cauchy_beta',
    'set_ga',
    'sw_max_its',
    'sw_max_succ',
    'sw_max_fail',
    'sw_rho',
    'sw_lb_rho',
    'ls_search_freq',
    'set_psw1',
    'unbound',  # set to user value
    'unbound_flag',  # 0 by default
    'unbound_model',  # set to "extended"
    'unbound_model_flag',  # 1 by default
    'epdb',  # NEW ???
    'epdb_flag',  # NEW ???
    # 'compute_unbound_extended',
    # 'compute_unbound_extended_flag',
    'ga_run',
    'write_all',
    'write_all_flag',
    'analysis'
]

local_search_list4 = [
    'autodock_parameter_version',
    'outlev',
    'custom_parameter_file',  # NEW
    'include_1_4_interactions',  # NEW
    'include_1_4_interactions_flag',  # NEW
    'intelec',
    'seed',
    'unbound_intnbp_coeffs',
    'unbound_intnbp_coeffs_flag',
    'ligand_types',  # NEW
    'fld',
    'map',
    'move',
    'flexres_flag',
    'flexres',
    'about',
    'reorient_flag',
    'tran0',
    # 'axisangle0',
    'quaternion0',
    'dihe0',
    'torsdof4',
    'rmstol',
    'rmsref',
    'rmsref_flag',
    'rmsatoms',
    'rmsatoms_flag',
    'extnrg',
    'e0max',
    'ga_pop_size',
    'sw_max_its',
    'sw_max_succ',
    'sw_max_fail',
    'sw_rho',
    'sw_lb_rho',
    'ls_search_freq',
    'set_psw1',
    'unbound',  # set to user value
    'unbound_flag',  # 0 by default
    'unbound_model',  # set to "extended"
    'unbound_model_flag',  # 1 by default
    'epdb',  # NEW
    'epdb_flag',  # NEW ???
    'do_local_only',
    'write_all',
    'write_all_flag',
    'analysis'
]

simulated_annealing_list4 = [
    'autodock_parameter_version',
    'outlev',
    'custom_parameter_file',  # NEW
    'include_1_4_interactions',  # NEW
    'include_1_4_interactions_flag',  # NEW
    'intelec',
    'seed',
    'unbound_intnbp_coeffs',
    'unbound_intnbp_coeffs_flag',
    'ligand_types',  # NEW
    'fld',
    'map',
    'move',
    'flexres_flag',
    'flexres',
    'about',
    'reorient_flag',
    'tran0',
    # 'axisangle0',
    'quaternion0',
    'dihe0',
    'tstep',
    'qstep',
    'dstep',
    'torsdof4',
    'rmsatoms',
    'rmsatoms_flag',
    'rmstol',
    'rmsref',
    'rmsref_flag',
    'extnrg',
    'e0max',
    'rt0',
    'trnrf',
    'quarf',
    'dihrf',
    'runs',
    'cycles',
    'accs',
    'rejs',
    'select',
    'linear_schedule',
    'unbound',  # set to user value
    'unbound_flag',  # 0 by default
    'unbound_model',  # set to "extended"
    'unbound_model_flag',  # 1 by default
    'epdb',  # NEW
    'epdb_flag',  # NEW ???
    'simanneal',
    'write_all',
    'write_all_flag',
    'analysis'
]

docking_parameter_list4 = [
    'autodock_parameter_version',
    'about',
    'accs',
    'analysis',
    'cluster',
    'compute_unbound_extended',
    'compute_unbound_extended_flag',  # 0 by default
    'cycles',
    'dihe0',
    'dihrf',
    'do_global_only',
    'do_local_only',
    'dstep',
    'e0max',
    'epdb',  # NEW
    'epdb_flag',  # NEW ???
    'extnrg',
    'fld',
    'fmap',
    'ga_cauchy_alpha',
    'ga_cauchy_beta',
    'ga_crossover_rate',
    'ga_crossover_mode_flag',
    'ga_elitism',
    'ga_mutation_rate',
    'ga_num_evals',
    'ga_num_generations',
    'ga_pop_size',
    'ga_run',
    'ga_window_size',
    'include_1_4_interactions',  # NEW
    'include_1_4_interactions_flag',  # NEW
    'intelec',
    'ligand_types',  # NEW
    'linear_schedule',
    'ls_search_freq',
    'map',
    'move',
    'reorient_flag',
    'flexres_flag',
    'flexres',
    'outlev',
    'parameter_library',  # NEW
    'qstep',
    'quarf',
    # 'axisangle0',
    'quaternion0',
    'rejs',
    'rmsatoms',
    'rmsatoms_flag',
    'rmsref',
    'rmstol',
    'rmsref',
    'rmsref_flag',
    'rt0',
    'rtrf',
    'runs',
    'seed',
    'select',
    'set_ga',
    'set_sw1',
    'set_psw1',
    'simanneal',
    'sw_lb_rho',
    'sw_max_fail',
    'sw_max_its',
    'sw_max_succ',
    'sw_rho',
    'torsdof4',
    'unbound',  # set to user value
    'unbound_flag',  # 0 by default
    'unbound_model',  # set to "extended"
    'unbound_model_flag',  # 1 by default
    'tran0',
    'trnrf',
    'tstep',
    'write_all',
    'write_all_flag']

#####AUTODOCK4_1####
##### and
#####AUTODOCK4_2####


###???parameter_file???
cluster_list4 = ['autodock_parameter_version', 'custom_parameter_file', 'ligand_types', 'rmstol', 'cluster', 'analysis']

genetic_algorithm_list4_1 = [
    'autodock_parameter_version',
    'outlev',
    # 'custom_parameter_file',       #REMOVE WHEN  autodock4.1 is released
    # 'parameter_file',              #REMOVE WHEN autodock4.1 is released
    'include_1_4_interactions',
    'include_1_4_interactions_flag',
    'intelec',
    'seed',
    'unbound_intnbp_coeffs',
    'unbound_intnbp_coeffs_flag',
    'ligand_types',
    'fld',
    'map',
    'move',
    'flexres_flag',
    'flexres',
    'about',
    'reorient_flag',
    'tran0',
    # 'axisangle0',
    'quaternion0',
    'dihe0',
    'torsdof4',
    'rmstol',
    'rmsref',
    'rmsref_flag',
    'rmsatoms',
    'rmsatoms_flag',
    'extnrg',
    'e0max',
    'ga_pop_size',
    'ga_num_evals',
    'ga_num_generations',
    'ga_elitism',
    'ga_mutation_rate',
    'ga_crossover_rate',
    'ga_crossover_mode',
    'ga_crossover_mode_flag',
    'ga_window_size',
    'ga_cauchy_alpha',
    'ga_cauchy_beta',
    'set_ga',
    'unbound',
    'unbound_flag',  # 0 by default
    'unbound_energy',
    'unbound_energy_flag',  # 0 by default
    'unbound_model',  # 'bound' by default
    'unbound_model_flag',
    'epdb',
    'epdb_flag',
    'do_global_only',  # 'do_global_only' sets 'DPF_GS' here
    'write_all',
    'write_all_flag',
    'analysis']

genetic_algorithm_list4_2 = genetic_algorithm_list4_1

genetic_algorithm_local_search_list4_1 = [
    'autodock_parameter_version',
    'outlev',
    # 'custom_parameter_file',   #REMOVE WHEN  autodock4.1 is released
    # 'parameter_file',          #REMOVE WHEN  autodock4.1 is released
    'include_1_4_interactions',
    'include_1_4_interactions_flag',
    'intelec',
    'seed',
    'unbound_intnbp_coeffs',
    'unbound_intnbp_coeffs_flag',
    'ligand_types',
    'fld',
    'map',
    'move',
    'flexres_flag',
    'flexres',
    'about',
    'reorient_flag',
    'tran0',
    'quaternion0',
    'dihe0',
    'torsdof4',
    'rmstol',
    'rmsref',
    'rmsref_flag',
    'rmsatoms',
    'rmsatoms_flag',
    'extnrg',
    'e0max',
    'ga_pop_size',
    'ga_num_evals',
    'ga_num_generations',
    'ga_elitism',
    'ga_mutation_rate',
    'ga_crossover_rate',
    'ga_crossover_mode',
    'ga_crossover_mode_flag',
    'ga_window_size',
    'ga_cauchy_alpha',
    'ga_cauchy_beta',
    'set_ga',
    'sw_max_its',
    'sw_max_succ',
    'sw_max_fail',
    'sw_rho',
    'sw_lb_rho',
    'ls_search_freq',
    'set_psw1',
    'epdb',
    'epdb_flag',
    'unbound',
    'unbound_flag',  # 0 by default
    'unbound_energy',
    'unbound_energy_flag',  # 0 by default
    'unbound_model',  # 'bound' by default
    'unbound_model_flag',
    'ga_run',
    'write_all',
    'write_all_flag',
    'analysis'
]

genetic_algorithm_local_search_list4_2 = genetic_algorithm_local_search_list4_1

genetic_algorithm_local_search_list4_1_with_parameter_file = [
    'autodock_parameter_version',
    'outlev',
    'custom_parameter_file',
    'parameter_file',
    'include_1_4_interactions',
    'include_1_4_interactions_flag',
    'intelec',
    'seed',
    'unbound_intnbp_coeffs',
    'unbound_intnbp_coeffs_flag',
    'ligand_types',
    'fld',
    'map',
    'move',
    'flexres_flag',
    'flexres',
    'about',
    'reorient_flag',
    'tran0',
    'quaternion0',
    'dihe0',
    'torsdof4',
    'rmstol',
    'rmsref',
    'rmsref_flag',
    'rmsatoms',
    'rmsatoms_flag',
    'extnrg',
    'e0max',
    'ga_pop_size',
    'ga_num_evals',
    'ga_num_generations',
    'ga_elitism',
    'ga_mutation_rate',
    'ga_crossover_rate',
    'ga_crossover_mode',
    'ga_crossover_mode_flag',
    'ga_window_size',
    'ga_cauchy_alpha',
    'ga_cauchy_beta',
    'set_ga',
    'sw_max_its',
    'sw_max_succ',
    'sw_max_fail',
    'sw_rho',
    'sw_lb_rho',
    'ls_search_freq',
    'set_psw1',
    'unbound',
    'unbound_flag',  # 0 by default
    'unbound_energy',
    'unbound_energy_flag',  # 0 by default
    'unbound_model',  # 'bound' by default
    'unbound_model_flag',
    'epdb',
    'epdb_flag',
    'ga_run',
    'write_all',
    'write_all_flag',
    'analysis'
]

genetic_algorithm_local_search_list4_2_with_parameter_file = genetic_algorithm_local_search_list4_1_with_parameter_file

local_search_list4_1 = [
    'autodock_parameter_version',
    'outlev',
    # 'custom_parameter_file',   #REMOVE WHEN  autodock4.1 is released
    # 'parameter_file',          #REMOVE WHEN  autodock4.1 is released
    'intelec',
    'ligand_types',
    'fld',
    'map',
    'move',
    'flexres_flag',
    'flexres',
    'about',
    'torsdof4',
    'sw_max_its',
    'sw_max_succ',
    'sw_max_fail',
    'sw_rho',
    'sw_lb_rho',
    'set_psw1',
    'do_local_only',
    'rmstol',
    'analysis'
]

local_search_list4_2 = local_search_list4_1

simulated_annealing_list4_1 = [
    'autodock_parameter_version',
    'outlev',
    'include_1_4_interactions',
    'include_1_4_interactions_flag',
    'intelec',
    'seed',
    'unbound_intnbp_coeffs',
    'unbound_intnbp_coeffs_flag',
    'ligand_types',
    'fld',
    'map',
    'move',
    'flexres_flag',
    'flexres',
    'about',
    'reorient_flag',
    'tran0',
    'quaternion0',
    'dihe0',
    'tstep',
    'qstep',
    'dstep',
    'torsdof4',
    'rmsatoms',
    'rmsatoms_flag',
    'rmstol',
    'rmsref',
    'rmsref_flag',
    'extnrg',
    'e0max',
    'rt0',
    # 'rtrf',
    'trnrf',
    'quarf',
    'dihrf',
    'runs',
    'cycles',
    'accs',
    'rejs',
    'select',
    'linear_schedule',
    'unbound',
    'unbound_flag',  # 0 by default
    'unbound_energy',
    'unbound_energy_flag',  # 0 by default
    'unbound_model',  # 'bound' by default
    'unbound_model_flag',
    'simanneal',
    'write_all',
    'write_all_flag',
    'analysis'
]

simulated_annealing_list4_2 = simulated_annealing_list4_1

docking_parameter_list4_1 = [
    'autodock_parameter_version',
    'about',
    'accs',
    'analysis',
    'cluster',
    'compute_unbound_extended',
    'compute_unbound_extended_flag',
    'cycles',
    'dihe0',
    'dihrf',
    'do_global_only',
    'do_local_only',
    'dstep',
    'e0max',
    'epdb',
    'epdb_flag',
    'extnrg',
    'fld',
    'fmap',
    'ga_cauchy_alpha',
    'ga_cauchy_beta',
    'ga_crossover_rate',
    'ga_crossover_mode',
    'ga_crossover_mode_flag',
    'ga_elitism',
    'ga_mutation_rate',
    'ga_num_evals',
    'ga_num_generations',
    'ga_pop_size',
    'ga_run',
    'ga_window_size',
    'include_1_4_interactions',
    'include_1_4_interactions_flag',
    'intelec',
    'ligand_types',
    'linear_schedule',
    'ls_search_freq',
    'map',
    'move',
    'reorient_flag',
    'flexres_flag',
    'flexres',
    'outlev',
    'parameter_library',
    'qstep',
    'quarf',
    'quaternion0',
    'rejs',
    'rmsatoms',
    'rmsatoms_flag',
    'rmsref',
    'rmstol',
    'rmsref',
    'rmsref_flag',
    'rt0',
    'rtrf',
    'runs',
    'seed',
    'select',
    'set_ga',
    'set_sw1',
    'set_psw1',
    'simanneal',
    'sw_lb_rho',
    'sw_max_fail',
    'sw_max_its',
    'sw_max_succ',
    'sw_rho',
    'torsdof4',
    'tran0',
    'trnrf',
    'tstep',
    'unbound',
    'unbound_flag',
    'unbound_energy',
    'unbound_energy_flag',  # 0 by default
    'unbound_model',
    'unbound_model_flag',
    'write_all',
    'write_all_flag']

docking_parameter_list4_2 = docking_parameter_list4_1

epdb_list4_2 = [
    'autodock_parameter_version',
    'outlev',
    'intelec',
    'ligand_types',
    'fld',
    'map',
    'move',
    'about',
    'epdb',
    'epdb_flag']

implemented = [
    'accs',
    'about',
    'analysis',
    'autodock_parameter_version',  # @@
    'axisangle0',
    'barrier',  #
    'bin_energies_by_rmsd',  #
    'charmap',  #
    'cluster',
    'compute_unbound_extended',
    'compute_unbound_extended_flag',
    'confsampler',  #
    'copyright',  #
    'custom_parameter_file',
    'cycles',
    'desolvmap',
    'dihe0',
    'dihrf',  #
    'do_gals',  # @@
    'do_global_only',
    'do_local_only',
    'do_pso',  # @@
    'dstep',
    'e0max',
    'elecmap',
    'epdb',
    'epdb_flag',
    'extnrg',
    'fld',
    'flex',  # @@
    'flexres',
    'flexres_flag',
    'flexible_residues',
    'ga_boltzman_selection',
    'ga_cauchy_alpha',
    'ga_cauchy_beta',
    'ga_crossover_rate',
    'ga_crossover_mode',
    'ga_crossover_mode_flag',
    'ga_elitism',
    'ga_high',
    'ga_linear_ranking_selection',
    'ga_low',
    'gals_run',
    'ga_mutation_rate',
    'ga_num_evals',
    'ga_num_generations',
    'ga_only_run',
    'ga_pop_size',
    'ga_proportional_selection',
    'ga_run',
    'ga_termination',
    'ga_termination_criterion',
    'ga_tournament_selection',
    'gausstorcon',
    'ga_window_size',
    'hardtorcon',
    'include_1_4_interactions',
    'include_1_4_interactions_flag',
    'intelec',
    'intelec4',  # @@
    'intnbp_coeffs',
    'intnbp_r_eps',
    'investigate',  # @@
    'ligand',  # @@
    'ligand_is_not_inhibitor',  # @@
    'ligand_types',
    'linear_schedule',
    'linsched',
    'ls_run',
    'ls_search_freq',
    'map',
    'fmap',
    'move',
    'ndihe',
    'outlev',
    'output_pop_file',
    'output_population_statistics',
    'output_resnum_as',
    'parameter_file',
    'parameter_library',
    'pso_adaptive_velocity',
    'pso_c1',
    'pso_c2',
    'pso_interpolate_as_scalars',
    'pso_k',
    'pso_neighbors',
    'pso_neighbors_dynamic',
    'pso_neighbors_symmetric',
    'pso_qvmax',
    'pso_random_by_dimension',
    'pso_regenerate_at_limit',
    'pso_rvmax',
    'pso_stage2constriction',
    'pso_tvmax',
    'pso_w_end',
    'pso_w_start',
    'psw_trans_scale',
    'psw_rot_scale',
    'psw_tors_scale',
    'qstep',
    'quarf',
    'quat0',
    'quaternion0',
    'receptor_types',
    'rejs',
    'reorient',
    'reorient_flag',
    'rmsatoms',
    'rmsatoms_flag',
    'rmsmode',
    'rmsnosym',
    'rmsref',
    'rmsref_flag',
    'rmstol',
    'rt0',
    'rtrf',
    'runs',
    'scale_eintermol',
    'schedlin',
    'schedule_linear',
    'trnrf',
    'quarf',
    'dihrf',
    'seed',
    'select',
    'set_ga',
    'set_pattern',
    'set_psw1',
    'set_sw1',
    'set_unbound_energy',
    'showtorpen',
    'simanneal',
    'smooth',
    'sw_lb_rho',
    'sw_max_fail',
    'sw_max_its',
    'sw_max_succ',
    'sw_rho',
    'torsdof',
    'torsdof4',
    'tran0',
    'trjbeg',
    'trjend',
    'trjfrq',
    'trjout',
    'trjsel',
    'trnrf',
    'tstep',
    'types',  # @@
    'unbound',
    'unbound_flag',
    'unbound_energy',
    'unbound_energy_flag',  # 0 by default
    'unbound_intnbp_coeffs',
    'unbound_intnbp_coeffs_flag',
    'unbound_model',
    'unbound_model_flag',
    'autodock_parameter_version',
    'warranty',
    'watch',
    'write_all',
    'write_all_flag',
    'write_all_cluster_members']


class DockingParameterFileMaker:
    """Accept a <ligand>.pdbq and <receptor>.pdbqs and create
    <ligand>_<receptor>.dpf
    """

    def __init__(self, verbose=None, autodock_parameter_version='3.05'):
        self.verbose = verbose
        self.dpo = DockingParameters()
        self.dpo.set_version(autodock_parameter_version)
        self.autodock_parameter_version = autodock_parameter_version

    def set_ligand(self, ligand_filename):
        verbose = self.verbose
        self.ligand_filename = os.path.basename(ligand_filename)
        if verbose: print("set ligand_filename to", self.ligand_filename)
        self.dpo.set_ligand(ligand_filename)
        # expect a filename like ind.out.pdbq: get 'ind' from it
        self.ligand_stem = self.ligand_filename.split('.')[0]
        if verbose: print("set ligand_stem to", self.ligand_stem)
        self.ligand = Read(ligand_filename)[0]
        if verbose: print("read ", self.ligand.name)
        # set dpo:
        # move
        self.dpo['move']['value'] = self.ligand_filename
        if verbose: print("set move to ", self.dpo['move']['value'])
        # ndihe
        # assumes ligand has torTree
        self.dpo['ndihe']['value'] = self.ligand.parser.keys.count("BRANCH")
        # self.dpo['ndihe']['value'] = len(self.ligand.torTree.torsionMap)
        if verbose: print("set ndihe to ", self.dpo['ndihe']['value'])
        # torsdof
        # caution dpo['torsdof']['value'] is a list [ndihe, 0.3113]
        self.dpo['torsdof']['value'][0] = self.ligand.TORSDOF
        if verbose: print("set torsdof to ", self.dpo['torsdof']['value'])
        # types
        d = {}
        for a in self.ligand.allAtoms:
            d[a.autodock_element] = 1
        sortKeyList = ['C', 'A', 'N', 'O', 'S', 'H', 'P', 'n', 'f', 'F', 'c', 'b', 'I', 'M']
        lig_types = ""
        for t in sortKeyList:
            if t in list(d.keys()):
                lig_types = lig_types + t
        self.ligand.types = lig_types
        self.dpo['types']['value'] = self.ligand.types
        if verbose: print("set types to ", self.dpo['types']['value'])
        # about
        self.ligand.getCenter()
        cen = self.ligand.center
        self.dpo['about']['value'] = [round(cen[0], 3), round(cen[1], 3), round(cen[2], 3)]
        # print "DP:line 2854: cen =", cen
        if verbose: print("set about to ", self.dpo['about']['value'])

    def set_receptor(self, receptor_filename):
        self.receptor_filename = os.path.basename(receptor_filename)
        self.receptor_stem = self.receptor_filename.split('.')[0]
        self.dpo.set_receptor(receptor_filename)
        self.dpo['types']['value'] = self.ligand.types

    def set_docking_parameters(self, **kw):
        """Any docking parameters should be set here
        """
        # like this:
        # newdict = {'ga_num_evals':1750000, 'ga_pop_size':150,
        #            'ga_run':20, 'rmstol':2.0}
        # self.mv.dpo['<parameter>']['value'] = <new value>
        for parm, newvalue in list(kw.items()):
            if self.verbose: print('set dpo for ', parm, ' to ', newvalue, ' check=', self.dpo[parm]['value'])
            self.dpo[parm]['value'] = newvalue
        # self.dpo['ga_num_evals']['value'] = 1750000
        # self.dpo['ga_run']['value'] = 20
        # self.dpo['ga_pop_size']['value'] = 150
        # self.dpo['rmstol']['value'] = 2.0

    def write_dpf(self, dpf_filename,
                  parm_list=genetic_algorithm_local_search_list):
        if not dpf_filename:
            dpf_filename = "%s%s%s%s" % \
                           (self.ligand_stem, "_",
                            self.receptor_stem, ".dpf")
        # now that we have a filename...
        if self.verbose: print("writing ", dpf_filename)
        self.dpo.write(dpf_filename, parm_list)


class DockingParameter4FileMaker:
    """Accept a <ligand>.pdbqt and <receptor>.pdbqt and create
    <ligand>_<receptor>4.dpf
    """

    def __init__(self, verbose=None, autodock_parameter_version='4.0'):
        self.verbose = verbose
        self.dpo = DockingParameters()
        # version in ['3.05', '4.0', '4.1', '4.2']:
        self.dpo.set_version(autodock_parameter_version)
        self.autodock_parameter_version = autodock_parameter_version

    def getTypes(self, molecule):
        if not len(molecule.allAtoms.bonds[0]):
            molecule.buildBondsByDistance()
        ad4_typer = AutoDock4_AtomTyper(verbose=self.verbose)
        ad4_typer.setAutoDockElements(molecule)
        dict = {}
        for a in molecule.allAtoms:
            dict[a.autodock_element] = 1
        d_types = list(dict.keys())
        d_types.sort()
        mol_types = d_types[0]
        for t in d_types[1:]:
            mol_types = mol_types + " " + t
        if self.verbose: print("end of getTypes: types=", mol_types, ' class=', mol_types.__class__)
        return mol_types

    def set_write_all(self, value):
        if value == 'True':
            value = 1
        if value == True:
            value = 1
        if value == 'False':
            value = 0
        if value == False:
            value = 0
        verbose = self.verbose
        self.dpo['write_all']['value'] = value
        if verbose: print("set write_all to", self.dpo['write_all']['value'])

    def set_ligand(self, ligand_filename):
        verbose = self.verbose
        self.ligand_filename = os.path.basename(ligand_filename)
        if verbose: print("set ligand_filename to", self.ligand_filename)
        self.dpo.set_ligand(ligand_filename)
        # expect a filename like ind.out.pdbq: get 'ind' from it
        self.ligand_stem = self.ligand_filename.split('.')[0]
        if verbose: print("set ligand_stem to", self.ligand_stem)
        self.ligand = Read(ligand_filename)[0]
        if self.ligand == None:
            print('ERROR reading: ', ligand_filename)
            return
        if verbose: print("read ", self.ligand.name)
        # set dpo:
        # move
        self.dpo['move']['value'] = self.ligand_filename
        if verbose: print("set move to ", self.dpo['move']['value'])
        # ndihe
        # assumes ligand has torTree
        self.dpo['ndihe']['value'] = self.ligand.parser.keys.count("BRANCH")
        # self.dpo['ndihe']['value'] = len(self.ligand.torTree.torsionMap)
        if verbose: print("set ndihe to ", self.dpo['ndihe']['value'])
        # torsdof
        # caution dpo['torsdof4']['value'] is a list [ndihe, 0.274]
        try:
            self.dpo['torsdof4']['value'][0] = self.ligand.TORSDOF
        except:
            print('!unable to use ligand.TORSDOF! torsdof always set to ligand.ndihe=', self.ligand.ndihe)
            self.dpo['torsdof4']['value'][0] = self.ligand.ndihe
        if verbose: print("set torsdof4 to ", self.dpo['torsdof4']['value'])
        # types
        self.ligand.types = self.getTypes(self.ligand)
        self.dpo['ligand_types']['value'] = self.ligand.types
        if verbose: print("set types to ", self.dpo['ligand_types']['value'])
        # about
        cen = self.ligand.getCenter()
        self.dpo['about']['value'] = [round(cen[0], 3), round(cen[1], 3), \
                                      round(cen[2], 3)]
        if verbose: print("set about to ", self.dpo['about']['value'])
        # print "DP:line 2976: cen =", cen

    def set_receptor(self, receptor_filename):
        self.receptor_filename = os.path.basename(receptor_filename)
        self.receptor_stem = self.receptor_filename.split('.')[0]
        self.dpo.set_receptor(receptor_filename)

    def set_flexres(self, flexres_filename):
        flexmol = Read(flexres_filename)[0]
        flexres_filename = os.path.basename(flexres_filename)
        self.dpo['flexres_flag']['value'] = True
        self.dpo['flexres']['value'] = flexres_filename
        # make sure each atom type in flexres molecule is in ligand_types
        d = {}
        current_types = self.dpo['ligand_types']['value'].split()
        for t in current_types:
            d[t] = 1
        for a in flexmol.allAtoms:
            d[a.autodock_element] = 1
        self.dpo['ligand_types']['value'] = ''.join(list(d.keys()))

    def set_docking_parameters(self, **kw):
        """Any docking parameters should be set here
        """
        # like this:
        # newdict = {'ga_num_evals':1750000, 'ga_pop_size':150,
        #            'ga_run':20, 'rmstol':2.0}
        # self.mv.dpo['<parameter>']['value'] = <new value>
        for parm, newvalue in list(kw.items()):
            self.dpo[parm]['value'] = newvalue
            if parm == 'set_sw1':
                self.dpo['set_psw1']['value'] = not newvalue
            if parm == 'set_psw1':
                self.dpo['set_sw1']['value'] = not newvalue
            if parm == 'flexres':
                self.set_flexres(newvalue)
            if parm == 'write_all':
                self.set_write_all(newvalue)

    def write_dpf(self, dpf_filename,
                  parm_list=genetic_algorithm_local_search_list4,
                  pop_seed=False):
        if not dpf_filename:
            dpf_filename = "%s%s%s%s" % \
                           (self.ligand_stem, "_",
                            self.receptor_stem, ".dpf")
        # set initial conformation
        if pop_seed:
            self.dpo['tran0']['value'] = self.dpo['about']['value']
            self.dpo['quaternion0']['value'] = '0.0 0. 0. 1.'
            # self.dpo['axisangle0']['value'] = '1.0 0. 0. 0.'
            ntors = self.dpo['ndihe']['value']
            if ntors > 0:
                dihe0 = '0. ' * self.dpo['ndihe']['value']
                dihe0.rstrip()
                self.dpo['dihe0']['value'] = dihe0
            else:
                parm_list.remove('dihe0')
        # now that we have a filename...
        if self.verbose:
            print("writing ", dpf_filename)
        self.dpo.write4(dpf_filename, parm_list)


class DockingParameter42FileMaker:
    """Accept a <ligand>.pdbqt and <receptor>.pdbqt and create
    <ligand>_<receptor>42.dpf
    """

    def __init__(self, verbose=None, autodock_parameter_version="4.2", pop_seed=False):
        self.verbose = verbose
        self.dpo = DockingParameters()
        self.dpo.set_version(autodock_parameter_version)
        self.autodock_parameter_version = autodock_parameter_version
        self.pop_seed = pop_seed

    def getTypes(self, molecule):
        if not len(molecule.allAtoms.bonds[0]):
            molecule.buildBondsByDistance()
        ad4_typer = AutoDock4_AtomTyper(verbose=self.verbose)
        ad4_typer.setAutoDockElements(molecule)
        mol_types = " ".join(list(set(molecule.allAtoms.autodock_element)))
        # print "DP42FM: end of getTypes on " + molecule.name + ": types=" + mol_types + ' class=', mol_types.__class__
        return mol_types

    def set_write_all(self, value):
        if value == 'True':
            value = 1
        if value == True:
            value = 1
        if value == 'False':
            value = 0
        if value == False:
            value = 0
        verbose = self.verbose
        self.dpo['write_all']['value'] = value
        if verbose: print("set write_all to", self.dpo['write_all']['value'])

    def set_ligand(self, ligand_filename):
        verbose = self.verbose
        self.ligand_filename = os.path.basename(ligand_filename)
        if verbose: print("set ligand_filename to", self.ligand_filename)
        self.dpo.set_ligand(ligand_filename)
        # expect a filename like ind.out.pdbq: get 'ind' from it
        self.ligand_stem = self.ligand_filename.split('.')[0]
        if verbose: print("set ligand_stem to", self.ligand_stem)
        self.ligand = Read(ligand_filename)[0]
        self.dpo.ligand = self.ligand
        if self.ligand == None:
            print('ERROR reading: ', ligand_filename)
            return
        if verbose: print("read ", self.ligand.name)
        # set dpo:
        # move
        self.dpo['move']['value'] = self.ligand_filename
        if verbose: print("set move to ", self.dpo['move']['value'])
        # ndihe
        # assumes ligand has torTree
        self.dpo['ndihe']['value'] = self.ligand.parser.keys.count("BRANCH")
        # self.dpo['ndihe']['value'] = len(self.ligand.torTree.torsionMap)
        if verbose: print("set ndihe to ", self.dpo['ndihe']['value'])
        # torsdof
        # caution dpo['torsdof4']['value'] is a list [ndihe, 0.274]
        try:
            self.dpo['torsdof4']['value'][0] = self.ligand.TORSDOF
        except:
            print('!unable to use ligand.TORSDOF! set torsdof to ligand.ndihe=', self.ligand.ndihe)
            self.dpo['torsdof4']['value'][0] = self.ligand.ndihe
        if verbose: print("set torsdof4 to ", self.dpo['torsdof4']['value'])
        # types
        self.ligand.types = self.getTypes(self.ligand)
        self.dpo['ligand_types']['value'] = self.ligand.types
        if verbose: print("set ligand types to ", self.dpo['ligand_types']['value'])
        # about
        cen = self.ligand.getCenter()
        self.dpo['about']['value'] = [round(cen[0], 3), round(cen[1], 3), \
                                      round(cen[2], 3)]
        if verbose: print("set about to ", self.dpo['about']['value'])
        # print "DP:line 3122: cen =", cen

    def set_receptor(self, receptor_filename):
        self.receptor_filename = os.path.basename(receptor_filename)
        self.receptor_stem = self.receptor_filename.split('.')[0]
        self.dpo.set_receptor(receptor_filename)

    def set_flexres(self, flexres_filename):
        flexmol = Read(flexres_filename)[0]
        flexres_filename = os.path.basename(flexres_filename)
        self.dpo['flexres_flag']['value'] = True
        self.dpo['flexres']['value'] = flexres_filename
        # make sure each atom type in flexres molecule is in ligand_types
        d = {}
        current_types = self.dpo['ligand_types']['value'].split()
        for t in current_types:
            d[t] = 1
        for a in flexmol.allAtoms:
            d[a.autodock_element] = 1
        self.dpo['ligand_types']['value'] = ''.join(list(d.keys()))

    def set_docking_parameters(self, **kw):
        """Any docking parameters should be set here
        """
        # like this:
        # newdict = {'ga_num_evals':1750000, 'ga_pop_size':150,
        #            'ga_run':20, 'rmstol':2.0}
        # self.mv.dpo['<parameter>']['value'] = <new value>
        for parm, newvalue in list(kw.items()):
            self.dpo[parm]['value'] = newvalue
            if parm == 'set_sw1':
                self.dpo['set_psw1']['value'] = not newvalue
            if parm == 'set_psw1':
                self.dpo['set_sw1']['value'] = not newvalue
            if parm == 'flexres':
                self.set_flexres(newvalue)
            if parm == 'write_all':
                self.set_write_all(newvalue)

    def write_dpf(self, dpf_filename,
                  parm_list=genetic_algorithm_local_search_list4_2,
                  pop_seed=False):
        if not dpf_filename:
            dpf_filename = "%s%s%s%s" % \
                           (self.ligand_stem, "_",
                            self.receptor_stem, ".dpf")
        # set initial conformation
        if pop_seed or self.pop_seed:
            self.dpo['tran0']['value'] = self.dpo['about']['value']
            self.dpo['quaternion0']['value'] = '0. 0. 0. 1.'
            ntors = self.dpo['ndihe']['value']
            if ntors > 0:
                dihe0 = '0. ' * self.dpo['ndihe']['value']
                dihe0.rstrip()
                self.dpo['dihe0']['value'] = dihe0
            else:
                parm_list.remove('dihe0')
        # now that we have a filename...
        if self.verbose:
            print("writing ", dpf_filename)
        self.dpo.write42(dpf_filename, parm_list)

        # now that we have a filename...


# ------------------------------------------------------------------
# ConfigFileMaker: added 9/26/2012 @@
# ------------------------------------------------------------------
class ConfigFileMaker():
    def __init__(self, receptor='', ligand='', **kw):
        self.receptor = receptor
        self.ligand = ligand
        # validate filenames:
        if not len(self.receptor):
            print("No input receptor specified!")
            raise IOError
        if not len(self.ligand):
            print("No input ligand specified!")
            raise IOError
        self.flexres = kw.get("flexres", "")
        self.center_x = kw.get("center_x", 0.0)
        self.center_y = kw.get("center_y", 0.0)
        self.center_z = kw.get("center_z", 0.0)
        self.size_x = kw.get("size_x", 0.0)
        self.size_y = kw.get("size_y", 0.0)
        self.size_z = kw.get("size_z", 0.0)
        self.out = kw.get("out", "")
        # optional parameters:
        self.log = kw.get("log", "")
        self.cpu = kw.get("cpu", "")
        self.seed = kw.get("seed", "")
        self.exhaustiveness = kw.get("exhaustiveness", 8)
        self.num_modes = kw.get("num_modes", 9)
        self.energy_range = kw.get("energy_range", 3)
        self.config = kw.get("config", "")
        self.score_only = kw.get("score_only", False)
        self.local_only = kw.get("local_only", False)
        self.randomize_only = kw.get("randomize_only", False)
        self.weight_gauss1 = kw.get("weight_gauss1", -0.035579)
        self.weight_gauss2 = kw.get("weight_gauss2", -0.005156)
        self.weight_repulsion = kw.get("weight_repulsion", 0.840245)
        self.weight_hydrophobic = kw.get("weight_hydrophobic", -0.035069)
        self.weight_hydrogen = kw.get("weight_hydrogen", -0.587439)
        self.weight_rot = kw.get("weight_rot", 0.058459999)
        self.verbose = kw.get('verbose', False)
        # validate input:
        if not self.score_only:
            if self.size_x <= 0:
                print("Invalid size_x specified: must be >0")
            if self.size_y <= 0:
                print("Invalid size_y specified: must be >0")
            if self.size_z <= 0:
                print("Invalid size_z specified: must be >0")
        if self.center_x == 0 and self.verbose:
            print("@@ setting center_x to 0.0")
        if self.center_y == 0 and self.verbose:
            print("@@ setting center_y to 0.0")
        if self.center_z == 0 and self.verbose:
            print("@@setting center_z to 0.0")

    #############################################################################
    # set methods
    #############################################################################
    # Molecule filenames:
    # input
    def set_receptor(self, receptor):
        # ?os.path.basename?
        self.receptor = receptor

    def set_ligand(self, ligand):
        # ?os.path.basename?
        self.ligand = ligand

    def set_flexres(self, flexres):
        # ?os.path.basename?
        self.flexres = flexres

    # output (will contain results):
    #  note: default vina output is: 'ligandstem'_out.pdbqt
    def set_out(self, out):
        # ?validate?
        self.out = out

    # Search space:
    # ?validate?@@
    # set center
    def set_center_x(self, center_x):
        self.center_x = center_x

    def set_center_y(self, center_y):
        self.center_y = center_y

    def set_center_z(self, center_z):
        self.center_z = center_z

    # set size (Angstroms)
    def set_size_x(self, size_x):
        self.size_x = size_x

    def set_size_y(self, size_y):
        self.size_y = size_y

    def set_size_z(self, size_z):
        self.size_z = size_z

    # Misc options
    def set_cpu(self, cpu):
        # ?validate? must be an integer
        self.cpu = cpu

    def set_seed(self, seed):
        # ?validate? an integer??
        self.seed = seed

    def set_exhaustiveness(self, exhaustiveness=8):
        # ?validate? an integer??
        self.exhaustiveness = exhaustiveness

    def set_num_modes(self, num_modes=9):
        # ?validate? an integer??
        self.num_modes = num_modes

    def set_energy_range(self, energy_range=3):
        # ?validate? must to be  positive??
        if energy_range < 0:
            print("energy_range less than zero! resetting to default 3")
            energy_range = 3
        self.energy_range = energy_range

    # Advanced options
    #  --score_only                                      score only - search space
    #                                                    can be omitted
    #  --local_only                                      do local search only
    #  --randomize_only                                  randomize input, attempting
    #                                                    to avoid clashes
    #  --weight_gauss1 arg (=-0.035579)                  gauss_1 weight
    #  --weight_gauss2 arg (=-0.005156)                  gauss_2 weight
    #  --weight_repulsion arg (=0.84024500000000002)     repulsion weight
    #  --weight_hydrophobic arg (=-0.035069000000000003) hydrophobic weight
    #  --weight_hydrogen arg (=-0.58743900000000004)     Hydrogen bond weight
    #  --weight_rot arg (=0.058459999999999998)          N_rot weight

    # ?validate these set methods? check for a specific type??
    def set_score_only(self, score_only):
        self.score_only = score_only

    def set_local_only(self, local_only):
        self.local_only = local_only

    def set_randomize_only(self, randomize_only):
        self.randomize_only = randomize_only

    def set_weight_gauss1(self, weight_gauss1):
        self.weight_gauss1 = weight_gauss1

    def set_weight_gauss2(self, weight_gauss2):
        self.weight_gauss2 = weight_gauss2

    def set_weight_repulsion(self, weight_repulsion):
        self.weight_repulsion = weight_repulsion

    def set_weight_hydrophobic(self, weight_hydrophobic):
        self.weight_hydrophobic = weight_hydrophobic

    def set_weight_hydrogen(self, weight_hydrogen):
        self.weight_hydrogen = weight_hydrogen

    def set_weight_rot(self, weight_rot):
        self.weight_rot = weight_rot

    #############################################################################
    # write method
    #############################################################################
    def write(self, filename=None):
        """Write the current state to a file
        file is a pointer to a writeable file
        param_list is a list of parameter strings.
        """
        if filename is None:
            filename = self.config
        fptr = open(filename, 'w')
        # parameter names
        mol_parm_list = ['receptor', 'flex', 'ligand']
        # variable names
        mol_var_list = [self.receptor, self.flexres, self.ligand]
        if self.flexres == "":
            mol_parm_list = ['receptor', 'ligand']
            mol_var_list = [self.receptor, self.ligand]
        for p, v in zip(mol_parm_list, mol_var_list):
            ostr = "%s = %s\n" % (p, str(v))
            fptr.write(ostr)
        if self.score_only != False:
            ostr = "score_only = True\n"
            fptr.write(ostr)
        else:
            searchspace_parms = ['center_x', 'center_y', 'center_z', 'size_x', 'size_y', 'size_z']
            searchspace_vars = [self.center_x, self.center_y, self.center_z, self.size_x, self.size_y, self.size_z]
            # VALIDATE SIZE
            if self.size_x == 0:
                print("invalid x-dim search space size: ", self.size_x)
                raise ValueError
            if self.size_y == 0:
                print("invalid y-dim search space size: ", self.size_y)
                raise ValueError
            if self.size_z == 0:
                print("invalid z-dim search space size: ", self.size_z)
                raise ValueError
            # searchspace required:
            for p, v in zip(searchspace_parms, searchspace_vars):
                ostr = "%s = %s\n" % (p, str(v))
                fptr.write(ostr)
            # optional search parameters:
            if self.log != "":
                ostr = "log = %s\n" % (self.log)
                fptr.write(ostr)
            if self.cpu != "":
                ostr = "cpu = %n\n" % (self.cpu)
                fptr.write(ostr)
            if self.seed != "":
                ostr = "seed = %s\n" % (self.seed)
                fptr.write(ostr)
            if self.exhaustiveness != 8:
                ostr = "exhaustiveness = %n\n" % (self.exhaustiveness)
                fptr.write(ostr)
            if self.num_modes != 9:
                ostr = "num_modes = %d\n" % (self.num_modes)
                fptr.write(ostr)
            if self.energy_range != 3.:  # @@ "" indicates var not set
                ostr = "energy_range = %f\n" % (self.energy_range)
                fptr.write(ostr)
            # OTHER optional parameters:
            if self.local_only != False:
                ostr = "local_only = True\n"
                fptr.write(ostr)
            # self.randomize_only = kw.get("randomize_only", False)
            if self.randomize_only != False:
                if self.out == "":
                    stem, ext = os.path.basename(self.ligand).split('.')
                    new_filename = stem + '_out.' + ext
                    if self.verbose:
                        print("@@ randomize_only is set to True@@ running vina using this config file ", filename,
                              " will write new coordinates for ", self.ligand, "in ", new_filename, "@@")
                ostr = "randomize_only = True\n"
                fptr.write(ostr)
            # self.weight_gauss1 = kw.get("weight_gauss1", -0.035579)
            if self.weight_gauss1 != -0.035579:
                ostr = "weight_gauss1 = %f\n" % (self.weight_gauss1)
                fptr.write(ostr)
            # self.weight_gauss2 = kw.get("weight_gauss2", -0.005156)
            if self.weight_gauss2 != -0.005156:
                ostr = "weight_gauss2 = %f\n" % (self.weight_gauss2)
                fptr.write(ostr)
            # self.weight_repulsion = kw.get("weight_repulsion", 0.840245)
            if self.weight_repulsion != 0.840245:
                ostr = "weight_repulsion = %f\n" % (self.weight_gauss2)
                fptr.write(ostr)
            # self.weight_hydrophobic = kw.get("weight_hydrophobic", -0.035069)
            if self.weight_hydrophobic != -0.035069:
                ostr = "weight_hydrophobic = %f\n" % (self.weight_gauss2)
                fptr.write(ostr)
            # self.weight_hydrogen = kw.get("weight_hydrogen", -0.587439)
            if self.weight_hydrogen != -0.587439:
                ostr = "weight_hydrogen = %f\n" % (self.weight_hydrogen)
                fptr.write(ostr)
            # self.weight_rot = kw.get("weight_rot", 0.058459999)
            if self.weight_rot != 0.058459999:
                ostr = "weight_rot = %f\n" % (self.weight_rot)
                fptr.write(ostr)
            if self.out != "":
                ostr = "out = %s\n" % (self.out)
                fptr.write(ostr)
            if self.log != "":
                ostr = "log = %s\n" % (self.log)
                fptr.write(ostr)
            if self.cpu != "":
                ostr = "cpu = %n\n" % (self.cpu)
                fptr.write(ostr)
            if self.seed != "":
                ostr = "seed = %s\n" % (self.seed)
                fptr.write(ostr)
            if self.num_modes != 9:
                ostr = "num_modes = %s\n" % (self.num_modes)
                fptr.write(ostr)
            if self.energy_range != 3:
                ostr = "energy_range = %s\n" % (self.energy_range)
                fptr.write(ostr)
            print("closing new config file: " + filename)
            fptr.close()


# ------------------------------------------------------------------
# END ConfigFileMaker
# ------------------------------------------------------------------

if __name__ == '__main__':
    import os

    if len(sys.argv) > 1:
        dpo = DockingParameters()
        try:
            dpo.read(sys.argv[1])
        except IOError as msg:
            print("IOError: %s" % (msg))
            exit(-1)

        tmp_filename = '/tmp/%d.dpf' % (os.getpid())
        dpo.write(tmp_filename, dpo.file_params)

        # instantiate a new DP object and read from file just written
        new_dpo = DockingParameters()
        new_dpo.read(tmp_filename)

        if dpo.file_params != new_dpo.file_params:
            for p in dpo.file_params:
                if not p in new_dpo.file_params:
                    print("%s missing from written dpf!" % p['keyword'])
            for p in new_dpo.file_params:
                if not p in dpo.file_params:
                    print("Unexpected %s in written dpf!" % p['keyword'])
        for p in new_dpo.file_params:
            if not dpo[p]['value'] == new_dpo[p]['value']:
                print("Parameter Discrepancy: %s" % (p))
                print("    I wrote this: %s" % (dpo[p]['value']))
                print("But, I read this: %s" % (new_dpo[p]['value']))
    else:
        dpo = DockingParameters('receptor.test.pdbqs', 'ligand.pdbq')

        # compare docking_parameter_list with keys of DockingParameter
        # dictionary and report discrepancies

        if docking_parameter_list != list(dpo.keys()):
            print("Comparing docking paramter dictionary and list...")
            for p in list(dpo.keys()):
                if not p in docking_parameter_list:
                    print(">> %s in dictionary but not docking_parameter_list" % (p))
            for p in docking_parameter_list:
                if not p in list(dpo.keys()):
                    print(">> %s in docking_parameter_list but not dictionary" % (p))
            print("done with comparison.")

        # This should also run through all the parameters, write default
        # values to a file and then read the file compare with the defaults
        print("Here are the parameters implemented so far (write):")
        dpo.write('', implemented)
        print()
        print("Here are the parameters not implemented yet (write):")
        for p in docking_parameter_list:
            if not p in implemented:
                print(p)

        tmp_filename = '/tmp/%d.dpf' % (os.getpid())
        dpo.write(tmp_filename, implemented)

        new_dpo = DockingParameters()
        new_dpo.read(tmp_filename)

        # write the new dpo to /dev/null to update the value fields
        for p in implemented:
            if not dpo[p]['value'] == new_dpo[p]['value']:
                print("%s: %s(%s) -> %s(%s)" % (dpo[p]['keyword'],
                                                dpo[p]['value'],
                                                type(dpo[p]['value']).__name__,
                                                new_dpo[p]['value'],
                                                type(new_dpo[p]['value']).__name__))
