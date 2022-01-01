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
#  Modification date: 1/1/22, 4:48 PM                                                              #
#                                                                                                  #
# ##################################################################################################

#############################################################################
#
# Author: Ruth HUEY, Michel SANNER
#
# Copyright: M. Sanner TSRI 2000
#
#############################################################################


# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/GridParameters.py,v 1.39.6.2 2016/02/11 09:24:07 annao Exp $
#
#
# $Id: GridParameters.py,v 1.39.6.2 2016/02/11 09:24:07 annao Exp $
#
#
#
#

import glob
import os.path
import sys
from collections import UserDict
from math import ceil

import numpy

from AutoDockTools.atomTypeTools import AutoDock4_AtomTyper
from MolKit import Read
from .energyConstants import Rij, epsij, SolVol, SolPar, SolCon

grid_parameter_list = [
    'receptor',
    'gridfld',
    'npts',
    'spacing',
    'gridcenter',
    'types',
    'smooth',
    'map',
    'elecmap',
    'dielectric',
    'fmap'
]

grid_parameter_list4 = [
    'npts',
    'custom_parameter_file',
    'gridfld',
    'spacing',
    'receptor_types',
    'ligand_types',
    'receptor',
    'gridcenter',
    'smooth',
    'map',
    'elecmap',
    'dsolvmap',
    'dielectric4',
]


class GridParameters(UserDict):
    def __init__(self, receptor_filename='', ligand_filename=''):
        UserDict.__init__(self)
        basename = os.path.basename(receptor_filename)
        self.receptor_filename = basename
        self.receptor_stem = os.path.splitext(basename)[0]
        # self.receptor_stem = basename[:string.rfind(basename, '.')]

        # if the grid parameters have been read from a file,
        # then the following instance variables will be set:
        self.gpf_filename = ''
        self.gpf_written_filename = ''
        self.file_params = []
        # begin dictionary
        self['constant'] = {
            'keyword': 'constant',
            'default': [],
            'comment': "grid map constant energy",
            'value': []
        }
        self['covalent_coords'] = {
            'keyword': 'covalent_coords',
            'default': [],
            'comment': "covalent_coords",
            'value': []
        }
        self['covalent_constant'] = {
            'keyword': 'covalent_constant',
            'default': -1.780,
            'comment': "covalent_constant",
            'value': -1.780
        }
        self['covalent_energy_barrier'] = {
            'keyword': 'covalent_energy_barrier',
            'default': 1000.,
            'comment': "covalent_energy barrier height",
            'value': 1000.
        }
        self['covalent_half_width'] = {
            'keyword': 'covalent_half_width',
            'default': 5.0,
            'comment': "covalent_half_width ",
            'value': 5.0
        }
        self['covalentmap'] = {
            'keyword': 'covalentmap',
            'default': 0,
            'comment': "covalent map",
            'value': 0
        }
        self['dielectric'] = {
            'keyword': 'dielectric',
            'default': -.1146,
            'comment': "<0, distance-dep.diel;>0, constant",
            'value': -.1146
        }
        self['dielectric4'] = {
            'keyword': 'dielectric',
            'default': -.1465,  # new sept/29/05:value from august recal
            'comment': "<0, AD4 distance-dep.diel;>0, constant",
            'value': -.1465  # new sept/29/05:value from august recal
        }
        self['dsolvmap'] = {
            'keyword': 'dsolvmap',
            'default': self.receptor_stem + '.d.map',
            'comment': "desolvation potential map",
            'value': self.receptor_stem + '.d.map'
        }
        self['elecmap'] = {
            'keyword': 'elecmap',
            'default': self.receptor_stem + '.e.map',
            'comment': "electrostatic potential map",
            'value': self.receptor_stem + '.e.map'
        }
        self['fmap'] = {
            'keyword': 'fmap',
            'default': 0,
            'comment': "floating point potential gridmap",
            'value': 0
        }
        self['gridcenter'] = {
            'keyword': 'gridcenter',
            'default': 'auto',
            'comment': "xyz-coordinates or auto",
            'value': 'auto'
        }
        self['gridcenterAuto'] = {
            'keyword': 'gridcenterAuto',
            'default': 1,
            'comment': "xyz-coordinates or auto",
            'value': 1
        }
        self['gridfld'] = {
            'keyword': 'gridfld',
            'default': self.receptor_stem + '.maps.fld',
            'comment': "grid_data_file",
            'value': self.receptor_stem + '.maps.fld'
        }
        self['ligand_types'] = {
            'keyword': 'ligand_types',
            'default': 'A C HD N NA OA SA',
            'comment': "ligand atom types",
            'value': 'A C HD N NA OA SA',
        }
        self['map'] = {
            'keyword': 'map',
            'default': "",
            'comment': "atom-specific affinity map",
            'value': ""
        }
        self['mset'] = {
            'keyword': 'mset',
            'default': "CNOSHHH",
            'comment': "atom-specific affinity map",
            'value': "CNOSHHH"
        }
        self['nbp_r_eps'] = {
            'keyword': 'nbp_r_eps',
            'default': [],
            'comment': "lj",
            'value': []
        }
        self['NHB'] = {
            'keyword': 'NHB',
            'default': 1,
            'comment': 'model N-H hydrogen bonds',
            'value': 1
        }
        self['npts'] = {
            'keyword': 'npts',
            'default': [40, 40, 40],
            'comment': "num.grid points in xyz",
            'value': [40, 40, 40]
        }
        self['OHB'] = {
            'keyword': 'OHB',
            'default': 1,
            'comment': 'model O-H hydrogen bonds',
            'value': 1
        }
        self['custom_parameter_file'] = {
            'keyword': 'custom_parameter_file',
            'default': 0,
            'comment': "use custom parameter library",
            'value': 0,
        }
        self['parameter_file'] = {
            'keyword': 'parameter_file',
            'default': 'AD4_parameters.dat',
            'comment': "force field default parameter file",
            'value': 'AD4_parameters.dat',
        }
        self['receptor'] = {
            'keyword': 'receptor',
            'default': self.receptor_stem + '.pdbqs',
            'comment': "macromolecule",
            'value': self.receptor_stem + '.pdbqs',
        }
        self['receptor_types'] = {
            'keyword': 'receptor_types',
            'default': 'A C HD N NA OA SA',
            'comment': "receptor atom types",
            'value': 'A C HD N NA OA SA',
        }
        self['SHB'] = {
            'keyword': 'SHB',
            'default': 1,
            'comment': 'model S-H hydrogen bonds',
            'value': 1
        }
        self['smooth'] = {
            'keyword': 'smooth',
            'default': 0.5,
            'comment': "store minimum energy w/in rad(A)",
            'value': 0.5
        }
        self['sol_par'] = {
            'keyword': 'sol_par',
            'default': [],
            'comment': "atomic fragmental volumen, solvation parm",
            'value': []
        }
        self['spacing'] = {
            'keyword': 'spacing',
            'default': 0.375,
            'comment': "spacing(A)",
            'value': 0.375
        }
        self['types'] = {
            'keyword': 'types',
            'default': 'CAONSH',
            'comment': "atom type names",
            'value': 'CAONSH',
        }
        # end dictionary

        self.set_receptor(receptor_filename)  # also sets self.receptor_stem
        self.set_ligand(ligand_filename)
        self.boolean_param_list = [
            'covalentmap',
            'fmap',
        ]
        # end __init__

    def set_ligand(self, ligand_filename):
        self.ligand_filename = os.path.basename(ligand_filename)
        # this should set types

    def set_ligand_types3(self, ligand_types4):
        d = {}
        for t in ligand_types4:
            if len(t) == 1:
                d[t] = 1
            elif t[1] in ['A', 'D']:  # NA,SA,OA,HD
                d[t[0]] = 1
            elif t in ['Cl', 'CL', 'cl']:  # special case: chlorine
                d['c'] = 1
            elif t in ['Br', 'BR', 'br']:  # special case: bromine
                d['b'] = 1
            elif t in ['Fe', 'FE', 'fe']:  # special case: iron
                d['f'] = 1
            else:
                print("unrecognized ligand_atom_type:", t)
        all_types = list(d.keys())
        all_types.sort()
        type_str = all_types[0]
        for t in all_types[1:]:
            type_str = type_str + t
        self['types']['value'] = type_str

    def set_receptor(self, receptor_filename):
        basename = os.path.basename(receptor_filename)
        self.receptor_filename = basename
        self.receptor_stem = os.path.splitext(basename)[0]
        # self.receptor_stem = basename[:string.rfind(basename, '.')]
        if receptor_filename != '':
            self['receptor']['value'] = basename
            self['gridfld']['value'] = self.receptor_stem + '.maps.fld'
            self['elecmap']['value'] = self.receptor_stem + '.e.map'

    #
    # read methods
    #
    def read(self, filename):
        """Read from and set the current state according to the file.
        """
        self.gpf_filename = filename
        gpf_ptr = open(filename)
        lines = gpf_ptr.readlines()
        gpf_ptr.close()

        self.file_params = []
        checkedTypes = []
        extraLigandTypes = []
        keys = list(self.keys())
        for line in lines:
            words = line.replace('\t', ' ').split()
            if words != [] and words[0][0] != '#':
                p = words[0]
                if p not in keys:
                    print("WARNING: unrecognized parameter in ", filename, ":\n", p)
                    continue
                # maintain a list of the parameters read from the file
                if self.file_params == [] or p != self.file_params[-1]:
                    self.file_params.append(p)

                # parse the line
                l = len(words)
                for i in range(l):
                    if words[i][0] == '#':
                        l = i
                        break
                values = words[1:l]
                if ((len(values) == 1) and
                        (type(self[p]['default']) != list)):
                    self[p]['value'] = self._get_val(values[0])
                    if words[0] == 'types':
                        # in this case have to set flags for possible new type
                        extraLigandTypes = self.checkLigTypes(values[0])
                elif words[0] == 'ligand_types':
                    self[p]['value'] = ''.join(words[1:l])
                elif words[0] == 'receptor_types':
                    self[p]['value'] = ''.join(words[1:l])
                elif words[0] == 'covalentmap':
                    # in this case set:
                    # covalent_ coords,constant,energy_barrier,half_width
                    self['covalentmap']['value'] = 1
                    self['covalent_half_width']['value'] = float(values[0])
                    self['covalent_energy_barrier']['value'] = float(values[1])
                    self['covalent_coords']['value'] = [float(values[2]), float(values[3]), float(values[4])]
                    self[p]['value'] = []
                elif words[0] == 'nbp_r_eps':
                    # in this case have to check for nhb,ohb,shb +mset
                    # in this case have to check for new type constants
                    ptype = words[-1]
                    if len(words[l]) == 1:
                        keyWord = words[l + 1]
                    else:
                        keyWord = words[l][1:]
                    mtype = keyWord.split('-')[0]
                    ntype = keyWord.split('-')[1]
                    if mtype in checkedTypes:
                        continue
                    if mtype in ['N', 'O', 'S'] and ntype == 'H':
                        # check for 12 6 vs 12 10 here
                        ind = mtype + 'HB'
                        if values[3] == '10':
                            self[ind]['value'] = 1
                        else:
                            self[ind]['value'] = 0
                        checkedTypes.append(mtype)
                    if mtype in extraLigandTypes:
                        i = ptype + mtype + ntype
                        Rij[i] = float(words[1])
                        epsij[i] = float(words[2])
                elif words[0] == 'sol_par':
                    if len(words[l]) == 1:
                        mtype = words[l + 1]
                    else:
                        mtype = words[l][1]
                    if mtype in extraLigandTypes:
                        SolVol[mtype] = float(values[0])
                        SolPar[mtype] = float(values[1])
                elif words[0] == 'constant':
                    if len(words[l]) == 1:
                        mtype = words[l + 1]
                    else:
                        mtype = words[l][1]
                    SolCon[mtype] = float(values[0])
                elif words[0] == 'gridcenter' and l > 1:
                    # need to convert to float
                    newvalue = [float(values[0]), float(values[1]), float(values[2])]
                    self['gridcenterAuto']['value'] = 0
                    self[p]['value'] = newvalue
                else:
                    self[p]['value'] = []
                    for v in values:
                        self[p]['value'].append(self._get_val(v))

    def checkLigTypes(self, typeStr):
        extraLigandTypes = []
        for t in typeStr:
            if t not in ['C', 'A', 'N', 'O', 'S', 'H', 'P', 'n', \
                         'f', 'F', 'c', 'b', 'I', 'M']:
                extraLigandTypes.append(t)
        return extraLigandTypes

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

    def read4(self, filename):
        """Read from and set the current state according to the AutoGrid4 file.
        """
        self.gpf_filename = filename
        gpf_ptr = open(filename)
        lines = gpf_ptr.readlines()
        gpf_ptr.close()

        keys = list(self.keys())
        self.file_params = []
        for line in lines:
            # print "reading ", line
            words = line.replace('\t', ' ').split()
            # print "words=", words
            if words != [] and words[0][0] != '#':
                p = words[0]
                if p not in keys:
                    print("WARNING: unrecognized parameter in ", filename, ":\n", p)
                    continue
                # print "p=", p
                # maintain a list of the parameters read from the file
                if self.file_params == [] or p != self.file_params[-1]:
                    self.file_params.append(p)

                # parse the line
                l = len(words)
                for i in range(l):
                    if words[i][0] == '#':
                        l = i
                        break
                values = words[1:l]
                if p == 'parameter_file':
                    self['custom_parameter_file']['value'] = 1
                    self['parameter_file']['value'] = values[0]
                elif ((len(values) == 1) and
                      (type(self[p]['default']) != list)):
                    self[p]['value'] = self._get_val(values[0])
                    # print "    value=", self[p]['value']
                    # if words[0]=='types':
                    #    #in this case have to set flags for possible new type
                    #    extraLigandTypes = self.checkLigTypes(values[0])
                # setting dielectric from a gpf is no longer supported
                # instead must be set in a parameter library file
                # elif p=='dielectric':
                #    self['dielectric4']['value'] = self._get_val(values[0])
                elif p == 'ligand_types':
                    self['ligand_types']['value'] = ''.join(words[1:l])
                elif p == 'receptor_types':
                    self['receptor_types']['value'] = ''.join(words[1:l])
                elif words[0] == 'covalentmap':
                    # in this case set:
                    # covalent_ coords,constant,energy_barrier,half_width
                    self['covalentmap']['value'] = 1
                    self['covalent_half_width']['value'] = float(values[1])
                    self['covalent_energy_barrier']['value'] = float(values[2])
                    self['covalent_coords']['value'] = [float(values[3]), float(values[4]), float(values[5])]
                    self[p]['value'] = []
                elif words[0] == 'gridcenter' and l > 1:
                    # need to convert to float
                    newvalue = [float(values[0]), float(values[1]), float(values[2])]
                    self['gridcenterAuto']['value'] = 0
                    self[p]['value'] = newvalue
                else:
                    # print "in else for ", p
                    self[p]['value'] = []
                    for v in values:
                        self[p]['value'].append(self._get_val(v))

    #
    # write methods
    #
    def write(self, filename, param_list):
        """Write the current state to a file

        file is a writeable file
        param_list is a list of parameter strings.
        For best results use the parameter_lists supplied by this class.
        """
        if filename == '':
            gpf_ptr = sys.stdout
        else:
            gpf_ptr = open(filename, 'w')

        types = self['types']['value']
        # FIX THIS:
        macroTypes = self['mset']['value']
        # macroTypes = ['C','N','O','S','H','H','H']

        for p in param_list:
            # maps are a special case
            if p == 'map':
                # hpos = 'H' in types
                for a in types:
                    gpf_ptr.write(self.make_map_string(p, a))
                    for t in macroTypes:
                        self.write_map_nbp(a, t, gpf_ptr)
                        # self.write_map_nbp(a, t, hpos, gpf_ptr)
                    self.write_constants(a, gpf_ptr)
            # all the other parameters handle themselves
            elif p == 'gridcenter' and self['gridcenterAuto']['value'] == 1:
                # if gridcenterAuto is true, reset p to 'auto' and write it
                self['gridcenter']['value'] = 'auto'
                gpf_ptr.write(self.make_param_string(p))
            elif p == 'fmap' and self['fmap']['value']:
                gpf_ptr.write(self.make_map_string(p, 'f'))
            elif p == 'covalentmap' and len(self['covalent_coords']['value']):
                gpf_ptr.write(self.make_covalentmap_string())
            else:
                gpf_ptr.write(self.make_param_string(p))

        if gpf_ptr != sys.stdout:
            gpf_ptr.close()
            self.gpf_filename = filename
            self.gpf_written_filename = filename

    def write_constants(self, a, gpf_ptr):
        try:
            outstring = 'sol_par  %5.2f %6.4f' % (SolVol[a], SolPar[a]) + \
                        '    # ' + a + ' atomic fragmental volume, solvation parameters\n'
        except KeyError:
            outstring = 'sol_par  0.000 0.000    #' \
                        + a + ' atomic fragmental volume, solvation parameters\n'
        gpf_ptr.write(outstring)
        try:
            outstring = 'constant  %5.3f  ' % SolCon[a] + \
                        '    # ' + a + ' grid map constant energy\n'
        except KeyError:
            outstring = 'constant  0.000          #' + a + ' grid map constant energy\n'
        gpf_ptr.write(outstring)

    def write_map_nbp(self, a, t, gpf_ptr):
        hbset = []
        for item in ['N', 'O', 'S']:
            ind = item + 'HB'
            if self[ind]['value']:
                hbset.append(item)
        if (a in hbset and t == 'H') or (a == 'H' and t in hbset):
            # if hpos and ((a in hbset and t=='H') or (a=='H' and t in hbset)):
            string_start = 'hb'
            string_nums = '12 10   # '
        else:
            string_start = 'lj'
            string_nums = '12  6   # '
        z = string_start + a + t
        try:
            outstring = 'nbp_r_eps %5.2f %9.7f ' % (Rij[z], epsij[z]) \
                        + string_nums + a + '-' + t + " " + string_start + '\n'
        except KeyError:
            outstring = 'nbp_r_eps 0.00 0.0000000 ' \
                        + string_nums + a + '-' + t + " " + string_start + '\n'
        gpf_ptr.write(outstring)

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

    def make_covalentmap_string(self):
        s = 'covalentmap ' + self['covalent_half_width']['value'] + ' '
        s = s + self['covalent_energy_barrier']['value'] + ' '
        s = s + self['covalent_coords']['value'] + '\n'
        return s

    def _make_string(self, p, val_str):
        # fix 1/2013 for bug report:
        # map bbbb_B99990001_mod_rigid.maps.fld# grid_data file
        return "%s %s%s # %s\n" % (p['keyword'],
                                   val_str,
                                   " " * (35 - (len(p['keyword']) + len(val_str))),
                                   p['comment'])

    # AD4

    def set_ligand4(self, ligand_filename, types=None):
        # this should set ligand_types
        # print "in set_ligand4: types=", types
        ftype = os.path.splitext(ligand_filename)[-1]
        if ftype != ".pdbqt":
            print("ligand_filename must be in pdbqt format")
            return "invalid input"
        self.ligand = Read(ligand_filename)[0]
        ligand = self.ligand
        ligand.buildBondsByDistance()
        if types is None:
            types = " ".join(list(set(ligand.allAtoms.autodock_element)))
        self['ligand_types']['value'] = types
        # print "set_ligand4: self['ligand_types']['value']=", self['ligand_types']['value']
        self.ligand_filename = os.path.basename(ligand_filename)
        self.ligand_stem = os.path.splitext(self.ligand_filename)[0]
        # print "GPO: set ligand_filename to ", self.ligand_filename

    def set_receptor4(self, receptor_filename, types=None):
        # this should set receptor_types
        ftype = os.path.splitext(receptor_filename)[-1]
        if ftype != ".pdbqt":
            print("receptor_filename must be in pdbqt format")
            return "invalid input"
        self.receptor = Read(receptor_filename)[0]
        receptor = self.receptor
        if types is None:
            types = " ".join(list(set(receptor.allAtoms.autodock_element)))
        self['receptor_types']['value'] = types
        basename = os.path.basename(receptor_filename)
        self.receptor_filename = basename
        self.receptor_stem = os.path.splitext(basename)[0]
        if receptor_filename != '':
            self['receptor']['value'] = basename
            self['gridfld']['value'] = self.receptor_stem + '.maps.fld'
            self['elecmap']['value'] = self.receptor_stem + '.e.map'
            self['dsolvmap']['value'] = self.receptor_stem + '.d.map'

    def write4(self, filename, param_list=grid_parameter_list4):
        """Write the current state to a file for AutoGrid4
        file is a writeable file
        param_list is a list of parameter strings.
        For best results use the parameter_lists supplied by this class.
        """
        if filename == '':
            gpf_ptr = sys.stdout
        else:
            gpf_ptr = open(filename, 'w')
        for p in param_list:
            if p == 'custom_parameter_file':
                if self['custom_parameter_file']['value']:
                    # self['parameter_file']['value'] = 'AD4_parameters.dat'
                    gpf_ptr.write(self.make_param_string('parameter_file'))
            elif p == 'map':
                # maps are a special case
                for s in self['ligand_types']['value'].split():
                    gpf_ptr.write(self.make_map_string(p, s))
            # all the other parameters handle themselves
            elif p == 'gridcenter' and self['gridcenterAuto']['value'] == 1:
                # if gridcenterAuto is true, reset p to 'auto' and write it
                self['gridcenter']['value'] = 'auto'
                gpf_ptr.write(self.make_param_string(p))
            elif p == 'dsolvmap':
                outstring = "dsolvmap %s              # desolvation potential map\n" % self['dsolvmap']['value']
                gpf_ptr.write(outstring)
            elif p == 'dielectric4':
                # now dielectric value can only be set in parameter file
                # val = self['dielectric4']['value']
                outstring = 'dielectric -0.1465                   # <0, AD4 distance-dep.diel;>0, constant\n'
                gpf_ptr.write(outstring)
            elif p == 'covalentmap' and len(self['covalent_coords']['value']):
                gpf_ptr.write(self.make_covalentmap_string())
            else:
                gpf_ptr.write(self.make_param_string(p))
        if gpf_ptr != sys.stdout:
            gpf_ptr.close()
            self.gpf_filename = filename
            self.gpf_written_filename = filename

    def write41(self, filename, param_list=grid_parameter_list4):
        """Write the current state to a file for AutoGrid41
        file is a writeable file
        param_list is a list of parameter strings.
        For best results use the parameter_lists supplied by this class.
        """
        if filename == '':
            gpf_ptr = sys.stdout
        else:
            gpf_ptr = open(filename, 'w')

        for p in param_list:
            if p == 'custom_parameter_file':
                # old_custom_parameter_file_value = self['custom_parameter_file']['value']
                # if old_parameter_file_value=='AD4_parameters.dat':
                #    self['parameter_file']['value'] = 'AD4.1_bound.dat'
                if self['custom_parameter_file']['value']:
                    old_parameter_file_value = self['parameter_file']['value']
                    gpf_ptr.write(self.make_param_string('parameter_file'))
                    self['parameter_file']['value'] = old_parameter_file_value
            elif p == 'map':
                # maps are a special case
                for s in self['ligand_types']['value'].split():
                    gpf_ptr.write(self.make_map_string(p, s))
            # all the other parameters handle themselves
            elif p == 'gridcenter' and self['gridcenterAuto']['value'] == 1:
                # if gridcenterAuto is true, reset p to 'auto' and write it
                self['gridcenter']['value'] = 'auto'
                gpf_ptr.write(self.make_param_string(p))
            elif p == 'dsolvmap':
                outstring = "dsolvmap %s              # desolvation potential map\n" % self['dsolvmap']['value']
                gpf_ptr.write(outstring)
            elif p == 'dielectric4':
                # now dielectric value can only be set in parameter file
                # val = self['dielectric4']['value']
                outstring = 'dielectric -0.1465                   # <0, AD4 distance-dep.diel;>0, constant\n'
                gpf_ptr.write(outstring)
            elif p == 'covalentmap' and len(self['covalent_coords']['value']):
                gpf_ptr.write(self.make_covalentmap_string())
            else:
                gpf_ptr.write(self.make_param_string(p))
        if gpf_ptr != sys.stdout:
            gpf_ptr.close()
            self.gpf_filename = filename
            self.gpf_written_filename = filename


class GridParameterFileMaker:
    """Accept a <ligand>.pdbq , <receptor>.pdbqs, reference.gpf and create
    <receptor>.gpf
    sets gridcenter to center of bounding box
    sets npts according to bounding box
    """

    def __init__(self, verbose=None, size_box_to_include_ligand=True):
        self.verbose = verbose
        self.gpo = GridParameters()
        self.size_box_to_include_ligand = size_box_to_include_ligand

    def read_reference(self, reference_filename):
        if self.verbose: print("reading ", reference_filename)
        self.gpo.read(reference_filename)

    def set_ligand(self, ligand_filename):
        self.ligand_filename = os.path.basename(ligand_filename)
        if self.verbose:
            print("set ligand_filename to", self.ligand_filename)
        self.gpo.set_ligand(ligand_filename)
        # expect a filename like ind.out.pdbq: get 'ind' from it
        self.ligand_stem = self.ligand_filename.split('.')[0]
        if self.verbose: print("set ligand_stem to", self.ligand_stem)
        self.ligand = Read(ligand_filename)[0]
        # IS THIS USEFUL???
        self.gpo.ligand = self.ligand
        if self.verbose: print("read ", self.ligand.name)
        # set gpo:
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
        self.gpo['types']['value'] = self.ligand.types
        if self.verbose: print("set types to ", self.gpo['types']['value'])
        # gridcenter
        self.ligand.center = self.ligand.getCenter()
        if self.size_box_to_include_ligand:
            self.getSideLengths(self.ligand)  # sets ligand.center
        cen = self.ligand.center
        self.gpo['gridcenter']['value'] = [round(cen[0], 4), round(cen[1], 4), \
                                           round(cen[2], 4)]
        self.gpo['gridcenterAuto']['value'] = 0
        if self.verbose: print("set gridcenter to ", self.gpo['gridcenter']['value'])
        # only make the box bigger from npts, do not make it smaller
        for ix, val in enumerate(self.gpo['npts']['value']):
            if hasattr(self.ligand, 'npts'):
                npts = self.ligand.npts
                if npts[ix] > val:
                    if self.verbose: print("increasing ", ix, " grid dimension to ", val)
                    self.gpo['npts']['value'][ix] = npts[ix]
        # if self.verbose: print "set npts to ", self.gpo['npts']['value']

    def getSideLengths(self, mol):
        c = mol.allAtoms.coords
        maxo = numpy.maximum.reduce(c)
        mino = numpy.minimum.reduce(c)
        sideLengths = maxo - mino
        mol.npts = list(map(int, list(map(ceil, sideLengths / (self.gpo['spacing']['value'])))))
        for ix, npts in enumerate(mol.npts):
            if npts > 126:
                mol.npts[ix] = 126
        # FIX THIS:
        # use this center instead of mol.getCenter which returns averaged
        # coords:
        # this should make sure the ligand fits inside the box
        # mino+(maxo-mino)/2.0
        mol.center = mino + (maxo - mino) / 2.0

    def set_receptor(self, receptor_filename, gpf_filename=None):
        self.receptor_filename = os.path.basename(receptor_filename)
        self.receptor_stem = self.receptor_filename.split('.')[0]
        self.gpo.set_receptor(receptor_filename)
        # FIX THIS
        # self.gpo['mset']['value'] = self.receptor.types
        self.gpo['types']['value'] = self.ligand.types

    def set_grid_parameters(self, **kw):
        """Any grid parameters should be set here
        """
        # like this:
        # should it be **kw
        # kw = {'spacing':1.0, 'mset':'CNOSHXM'}
        # self.mv.gpo['parm']['value'] = <new value>
        # EXCEPT for 'npts' for which value must be 60,60,60
        for parm, newvalue in list(kw.items()):
            self.gpo[parm]['value'] = newvalue
            if parm == 'npts':
                self.gpo['npts']['value'] = list(map(int, newvalue.split(',')))

    def write_gpf(self, gpf_filename=None,
                  parm_list=grid_parameter_list):
        if not gpf_filename:
            gpf_filename = self.receptor_stem + ".gpf"
        # now that we have a filename...
        if self.verbose:
            print("writing ", gpf_filename)
        self.gpo.write(gpf_filename, parm_list)


class GridParameter4FileMaker:
    """Accept a <ligand>.pdbqt, <receptor>.pdbqt, reference4.gpf and create
    <receptor>4.gpf with help of its "gpo" an instance of a GridParameters
    sets gridcenter to center of bounding box
    sets npts according to bounding box
    """

    def __init__(self, verbose=None, size_box_to_include_ligand=True):
        self.verbose = verbose
        self.gpo = GridParameters()
        self.size_box_to_include_ligand = size_box_to_include_ligand

    def read_reference(self, reference_filename):
        if self.verbose: print("reading ", reference_filename)
        self.gpo.read4(reference_filename)

    def set_types_from_directory(self, directory):
        if self.verbose:
            print("reading directory ", directory)
        filelist = glob.glob(directory + "/*.pdb*")
        if self.verbose:
            print("len(filelist)=", len(filelist))
        ad4_typer = AutoDock4_AtomTyper()
        type_dict = {}
        all_types = ""
        for f in filelist:
            ftype = os.path.splitext(f)[-1]
            if ftype != ".pdbqt":
                print("skipping ", f, " not in PDBQT format!")
                continue
            m = Read(f)[0]
            m_types = ""
            m_types = " ".join(list(set(m.allAtoms.autodock_element)))
            self.getSideLengths(m)  # sets ligand.center
            npts = m.npts
            # only make the box bigger, do NOT make it smaller
            for ix, val in enumerate(self.gpo['npts']['value']):
                if npts[ix] > val:
                    self.gpo['npts']['value'][ix] = npts[ix]
                    if self.verbose:
                        print(m.name, " increased grid dimension ", ix, " to ", npts[ix])
            all_types = all_types + m_types
            if self.verbose:
                print("added ", m_types, " atom types in directory ", directory)
        print("end: all_types = ", all_types)
        self.gpo['ligand_types']['value'] = all_types
        if self.verbose:
            print("all ligand_types for ", directory, "= ", self.gpo['ligand_types']['value'])

    def set_ligand(self, ligand_filename, center_on_ligand=False):
        ftype = os.path.splitext(ligand_filename)[-1]
        if ftype != ".pdbqt":
            print("set_ligand:only ligands in 'pdbqt' files are valid.  ", ftype, " files are not supported!")
            return "ERROR"
        self.ligand = Read(ligand_filename)[0]
        if self.ligand == None:
            print('ERROR reading: ', ligand_filename)
            return
        if self.verbose:
            print("read ", self.ligand.name)
        ligand_types = self.getTypes(self.ligand)
        self.gpo.set_ligand4(ligand_filename, types=ligand_types)
        # this sets ligand_types, gpo.ligand_stem and gpo.ligand_filename
        if self.verbose:
            print("set gpo.ligand_stem to", self.gpo.ligand_stem)
            print("set gpo.ligand_filename to", self.gpo.ligand_filename)
            print("set gpo.ligand_types to", self.gpo['ligand_types']['value'].__class__)
        # need to get npts
        if self.size_box_to_include_ligand:
            self.getSideLengths(self.ligand)  # sets ligand.center
        # gridcenter IS NOT SET BY THIS!!!
        if center_on_ligand:
            # cen = self.ligand.getCenter()
            self.getSideLengths(self.ligand)
            cen = self.ligand.center  # set by call to getSideLengths NOT self.ligand.getCenter
            self.gpo['gridcenter']['value'] = [round(cen[0], 4), round(cen[1], 4), \
                                               round(cen[2], 4)]
            self.gpo['gridcenterAuto']['value'] = 0
            if self.verbose: print("set gridcenter to ", self.gpo['gridcenter']['value'])
        # only make the box bigger, do NOT make it smaller
        for ix, val in enumerate(self.gpo['npts']['value']):
            # npts
            if hasattr(self.ligand, 'npts'):
                npts = self.ligand.npts
                if npts[ix] > val:
                    self.gpo['npts']['value'][ix] = npts[ix]
        if self.verbose: print("set npts to ", self.gpo['npts']['value'])

    def getTypes(self, molecule):
        mol_types = ""
        mol_types = " ".join(list(set(molecule.allAtoms.autodock_element)))
        if self.verbose:
            print("end of getTypes: mol_types=", mol_types, ' class=', mol_types.__class__)
        return mol_types

    def getSideLengths(self, mol):
        c = mol.allAtoms.coords
        maxo = numpy.maximum.reduce(c)
        mino = numpy.minimum.reduce(c)
        sideLengths = maxo - mino
        mol.npts = list(map(int, list(map(ceil, sideLengths / (self.gpo['spacing']['value'])))))
        for ix, npts in enumerate(mol.npts):
            if npts > 126:
                mol.npts[ix] = 126
        # FIX THIS:
        # use this center instead of mol.getCenter which returns averaged
        # coords:
        # this should make sure the ligand fits inside the box
        # mino+(maxo-mino)/2.0
        mol.center = mino + (maxo - mino) / 2.0

    def set_receptor(self, receptor_filename, gpf_filename=None):
        ftype = os.path.splitext(receptor_filename)[-1]
        if ftype != ".pdbqt":
            print("set_receptor:only pdbqt files valid.  ", ftype, " files are not supported.")
            return "ERROR:"
        self.receptor = Read(receptor_filename)[0]
        receptor_filename = os.path.basename(receptor_filename)
        if self.receptor == None:
            print('ERROR reading: ', receptor_filename)
            return
        if self.verbose: print("set_receptor filename to ", receptor_filename)
        receptor_types = self.getTypes(self.receptor)
        self.gpo.set_receptor4(receptor_filename, types=receptor_types)
        self.receptor_filename = os.path.basename(receptor_filename)
        if hasattr(self, 'receptor'):
            self.receptor_stem = self.receptor.name
        else:
            self.receptor_stem = os.path.splitext(self.receptor_filename)[0]
        # all of this is handled by set_receptor4
        # self.gpo['gridfld']['value'] = self.receptor_stem + '.maps.fld'
        # self.gpo['elecmap']['value'] = self.receptor_stem + '.e.map'
        # self.gpo['dsolvmap']['value'] = self.receptor_stem + '.d.map'
        # this sets gpo.receptor_types, gpo.receptor_stem and gpo.receptor_filename

    def set_grid_parameters(self, **kw):
        """Any grid parameters should be set here
        """
        # like this:
        # should it be **kw
        # kw = {'spacing':1.0, 'receptor_types':'C A NA OA N SA HD MG'}
        # self.mv.gpo['parm']['value'] = <new value>
        # EXCEPT for 'npts' for which value must be 60,60,60
        for parm, newvalue in list(kw.items()):
            if self.verbose:
                print("parm=", parm)
                print("newvalue=", newvalue)
            if parm == 'gridcenter':
                self.gpo['gridcenterAuto']['value'] = newvalue == 'auto'
            self.gpo[parm]['value'] = newvalue
            if parm == 'npts':
                self.gpo['npts']['value'] = list(map(int, newvalue.split(',')))
            if parm == 'ligand_types':
                if newvalue.find(',') > -1:
                    newvalue = newvalue.replace(',', ' ')
                print("setting ligand_types: newvalue=", newvalue)
                self.gpo[parm]['value'] = newvalue

    def write_gpf(self, gpf_filename=None,
                  parm_list=grid_parameter_list4):
        if not gpf_filename:
            gpf_filename = self.receptor_stem + ".gpf"
        # now that we have a filename...
        if self.verbose:
            print("writing ", gpf_filename)
            for item in parm_list:
                print(item, end=' ')
            print()
        self.gpo.write4(gpf_filename, parm_list)
