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
#  Modification date: 10/5/20 18:51                                                                #
#                                                                                                  #
# ##################################################################################################

#############################################################################
#
# Author: Ruth HUEY, Michel F. SANNER
#
# Copyright: M. Sanner TSRI 2002
#
#############################################################################


# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/EpdbParser.py,v 1.7 2010/02/03 00:06:39 rhuey Exp $
#
# $Id: EpdbParser.py,v 1.7 2010/02/03 00:06:39 rhuey Exp $
#
#
#
#
#
#
#

"""
This Object parses the result of an AutoDock command mode epdb operation. It builds a dictionary. 

"""
import os

from AutoDockTools.ResultParser import ResultParser


class EpdbParser(ResultParser):
    """ reads log from a AutoDock docking and return structured data"""

    keywords = ResultParser.keywords + [
        'coords',
        'vdw_energies',
        'estat_energies',
        'inhib_constant',
        'intermol_energy',  # (1) 
        'internal_energy',  # (2) NB: 1+2->final docked energy
        'torsional_energy',  # (3) NB: 1+3->free energy of binding

    ]

    def __init__(self, dlgFile=None):
        """selected dlgFile,ok sets which docked conformations to show"""
        ResultParser.__init__(self)
        self.filename = dlgFile
        self.version = 4.2
        self.ntors = 0
        self.found_ntors = 0
        if dlgFile:
            self.filename = os.path.basename(dlgFile)
            self.parse(dlgFile)

    def parse(self, filename):
        """
        uses key 'NOW IN COMMAND MODE.' to start matching:
            'ATOM',
        next uses 'Intermolecular Energy Analysis' to start
            capturing individual energy breakdowns
        finally captures '^epdb: USER' lines
        after parsing: 
        """
        self.filename = filename
        # reset
        dlgptr = open(filename, 'r')
        allLines = dlgptr.readlines()
        self.clusterRecord = None
        self._parse(allLines)

    def _parse(self, allLines):
        if not len(allLines):
            return 'ERROR'
        for item in ['ligLines', 'dpfLines', 'energyLines', 'epdbLines', \
                     'atTypes', 'vdw_energies', 'estat_energies', 'clist', 'clusterlines', \
                     'histogramlines', 'modelList', 'total_energies']:
            setattr(self, item, [])
        lineLen = len(allLines)
        # print 'lineLen=', lineLen
        atmCtr = self.atmCtr = 0
        ligLines = self.ligLines
        dpfLines = self.dpfLines
        energyLines = self.energyLines
        epdbLines = self.epdbLines
        for i in range(lineLen):
            l = allLines[i]
            # while l.find( 'NOW IN COMMAND MODE')<0 and i<lineLen-1:
            #    continue
            if l.find('INPUT-PDBQ: ATOM') == 0:
                ligLines.append(l[12:])
                atmCtr = atmCtr + 1
            elif l.find('INPUT-PDBQ: HETA') == 0:
                ligLines.append(l[12:])
                atmCtr = atmCtr + 1
            elif not self.found_ntors and l.find('active torsions:') > -1:
                self.found_ntors = 1
                self.ntors = int(l.split()[2])
            elif l.find('INPUT-LIGAND-PDBQT: HETA') == 0:
                ligLines.append(l[12:])
                atmCtr = atmCtr + 1
            elif l.find('INPUT-LIGAND-PDBQT: ATOM') == 0:
                ligLines.append(l[12:])
                atmCtr = atmCtr + 1
            elif l.find('DPF>') == 0:
                dpfLines.append(l[5:-1])
            elif l.find('Intermolecular Energy Analysis') > 0:
                # elif l.find( 'Intermolecular Energy Analysis')>-1:
                # print 'found Intermolecular Energy Analysis at line ', i
                break
        self.atmCtr = atmCtr
        i = i + 5
        for x in range(i, lineLen):
            # process energy lines
            i = i + 1
            l = allLines[i]
            if l.find('Total') == 0:
                self.energyLines = self.energyLines[:-1]
                break
            self.energyLines.append(l)
        for l in self.energyLines:
            ll = l.split()
            self.atTypes.append(int(ll[0]))
            self.vdw_energies.append(float(ll[2]))
            self.estat_energies.append(float(ll[3]))
            self.total_energies.append(float(ll[2]) + float(ll[3]))
        while l.find('epdb') < 0:
            # skip some stuff
            l = allLines[i]
            ll = l.split()
            i = i + 1
        for x in range(i - 1, lineLen):
            # process epdb lines
            l = allLines[x]
            ll = l.split()
            if l.find('Estimated Free Energy') > 0:
                self.estFreeEnergy = float(ll[8])
            elif l.find('Final Docked Energy') > 0:
                self.finalDockedEnergy = float(ll[6])
            elif l.find('Final Intermolecular Energy') > 0:
                self.finalIntermolEnergy = float(ll[7])
            elif l.find('Final Internal Energy') > 0:
                self.finalInternalEnergy = float(ll[9])
            elif l.find('Final Total Internal Energy') > 0:
                self.finalTotalInternalEnergy = float(ll[8])
            elif l.find('Torsional Free Energy') > 0:
                self.torsionalFreeEnergy = float(ll[7])
