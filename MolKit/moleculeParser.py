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
#  Modification date: 4/7/20 5:48                                                                  #
#                                                                                                  #
# ##################################################################################################

#############################################################################
#
# Author: Michel F. SANNER
#
# Copyright: M. Sanner TSRI 2000
#
#############################################################################

import warnings

# dict used to guess parser based on file extension
parserToExt = {'PDB': '.pdb', 'PDBQ': '.pdbq',
               'PDBQS': '.pdbqs', 'PDBQT': '.pdbqt',
               'MOL2': '.mol2',
               'PQR': '.pqr',
               'GRO': '.gro',
               'F2D': '.f2d',
               'SDF': '.sdf',
               'CIF': '.cif',
               }


class MoleculeParser:
    def __init__(self, filename=None, allLines=None):
        """Supply the filename for reading by readFile, or
        supply the lines directly via allLines
        """
        self.filename = str(filename)
        self.allLines = allLines  # stores all lines from file

    def readFile(self):
        f = open(self.filename)
        self.allLines = f.readlines()
        if len(self.allLines) == 1:
            # this file probably has \r instead or \n
            self.allLines = self.allLines[0].split('\r')
            warnings.warn('Only 1 line read from PDB file, splitting on \r')
        f.close()
        self.allLines = list(filter(lambda x, s=str.strip: len(s(x)),
                                    self.allLines))
