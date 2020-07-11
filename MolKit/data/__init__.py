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
# Author: Michel F. SANNER, Sophie COON
#
# Copyright: M. Sanner TSRI 2000
#
#############################################################################

def Read(filename):
    from MolKit.pdbParser import PdbParser, PdbqParser, PdbqsParser, PQRParser
    from MolKit.mol2Parser import Mol2Parser
    ext = filename.split('.')
    if ext[-1] == 'pdb':
        parser = PdbParser(filename)

    elif ext[-1] == 'pdbq':
        parser = PdbqParser(filename)

    elif ext[-1] == 'pdbqs':
        parser = PdbqsParser(filename)

    elif ext[-1] == 'pqr':
        parser = PQRParser(filename)

    elif ext[-1] == 'mol2':
        parser = Mol2Parser(filename)

    else:
        print("File Format unknown can't parse it")
        return []
    molecules = parser.parse()
    return molecules


def WritePDB(filename, node):
    from MolKit.pdbWriter import PdbWriter
    writer = PdbWriter()
    writer.write(filename, node)
