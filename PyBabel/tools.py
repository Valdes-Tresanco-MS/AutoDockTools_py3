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
#  Modification date: 10/5/20 22:13                                                                #
#                                                                                                  #
# ##################################################################################################

#############################################################################
#
# Author: Michel F. SANNER
# Reimplemented from Babel v1.6 from Pat Walters and Math Stahl
#
# Copyright: M. Sanner TSRI 2000
#
#############################################################################


def read_element_table(filename):
    """void <- read_element_table(filename)
    populates the elementsTable dictionary from the a given file.
    the file provides:
    line number, element string, cov_rad, bond_ord_rad, vdw_rad, bs_rad,
    max_bonds, red, green, blue
    """
    f = open(filename)
    lines = f.readlines()
    f.close()
    elemTable = {}
    for i in range(len(lines)):
        dd = lines[i].split()
        elemTable[dd[1]] = {'num': i,
                            'cov_rad': float(dd[2]),
                            'bond_ord_rad': float(dd[3]),
                            'vdw_rad': float(dd[4]),
                            'bs_rad': float(dd[5]),
                            'max_bonds': int(dd[6]),
                            'rgb': (float(dd[7]), float(dd[8]), float(dd[9]))
                            }
    return elemTable


def writeElementTableAsPythonCode(elemTab, inFileName, outFileName):
    """write elemTable as a python dictionary that can be imported"""

    f = open(outFileName, 'w')
    f.write("# File generated from %s\n#\n" % inFileName)
    f.write("babel_elements = {\n")
    for k, v in list(elemTab.items()):
        f.write("  '%s': %s, \n" % (k, str(v)))
    f.write('}\n#END\n')
    f.close()


def read_types_table(filename):
    f = open(filename)
    typestab = {}
    nrow, ncol = list(map(int, f.readline().split()))
    typeFormats = f.readline().split()
    for t in typeFormats:
        typestab[t] = []
    for i in range(nrow - 1):
        typeNames = f.readline().split()
        for j in range(ncol):
            typestab[typeFormats[j]].append(typeNames[j])
    f.close()
    return typestab


def writeTypesTableAsPythonCode(typestab, inFileName, outFileName):
    """write typestab as a python dictionary that can be imported"""

    f = open(outFileName, 'w')
    f.write("# File generated from %s\n#\n" % inFileName)
    f.write("babel_types = {\n")
    for k, v in list(typestab.items()):
        f.write("  '%s': %s, \n" % (k, str(v)))
    f.write('}\n#END\n')
    f.close()


if __name__ == '__main__':
    # write tables
    et = read_element_table('element.lis')
    writeElementTableAsPythonCode(et, 'element.lis', 'babelElements.py')

    tt = read_types_table('types.lis')
    writeTypesTableAsPythonCode(tt, 'types.lis', 'babelAtomTypes.py')
