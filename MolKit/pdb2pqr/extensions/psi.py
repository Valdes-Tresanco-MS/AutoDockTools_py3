"""
    Psi extension

    Print the psi backbone angle for each residue in the structure.
    Psi angle is determined by the coordinates of the N(i), CA(i), C(i), N(i+1)
    atoms.

    Author:  Mike Bradley and Todd Dolinsky
"""

__date__ = "17 February 2006"
__author__ = "Mike Bradley, Todd Dolinsky"

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

from ..src.utilities import *


def usage():
    text = "        --psi         :  Print the per-residue backbone psi\n"
    text += "                         angle to {output-path}.psi\n"
    return text


def psi(routines, outroot):
    """
        Print the list of psi angles

        Parameters
            routines:  A link to the routines object
            outroot:   The root of the output name
    """

    outname = outroot + ".psi"
    file = open(outname, "w")

    routines.write("\nPrinting psi angles for each residue...\n")
    routines.write("Residue     Psi\n")
    routines.write("----------------\n")

    # Initialize some variables

    protein = routines.protein

    for residue in protein.getResidues():
        if residue.hasAtom("N"):
            ncoords = residue.getAtom("N").getCoords()
        else:
            continue

        if residue.hasAtom("CA"):
            cacoords = residue.getAtom("CA").getCoords()
        else:
            continue

        if residue.hasAtom("C"):
            ccoords = residue.getAtom("C").getCoords()
        else:
            continue

        try:
            if residue.peptideN is not None:
                pepcoords = residue.peptideN.getCoords()
            else:
                continue
        except AttributeError:  # Non amino acids
            continue

        psi = getDihedral(ncoords, cacoords, ccoords, pepcoords)
        routines.write("%s\t%.4f\n" % (residue, psi))
        file.write("%s\t%.4f\n" % (residue, psi))

    routines.write("\n")
    file.close()
