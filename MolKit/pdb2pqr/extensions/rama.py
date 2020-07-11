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

"""
    Ramachandran extension

    Print both the phi and psi angles to standard out.  See the individual
    functions for more info.

    Author:  Mike Bradley and Todd Dolinsky
"""

__date__ = "17 February 2006"
__author__ = "Mike Bradley, Todd Dolinsky"

from ..src.utilities import *


def usage():
    text = "        --rama        :  Print the per-residue phi and psi\n"
    text += "                         angles to {output-path}.rama for\n"
    text += "                         Ramachandran plots\n"
    return text


def rama(routines, outroot):
    """
        Print the list of phi and psi angles for use in a Ramachandran plot.

        Parameters
            routines:  A link to the routines object
            outroot:   The root of the output name
    """

    outname = outroot + ".rama"
    file = open(outname, "w")

    routines.write("\nPrinting phi and psi angles for each residue...\n")
    routines.write("Residue        Phi          Psi\n")
    routines.write("-------------------------------\n")

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
                pepncoords = residue.peptideN.getCoords()
            else:
                continue

            if residue.peptideC is not None:
                pepccoords = residue.peptideC.getCoords()
            else:
                continue
        except AttributeError:  # Non amino acids
            continue

        phi = getDihedral(pepccoords, ncoords, cacoords, ccoords)
        psi = getDihedral(ncoords, cacoords, ccoords, pepncoords)
        routines.write("%s\t%.4f\t%.4f\n" % (residue, phi, psi))
        file.write("%s\t%.4f\t%.4f\n" % (residue, phi, psi))

    routines.write("\n")
    file.close()
