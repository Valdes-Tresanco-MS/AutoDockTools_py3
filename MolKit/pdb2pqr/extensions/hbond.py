"""
    Hbond extension

    Find all hydrogen bonds as determined by the cutoffs below.
    Uses PDB2PQR to determine donors and acceptors, and displays
    all available bonds to stdout.

    Author:  Todd Dolinsky
"""

__date__ = "17 February 2006"
__author__ = "Todd Dolinsky"

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

from ..src.routines import *
from ..src.utilities import *

ANGLE_CUTOFF = 20.0  # A - D - H(D) angle
DIST_CUTOFF = 3.3  # H(D) to A distance


def usage():
    text = "        --hbond       :  Print a list of hydrogen bonds to\n"
    text += "                         {output-path}.hbond\n"
    return text


def hbond(routines, outroot):
    """
        Print a list of hydrogen bonds.

        Parameters
            routines:  A link to the routines object
            outroot:   The root of the output name
    """
    outname = outroot + ".hbond"
    file = open(outname, "w")

    routines.write("Printing hydrogen bond list...\n")

    # Initialize - set nearby cells, donors/acceptors
    # The cell size adds one for the D-H(D) bond, and rounds up

    cellsize = int(DIST_CUTOFF + 1.0 + 1.0)
    protein = routines.protein
    routines.setDonorsAndAcceptors()
    routines.cells = Cells(cellsize)
    routines.cells.assignCells(protein)

    for donor in protein.getAtoms():

        # Grab the list of donors
        if not donor.hdonor:
            continue
        donorhs = []
        for bond in donor.bonds:
            if bond.isHydrogen():
                donorhs.append(bond)
        if not donorhs:
            continue

        # For each donor, grab all acceptors

        closeatoms = routines.cells.getNearCells(donor)
        for acc in closeatoms:
            if not acc.hacceptor:
                continue
            if donor.residue == acc.residue:
                continue
            for donorh in donorhs:

                # Do distance and angle checks

                dist = distance(donorh.getCoords(), acc.getCoords())
                if dist > DIST_CUTOFF:
                    continue
                angle = getAngle(acc.getCoords(), donor.getCoords(), donorh.getCoords())
                if angle > ANGLE_CUTOFF:
                    continue
                routines.write("Donor: %s %s\tAcceptor: %s %s\tHdist: %.2f\tAngle: %.2f\n" %
                               (donor.residue, donor.name, acc.residue, acc.name, dist, angle))
                file.write("Donor: %s %s\tAcceptor: %s %s\tHdist: %.2f\tAngle: %.2f\n" %
                           (donor.residue, donor.name, acc.residue, acc.name, dist, angle))

    routines.write("\n")
    file.close()
