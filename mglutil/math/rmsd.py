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
#  Modification date: 10/5/20 18:04                                                                #
#                                                                                                  #
# ##################################################################################################

#############################################################################
#
# Author: Sophie I. COON, William LINDSTROM, Michel F. SANNER
#
# Copyright: M. Sanner TSRI 2000
#
#############################################################################
#


import math

import numpy

from mglutil.math.munkres import Munkres


def getAtomIndicesPerType(atoms, typename='autodock_element'):
    d1 = {}
    for i, a in enumerate(atoms):
        try:
            d1[getattr(a, typename)].append(i)
        except KeyError:
            d1[getattr(a, typename)] = [i]
    return d1


class HungarianMatchingRMSD:
    """
    class to compute RMSD between 2 poses of the same molecule with pairing
    calculated using the Hungarian matching algorithm.

    typeIndicesRef are dictionary of where the key is an atom type and the value
    is a 0-based list of indices for atoms of that type in the list of atoms provided
    to the constructor (i.e. the reference atoms).

    the
    """

    def __init__(self, atoms, typeIndicesRef, typeIndicesMoving,
                 ignoreTypes=['HD']):
        # create rmsd calculator
        self.sortedRefAts = atoms
        self.typeIndicesRef = typeIndicesRef
        self.typeIndicesMoving = typeIndicesMoving
        self.atypes = list(typeIndicesRef.keys())
        for typeName in ignoreTypes:
            if typeName in self.atypes:
                self.atypes.remove(typeName)
        self.matching = None  # will hold a list of computed pairs after matching

    def setRefCoords(self, coords):
        """
        set the reference atoms
        """
        self.sortedRefAts.updateCoords(coords)

    def computeRMSD(self, coords):
        """
        compute RMSD with reference atoms. coords are assumed to be in the same order
        as self.sortedRefAts
        """

        # use the Hungarian matching algorithm to pair up atoms
        # of the same time while minimizing the sum of the distances squared
        matching = []
        total = 0  # sum up square of distances

        # loop over atoms types
        for atype in self.atypes:
            inds1 = self.typeIndicesRef[atype]
            inds2 = self.typeIndicesMoving.get(atype, None)
            if inds2 is None:
                continue

            if len(inds1) == 1 and len(inds2) == 1:  # only one atom of this type, matching is obvious
                matching.append((inds1[0], inds2[0]))
                x1, y1, z1 = self.sortedRefAts[inds1[0]].coords
                x2, y2, z2 = coords[inds2[0]]
                total += (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2)

            elif len(inds1) == 2 and len(inds2) == 2:  # only two atoms of this type, matching is obvious
                x1i, y1i, z1i = self.sortedRefAts[inds1[0]].coords
                x1j, y1j, z1j = self.sortedRefAts[inds1[1]].coords
                x2i, y2i, z2i = coords[inds2[0]]
                x2j, y2j, z2j = coords[inds2[1]]
                # compute dist(i,i) + dist(j,j)
                sum1 = ((x1i - x2i) * (x1i - x2i) + (y1i - y2i) * (y1i - y2i) + (z1i - z2i) * (z1i - z2i) +
                        (x1j - x2j) * (x1j - x2j) + (y1j - y2j) * (y1j - y2j) + (z1j - z2j) * (z1j - z2j))
                sum2 = ((x1i - x2j) * (x1i - x2j) + (y1i - y2j) * (y1i - y2j) + (z1i - z2j) * (z1i - z2j) +
                        (x1j - x2i) * (x1j - x2i) + (y1j - y2i) * (y1j - y2i) + (z1j - z2i) * (z1j - z2i))
                if sum1 < sum2:
                    matching.append((inds1[0], inds1[0]))
                    matching.append((inds1[1], inds1[1]))
                    total += sum1
                else:
                    matching.append((inds1[0], inds1[1]))
                    matching.append((inds1[1], inds1[0]))
                    total += sum2

            else:  # use Hungarian matching algorithm to find best assignment
                l1 = len(inds1)
                l2 = len(inds2)
                # print atype, inds
                matrix = numpy.zeros((l1, l2), 'f')
                for i, n1 in enumerate(inds1):
                    x, y, z = self.sortedRefAts[n1].coords
                    for j, n2 in enumerate(inds2):
                        x1, y1, z1 = coords[n2]
                        matrix[i][j] = (x - x1) * (x - x1) + (y - y1) * (y - y1) + (z - z1) * (z - z1)

                # compute best assignment
                m = Munkres()
                indexes = m.compute(matrix.tolist())
                ltotal = 0.
                for row, column in indexes:
                    value = matrix[row][column]
                    ltotal += value
                    matching.append((inds1[row], inds2[column]))
                # print matrix
                # print 'total cost: %s %f %d' % (atype, ltotal, len(inds))
                total += ltotal

        self.matching = matching

        from math import sqrt
        rmsd = sqrt(total / len(matching))


        return rmsd


class RMSDCalculator:
    """
    This class implements method to compute RMSD and distance vector
    between two given lists of coordinates.
    """

    def __init__(self, refCoords=None):
        self.refCoords = refCoords

    def setRefCoords(self, refCoords):
        self.refCoords = refCoords

    def computeRMSD(self, listCoords):
        """rmsd <- computRMSD(listCoords)
        rmsd returns the overall root mean square distance (rmsd) and
        also sets self.distVect as the vector of distances between each
        pair of points.
        """
        if self.refCoords is None:
            raise ValueError("no reference coordinates set")
        if len(self.refCoords) != len(listCoords):
            raise ValueError("input vector length mismatch")

        deltaVect = numpy.array(self.refCoords) - numpy.array(listCoords)
        distSquaredVect = numpy.sum(numpy.transpose(deltaVect * deltaVect))
        self.distVect = numpy.sqrt(distSquaredVect)
        self.rmsd = math.sqrt(numpy.sum(distSquaredVect) / len(self.refCoords))
        return self.rmsd
