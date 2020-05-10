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

#
# Last modified on Tue Sep  4 16:32:29 PDT 2001 by lindy
#

"""ncoords.py - Numeric coordinates

This class is intented to be the base class of a number
of classes which transform and generally operate on lists
of homogeneous coordinates.
"""

import numpy


class Ncoords:
    def __init__(self, refCoords, tolist=1):
        """refCoords is an nx3 list of n points
        
        resultCoords is set up and maintained as homogeneous coords
        if tolist then return the result coords as a python list
        """
        try:
            self.refCoords = numpy.array(numpy.concatenate(
                (refCoords, numpy.ones((len(refCoords), 1), 'f')), 1))
        except TypeError:
            raise ValueError("invalid input array")

        self.resultCoords = self.refCoords
        self.tolist = tolist

    def reset(self):
        self.resultCoords = self.refCoords

    def getResultCoords(self):
        """Return the list of result coordinates

        if tolist is set, return an nx3 Python ListType.
        if tolist is not set, return an nx4 numpy array.
        """
        if self.tolist:
            return self.resultCoords[:, :3].tolist()
        else:
            return self.resultCoords
