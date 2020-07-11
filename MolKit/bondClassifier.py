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
# Author: Ruth Huey, William M. Lindstrom
#
# Copyright: R. Huey, W. M. Lindstrom TRSI 2003
#
#############################################################################

"""
This module implements a classifier which select bonds based on a 
dictionary of key, bondSelector.
It returns  a dictionary with keys the specified bond types and 
values the bonds which have been classified.
"""

from MolKit.molecule import BondSet


class BondClassifier:
    """ Base class that sorts bonds based on an input dictionary with keys
    and bondSelector values
    """

    def __init__(self, d={}):
        self.dict = d

    def classify(self, bonds=None):
        """ 
        select using each bondselector (the values of the dict); store result
        in resultDict and return the result dict when finished...
        """
        # make sure that classify is called with some bonds
        assert isinstance(bonds, BondSet)
        resultDict = {}
        for k, v in list(self.dict.items()):
            resultDict[k] = v.select(bonds)
        return resultDict
