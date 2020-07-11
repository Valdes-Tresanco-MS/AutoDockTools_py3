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
# Authors: Michel F. SANNER, Ruth Huey
#
# Copyright: M. Sanner TSRI 2005
#
#############################################################################

from .tree import TreeNodeSet


class Sets(dict):
    """
Object used to manage a collection of explicit sets of TreeNodes
"""

    def add(self, name, set, overwrite=True):
        assert isinstance(set, TreeNodeSet)
        assert type(name) in (str,)

        if not overwrite:
            assert name not in list(self.keys())
        self[name] = set

    def remove(self, name):
        # remove a set by name. Silently ignore non existing sets but returns
        # true when a set gets deleted, else returns False
        if name in list(self.keys()):
            del self[name]
            return True
        return False

    def removeByInstance(self, set):
        # remove a set that is specified by a TreeNodeSet.
        # Silently ignore non existing sets but returns
        # true when a set gets deleted, else returns False
        for n, s in list(self.items()):
            if s == set:
                del self[n]
                return True
        return False

    def get(self, stype=None):
        # return a dict of sets optionally restricted to a user specified type
        # if stype is specified it has to be a subclass of TreeNodeSet
        if stype is None:
            return self
        else:  # select the sets of a given type
            assert issubclass(stype, TreeNodeSet)
            result = {}
            for name, set in list(self.items()):
                if isinstance(set, stype):
                    result[name] = set
            return result
