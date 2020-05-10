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
# Last modified on Tue Apr 23 09:20:22 PDT 2002 by lindy
#
# $Header: /opt/cvs/python/packages/share1.5/mglutil/math/statetocoords.py,v 1.11.10.1 2016/02/11 23:15:05 annao Exp $
#

"""statetocoords.py - state to coordinates

The StateToCoords class inherits from Kinematics and Ncoords.
The StateToCoords class handles transformations that apply to
the rootNode of the torTree, changing the coordinates in world
space. The Kinematics class handles those transformations that
apply to the internal nodes to the torTree (ie. torsions) changing
the coordinates in the molecules local coordinate system.
"""
import numpy

from mglutil.math.kinematics import Kinematics
from mglutil.math.transformation import Transformation


class StateToCoords(Kinematics):
    def __init__(self, mol, origin, confIndex):
        Kinematics.__init__(self, mol.allAtoms.coords, mol.torTree, tolist=1)
        # this stoc object will leave always deposite it's coords
        # in the given confIndex slot.
        self.confIndex = confIndex

        mol.allAtoms.setConformation(confIndex)

        def __prepareNode(node, allAtoms, o):
            """Supply each node with atomSet, coords, and atomRange,
            Pre-compute and save the torsionUnitVector, 
            Transform the coords to their local space by subtracting
            the origin
            """
            atomSet = []
            coords = []
            for i in node.atomList:
                atom = allAtoms[i]
                atomSet.append(atom)
                # start with the original coordinates
                c = atom.coords
                # subract the origin
                coords.append((c[0] - o[0], c[1] - o[1], c[2] - o[2], 1.0))
            node.atomSet = atomSet
            node.coords = coords
            node.atomRange = list(range(len(atomSet)))
            if node.bond[0] != None:  # skip the root node
                node.a = allAtoms[node.bond[0]]
                node.b = allAtoms[node.bond[1]]

        # add atomSets to each node
        root = mol.torTree.rootNode
        root.pre_traverse(__prepareNode, root, mol.allAtoms, origin)

    def applyState(self, state):
        """
        """
        q = state.quaternion
        t = numpy.array(state.translation)
        o = numpy.array(state.origin)

        # construct rootNode transformation matrix
        # mtx = Transformation(t+o, q).getMatrix(transpose=1)
        # Corrected by AG 08/28/2008
        mtx = Transformation(t, q).getMatrix().transpose()

        # apply the torsions
        self.applyAngList(state.torsions, mtx)

    def applyStateOld(self, state):
        """
        """
        q = state.quaternion
        t = numpy.array(state.translation)
        o = numpy.array(state.origin)

        # center the coordinates
        ##          self.resultCoords = (self.resultCoords -
        ##                               numpy.array([o[0], o[1], o[2], 0.0]))

        # center the coordinates (node-by-node)
        def __center(node, o): node.coords = node.coords - o

        root = self.torTree.rootNode
        root.pre_traverse(__center, root, numpy.array([o[0], o[1], o[2], 0.0]))

        # construct rootNode transformation matrix
        mtx = Transformation(t + o, q).getMatrix(transpose=1)

        # apply the torsions
        coords = self.applyAngList(state.torsions, mtx)

        # must "reset" each nodes coords
        def __uncenter(node, o): node.coords = node.coords + o

        root.pre_traverse(__uncenter, root, numpy.array([o[0], o[1], o[2], 0.0]))

        return coords

    def applyOrientation(self, q=(0., 0., 0., 0.), t=(0., 0., 0.), o=(0., 0., 0.)):
        """origin specifies where the local origin is in world coordinates
        (i.e., where is this object's origin in the world)
        """
        # center the coordinates
        self.resultCoords = (self.resultCoords -
                             numpy.array([o[0], o[1], o[2], 0.0]))
        sum = numpy.array(t) + numpy.array(o)
        self.resultCoords = Transformation(sum, q).apply(self.resultCoords)
        return self.getResultCoords()

    def applyQuaternion(self, q, o=(0., 0., 0.)):
        """Apply the given quaterion.
        """
        # center the coordinates
        self.resultCoords = (self.resultCoords -
                             numpy.array([o[0], o[1], o[2], 0.0]))
        self.resultCoords = Transformation(o, q).apply(self.resultCoords)
        return self.getResultCoords()

    def applyTranslation(self, t=(0., 0., 0.)):
        """Translate by (x, y, z)
        """
        translation = numpy.array([t[0], t[1], t[2], 0.0])
        self.resultCoords = self.resultCoords + translation
        return self.getResultCoords()
