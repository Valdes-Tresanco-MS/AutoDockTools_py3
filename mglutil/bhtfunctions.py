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
#  Modification date: 10/5/20 22:22                                                                #
#                                                                                                  #
# ##################################################################################################

import numpy
from bhtree import bhtreelib


# ClosePointsDist2: result and dist are empty arrays large enough to contain
# all the points expected to be found. To be safe, they should be the same
# size as the list of coordinates. This function then puts in the results
# array the indices of the close points in the list (as supplied in the array
# ids): the dist array contains the corresponding distances.

def findNearestAtoms(mol, vertices, **kw):
    """None <- color(mol,vertices2,**kw)
    mol:           reference molecule
    vertices:      list of lists(coordinates): the first three items in each list
                   must be coordinates x,y,z of a point.
                   
    atomIndices is the index of the nearest atom to the vertex, such that
    mol.allAtoms[atomIndices[x]] is the nearest atom to vertices[x]
    vertexIndices is the list of nearest vertices to an atom, such that
    vertexIndices[x] = [vertex1,vertex2,...] are the vertices associated with
    mol.allAtoms[x]
    """

    coords = mol.allAtoms.coords
    if not hasattr(mol, 'bhtree'):
        print("Building bhtree for ", mol)
        ids = numpy.arange(len(coords)).astype('i')
        bhtree = bhtreelib.TBHTree(coords, ids, 10, 10, 9999.0)
        mol.bhtree = bhtree

    vertexIndices = {}
    atomIndices = {}
    for x in range(len(coords)):
        vertexIndices[x + 1] = []

    cutoff = 5.
    for x in range(len(vertices)):
        xyz = vertices[x]
        result = numpy.zeros((len(vertices),)).astype('i')
        dist = numpy.zeros((len(vertices),)).astype('f')
        nb2 = mol.bhtree.ClosePointsDist2(tuple(xyz[:3]), cutoff, result, dist)
        while nb2 == 0:
            cutoff = cutoff + 5.
            nb2 = mol.bhtree.ClosePointsDist2(tuple(xyz[:3]), cutoff, result, dist)
        result = result[:nb2]
        dist = dist[:nb2]
        idx = dist.tolist().index(min(dist))
        atnum = result[idx] + 1
        atomIndices[x] = atnum
        vertexIndices[atnum].append(x)

    return atomIndices, vertexIndices
