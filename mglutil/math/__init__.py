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

import numpy


def crossProduct(A, B, normal=True):
    """     Return cross product of two vectors A and B
normal: return normalized vector
"""
    res = [A[1] * B[2] - A[2] * B[1],
           A[2] * B[0] - A[0] * B[2],
           A[0] * B[1] - A[1] * B[0]]
    if normal:
        return norm(res)
    else:
        return res


def norm(A):
    """     Return normalized vector A.
"""
    if type(A) == list:
        A = numpy.array(A, 'f')
        res = A / numpy.sqrt(numpy.dot(A, A))
        return res.tolist()
    elif type(A) == numpy.ndarray:
        return A / numpy.sqrt(numpy.dot(A, A))
    else:
        print("Need a list or numpy array")
        return None


def getCenter(coords):
    """ get center of all the coords """
    coords = numpy.array(coords, 'f')
    return (numpy.sum(coords, 0) / len(coords)).tolist()
