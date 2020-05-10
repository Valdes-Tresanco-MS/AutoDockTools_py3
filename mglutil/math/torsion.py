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

# taken from Pmv/measureCommands.py
import numpy


def torsion(x1, x2, x3, x4):
    """
    Compute the torsion angle between x1, x2, x3, x4.
    All coordinates are cartesian; result is in degrees.
    Raises a ValueError if angle is not defined.
    """
    from math import sqrt, acos

    tang = 0.0
    x1 = numpy.array(x1, 'f')
    x2 = numpy.array(x2, 'f')
    x3 = numpy.array(x3, 'f')
    x4 = numpy.array(x4, 'f')

    assert x1.shape == (3,)
    assert x2.shape == (3,)
    assert x3.shape == (3,)
    assert x4.shape == (3,)

    a = x1 - x2
    b = x3 - x2
    c = vvmult(a, b)

    a = x2 - x3
    b = x4 - x3
    d = vvmult(a, b)

    dd = sqrt(numpy.sum(c * c))
    de = sqrt(numpy.sum(d * d))

    if dd < 0.001 or de < 0.001:
        raise ValueError('Torsion angle undefined, degenerate points')

    vv = numpy.dot(c, d) / (dd * de);
    if vv < 1.0:
        tang = vv
    else:
        tang = 1.0
    if tang < -1.0: tang = -1.0
    tang = acos(tang)
    tang = tang * 57.296

    b = vvmult(c, d)
    if numpy.dot(a, b) > 0.0: tang = -tang
    return tang


def vvmult(a, b):
    """
    Compute a vector product for 3D vectors
    """
    res = numpy.zeros(3, 'f')
    res[0] = a[1] * b[2] - a[2] * b[1]
    res[1] = a[2] * b[0] - a[0] * b[2]
    res[2] = a[0] * b[1] - a[1] * b[0]
    return res
