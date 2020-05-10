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

##  Copyright 1997-1999 by Konrad Hinsen, except as noted below.

##  Permission to use, copy, modify, and distribute this software and its
##  documentation for any purpose and without fee is hereby granted,
##  provided that the above copyright notice appear in all copies and that
##  both that copyright notice and this permission notice appear in
##  supporting documentation.

##  THE AUTHOR DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
##  INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO
##  EVENT SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, INDIRECT OR
##  CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF
##  USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR
##  OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
##  PERFORMANCE OF THIS SOFTWARE.

# This module defines 3d geometrical vectors with the standard
# operations on them. The elements are stored in an
# array.
#
# Written by: Konrad Hinsen <hinsen@cnrs-orleans.fr>
# Last revision: 1999-7-23
#

_undocumented = 1

import numpy

from . import TensorModule


class Vector:
    """Vector in 3D space

    Constructor:

    - Vector(|x|, |y|, |z|)   (from three coordinates)
    - Vector(|coordinates|)   (from any sequence containing three coordinates)

    Vectors support the usual arithmetic operations
    ('v1', 'v2': vectors, 's': scalar):

    -  'v1+v2'           (addition)
    -  'v1-v2'           (subtraction)
    -  'v1*v2'           (scalar product)
    -  's*v1', 'v1*s'    (multiplication with a scalar)
    -  'v1/s'            (division by a scalar)

    The three coordinates can be extracted by indexing.

    Vectors are *immutable*, i.e. their elements cannot be changed.

    Vector elements can be any objects on which the standard
    arithmetic operations plus the functions sqrt and arccos are defined.
    """

    is_vector = 1

    def __init__(self, x=None, y=None, z=None):
        if x is None:
            self.array = [0., 0., 0.]
        elif y is None and z is None:
            self.array = x
        else:
            self.array = [x, y, z]
        self.array = numpy.array(self.array)

    def __getstate__(self):
        return list(self.array)

    def __setstate__(self, state):
        self.array = numpy.array(state)

    def __copy__(self, memo=None):
        return self

    __deepcopy__ = __copy__

    def __repr__(self):
        return 'Vector(%s,%s,%s)' % (repr(self.array[0]), \
                                     repr(self.array[1]), repr(self.array[2]))

    def __str__(self):
        return repr(list(self.array))

    def __add__(self, other):
        return Vector(self.array + other.array)

    __radd__ = __add__

    def __neg__(self):
        return Vector(-self.array)

    def __sub__(self, other):
        return Vector(self.array - other.array)

    def __rsub__(self, other):
        return Vector(other.array - self.array)

    def __mul__(self, other):
        if isVector(other):
            return numpy.add.reduce(self.array * other.array)
        elif TensorModule.isTensor(other):
            product = TensorModule.Tensor(self.array).dot(other)
            if product.rank == 1:
                return Vector(product.array)
            else:
                return product
        elif hasattr(other, "_product_with_vector"):
            return other._product_with_vector(self)
        else:
            return Vector(numpy.multiply(self.array, other))

    def __rmul__(self, other):
        if TensorModule.isTensor(other):
            product = other.dot(TensorModule.Tensor(self.array))
            if product.rank == 1:
                return Vector(product.array)
            else:
                return product
        else:
            return Vector(numpy.multiply(self.array, other))

    def __div__(self, other):
        if isVector(other):
            raise TypeError("Can't divide by a vector")
        else:
            return Vector(numpy.divide(self.array, 1. * other))

    def __rdiv__(self, other):
        raise TypeError("Can't divide by a vector")

    def __cmp__(self, other):
        if numpy.add.reduce(abs(self.array - other.array)) < 0:
            return -1
        elif numpy.add.reduce(abs(self.array - other.array)) == 0:
            return 0
        else:
            return 1

    def __len__(self):
        return 3

    def __getitem__(self, index):
        return self.array[index]

    def x(self):
        "Returns the x coordinate."
        return self.array[0]

    def y(self):
        "Returns the y coordinate."
        return self.array[1]

    def z(self):
        "Returns the z coordinate."
        return self.array[2]

    def length(self):
        "Returns the length (norm)."
        return numpy.sqrt(numpy.add.reduce(self.array * self.array))

    def normal(self):
        "Returns a normalized copy."
        len = numpy.sqrt(numpy.add.reduce(self.array * self.array))
        if len == 0:
            raise ZeroDivisionError("Can't normalize a zero-length vector")
        return Vector(numpy.divide(self.array, len))

    def cross(self, other):
        "Returns the cross product with vector |other|."
        if not isVector(other):
            raise TypeError("Cross product with non-vector")
        return Vector(self.array[1] * other.array[2]
                      - self.array[2] * other.array[1],
                      self.array[2] * other.array[0]
                      - self.array[0] * other.array[2],
                      self.array[0] * other.array[1]
                      - self.array[1] * other.array[0])

    def asTensor(self):
        "Returns an equivalent tensor object of rank 1."
        return TensorModule.Tensor(self.array, 1)

    def dyadicProduct(self, other):
        "Returns the dyadic product with vector or tensor |other|."
        if isVector(other):
            return TensorModule.Tensor(self.array, 1) * \
                   TensorModule.Tensor(other.array, 1)
        elif TensorModule.isTensor(other):
            return TensorModule.Tensor(self.array, 1) * other
        else:
            raise TypeError("Dyadic product with non-vector")

    def angle(self, other):
        "Returns the angle to vector |other|."
        if not isVector(other):
            raise TypeError("Angle between vector and non-vector")
        cosa = numpy.add.reduce(self.array * other.array) / \
               numpy.sqrt(numpy.add.reduce(self.array * self.array) * \
                          numpy.add.reduce(other.array * other.array))
        cosa = max(-1., min(1., cosa))
        return numpy.arccos(cosa)


# Type check

def isVector(x):
    "Return 1 if |x| is a vector."
    return hasattr(x, 'is_vector')


# Some useful constants

ex = Vector(1, 0, 0)
ey = Vector(0, 1, 0)
ez = Vector(0, 0, 1)
null = Vector(0., 0., 0.)
