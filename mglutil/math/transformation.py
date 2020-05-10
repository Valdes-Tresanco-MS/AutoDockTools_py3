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
# $Header: /opt/cvs/python/packages/share1.5/mglutil/math/transformation.py,v 1.45.12.1 2016/02/11 23:15:05 annao Exp $
#

import numpy

from mglutil.math import VectorModule
from mglutil.math import rotax

Vector = VectorModule.Vector


class Quaternion:
    """ Base Quaternion class
    """

    def __init__(self, data=(1.0, numpy.array((0., 0., 0.), 'f'))):
        """data is in the form ( c, (x y, z)), where c is the 
        real part (float) and (x,y,z) is the pure part (
        array of floats)
        """
        try:
            self.real = float(data[0])
            self.pure = numpy.array((data[1][0], data[1][1], data[1][2]), 'f')
        except:
            raise ValueError("1Arguments must be (c,(x,y,z))")
        if len(self.pure) != 3:
            raise ValueError("2Arguments must be (c,(x,y,z))")

    def __repr__(self):
        """ representation of a general quaternion must be (real,pure),
        since not all quaternions are rotations
        """
        result = "Quaternion (%g (%g %g %g))" % \
                 (self.real, self.pure[0],
                  self.pure[1], self.pure[2])
        return result

    def __add__(self, other):
        """ Get the sum of two quaternions.
        """
        real = self.real + other.real
        pure = self.pure + other.pure
        return Quaternion((real, pure))

    def __mul__(self, other):
        """ Multiply two quaternions together.
        For unit Quaternons, this is equivalent to concatenating rotations"""
        real = self.real * other.real - \
               numpy.inner(self.pure, other.pure)
        v1 = self.pure
        v2 = other.pure
        cofactor1 = v1[1] * v2[2] - v1[2] * v2[1]
        cofactor2 = v1[2] * v2[0] - v1[0] * v2[2]
        cofactor3 = v1[0] * v2[1] - v1[1] * v2[0]
        pure = numpy.array([cofactor1, cofactor2, cofactor3]) \
               + self.real * other.pure \
               + other.real * self.pure
        return Quaternion((real, pure))

    def conjugate(self):
        """ The conjugate of quaternion (c,(x,y,z)) is (c,(-x,-y,-z))
        So the product of a quaternion and its conjugate is its
        magnitude
        """
        pure = -self.pure
        real = self.real
        return Quaternion((real, pure))

    def magnitude(self):
        """ Quicker than multiplying conjugates"""
        return self.real ** 2 + numpy.inner(self.pure, self.pure)

    def inverse(self):
        """Get the multiplicative inverse of a quaternion"""
        real = self.real / self.magnitude()
        pure = -self.pure / self.magnitude()
        return Quaternion((real, pure))

    def normal(self):
        """Normalise a quaternion by dividing throughout by the
        magnitude
        """
        M = numpy.sqrt(self.magnitude())
        self.pure = self.pure / M
        self.real = self.real / M


class UnitQuaternion(Quaternion):
    """ Special subclass of Quaternions with magnitude 1.0
    Can be used to represent rotations, in which case real =
    cos(theta/2) and pure = sin(theta/2)*(unit direction vector)
    Input can also be given in the form (x,y,z,theta), where (x,y,z)
    is the rotation axis (not necessarily normalized) and theta is the
    rotation angle in degrees.
    """

    def __init__(self, data=(1.0, numpy.array((0., 0., 0.), 'f'))):
        """ (real,(pure x,pure y,pure z)) or (x,y,z,theta) (theta in degrees)
        """
        if len(data) == 2:
            self.real = data[0]
            try:
                theta = numpy.arccos(self.real)
                self.pure = numpy.array((data[1][0], data[1][1], data[1][2]), 'f')
            except:
                raise ValueError("The real part must be between -1.0 and 1.0")
        elif len(data) == 4:
            theta = numpy.pi * data[3] / 360.
            self.real = numpy.cos(theta)
            self.pure = numpy.sin(theta) * numpy.array((data[0], data[1],
                                                        data[2]), 'f')
        else:
            raise ValueError("Args must be (x,y,z,theta) or (real,pure)")
        self.normal()

    def normal(self):
        if self.real != 1.:
            theta = numpy.arccos(self.real)
            vector = self.pure / numpy.sin(theta)
            vector = vector / numpy.sqrt(numpy.inner(vector, vector))
            self.pure = numpy.sin(theta) * vector
        else:
            self.pure = numpy.zeros(3, 'f')

    def __repr__(self):
        """Representation of a unit quaternion is as rx,ry,rz,theta,
        so we can see what it does
        """
        if self.real != 1.:
            # if it is not the identity
            theta = numpy.arccos(self.real)
            angle = 360 * theta / numpy.pi
            xyz = self.pure / numpy.sin(theta)
        else:
            # if it is the identity
            angle = 0.
            xyz = self.pure
        return "Unit Quaternion %7.4f %7.4f %7.4f %7.3f" % \
               (xyz[0], xyz[1], xyz[2], angle)

    def __mul__(self, other):
        # same as Quaternion, except return another UnitQuaternion
        result = Quaternion.__mul__(self, other)
        return UnitQuaternion((result.real, result.pure))

    def conjugate(self):
        result = Quaternion.conjugate(self)
        return UnitQuaternion((result.real, result.pure))

    def inverse(self):
        return self.conjugate()

    def getAxisAndAngleDegres(self):
        """Given a quaternion, compute axis and angle.
"""
        theta = numpy.arccos(self.real)
        angle = 360 * theta / numpy.pi
        xyz = self.pure / numpy.sin(theta)
        return xyz, angle

    def getRotMatrix(self, shape=(4, 4), transpose=None):
        """return the rotation matrix as a  array of shape shape.
        """
        try:
            assert (shape in [(3, 3), (4, 4), (9,), (16,)])
        except:
            raise ValueError('shape must be (3,3), (4,4), (9,) or (16,)')

        # get the inverse 4x4 from rotax
        mtx = rotax.rotax(numpy.array([0., 0., 0.], 'f'), self.pure, 2 * numpy.arccos(self.real))

        # strip if necessary
        if shape in ((3, 3), (9,)):
            mtx = [x[:3] for x in mtx]
            mtx = mtx[:3]

        if not transpose:
            return numpy.reshape(numpy.transpose(mtx), shape)
        else:
            return numpy.reshape(mtx, shape)

    def apply(self, points):
        # apply the rotational part alone to a point or list of points
        # can be homogeneous coordinates or not.
        pshape = numpy.shape(points)
        homogeneous = 1
        if len(pshape) == 1:
            if pshape[0] == 3:
                points = numpy.array(numpy.concatenate((points, numpy.ones(1, 'f')), 1))
                homogeneous = 0
        elif len(pshape) == 2:
            if pshape[1] == 3:
                points = numpy.array(numpy.concatenate(
                    (numpy.array(points), numpy.ones((pshape[0], 1), 'f')), 1))
                homogeneous = 0
        mtx = self.getRotMatrix((4, 4), transpose=1)
        newpoints = numpy.dot(points, mtx)
        if homogeneous:
            return newpoints
        else:
            # strip the final zero off the coordinates
            if len(pshape) == 1:
                return newpoints[:3]
            else:
                newpoints = [x[:3] for x in newpoints]
                return newpoints


class Transformation(UnitQuaternion):
    """ Base class for manipulating transformations.
    """

    def __init__(self, trans=numpy.array([0., 0., 0., 1.], 'f'),
                 quaternion=numpy.array([0., 0., 0., 0.], 'f'),
                 scale=numpy.array([1., 1., 1., 1.], 'f')):
        UnitQuaternion.__init__(self, quaternion)
        # make the translation homogeneous if it isn't
        if len(trans) == 3:
            trans = list(trans)
            trans.append(1.)
        self.trans = numpy.array((trans[0], trans[1], trans[2], trans[3]), 'f')

    def __repr__(self):
        """ Representation is of the form tx,ty,tz,qx,qy,qz,theta
        """
        #  first check for identity quaternion to avoid nans
        if self.real != 1:
            theta = numpy.arccos(self.real)
            angle = 360 * theta / numpy.pi
            xyz = self.pure / numpy.sin(theta)
        else:
            angle = 0.
            xyz = self.pure
        result = "Transformation: tx ty tz rx ry rz angle\n %g %g %g %g %g %g %g" \
                 % (self.trans[0], self.trans[1], self.trans[2],
                    xyz[0], xyz[1], xyz[2], angle)
        return result

    def output(self):
        """ As __repr__ but without the explanation. For getting the numbers only
        """
        if self.real != 1:
            theta = numpy.arccos(self.real)
            angle = 360 * theta / numpy.pi
            xyz = self.pure / numpy.sin(theta)
        else:
            angle = 0.
            xyz = self.pure
        result = "%g %g %g %g %g %g %g" % (self.trans[0], self.trans[1], self.trans[2],
                                           xyz[0], xyz[1], xyz[2], angle)
        return result

    def __mul__(self, other):
        """ concatenate two transformations. self*other (other performed first).
        """
        # combined rotation is the product of the two rotations (Rself*Rother):
        v1 = self.pure
        v2 = other.pure
        real = self.real * other.real - \
               numpy.inner(v1, v2)
        cofactor1 = v1[1] * v2[2] - v1[2] * v2[1]
        cofactor2 = v1[2] * v2[0] - v1[0] * v2[2]
        cofactor3 = v1[0] * v2[1] - v1[1] * v2[0]
        pure = numpy.array([cofactor1, cofactor2, cofactor3]) \
               + self.real * other.pure \
               + other.real * self.pure
        # combined translation
        trans = self.getQuaternion().apply(other.trans) + self.trans
        trans[3] = 1.
        return Transformation(trans=trans, quaternion=(real, pure))

    def reset(self):
        self.real = 1.0
        self.pure = numpy.array((0., 0., 0.))
        self.trans = numpy.array([0., 0., 0., 1.])

    def getQuaternion(self):
        return UnitQuaternion((self.real, self.pure))

    def getTranslation(self, shape=(4,)):
        """ get the translation vector with shape = (3,) or (4,)
        (default is (4,))
        """
        if shape == (3,):
            return self.trans[:3]
        elif shape == (4,):
            return self.trans
        else:
            raise ValueError("Shape must be (3,) or (4,)")

    def getMatrix(self, shape=(4, 4), transpose=None):
        mtx = self.getRotMatrix((4, 4), transpose=transpose)  # from Quaternion
        mtx[3] = self.getTranslation()
        if transpose:
            return numpy.reshape(mtx, shape)
        else:
            return numpy.reshape(numpy.transpose(mtx), shape)

    def getDejaVuMatrix(self):
        """returns a 4x matrix usable as an instance matrix"""
        mtx = self.getRotMatrix((4, 4), transpose=None)  # from Quaternion
        mtx[3] = self.getTranslation()
        mtx[:3, 3] = mtx[3, :3]
        mtx[3, :3] = [0, 0, 0]
        return mtx

    def apply(self, points):
        """ Apply the entire transformation to a list of points
        """
        pshape = numpy.shape(points)
        homogeneous = 1
        if len(pshape) == 1:
            if pshape[0] == 3:
                points = numpy.array(numpy.concatenate((points, numpy.ones(1, 'f')), 1))
                homogeneous = 0
        elif len(pshape) == 2:
            if pshape[1] == 3:
                points = numpy.array(numpy.concatenate(
                    (numpy.array(points), numpy.ones((pshape[0], 1), 'f')), 1))
                homogeneous = 0
        mtx = self.getMatrix((4, 4), transpose=1)
        newpoints = numpy.dot(points, mtx)
        if homogeneous:
            return newpoints
        else:
            # strip the final one off the coordinates
            if len(pshape) == 1:
                return newpoints[:3]
            else:
                newpoints = [x[:3] for x in newpoints]
                return newpoints

    def inverse(self):
        # inverse rotation is the same as for a pure rotation
        real = self.real
        pure = -self.pure
        # inverse translation is application of inverse rotation
        transl = -numpy.dot(self.getRotMatrix(transpose=1, shape=(3, 3)),
                            self.trans[:3])
        return Transformation(trans=transl, quaternion=(real, pure))

    def getScrewAxis(self, center=None, linelength=None):
        """ Get the representation of a transformation in screw
        format. Returns two points on the axis and the translation
        component along the axis.
        Takes an optional center argument. The first point returned is
        then the point on the axis nearest to the center.
        The optional linelength argument defines the distance between the
        two points returned. The default is the translation component.    
        """

        # first check that there is a rotational component. If not, abort
        # if there is a rotation, self.real != 1.0
        if self.real <= 0.99999999:
            # need the direction to determine which way to draw the line
            trans = Vector(self.trans[:3])
            theta = numpy.arccos(self.real)
            axis = self.pure / numpy.sin(theta)
            axis = Vector(axis)
            screw = (trans * axis)
            tpar = screw * axis
            tper = trans - tpar
            cpt1 = tper / 2.
            length = tper.length()
            height = length / (2 * numpy.tan(theta))
            cpt2 = height * (axis.cross(tper)).normal()
            point = cpt1 + cpt2
            if center:
                try:
                    center = Vector(center)
                except:
                    raise ValueError('center must be a numpy array of shape (3,)')
                m = (center - point) * axis
                point = point + m * axis
            if not linelength:
                return point, point + axis * screw, screw
            else:
                return point, point + linelength * numpy.sign(screw) * axis, screw

        else:
            return None
