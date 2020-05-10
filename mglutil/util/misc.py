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
#  Modification date: 10/5/20 17:36                                                                #
#                                                                                                  #
# ##################################################################################################

#############################################################################
#
# Author: Michel F. SANNER
#
# Copyright: M. Sanner TSRI 2000
#
#############################################################################


import os
import sys

import numpy

_proc_status = '/proc/%d/status' % os.getpid()

_scale = {'kB': 1024.0, 'mB': 1024.0 * 1024.0,
          'KB': 1024.0, 'MB': 1024.0 * 1024.0}


def _VmB(VmKey):
    """Private.
    """
    global _proc_status, _scale
    # get pseudo file  /proc/<pid>/status
    try:
        t = open(_proc_status)
        v = t.read()
        t.close()
    except:
        return 0.0  # non-Linux?
    # get VmKey line e.g. 'VmRSS:  9999  kB\n ...'
    i = v.index(VmKey)
    v = v[i:].split(None, 3)  # whitespace
    if len(v) < 3:
        return 0.0  # invalid format?
    # convert Vm value to bytes
    return float(v[1]) * _scale[v[2]]


def memory(since=0.0):
    """Return memory usage in bytes.
    """
    return _VmB('VmSize:') - since


def resident(since=0.0):
    """Return resident memory usage in bytes.
    """
    return _VmB('VmRSS:') - since


def stacksize(since=0.0):
    """Return stack size in bytes.
    """
    return _VmB('VmStk:') - since


def issequence(a):
    return type(a) is tuple or \
           type(a) is list or \
           isinstance(a, numpy.ndarray)


def isnumericstring(a):
    try:
        float(a)
        return 1
    except:
        return 0


def uniq(objectSequence):
    """Remove the duplicates from a list while keeping the original
    list order """
    l = []
    d = {}
    for o in objectSequence:
        if o not in d:
            d[o] = None
            l.append(o)
    return l


def deepCopySeq(sequence):
    """ returns the deep copy of the given sequence """
    assert type(sequence) in (tuple, list, type(numpy.array([1, 2, 3])))
    if hasattr(sequence, 'copy'):
        dcSeq = sequence.copy()
    else:
        dcSeq = sequence[:]

    return dcSeq


def ensureFontCase(font):
    return font


#    from Tkinter import TkVersion
#    lFont = font[0].upper() + font[1:].lower()
#    if TkVersion == '8.4' and sys.platform != "win32":
#        lFont = font.lower()
#    return lFont


def isInstance(lObject):
    if sys.version.startswith('2.5'):  # detect python25
        if type(lObject) == list:
            return True
        else:
            return False
    else:
        import inspect
        ltype = type(lObject)
        if ltype == type:
            return True
        elif inspect.isclass(lObject) is False and isinstance(lObject, ltype) is True:
            from abc import ABCMeta
            if ltype == type is True:
                return True
            elif type(ltype) == ABCMeta:
                return True
            else:
                return False
        else:
            return False


def importMainOrIPythonMain():
    try:
        from IPython import ipapi
        mainDict = ipapi.get().user_ns
    except:
        mainDict = __import__('__main__').__dict__
    return mainDict


def suppressMultipleQuotes(aString):
    lStringToSimplify = aString
    lSimplifiedString = lStringToSimplify
    while type(lStringToSimplify) == bytes:
        lSimplifiedString = lStringToSimplify
        try:
            lStringToSimplify = eval(lSimplifiedString)
        except:
            break
    return lSimplifiedString


class IntVar:
    def __init__(self, val=0):
        self.set(val)

    def get(self):
        return self.val

    def set(self, val):
        self.val = int(val)


class StringVar:
    def __init__(self, val=""):
        self.set(val)

    def get(self):
        return self.val

    def set(self, val):
        self.val = str(val)


class BooleanVar:
    def __init__(self, val=False):
        self.set(val)

    def get(self):
        return self.val

    def set(self, val):
        self.val = (val == True)
