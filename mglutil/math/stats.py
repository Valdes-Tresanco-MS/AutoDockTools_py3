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

import warnings


def stats(values):
    """returns the mimn, max, mean and standard deviation of a list of values"""

    warnings.warn(
        "\n\nWARNING!! This function has been deprecated!!\n \
        Use the stats in Volume/Grid3D.py\n", DeprecationWarning, 2)

    npts = len(values)
    if npts:
        from math import sqrt
        sum = 0.0
        sumsq = 0.0
        mini = maxi = values[0]
        for v in values:
            sum += v
            sumsq += float(v) * float(v)
            if v < mini:
                mini = v
            if v > maxi:
                maxi = v
        mean = float(sum) / npts
        stdev = sqrt((sumsq - (sum * sum / float(npts))) / (npts - 1))
        return mini, maxi, mean, stdev

    else:
        return (0., 0., 1., 1.)
