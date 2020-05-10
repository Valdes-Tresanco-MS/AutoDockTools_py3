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
#  Modification date: 10/5/20 17:47                                                                #
#                                                                                                  #
# ##################################################################################################

def uniq(alist):  # Fastest order preserving
    set = {}
    return [set.setdefault(e, e) for e in alist if e not in set]


def uniq3(alist):  # Fastest without order preserving
    set = {}
    list(map(set.__setitem__, alist, []))
    return list(set.keys())


"""
from mglutil.util.uniq import uniq, uniq2, uniq3
import time
a=range(100)
b=range(10)
c=a+b

t1=time.time()
for i in range(5000):  x=uniq(c)
print time.time()-t1

t1=time.time()
for i in range(5000):  x=uniq2(c)
print time.time()-t1

t1=time.time()
for i in range(5000):  x=uniq3(c)
print time.time()-t1


>>>
0.865363121033
0.463307857513
0.260641098022
"""
