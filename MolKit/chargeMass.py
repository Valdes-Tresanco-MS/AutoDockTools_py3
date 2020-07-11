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


def getChargeMass(atoms):
    elist = []
    totmass = 0.0

    for a in atoms:
        if a.element == 'H':
            e = 1.0
            totmass += 1.00794
        elif a.element == 'C':
            e = 6.0
            totmass += 12.0107
        elif a.element == 'AU':
            e = 79.0
            totmass += 196.96655
        elif a.element == 'N':
            e = 7.0
            totmass += 14.00674
        elif a.element == 'O':
            e = 8.0
            totmass += 15.9994
        elif a.element == 'P':
            e = 15.0
            totmass += 30.973761
        elif a.element == 'S':
            e = 16.0
            totmass += 32.066
        elif a.element == 'W':
            e = 18.0
            totmass += 1.00794 * 2 + 15.9994  # ficticious water 'atom'
        else:
            print("skipping unknown atom %s" % a.element)
            e = -1
        elist.append(e)

    return elist, totmass
