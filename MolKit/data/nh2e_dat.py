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

nh2e_dat = {
    "filename": 'nh2e.in',
    "NAMDBF": 'db4.dat',
    "IDBGEN,IREST,ITYPF": ['1', '1', '2'],
    "NHE": {"N": {'torsion': 180.0, 'tree': 'M', 'NC': 1, 'NB': 2, 'NA': 3, 'I': 4, 'angle': 116.6, 'blen': 1.335,
                  'charge': -0.463, 'type': 'N'},
            "INTX,KFORM": ['INT', '1'],
            "atNameList": ['N', 'HN1', 'HN2'],
            "DUMM": [['1', 'DUMM', 'DU', 'M', '0', '-1', '-2', '0.0000', '0.0000', '0.0000'],
                     ['2', 'DUMM', 'DU', 'M', '1', '0', '-1', '1.0000', '0.0000', '0.0000'],
                     ['3', 'DUMM', 'DU', 'M', '2', '1', '0', '1.0000', '90.0000', '0.0000']],
            "IFIXC,IOMIT,ISYMDU,IPOS": ['CORRECT', 'OMIT', 'DU', 'BEG'],
            "HN2": {'torsion': 180.0, 'tree': 'E', 'NC': 2, 'NB': 3, 'NA': 4, 'I': 6, 'angle': 119.8, 'blen': 1.01,
                    'charge': 0.2315, 'type': 'H'},
            "impropTors": [['-M', 'HN1', 'N', 'HN2']],
            "CUT": ['0.00000'],
            "HN1": {'torsion': 0.0, 'tree': 'E', 'NC': 2, 'NB': 3, 'NA': 4, 'I': 5, 'angle': 119.8, 'blen': 1.01,
                    'charge': 0.2315, 'type': 'H'},
            "NAMRES": 'NH2 ENDING GROUP',
            },
}
