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

opls_nacl_dat = {
    "CL": {"INTX,KFORM": ['INT', '1'],
           "atNameList": ['CL'],
           "DUMM": [['1', 'DUMM', 'DU', 'M', '0', '-1', '-2', '0.0000', '0.0000', '0.0000', '0.000'],
                    ['2', 'DUMM', 'DU', 'M', '1', '0', '-1', '1.0000', '0.0000', '0.0000', '0.000'],
                    ['3', 'DUMM', 'DU', 'M', '2', '1', '0', '1.0000', '90.0000', '0.0000', '0.000']],
           "IFIXC,IOMIT,ISYMDU,IPOS": ['CORR', 'OMIT', 'DU', 'BEG'],
           "CL": {'torsion': 180.0, 'tree': 'M', 'NC': 1, 'NB': 2, 'NA': 3, 'I': 4, 'angle': 90.0, 'blen': 1.0,
                  'charge': -1.0, 'type': 'CL'},
           "CUT": ['0.000000'],
           "NAMRES": 'Chloride Ion',
           },
    "NAMDBF": 'db4.dat',
    "filename": 'opls_nacl.in',
    "NA": {"INTX,KFORM": ['INT', '1'],
           "atNameList": ['NA'],
           "DUMM": [['1', 'DUMM', 'DU', 'M', '0', '-1', '-2', '0.0000', '0.0000', '0.0000', '0.000'],
                    ['2', 'DUMM', 'DU', 'M', '1', '0', '-1', '1.0000', '0.0000', '0.0000', '0.000'],
                    ['3', 'DUMM', 'DU', 'M', '2', '1', '0', '1.0000', '90.0000', '0.0000', '0.000']],
           "NA": {'torsion': 180.0, 'tree': 'M', 'NC': 1, 'NB': 2, 'NA': 3, 'I': 4, 'angle': 90.0, 'blen': 1.0,
                  'charge': 1.0, 'type': 'SO'},
           "IFIXC,IOMIT,ISYMDU,IPOS": ['CORR', 'OMIT', 'DU', 'BEG'],
           "CUT": ['0.000000'],
           "NAMRES": 'Sodium Ion',
           },
    "IDBGEN,IREST,ITYPF": ['1', '1', '3'],
}
