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

templates_attributes = {}
templates_attributes['IDD594'] = {'type': 'Carboxylic_Acid', 'modelpka': 3.90, 'titratableatoms': ['O1', 'O2', 'C']}
templates_attributes['CarboxyGroup'] = {'type': 'Carboxylic_Acid', 'modelpka': 4.75,
                                        'titratableatoms': ['O10', 'O20', 'C3']}
templates_attributes['AceticAcid'] = {'type': 'Carboxylic_Acid', 'modelpka': 4.80, 'titratableatoms': ['O1', 'O2', 'C']}
templates_attributes['imidazole'] = {'type': 'Base', 'modelpka': 6.30,
                                     'titratableatoms': ['N10', 'C10', 'N20', 'C30', 'C40']}
templates_attributes['piperidine'] = {'type': 'Base', 'modelpka': 6.33, 'titratableatoms': ['XXX1']}
templates_attributes['PropanoicAcid'] = {'type': 'Carboxylic_Acid', 'modelpka': 4.20,
                                         'titratableatoms': ['O10', 'O20', 'C30']}
templates_attributes['acetylsalicylicacid'] = {'type': 'Carboxylic_Acid', 'modelpka': 4.99,
                                               'titratableatoms': ['O1', 'O2', 'C', 'C30', 'C50', 'C80', 'C60', 'C40']}

# templates_attributes['CRAP']         ={'type': 'Acid','modelpka':4.555, 'titratableatoms': ['O10', 'O20','C30']}
