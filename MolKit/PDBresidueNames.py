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

#############################################################################
#
# Author: Michel F. SANNER
#
# Copyright: M. Sanner TSRI 2010
#
#
#############################################################################

DNAnames = {
    'DC': 'D', 'DG': 'G', 'DA': 'A', 'DT': 'T', 'T': 'T', 'DI': 'I', 'N': 'N',
}

RNAnames = {
    'C': 'C', 'G': 'G', 'A': 'A', 'U': 'U', 'I': 'I', 'N': 'N',
}

Nucleotides = DNAnames.copy()
Nucleotides.update(RNAnames)

AAnames = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', 'GLY': 'G',
    'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
    'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S', 'THR': 'T', 'VAL': 'V',
    'TRP': 'W', 'TYR': 'Y',
    # the follwould be added automatically if the
    # MODRES ris present in the pdb file but we put
    # them herays
    'HID': '?', 'HSP': '?', 'HIE': '?', 'HIP': '?', 'CYX': '?',
    'CSS': '?', 'ACE': '?', 'MSE': '?', '5HP': '?', 'SOC': '?',
}

#
# list of resames for ions taken from
# http://decogroup.org/ion_list.txt
#
ionNames = {
    '1CU': '?',
    '2HP': '?',
    '2MO': '?',
    '2OF': '?',
    '3CO': '?',
    '3MT': '?',
    '3NI': '?',
    '4MO': '?',
    '543': '?',
    '6MO': '?',
    'ACT': '?',
    'AG': '?',
    'AL': '?',
    'ALF': '?',
    'ATH': '?',
    'AU': '?',
    'AU3': '?',
    'AUC': '?',
    'AZI': '?',
    'BA': '?',
    'BCT': '?',
    'BEF': '?',
    'BF4': '?',
    'BO4': '?',
    'BR': '?',
    'CA': '?',
    'CAC': '?',
    'CD': '?',
    'CD1': '?',
    'CD3': '?',
    'CD5': '?',
    'CE': '?',
    'CHT': '?',
    'CL': '?',
    'CO': '?',
    'CO3': '?',
    'CO5': '?',
    'CON': '?',
    'CR': '?',
    'CS': '?',
    'CU': '?',
    'CU1': '?',
    'CUA': '?',
    'CUZ': '?',
    'CYN': '?',
    'DMI': '?',
    'E4N': '?',
    'EMC': '?',
    'EU': '?',
    'EU3': '?',
    'F': '?',
    'FE': '?',
    'FE2': '?',
    'FPO': '?',
    'GA': '?',
    'GD3': '?',
    'HAI': '?',
    'HG': '?',
    'HGC': '?',
    'HO': '?',
    'IN': '?',
    'IOD': '?',
    'IR': '?',
    'IR3': '?',
    'IRI': '?',
    'IUM': '?',
    'K': '?',
    'KO4': '?',
    'LA': '?',
    'LCO': '?',
    'LCP': '?',
    'LI': '?',
    'LU': '?',
    'MAC': '?',
    'MG': '?',
    'MH2': '?',
    'MH3': '?',
    'MLI': '?',
    'MLT': '?',
    'MMC': '?',
    'MN': '?',
    'MN3': '?',
    'MN5': '?',
    'MO1': '?',
    'MO2': '?',
    'MO3': '?',
    'MO4': '?',
    'MO5': '?',
    'MO6': '?',
    'MOO': '?',
    'MOS': '?',
    'MW1': '?',
    'MW2': '?',
    'MW3': '?',
    'NA': '?',
    'NA2': '?',
    'NA5': '?',
    'NA6': '?',
    'NAO': '?',
    'NAW': '?',
    'NC': '?',
    'NET': '?',
    'NH4': '?',
    'NI': '?',
    'NI1': '?',
    'NI2': '?',
    'NI3': '?',
    'NO2': '?',
    'NO3': '?',
    'O4M': '?',
    'OAA': '?',
    'OC1': '?',
    'OC2': '?',
    'OC3': '?',
    'OC4': '?',
    'OC5': '?',
    'OC6': '?',
    'OC7': '?',
    'OCL': '?',
    'OCM': '?',
    'OCN': '?',
    'OCO': '?',
    'OF1': '?',
    'OF2': '?',
    'OF3': '?',
    'OH': '?',
    'OS': '?',
    'OXL': '?',
    'PB': '?',
    'PBM': '?',
    'PD': '?',
    'PER': '?',
    'PI': '?',
    'PO3': '?',
    'PO4': '?',
    'PR': '?',
    'PT': '?',
    'PT4': '?',
    'PTN': '?',
    'RB': '?',
    'RHD': '?',
    'RU': '?',
    'SB': '?',
    'SCN': '?',
    'SE4': '?',
    'SM': '?',
    'SMO': '?',
    'SO3': '?',
    'SO4': '?',
    'SOH': '?',
    'SR': '?',
    'TB': '?',
    'TCN': '?',
    'TEA': '?',
    'THE': '?',
    'TL': '?',
    'TMA': '?',
    'TRA': '?',
    'UNX': '?',
    'V': '?',
    'VO4': '?',
    'W': '?',
    'WO5': '?',
    'Y1': '?',
    'YB': '?',
    'YT3': '?',
    'ZN': '?',
    'ZN2': '?',
    'ZN3': '?',
    'ZNO': '?',
    'ZO3': '?',
}

waterNames = {'HOH': '?', 'WAT': '?'}

allResidueNames = {}
allResidueNames.update(waterNames)
allResidueNames.update(RNAnames)
allResidueNames.update(AAnames)
allResidueNames.update(DNAnames)
allResidueNames.update(ionNames)
