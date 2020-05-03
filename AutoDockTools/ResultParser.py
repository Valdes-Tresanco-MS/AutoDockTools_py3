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
#  Modification date: 2/5/20 19:51                                                                 #
#                                                                                                  #
# ##################################################################################################

#
# Last modified on Mon Mar  4 14:35:36 PST 2002 by lindy
#
# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/ResultParser.py,v 1.3 2009/01/09 19:53:49 rhuey Exp $
#




class ResultParser:
    """Base class of result parsers. Supply your own parseline
    function which appends a conformation dictionary to self.clist
    """
    def __init__(self):
        self.clist = []

    # add more keywords in subclass if you like 
    keywords = [
        'cluster',            # number of cluster
        'cluster_rank',       # rank within cluster
        'rmsd_ref',           # distance to reference structure
        'rmsd_seed',          # distance to lowest energy conf. in cluster
        'binding_energy',     # estimated free energy of binding
        'docking_energy',     # final docked energy
        'internal_energy',
        'trn_x',              # translation x, y, z
        'trn_y',
        'trn_z',
        'qtn_nx',             # quaternion unit vector x, y, z
        'qtn_ny',
        'qtn_nz',
        'qtn_ang_deg',        # quaternion rotation angle
        'num_torsions',
        'torsion_values',
        'rseed1',             # the random number seeds for this conformation
        'rseed2',
        ]

    def parse(self, filename):
        """
        """
        file_ptr = open(filename)
        
        self.clist = []
        for line in file_ptr.readlines():
            self.parseline(line)
        file_ptr.close()

        self.filename = filename
        
        return self.clist


    def parseline(self, line):
        """over ride me"""
        pass
    
