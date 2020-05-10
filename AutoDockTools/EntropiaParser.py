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
#  Modification date: 10/5/20 18:51                                                                #
#                                                                                                  #
# ##################################################################################################

#
# Last modified on Thu Mar 11 10:50:58 PST 2004 by lindy
#
# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/EntropiaParser.py,v 1.7 2004/03/11 19:17:04 lindy Exp $
#
#


from AutoDockTools.ResultParser import ResultParser


class EntropiaParser(ResultParser):
    """<class docstring>
    """
    # keywords are the keys of the dictionary used by the
    # ConformationHandler to instantiate a Conformation.
    keywords = ResultParser.keywords + [
        # add keywords here
    ]

    def __init__(self):
        ResultParser.__init__(self)

    def parseline(self, line):
        """Parse the given line and return append a dictionary
        onto self.clist, the superclass-declared list of conformations.
        """
        line_list = line.split()
        # make a heuristic decision about which keys to use!!
        if line_list[2] == repr(1.01) or line_list[2][:3] == repr(1.0):
            # version 1.01 Entropia result file
            dict = self._parseStateLineList(line_list, self.result_keys_v101)
        elif float(line_list[9]) == 1.0:
            # version 1.0 Entropia result file
            dict = self._parseStateLineList(line_list)
        else:
            print("Unparsable result file line: %s" % (line))
            return  # without doing anything
        self.clist.append(dict)

    result_keys = [  # original version 1.00
        'output_id',
        'data_run_id',
        'dpf_id',
        'creation_dtime',
        'last_update_dtime',
        'ei_version',  # Entropia Interface version
        'ag_version',  # AutoGrid version
        'ad_version',  # AutoDock version
        'run_rank',  # rank of this run in cluster
        'run_number',  # number of this run in dpf
        'cluster_rank',  # rank of this cluster
        'cluster_size',  # number of conformations in this cluster
        'run_size',  # number of runs specified in dpf
        'rseed1',
        'rseed2',
        'rmsd',  # RMSD from reference structure
        'binding_energy',  # estimated free energy of binding
        'docking_energy',  # final docked energy
        'trn_x',  # translation x, y, z
        'trn_y',
        'trn_z',
        'qtn_nx',  # quaternion unit vector x, y, z
        'qtn_ny',
        'qtn_nz',
        'qtn_ang_deg',  # quaternion rotation angle
        'num_torsions',
        'torsion_values']

    result_keys_v101 = [  # for the 1.01 version of the Entropia Interface
        'data_run_id',
        'dpf_id',
        'ei_version',  # Entropia Interface version
        'ag_version',  # AutoGrid version
        'ad_version',  # AutoDock version
        'run_rank',  # rank of this run in cluster
        'run_number',  # number of this run in dpf
        'cluster_rank',  # rank of this cluster
        'cluster_size',  # number of conformations in this cluster
        'run_size',  # number of runs specified in dpf
        'rseed1',
        'rseed2',
        'rmsd',  # RMSD from reference structure
        'binding_energy',  # estimated free energy of binding
        'docking_energy',  # final docked energy
        'trn_x',  # translation x, y, z
        'trn_y',
        'trn_z',
        'qtn_nx',  # quaternion unit vector x, y, z
        'qtn_ny',
        'qtn_nz',
        'qtn_ang_deg',  # quaternion rotation angle
        'num_torsions',
        'torsion_values']

    type_conv = {
        int: lambda x: int(x),
        float: lambda x: float(x),
        bytes: lambda x: x,
        'PercentType': lambda x: float(x[:-1]),
        'TimeType': lambda x: x
    }

    field_defn = {
        #    key                : (type, whitespace delimited sub-fields)
        'output_id': (int, 1),
        'data_run_id': (int, 1),
        'dpf_id': (int, 1),
        'creation_dtime': ('TimeType', 3),
        'last_update_dtime': ('TimeType', 3),
        'ei_version': (float, 1),
        'ag_version': (float, 1),
        'ad_version': (float, 1),
        'run_rank': (int, 1),
        'run_number': (int, 1),
        'cluster_rank': (int, 1),
        'cluster_size': (int, 1),
        'run_size': (int, 1),
        'rseed1': (int, 1),
        'rseed2': (int, 1),
        'rmsd': (float, 1),
        'binding_energy': (float, 1),
        'docking_energy': (float, 1),
        'trn_x': (float, 1),
        'trn_y': (float, 1),
        'trn_z': (float, 1),
        'qtn_nx': (float, 1),
        'qtn_ny': (float, 1),
        'qtn_nz': (float, 1),
        'qtn_ang_deg': (float, 1),
        'num_torsions': (int, 1),
        'torsion_values': (float, None)
    }

    def _parseStateLineList(self, lineList=None, keys=None):
        """Return a dictionary of values.
        """
        if not keys:
            keys = self.result_keys

        # initalize the stateDict with all fields
        stateDict = {}
        for field in self.result_keys:
            stateDict[field] = None

        # now run through the specified fields
        subFldIx = 0
        for field in keys:

            # most fields just get converted in the first case
            if self.field_defn[field][1] == 1:
                stateDict[field] = self.type_conv[
                    self.field_defn[field][0]](lineList[subFldIx])
                subFldIx = subFldIx + 1

            # the time fields have spaces in them; handled here
            elif self.field_defn[field][0] == 'TimeType':
                timeList = []
                for subIx in range(self.field_defn[field][1]):
                    # just append the strings to the timeList
                    timeList.append(lineList[subFldIx])
                    subFldIx = subFldIx + 1
                stateDict[field] = timeList

            # this case handles the variable number of torsions
            elif field == 'torsion_values':
                torsionList = []
                for subIx in range(stateDict['num_torsions']):
                    # loop through the torsion_values
                    torsionList.append(self.type_conv[
                                           self.field_defn[field][0]](lineList[subFldIx]))
                    subFldIx = subFldIx + 1
                stateDict[field] = torsionList

        # make sure there are no left over fields
        assert len(lineList) == subFldIx
        return stateDict
