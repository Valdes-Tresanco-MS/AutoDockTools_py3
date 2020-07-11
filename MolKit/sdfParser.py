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

#
# parse an SDF file V2000 and build a MolKit molecule
#
# Format source http://en.wikipedia.org/wiki/Chemical_table_file
#
# Lines	Section	Description
#
# 1-3	Header
# 1		Molecule name ("benzene") unformated 80 characters
# 2		User/Program/Date/etc information in the following format
#              IIPPPPPPPPMMDDYYHHmmddSSssssssssssEEEEEEEEEEEERRRRRR
#              A2<--A8--><---A10-->A2I2<--F10.5-><---F12.5--><-I6->
#               A2 User's first and last initials (l), \
#               A8 program name (P),
#               A10 date/time (M/D/Y,H:m),
#               A2 dimensional codes (d),
#               I2 F10.5 scaling factors (S, s),
#               F12.5 energy (E) if modeling program input,
#               I6 internal registry number (R) if input through MDL form.
# 3		Comment (blank if no comment)
#
# 4-17	Connection table (Ctab)
# 4		Counts line: 6 atoms, 6 bonds, ..., V2000 standard
#              aaabbblllfffcccsssxxxrrrpppiiimmmvvvvvv
#              Where::
#              aaa       = number of atoms (current max 255)*        [Generic]
#              bbb       = number of bonds (current max 255)*        [Generic]
#              lll       = number of atom lists (max 30)*            [Query]
#              fff       = (obsolete)
#              ccc       = chiral flag: 0=not chiral, 1=chiral       [Generic]
#              sss       = number of stext entries                   [ISIS/Desktop]
#              xxx       = (obsolete)
#              rrr       = (obsolete)
#              ppp       = (obsolete)
#              iii       = (obsolete)
#              mmm       = number of lines of additional properties, [Generic]
#                        including the M END line. No longer
#                        supported, the default is set to 999.

# 5-10		Atom block (1 line for each atom): x, y, z, element, etc.
# 11-16	Bond block (1 line for each bond): 1st atom, 2nd atom, type, etc.
# 17		Properties block (empty)
# 18	$$$$	See note
#

from MolKit.molecule import Atom, Bond
# About blank lines
#
# - Only one blank line should terminate a data item.
# - There should only be one blank line between the last data item and the
#   $$$$ delimiter line.
# - If the SDfile only contains structures, there can be no blank line between
#   the last "M END" and the $$$$ delimiter line.
from MolKit.moleculeParser import MoleculeParser
from MolKit.protein import Protein, Residue, Chain


class SDFParser(MoleculeParser):

    def __init__(self, filename=None, allLines=None, modelsAs='molecules'):
        """Constructor for sdfParser"""
        MoleculeParser.__init__(self, filename, allLines)
        self.model = False
        self.modelsAs = modelsAs

    def readFile(self):
        f = open(self.filename)
        self.allLines = f.readlines()
        f.close()

    def parse(self, objClass=Protein):
        if self.allLines is None and self.filename:
            self.readFile()
            if self.allLines is None or len(self.allLines) == 0:
                return
        lines = self.allLines

        # index molecules in file
        # molIndex containd the line number of first line for each molecule
        molIndex = self.molIndex = [0]
        molNames = self.molNames = [lines[0].strip()]
        lind = 1
        nbLines = len(lines)
        while True:
            if lind >= nbLines - 1:
                break
            if lines[lind][0:4] == "$$$$":
                lind += 1
                molIndex.append(lind)
                molNames.append(lines[lind].strip())
                lind += 3
            lind += 1
        # print "  molIndex:", molIndex, "molNames;", molNames
        return self.getMolecule(0)

    def getMolIndex(self, molname):
        return self.molNames.index(molname)

    def getMolName(self, molind):
        return self.molNames[molind]

    def numOfMolecules(self):
        return len(self.molNames)

    def getMolecule(self, molInd):

        molecules = []
        if molInd == len(self.molIndex) - 1:
            lastLine = -1
        else:
            lastLine = self.molIndex[molInd + 1]
        # lines fotr that molecule
        lines = self.allLines[self.molIndex[molInd]:lastLine]
        lineIndex = 0
        atomsSeen = {}  # dict of atom types and number of atoms seen

        # parser header
        molName = lines[lineIndex].strip()
        lineIndex += 3

        # create molecule
        mol = Protein(name=molName)
        mol.info = lines[lineIndex + 1]
        mol.comment = lines[lineIndex + 1]
        # self.mol.parser = self
        chain = Chain(id='1', parent=mol, top=mol)
        res = Residue(type='UNK', number='1', parent=chain, top=mol)
        mol.levels = [Protein, Chain, Residue, Atom]

        # parse count line
        line = lines[lineIndex]
        assert line[33:39] == " V2000", "Format error: only V2000 is suported, got %s" % line[33:39]
        nba = int(line[0:3])  # number of atoms
        nbb = int(line[3:6])  # number of bonds
        nbal = int(line[6:9])  # number of atom lists
        ccc = int(line[12:15])  # chiral flag: 0=not chiral, 1=chiral
        sss = int(line[15:18])  # number of stext entries
        lineIndex += 1

        # parse atoms
        for anum in range(nba):
            line = lines[lineIndex]
            element = line[31:34].strip()
            if element in atomsSeen:
                atomsSeen[element] += 1
            else:
                atomsSeen[element] = 1
            atom = Atom(name='%s_%s' % (element, atomsSeen[element]), parent=res,
                        chemicalElement=element, top=mol)

            atom._coords = [[float(line[0:10]), float(line[10:20]),
                             float(line[20:30])]]
            atom._charges['sdf'] = int(line[35:38])
            atom.chargeSet = 'sdf'
            mol.allAtoms.append(atom)

            atom.massDiff = int(line[34:36])
            atom.stereo = int(line[38:41])
            atom.hcount = line[41:44]
            atom.valence = int(line[47:50])
            atom.hetatm = 1
            atom.occupancy = 0.0
            atom.temperatureFactor = 0.0
            lineIndex += 1

        # parse bonds
        for bnum in range(nba):
            line = lines[lineIndex]
            at1 = mol.allAtoms[int(line[0:3]) - 1]
            at2 = mol.allAtoms[int(line[3:6]) - 1]
            if at1.isBonded(at2):
                continue
            bond = Bond(at1, at2, check=0)

            bond.bondOrder = int(line[6:9])
            # 1 = Single, 2 = Double,
            # 3 = Triple, 4 = Aromatic,
            # 5 = Single or Double,
            # 6 = Single or Aromatic,
            # 7 = Double or Aromatic, 8 = Any

            bond.stereo = int(line[9:12])
            # Single bonds: 0 = not stereo,
            # 1 = Up, 4 = Either,
            # 6 = Down, Double bonds: 0 = Use x-, y-, z-coords
            # from atom block to determine cis or trans,
            # 3 = Cis or trans (either) double bond

            bond.topo = int(line[15:18])
            # 0 = Either, 1 = Ring, 2 = Chain

            try:
                bond.ReactionCenter = int(line[18:21])
            except ValueError:
                bond.ReactionCenter = 0
            # 0 = unmarked, 1 = a center, -1 = not a center,
            # Additional: 2 = no change,
            # 4 = bond made/broken,
            # 8 = bond order changes
            # 12 = 4+8 (both made/broken and changes);
            # 5 = (4 + 1), 9 = (8 + 1), and 13 = (12 + 1)

        # "M END" and properties are not parsed at this point
        self.mol = mol
        mname = mol.name
        strRpr = mname + ':::'
        mol.allAtoms.setStringRepr(strRpr)
        strRpr = mname + ':'
        mol.chains.setStringRepr(strRpr)
        for c in mol.chains:
            cname = c.id
            strRpr = mname + ':' + cname + ':'
            c.residues.setStringRepr(strRpr)
            for r in c.residues:
                rname = r.name
                strRpr = mname + ':' + cname + ':' + rname + ':'
                r.atoms.setStringRepr(strRpr)
        molList = mol.setClass()
        molList.append(mol)
        mol.parser = self
        for n in molList.name:
            name = n + ','
        name = name[:-1]
        molList.setStringRepr(name)
        strRpr = name + ':::'
        molList.allAtoms.setStringRepr(strRpr)

        return molList

    def getMoleculeInformation(self):
        """ Function to retrieve the general informations on the molecule.
        This information is used by the molecule chooser to provide
        informations on the molecule selected.
        """
        molStr = ''
        return molStr

    def configureProgressBar(self, **kw):
        # this method is to be implemented by the user from outside
        pass

    def hasSsDataInFile(self):
        """ Function testing if the informations on the secondary structure
        are in the file"""
        return 0

    # from MolKit.sdfParser import SDFParser
# pr = SDFParser("ZINC_results2.sdf")
# mol = pr.parse()
