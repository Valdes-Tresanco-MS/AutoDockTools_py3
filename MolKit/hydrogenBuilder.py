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

############################################################################
#
# Author:  Ruth Huey
#
# Copyright: M. Sanner TSRI 2004
#
#############################################################################

"""
This module implements the HydrogenBuilder classes which add hydrogens to AtomSets.

"""

from MolKit.molecule import Atom, AtomSet, Bond
from PyBabel.addh import AddHydrogens
from PyBabel.aromatic import Aromatic
from PyBabel.atomTypes import AtomHybridization
from PyBabel.bo import BondOrder
from PyBabel.cycle import RingFinder


# this is for future use
# from AutoDockTools.observer import Subject


class HydrogenBuilder:
    """Base Class for adding hydrogen atoms to a molecule.

NB: test using 'withBondOrder' put hydrogens on aromatic ring of indinavir in
grossly incorrect positions AND added hydrogen to N1 which is also incorrect.
it is included in order to duplicate functionality of addHydrogens in Pmv/editCommands.py
but be warned that it apparently fails to add hydrogens correctly to cyclic
carbons.....
    """

    def __init__(self, htype='all', renumber=1, method='noBondOrder'):
        # NB: noBondOrder is for pdb files
        self.htype = 'all'
        self.renumber = renumber
        self.method = method

    def addHydrogens(self, mol):
        # check for bonds
        if len(mol.allAtoms.bonds[0]) == 0:
            mol.buildBondsByDistance()
        bonds = mol.allAtoms.bonds[0]
        # could have preset babel_types
        # so check if allAtoms are already typed
        try:
            t = mol.allAtoms.babel_type
        except:
            # if all are not pretyped, type them
            babel = AtomHybridization()
            babel.assignHybridization(mol.allAtoms)

        if self.method == 'withBondOrder':
            mol.rings = RingFinder()
            mol.rings.findRings2(mol.allAtoms, mol.allAtoms.bonds[0])
            mol.rings.bondRings = {}
            for ind in range(len(mol.rings.rings)):
                r = mol.rings.rings[ind]
                for b in r['bonds']:
                    if b not in mol.rings.bondRings:
                        mol.rings.bondRings[b] = [ind, ]
                    else:
                        mol.rings.bondRings[b].append(ind)
            bo = BondOrder()
            bo.assignBondOrder(mol.allAtoms, bonds, mol.rings)
            mol.allAtoms._bndtyped = 1
            # do aromatic here
            arom = Aromatic(mol.rings)
            arom.find_aromatic_atoms(mol.allAtoms)

        hat = AddHydrogens().addHydrogens(mol.allAtoms, method=self.method)
        bondedAtomDict = {}  # key is heavy atom
        for a in hat:
            if a[1] in bondedAtomDict:
                bondedAtomDict[a[1]].append(a)
            else:
                bondedAtomDict[a[1]] = [a]

        # now create Atom object for hydrogens
        # and add the to the residues's atom list
        molNewHs = AtomSet([])  # list of created H atoms for this molecule
        heavyAtoms = AtomSet([])  # list of atoms that need new radii

        for heavyAtom, HatmsDscr in list(bondedAtomDict.items()):
            # don't add hydrogens to carbons: polar Only!!!
            if self.htype != 'all' and heavyAtom.element == 'C':
                continue
            res = heavyAtom.parent
            # find where to insert H atom
            childIndex = res.children.index(heavyAtom) + 1

            # loop over H atoms description to be added
            # start at the end to number correctly
            l = len(HatmsDscr)
            for i in range(l - 1, -1, -1):
                a = HatmsDscr[i]
                # build H atom's name
                if len(heavyAtom.name) == 1:
                    name = 'H' + heavyAtom.name
                else:
                    name = 'H' + heavyAtom.name[1:]

                # if more than 1 H atom, add H atom index
                # for instance HD11, HD12, Hd13 (index is 1,2,3)
                if l > 1:
                    name = name + str(i + 1)

                # create the H atom object
                atom = Atom(name, res, top=heavyAtom.top,
                            chemicalElement='H',
                            childIndex=childIndex, assignUniqIndex=0)

                # set atoms attributes
                atom._coords = [a[0]]
                if hasattr(a[1], 'segID'):
                    atom.segID = a[1].segID
                atom.hetatm = 0
                atom.alternate = []
                # atom.element = 'H'
                atom.occupancy = 1.0
                atom.conformation = 0
                atom.temperatureFactor = 0.0
                atom.babel_atomic_number = a[2]
                atom.babel_type = a[3]
                atom.babel_organic = 1
                atom.radius = 1.2

                # create the Bond object bonding Hatom to heavyAtom
                bond = Bond(a[1], atom, bondOrder=1)

                # add the created atom the the list
                molNewHs.append(atom)
                # in case this new hydrogen atom ever ends up in pmv 
                # HAVE TO CREATE THESE ENTRIES 
                # create the color entries for all geoemtries
                # available for the heavyAtom
                for key, value in list(heavyAtom.colors.items()):
                    atom.colors[key] = (0.0, 1.0, 1.0)
                    atom.opacities[key] = 1.0

        mol.allAtoms = mol.chains.residues.atoms
        if self.renumber:
            mol.allAtoms.number = list(range(1, len(mol.allAtoms) + 1))
        return len(molNewHs)


class PolarHydrogenBuilder(HydrogenBuilder):
    """Base Class for adding hydrogen atoms to a molecule.
    NB: molecule must have bonds built
    """

    def __init__(self, htype='polarOnly', renumber=1, method='noBondOrder'):
        HydrogenBuilder.__init__(self, htype=htype, renumber=renumber,
                                 method=method)
        self.htype = 'polarOnly'
