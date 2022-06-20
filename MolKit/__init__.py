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
#  Modification date: 8/7/20 15:11                                                                 #
#                                                                                                  #
# ##################################################################################################

#############################################################################
#
# Author: Michel F. SANNER, Sophie COON
#
# Copyright: M. Sanner TSRI 2000
#
#############################################################################

import os

from MolKit.mmcifParser import MMCIFParser
from MolKit.mol2Parser import Mol2Parser
from MolKit.pdbParser import PdbParser, PdbqParser, PdbqsParser, PdbqtParser, PQRParser, F2DParser


def Read(filename=None, alllines=None, dataformat='pdb', modelsAs='molecules'):
    if not alllines and not filename:
        raise AssertionError("%s invalid molecule" % filename)
    elif alllines and filename:
        alllines = None
    elif filename:
        if not os.path.exists(filename):
            raise AssertionError("%s does't exist" % filename)

    if alllines:
        ext = dataformat
        args = [None, alllines]
    else:
        ext = filename.split('.')[-1]
        args = [filename, None]
    if ext == 'pdb':
        parser = PdbParser(*args, modelsAs=modelsAs)
    elif ext == 'pdbq':
        parser = PdbqParser(*args, modelsAs=modelsAs)
    elif ext == 'pdbqt':
        parser = PdbqtParser(*args, modelsAs=modelsAs)
    elif ext == 'pdbqs':
        parser = PdbqsParser(*args, modelsAs=modelsAs)
    elif ext == 'pqr':
        parser = PQRParser(*args, modelsAs=modelsAs)
    # FIXME: pass all lines???
    elif ext == 'mol2':
        parser = Mol2Parser(filename)  # ??should modelsAs be available for mol2 format??
    elif ext == 'cif':
        parser = MMCIFParser(filename)
    elif ext == 'f2d':
        parser = F2DParser(filename)
    else:
        print("File Format unknown can't parse it")
        return []
    molecules = parser.parse()
    return molecules


def WritePDB(filename, node):
    from MolKit.pdbWriter import PdbWriter
    writer = PdbWriter()
    writer.write(filename, node)


def makeMoleculeFromAtoms(molname, atomSet):
    """
    create a new molecule from a list of atoms

    mol <- makeMoleculeFromAtoms(molname, atomSet)
    """
    from MolKit.molecule import Atom, AtomSet
    from MolKit.protein import Protein, Chain, Residue

    # create the top object
    mol = Protein(name=molname)

    # find out all residues
    residues = atomSet.parent.uniq()

    # find out all chains
    chains = residues.parent.uniq()

    # create all chains
    chainsd = {}
    for c in chains:
        newchain = Chain(c.id, mol, top=mol)
        chainsd[c] = newchain

    # create all residues
    resd = {}
    for res in residues:
        if res.icode == '':
            rnum = res.name[3:]
        else:
            rnum = res.name[3:-1]
        newres = Residue(res.name[:3], rnum, res.icode,
                         chainsd[res.parent], top=mol)
        resd[res] = newres
        newres.hasCA = 0
        newres.hasO = 0
        newres.CAatom = None
        newres.Oatom = None
        newres.C1atom = None

    # create all the atoms
    newats = []
    for num, at in enumerate(atomSet):
        name = at.name
        res = resd[at.parent]
        name1 = name
        if hasattr(at, "altname") and at.altname:
            name = at.name.split("@")[0]

        newat = Atom(name, res, at.element, top=mol)
        if name == 'CA' or name[:3] == 'CA@':
            res.hasCA = 1
            res.CAatom = newat
        elif name == 'O' or name == 'OXT' or (len(name) > 3 and name[:3] == 'OCT'):
            res.hasO = 2
            res.Oatom = newat
        elif name == 'C1*':
            res.C1atom = newat
        if name != name1:
            newat.name = name1
            newat.altname = at.altname
            # Check this: MolKit/pdbParser.py/parse_PDB_ATOM_record()
            # has code for atom.alternate list. Not sure if we need
            # to do it here:
            # newatom.alternate --- ???
        newats.append(newat)
        # set constructotr attributes
        newat._coords = []
        for coords in at._coords:
            newat._coords.append(coords[:])
        newat.conformation = at.conformation
        newat.chemElem = at.chemElem
        newat.atomicNumber = at.atomicNumber
        newat.bondOrderRadius = at.bondOrderRadius
        newat.covalentRadius = at.covalentRadius
        newat.vdwRadius = at.vdwRadius
        newat.maxBonds = at.maxBonds
        newat.organic = at.organic
        newat.colors = at.colors.copy()
        newat.opacities = at.opacities.copy()
        newat._charges = at._charges.copy()
        newat.chargeSet = at.chargeSet

        # set attributes from PDB parser
        try:  # pdbqs do not have this
            newat.segID = at.segID
        except AttributeError:
            pass
        newat.hetatm = at.hetatm
        try:  # pdbqs do not have this
            newat.normalname = at.normalname
        except AttributeError:
            pass
        newat.number = num  # at.number
        newat.occupancy = at.occupancy
        newat.temperatureFactor = at.temperatureFactor
        newat.altname = at.altname

        # attribute created by PQR parser
        if hasattr(at, 'pqrRadius'):
            newat.pqrRadius = at.pqrRadius

        # attribute created by F2D parser
        if hasattr(at, 'hbstatus'):
            newat.hbstatus = at.hbstatus

        # attribute created by PDBQ parser
        if hasattr(at, 'autodock_element'):
            newat.autodock_element = at.autodock_element

        # attribute created by PDBQT parser
        # if hasattr(at, ''):
        #    newat. = at.

        # attribute created by PDBQS parser
        if hasattr(at, 'AtVol'):
            newat.AtVol = at.AtVol
            newat.AtSolPar = at.AtSolPar

    mol.allAtoms = AtomSet(newats)
    return mol
