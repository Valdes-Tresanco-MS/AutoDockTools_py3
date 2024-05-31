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
#  Modification date: 15/06/20, 8:27 p. m.                                                         #
#                                                                                                  #
# ##################################################################################################

from setuptools import setup
import versioneer

setup(
    name='AutoDockTools_py3',
    version = versioneer.get_version(),
    cmdclass = versioneer.get_cmdclass(),
    packages=['MolKit', 'MolKit.data', 'MolKit.pdb2pqr', 'MolKit.pdb2pqr.src', 'MolKit.pdb2pqr.propka',
              'MolKit.pdb2pqr.pdb2pka', 'MolKit.pdb2pqr.pdb2pka.substruct', 'MolKit.pdb2pqr.pdb2pka.ligandclean',
              'MolKit.pdb2pqr.extensions', 'PyBabel', 'mglutil', 'mglutil.math', 'mglutil.util',
              'AutoDockTools', 'AutoDockTools.Utilities24'],
    url='https://github.com/Valdes-Tresanco-MS/AutoDockTools_py3',
    license='MGLTools LICENSE',
    author='MS. Valdes-Trasanco and ME. Valdes-Tresanco ',
    author_email='bioinfobrothers@gmail.com',
    description='Translation of ADT to python3.x',
    entry_points = {
    'console_scripts': [
        'prepare_ligand4 = AutoDockTools.Utilities24.prepare_ligand4:main',
        'prepare_receptor4 = AutoDockTools.Utilities24.prepare_receptor4:main',
        'prepare_gpf4 = AutoDockTools.Utilities24.prepare_gpf4:main',
        'prepare_dpf42 = AutoDockTools.Utilities24.prepare_dpf42:main',
        'prepare_flexreceptor4 = AutoDockTools.Utilities24.prepare_flexreceptor4:main',
        'prepare_covalent_flexres = AutoDockTools.Utilities24.prepare_covalent_flexres:main',
        'AutoLigand = AutoDockTools.AutoLigand:main',
    ],
},
)
