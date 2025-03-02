import ase
from random import randint
from ase.io import *
from ase import Atom
from ase.db import connect
from ase.build import bulk
from icet.tools import enumerate_structures
from icet import ClusterSpace
from ase.build import bulk
from icet.tools.structure_generation import (generate_sqs,
                                             generate_sqs_by_enumeration,
                                             generate_target_structure)

db = connect('db/TiZrHfTaNb.db')
prim = bulk('Ta', 'bcc', a = 3.3087917837241765)
prim.write('Ta_bcc_prim.cif')
# prim = ase.io.read("rerepaperonheas/bcc.vasp",format='vasp')
# print(prim)

c=0
for structure in enumerate_structures(structure=prim,
                                      sizes=range(1, 6),
                                      chemical_symbols=['Ti', 'Zr', 'Hf', 'Ta', 'Nb']):
    print(c, structure.get_number_of_atoms(), structure)
    db.write(structure)
    structure.write('cif/' + str(c).zfill(3) + '.cif')
    c+=1