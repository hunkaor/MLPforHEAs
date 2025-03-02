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

db = connect('db/PdNi.db')
prim = bulk('Pd')
print(prim)

c=0
for structure in enumerate_structures(structure=prim,
                                      sizes=range(1, 6),
                                      chemical_symbols=['Pd', 'Ni', 'Al', 'Si', 'Zr']):
    # print(c, structure.get_number_of_atoms(), structure)
    db.write(structure)
    # structure.write('cif/' + str(c).zfill(3) + '.cif')
    c+=1