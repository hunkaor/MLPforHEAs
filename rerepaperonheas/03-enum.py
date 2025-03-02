import ase
from random import randint
from ase.io import *
from ase import Atom
from ase.db import connect
from ase.build import bulk
from icet.tools import enumerate_structures
from icet import ClusterSpace
from icet.tools.structure_generation import (generate_sqs,
                                             generate_sqs_by_enumeration,
                                             generate_target_structure)
#from icet.io.logging import set_log_config
#set_log_config(level='INFO')

#prim = ase.io.read("Zr2r.vasp",format='vasp')
prim = ase.io.read("bcc.vasp",format='vasp')


conc = {'Zr':(0.32,0.34),'Ti':(0.32,0.34)}
db = connect('DB/TiZrNb_enum0.db')
for structure in enumerate_structures(prim, range(6, 7), ['Ti','Zr', 'Nb'], concentration_restrictions=conc):
    db.write(structure)

conc = {'Hf':(0.32,0.34),'Ti':(0.32,0.34)}
db = connect('DB/TiHfNb_enum0.db')
for structure in enumerate_structures(prim, range(6, 7), ['Ti','Hf', 'Nb'], concentration_restrictions=conc):
    db.write(structure)

conc = {'Ta':(0.32,0.34),'Ti':(0.32,0.34)}
db = connect('DB/TiTaNb_enum0.db')
for structure in enumerate_structures(prim, range(6, 7), ['Ti','Ta', 'Nb'], concentration_restrictions=conc):
    db.write(structure)

conc = {'Zr':(0.32,0.34),'Hf':(0.32,0.34)}
db = connect('DB/ZrHfNb_enum0.db')
for structure in enumerate_structures(prim, range(6, 7), ['Hf','Zr', 'Nb'], concentration_restrictions=conc):
    db.write(structure)

conc = {'Zr':(0.32,0.34),'Ta':(0.32,0.34)}
db = connect('DB/ZrTaNb_enum0.db')
for structure in enumerate_structures(prim, range(6, 7), ['Ta','Zr', 'Nb'], concentration_restrictions=conc):
    db.write(structure)

conc = {'Ta':(0.32,0.34),'Hf':(0.32,0.34)}
db = connect('DB/HfTaNb_enum0.db')
for structure in enumerate_structures(prim, range(6, 7), ['Hf','Ta', 'Nb'], concentration_restrictions=conc):
    db.write(structure)



conc = {'Zr':(0.25,0.26),'Hf':(0.25,0.26),'Ta':(0.25,0.26)}
db = connect('DB/ZrHfTaNb_enum0.db')
for structure in enumerate_structures(prim, [2,4], ['Ta','Zr', 'Hf', 'Nb'], concentration_restrictions=conc):
    db.write(structure)

conc = {'Ti':(0.25,0.26),'Hf':(0.25,0.26),'Ta':(0.25,0.26)}
db = connect('DB/TiHfTaNb_enum0.db')
for structure in enumerate_structures(prim, [2,4], ['Ti','Hf', 'Ta', 'Nb'], concentration_restrictions=conc):
    db.write(structure)

conc = {'Ti':(0.25,0.26),'Zr':(0.25,0.26),'Ta':(0.25,0.26)}
db = connect('DB/TiZrTaNb_enum0.db')
for structure in enumerate_structures(prim, [2,4], ['Ti','Zr', 'Ta', 'Nb'], concentration_restrictions=conc):
    db.write(structure)

conc = {'Ti':(0.25,0.26),'Zr':(0.25,0.26),'Hf':(0.25,0.26)}
db = connect('DB/TiZrHfNb_enum0.db')
for structure in enumerate_structures(prim, [2,4], ['Ti','Zr', 'Hf', 'Nb'], concentration_restrictions=conc):
    db.write(structure)


conc = {'Ti':(0.2,0.21),'Zr':(0.2,0.21),'Hf':(0.2,0.21),'Ta':(0.2,0.21)}
db = connect('DB/TiZrHfTaNb_enum10.db')
for structure in enumerate_structures(prim, [5], ['Ti', 'Ta','Zr', 'Hf', 'Nb'], concentration_restrictions=conc):
    db.write(structure)



