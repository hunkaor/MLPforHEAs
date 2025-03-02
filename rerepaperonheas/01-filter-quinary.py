import ase
from ase.db import connect
from ase.io import read
import numpy as np
from icet import (ClusterSpace,StructureContainer)
from icet.tools.structure_mapping import *

###################  method 3 to add sc  #########################
prim = ase.io.read("../MB2.vasp",format='vasp')
species = [['Ti','Zr', 'Hf', 'Ta', 'Nb'], ['B'], ['B']]
cs = ClusterSpace(prim, [7.5,5.8,4.6], species)    #N=10



#ll *C*/DB/*.db | awk '{print "+read("$9"@calculator=vasp)"}'
a=read('TiZrHfTaNb_enum0.db@Nb,Ti,Zr,Hf,Ta,B>10')

db0=connect('TiZrHfTaNb_n24+_enum00.db')
db=connect('TiZrHfTaNb_n24+_enum0.db')

for i in range(len(a)):
    db0.write(a[i])



id2rm=[]
sc = StructureContainer(cluster_space=cs)
for row in db0.select():
    try: atoms,info=map_structure_to_reference(row.toatoms(), prim)
    except ValueError: pass
    try: sc.add_structure(atoms, allow_duplicate=False,user_tag=str(row.id))
    except ValueError:
        pass  #pass/skip the duplicate structures
        id2rm.append(row.id)


for row in db0.select():
    if row.id not in id2rm:
       db.write(row.toatoms())

