import os,sys, re
from ase.db import connect
from ase.io import *
from ase.io.trajectory import Trajectory as traj
from ase import Atoms
from ase.constraints import StrainFilter
from ase.constraints import UnitCellFilter
from ase.calculators.vasp import Vasp
from ase.optimize import BFGS
import numpy as np
from ase.units import GPa

db=connect('../TiZrHf_bcc.db')

# -------------------------------------------------------------
# Parameters in VASP
# -------------------------------------------------------------
calc = Vasp(prec   = 'high',
            xc     = 'PBE',
            setups = {'Ti':'_sv','Zr':'_sv','Hf':'_sv'},
            lreal  = 'False',
            encut  = 560.0,
            istart = 0,
            icharg = 2,
            ncore  = 28,
            nbands = 98,  # 8*ncores
            nelm   = 200,                  #without exit if no convergence reached. continue with wrong WF
            nelmin = 5,
            kspacing  = 0.1,
            gamma  = True,
            ismear = 1,
            sigma  = 0.16,
#            algo   = "Very_Fast",
            algo   = "Normal",
            lwave  = False,                 #bool_key
            lcharg = False,                 #bool_key
#           amix   = 0.1,
#           bmix   = 0.01,
            symprec = 1e-4,
            ediff  = 1e-6)
#calc.set(prec='Accurate', ediff=1E-5)

for row in db.select():
#    if row.id >= start and row.id <= end:
    if row.id == idd:
        atoms = row.toatoms()
        atoms.set_calculator(calc)
#        write('POSCAR'+str(row.id), read('POSCAR'),format='vasp')
#only relax unit cell until the stress(1-true) is zero,
        sf = StrainFilter(atoms,mask=[1,1,1,1,1,1])
        opt = BFGS(sf,trajectory=traj(str(row.id)+'_cell.traj','w',atoms))
        opt.run(fmax=0.10)
#simultaneously minimize the atomic forces and unit cell stresses
        ucf = UnitCellFilter(atoms,mask=[1,1,1,1,1,1]) # 1-stress tensor relax to zero, 
        opt = BFGS(ucf, trajectory=traj(str(row.id)+'_atom.traj','w',atoms))
        opt.run(fmax=0.01)
        db.update(row.id, atoms)
        write('CONTCAR'+str(row.id), read('CONTCAR'),format='vasp')


