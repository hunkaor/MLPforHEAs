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
from ase.build import bulk
import matplotlib.pyplot as plt
import json

output_dir = 'vasp_outputs/'
calc = Vasp(prec   = 'Accurate',
            xc     = 'PBE',
            # setups = {'Ti':'_sv','Zr':'_sv','Hf':'_sv'},
            lreal  = 'False',
            # encut  = 560.0,
            encut = 250.0, # ! dont forget to change it back
            istart = 0,
            icharg = 2,
            # ncore  = 28,
            # nbands = 98,  # 8*ncores
            nelm  = 80,                  #without exit if no convergence reached. continue with wrong WF
            nelmin = 5,
            kspacing  = 0.5, # ! dont forget to change it back
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
            ediff  = 1e-6, 
            # command="mpirun -np 1 vasp_std"
            directory= output_dir,
            )

output_name = 'Pd'
atoms = bulk('Pd')
atoms = atoms * (2,1,1)
print(atoms.positions)

atoms.positions += [[0.1, 0.2, 0.3], [-0.1, 0, 0.2]]
atoms.set_calculator(calc)

# Apply random purturbation on atomic positions
# perturbation_magnitude = 0.2  # Adjust as needed
# atoms.positions += np.random.uniform(-perturbation_magnitude, 
#                                      perturbation_magnitude, 
#                                      atoms.positions.shape)
# atoms.positions += np.array(-perturbation_magnitude, 
#                                      perturbation_magnitude, 
#                                      atoms.positions.shape)

# # Apply small random strain to the cell
# cell_perturbation = 0.2 # Fractional change (2%)
# cell = atoms.cell.copy()  # Get current cell
# strain_matrix = np.eye(3) + np.random.uniform(-cell_perturbation, 
#                                               cell_perturbation, 
#                                               (3, 3))  # Small random strain
# atoms.set_cell(np.dot(cell, strain_matrix), scale_atoms=True)  # Apply strai

# atoms.set_calculator(calc)
# atoms.get_potential_energy()

# -------------------------------------------------------------------------------------------------------

def check_convergence(outcar_content):
    """
    Check whether the VASP converged within EDIFF.
    """
    EDIFF_keyword = "EDIFF  ="
    dE_keyword = "total energy-change (2. order) :"

    EDIFF_val = 0
    dE_val = 100000
    for line in outcar_content.splitlines():
        if EDIFF_keyword in line:
            EDIFF_val = float(line.split("=")[-1].split("  ")[0])

        if dE_keyword in line:
            dE_val = abs(float(line.split(":")[-1].split("  ")[0]))
            if dE_val < EDIFF_val:
                return True
    return False

def check_vasp_status(outcar_content):
    """
    Check whether the VASP calculation in OUTCAR was successful.
    """
    error_keywords = ["ERROR", "ZBRENT: fatal error"]
    success_keywords = [
        "aborting loop because EDIFF is reached",
        "General timing and accounting informations for this job"
    ]

    failed = any(keyword in outcar_content for keyword in error_keywords)
    if failed:
        return "FAILED"
    
    if check_convergence(outcar_content) == False:
        return "NOT_CONVERGED"

    success = any(keyword in outcar_content for keyword in success_keywords)
    if success:
        return "SUCCESS"
    else:
        return "INCOMPLETE"
    
def write_to_db(database, atoms, output_dir, calculator, configuration_index, 
                relax_index, relax_type, relax_iteration, log_msg):
    
    with open(output_dir + 'OUTCAR', "r") as f:
        outcar_content = f.read()

    calculation_status = check_vasp_status(outcar_content)

    database.write(atoms, 
            parameters = json.dumps(calculator.todict()),
            calculation_status = calculation_status,
            configuration_index = configuration_index,
            relax_index = relax_index,
            relax_type = relax_type,
            relax_iteration = relax_iteration,
            log_msg = log_msg,
            OUTCAR = outcar_content)
    
def run(optimizer, atoms, configuration_index, calculator, output_dir, database, relax_type, relax_index):
    while True:
        msg = optimizer.log()
        if optimizer.nsteps == 0 and relax_index > 0: # skip the last configuration from previous relaxation
            pass
        else:
            write_to_db(database = database,
                        atoms = atoms,
                        output_dir = output_dir,
                        calculator = calculator,
                        configuration_index = configuration_index,
                        relax_index = relax_index,
                        relax_type = relax_type,
                        relax_iteration = optimizer.nsteps,
                        log_msg = msg
                        )
        print(atoms.get_forces())

        if optimizer.converged():
            break
        
        if optimizer.nsteps == optimizer.max_steps:
            return False # not converged

        optimizer.step()
        optimizer.nsteps += 1

    return True

db_path = 'results/test.db'

# Remove old database file if it exists
if os.path.exists(db_path):
    os.remove(db_path)

i = 0
db = connect('results/test.db')


# with open(output_dir + 'OUTCAR', "r") as f:
#     outcar_content = f.read()

# calculation_status = check_vasp_status(outcar_content)

# db.write(atoms, 
#          parameters = json.dumps(calc.todict()),
#          calculation_status = True,
#          configuration_index = i,
#          relax_index = 0,
#          relax_type = 'position',
#          relax_iteration = 5,
#          fully_relax = False,
#          data={"OUTCAR": outcar_content})  

# # For how to set/not-set ibrion
# # https://matsci.org/t/problem-with-using-vasp-with-ase-to-predict-the-properties-of-materials/54989/2

# --------------------------------------------------------------------------------------------------------

# atoms.write('vasp_outputs/atom.cif')
# atoms.write('vasp_outputs/atom.xyz')

# print(atoms.get_forces())

print('\n\n######## Start relaxing positions ##########')
opt = BFGS(atoms)
opt.fmax = 0.01
opt.max_steps = 100
# opt.run(fmax=0.01, steps=100)

relax_index = 0
run(optimizer = opt,
    atoms = atoms,
    configuration_index = 0,
    calculator = calc,
    output_dir = output_dir,
    database = db,
    relax_type = 'positions',
    relax_index = relax_index 
    )

relax_index +=1
# ---------------------------------------------------------------------------------------------------------------

print('\n\n######## Start relaxing cell ##########')
mask = [True, True, True, True, True, True]
ecf = StrainFilter(atoms, mask=mask)  # keep scaled positions fixed

opt = BFGS(ecf)
opt.fmax = 0.01
opt.max_steps = 100
# dyn.run(fmax=0.01, steps=100)


run(optimizer = opt,
    atoms = atoms,
    configuration_index = 0,
    calculator = calc,
    output_dir = output_dir,
    database = db,
    relax_type = 'cell',
    relax_index = relax_index
    )
relax_index +=1

# ---------------------------------------------------------------------------------------------------------------


# create class - position_optimizer, cell_optimizer