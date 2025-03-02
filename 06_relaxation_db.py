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
        # print(atoms.get_forces())

        if optimizer.converged():
            break
        
        if optimizer.nsteps == optimizer.max_steps:
            return False # not converged

        optimizer.step()
        optimizer.nsteps += 1

    return True

output_dir = 'vasp_outputs/'
calc = Vasp(prec   = 'high',
            xc     = 'PBE',
            setups = {'Ti':'_sv','Zr':'_sv','Hf':'_sv', 'Nb':'_sv'},
            lreal  = 'False',
            encut  = 560.0,
            istart = 0,
            icharg = 2,
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
            ediff  = 1e-6,
            directory= output_dir,
            ncore=1,
            kpar=11
            )

input_db_path = 'db/TiZrHfTaNb.db'
output_db_path = input_db_path.split('.')[0] + '_relax.db'
print('output_path: %s' %output_db_path) 

# Remove old database file if it exists
if os.path.exists(output_db_path):
    os.remove(output_db_path)

input_db = connect(input_db_path)
output_db = connect(output_db_path)

input_index = 0
for row in input_db.select():
    atoms = row.toatoms()
    atoms.set_calculator(calc)

    print(input_index, atoms)

    opt = BFGS(atoms)
    opt.fmax = 0.01
    opt.max_steps = 100

    # initialize counters
    relax_index = 0

    run(optimizer = opt,
        atoms = atoms,
        configuration_index = input_index,
        calculator = calc,
        output_dir = output_dir,
        database = output_db,
        relax_type = 'positions',
        relax_index = relax_index 
        )

    relax_index +=1

    # -----------------------------------------

    mask = [True, True, True, True, True, True]
    ecf = StrainFilter(atoms, mask=mask)  # keep scaled positions fixed

    opt = BFGS(ecf)
    opt.fmax = 0.01
    opt.max_steps = 100

    run(optimizer = opt,
        atoms = atoms,
        configuration_index = input_index,
        calculator = calc,
        output_dir = output_dir,
        database = output_db,
        relax_type = 'cell',
        relax_index = relax_index
        )
    relax_index +=1

    input_index+=1