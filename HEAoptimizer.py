import numpy as np
from typing import List, Dict, Optional, TextIO
from ase import Atoms
from ase.db import connect
from ase.calculators.calculator import Calculator
from ase.optimize.optimize import Optimizer
from ase.optimize import BFGS
from ase.constraints import StrainFilter
import json
import os
from ase.calculators.vasp import Vasp

class HEAoptimizer:
    """
    A class to perform structure optimizations (from ase-atoms database) using VASP, and log calculation detail and results to a new database.

    Attributes:
        input_path (str): Path to the input database containing ase-atoms.
        vasp_params: (Dict): Dictionary containing vasp settings, for example {'prec':'high', 'xc':'PBE'}.
        output_path (str): Path of the ouput database, where relaxed trajectories are logged (default = input_path + "_relax.db")
        vasp_output_dir (str): Directory where vasp files will be stored (default = "/vasp_outputs")
        start_index (int): Starting index of the ase-atoms (first row is 1) for the optimization (default = 1)
        end_index (int): End index of the ase-atoms for the optimization (default = -1)
        fmax: minimum forces for the relaxation (default=0.01)
        max_steps: maximum steps for the relaxation (default=100)
        start_from_scratch (bool): Remove an exsiting output database and create an emtpy one (default = True)
    """
    def __init__(self, 
                 input_path: str, 
                 vasp_params: Dict,
                 output_path: str = None,
                 vasp_output_dir: str = None,
                 start_index: int = 1, 
                 end_index: int = - 1,
                 fmax: float = 0.01,
                 max_steps: int = 100,
                 start_from_scratch: bool = True
                 ):
        
        self.input_db = connect(input_path)

        if output_path == None:
            self.output_path = input_path.split('.')[0] + '_relax.db'
        else:
            self.output_path = output_path
        print("The result will be written to: %s" %self.output_path)

        if start_from_scratch:
            if os.path.exists(self.output_path):
                os.remove(self.output_path)

        self.output_db = connect(self.output_path)
        # -------------
        if vasp_output_dir == None:
            self.vasp_output_dir = 'vasp_outputs/'

        os.makedirs(self.vasp_output_dir, exist_ok=True)
        vasp_params['directory'] = self.vasp_output_dir
        # -------------
        self.calculator = Vasp(**vasp_params)

        self.start_index = start_index
        if end_index == -1:
            self.end_index = len(self.input_db)+1
        else:
            self.end_index = end_index

        print("Perform relaxation on configurations %d:%d \n\n" %(start_index, end_index))

        self.fmax = fmax
        self.max_steps = max_steps

    def check_convergence(self, outcar_content: TextIO):
        """
        Check whether the VASP converged within EDIFF.

        Args:
        outcar_content: A file-like object loaded from OUTCAR.

        Returns:
            bool: True if converged, False otherwise
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

    def check_vasp_status(self, outcar_content: TextIO):
        """
        Check whether the VASP calculation in the OUTCAR was successful.

        Args:
        outcar_content: A file-like object loaded from OUTCAR.

        Returns:
            str: The status of the OUTCAR file: 
                - "FAILED" if any error keyword is found,
                - "NOT_CONVERGED" if convergence is not reached,
                - "SUCCESS" if any success keyword is found,
                - "INCOMPLETE" if none of the above conditions are met.
        """
        error_keywords = ["ERROR", "ZBRENT: fatal error"] # add more if needed
        success_keywords = [
            "aborting loop because EDIFF is reached",
            "General timing and accounting informations for this job"
        ]

        failed = any(keyword in outcar_content for keyword in error_keywords)
        if failed:
            return "FAILED"
        
        if self.check_convergence(outcar_content) == False:
            return "NOT_CONVERGED"

        success = any(keyword in outcar_content for keyword in success_keywords)
        if success:
            return "SUCCESS"
        else:
            return "INCOMPLETE"
        
    def write_to_db(self, info_dict: Dict, calculation_status: str):
        """
        Write updated atoms and related metadata into the database.

        Args:
        info_dict: dictionary containing key/value needed to be written to a database.
        calculation_status: check "check_vasp_status" module

        Return: None
        """
        outcar_content = ""
        if calculation_status != "FAILED":
            try:
                with open(self.vasp_output_dir + '/OUTCAR', "r") as f:
                    outcar_content = f.read()
            except:
                calculation_status = "INCOMPLETE"
            
            if outcar_content != "":
                calculation_status = self.check_vasp_status(outcar_content)

        info_dict['OURCAR'] = outcar_content
        info_dict['calculation_status'] = calculation_status
        info_dict['vasp_params'] = json.dumps(self.calculator.todict())

        self.output_db.write(**info_dict)
        
    def optimize(self, 
                 atoms: Atoms,
                 configuration_index: int,
                 relax_type: str,
                 relax_index: int,
                 ):
        '''
            Optimize ase atoms and log all iterations into the database.

        Args:
        atoms: Atom object.
        configuration_index: current index of the atoms
        relax_type: "positions" or "cell"
        relax_index: index of the optimization starting from 0

        Returns: True if success, False otherwise
        '''
        if relax_type == "positions":
            optimizer = BFGS(atoms)
        elif relax_type == "cell":
            mask = [True, True, True, True, True, True]
            ecf = StrainFilter(atoms, mask=mask)  # keep scaled positions fixed
            optimizer = BFGS(ecf)
        else:
            raise ValueError(f"Unknown relax_type: {relax_type}. Expected 'positions' or 'cell'.")

        optimizer.fmax = self.fmax
        optimizer.max_steps = self.max_steps
        terminate_relaxation = False

        print('----- Relax %s, configuration: %d, relax_index: %d -----' %(relax_type, configuration_index, relax_index))
        while True:
            try:
                msg = optimizer.log() # ase start vasp calculation here
                calculation_status = "SUCCESS"
            except:
                msg = ""
                calculation_status = "FAILED"
                terminate_relaxation = True
                print('DFT calculation failed')

            info_dict = {
                "atoms": atoms,
                "configuration_index": configuration_index,
                "relax_type": relax_type,
                "relax_index": relax_index,
                "relax_iteration": optimizer.nsteps,
                "log_msg": msg
            }

            if relax_index > 0 and optimizer.nsteps == 0: # skip the last configuration from previous relaxation
                pass
            else:
                self.write_to_db(info_dict, calculation_status)

            if terminate_relaxation:
                return False # DFT failed, log first and terminate here

            if optimizer.converged():
                return True
            
            if optimizer.nsteps == optimizer.max_steps:
                return False # not converged

            optimizer.step()
            optimizer.nsteps += 1
        
    def check_relax_convergence(self, atoms: Atoms, relax_type: str):
        '''
            Check whether the relaxation is converged.

        Args:
        atoms: Atom object.
        relax_type: "positions" or "cell"

        Returns: True if converged, False otherwise
        '''
        if relax_type == "positions":
            optimizer = BFGS(atoms)
        elif relax_type == "cell":
            mask = [True, True, True, True, True, True]
            ecf = StrainFilter(atoms, mask=mask)  # keep scaled positions fixed
            optimizer = BFGS(ecf)
        else:
            raise ValueError(f"Unknown relax_type: {relax_type}. Expected 'positions' or 'cell'.")

        optimizer.fmax = self.fmax
        optimizer.max_steps = self.max_steps

        if optimizer.converged():
            return True
        else:
            return False
    
    # def get_optimizer_log_msg(self, optimizer):
        # to do : read log_msg without modifying the original ase code

    def run(self):
        """
            Run the main loops.
        """
        for i in range(self.start_index, self.end_index): # ase db has 1-based row indices
            row = self.input_db.get(id=i)
            atoms = row.toatoms()
            atoms.set_calculator(self.calculator)

            relax_index = 0
            while True:
                if self.optimize(atoms = atoms,
                                configuration_index = i,
                                relax_type = "positions",
                                relax_index = relax_index) == False:
                    break  # skip this configuration if encountered any error
                else:
                    relax_index+=1

                if self.check_relax_convergence(atoms = atoms, 
                                                relax_type = "cell") == True:
                    print('Cell is already converged')
                    break  # this calculation is done

                if self.optimize(atoms = atoms,
                            configuration_index = i,
                            relax_type = "cell",
                            relax_index = relax_index) == False:
                    break # skip this configuration if encountered any error
                else:
                    relax_index+=1

                if self.check_relax_convergence(atoms = atoms, 
                                                relax_type = "positions") == True:
                    print('Atomic positions are already converged')
                    break  # this calculation is done

            

