from ase.db import connect
import json

db = connect('results/test.db')

# row = db.get(1)  # Retrieve the first entry
rows = db.select(configuration_index=0)  # retrieve by key

for row in rows:
    calc_params_loaded = json.loads(row.parameters)  # Convert JSON back to dictionary
    atoms = row.toatoms()

    # print(atoms)
    # print(calc_params_loaded)
    print(atoms.get_potential_energy())
    print(atoms.get_forces())
    print(atoms.get_stress())

    # print(row['parameters'])
    print(row['relax_type'])
    print(row['relax_iteration'])
    # print(row['log_msg'])
    # print(row['calculation_status'])

    # print(row.data['OUTCAR'])