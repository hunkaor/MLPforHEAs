import numpy as np

output_dir = 'vasp_outputs/'
with open(output_dir + 'OUTCAR', "r") as f:
    outcar_content = f.read()

EDIFF_keyword = "EDIFF  ="
dE_keyword = "total energy-change (2. order) :"

EDIFF_val = 0
dE_val = 100000
for line in outcar_content.splitlines():
    if EDIFF_keyword in line:
        EDIFF_val = float(line.split("=")[-1].split("  ")[0])
        print(EDIFF_val)

    if dE_keyword in line:
        dE_val = abs(float(line.split(":")[-1].split("  ")[0]))
        print(dE_val)

        