from HEAoptimizer import HEAoptimizer

vasp_params = {
    'prec': 'high',
    'xc': 'PBE',
    'setups': {'Ti': '_sv', 'Zr': '_sv', 'Hf': '_sv', 'Nb':'_sv'},
    'lreal': 'False',
    'encut': 300.0,  # 560.0
    'istart': 0,
    'icharg': 2,
    'nelm': 200,
    'nelmin': 5,
    'kspacing': 0.5,  # 0.1
    'gamma': True,
    'ismear': 1,
    'sigma': 0.16,
    'algo': 'Normal',
    'lwave': False,
    'lcharg': False,
    'symprec': 1e-4,
    'ediff': 1e-6,
    'ncore': 1,
    'kpar': 11
}

opt = HEAoptimizer(input_path = 'db/TiZrHfTaNb.db',
                   vasp_params = vasp_params,
                   start_index = 20,
                   end_index = 30)
opt.run()
