import sys,os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
from dpgen.generator.run import *
from dpgen.generator.lib.gaussian import detect_multiplicity

param_file = 'param-mg-vasp.json'
param_old_file = 'param-mg-vasp-old.json'
param_pwscf_file = 'param-pyridine-pwscf.json'
param_pwscf_old_file = 'param-pyridine-pwscf-old.json'
param_gaussian_file = 'param-pyridine-gaussian.json'
machine_file = 'machine-local.json'
param_diy_file = 'param-mg-vasp-diy.json'
