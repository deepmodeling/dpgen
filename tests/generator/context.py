import sys,os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
from dpgen.generator.run import *
from dpgen.generator.lib.gaussian import detect_multiplicity

param_file = 'param-mg-vasp.json'
param_old_file = 'param-mg-vasp-old.json'
param_pwscf_file = 'param-pyridine-pwscf.json'
param_pwscf_old_file = 'param-pyridine-pwscf-old.json'
param_gaussian_file = 'param-pyridine-gaussian.json'
param_siesta_file = 'param-pyridine-siesta.json'
param_cp2k_file = 'param-pyridine-cp2k.json'
machine_file = 'machine-local.json'
param_diy_file = 'param-mg-vasp-diy.json'

def my_file_cmp(test, f0, f1):
    with open(f0) as fp0 :
        with open(f1) as fp1:
            test.assertTrue(fp0.read() == fp1.read())

def setUpModule():
    os.chdir(os.path.abspath(os.path.dirname(__file__)))
