import sys,os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
from dpgen.data.gen import *

param_file = 'al.json'
abacus_param_file = 'ch4.json'
abacus_stru_file = 'STRU.hcp'
def setUpModule():
    os.chdir(os.path.abspath(os.path.dirname(__file__)))
