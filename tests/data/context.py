import sys,os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
from dpgen.data.gen import *

param_file = 'al.json'

def setUpModule():
    os.chdir(os.path.abspath(os.path.dirname(__file__)))
