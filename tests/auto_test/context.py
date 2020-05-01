import sys,os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
from dpgen.auto_test.lib.vasp import *

def setUpModule():
    os.chdir(os.path.abspath(os.path.dirname(__file__)))
