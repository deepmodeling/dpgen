import os,sys,json,glob,shutil,textwrap
import dpdata
import numpy as np
import unittest
from pathlib import Path

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
__package__ = 'generator'

from .context import make_calypso_input
from .context import write_model_devi_out
from .context import _parse_calypso_input
from .context import _parse_calypso_dis_mtx

# temp dir
test_path = Path('.').joinpath('calypso_test_path')
test_path.mkdir(parents=True,exist_ok=True)
os.system('rm calypso_test_path/*')
fmax = 0.01
cwd = os.getcwd()

model_devi = np.array([[0.000000e+00, 2.328491e-02, 5.476687e-09, 1.009454e-02,
        3.279617e-02, 4.053224e-03, 1.869795e-02, 2.184905e+00],
       [1.000000e+00, 3.668334e-02, 8.200870e-09, 1.706517e-02,
        2.844074e-02, 7.093109e-03, 1.623275e-02, 2.424708e+00],
       [2.000000e+00, 2.832296e-02, 4.828951e-08, 1.573961e-02,
        2.443331e-02, 2.871548e-03, 1.489787e-02, 2.564113e+00]])

model_devi_jobs = {"model_devi_jobs": {"times":[4],"NameOfAtoms":["Mg","Al","Cu"],"NumberOfAtoms":[1,1,1],"NumberOfFormula":[1,4],"Volume":[30],"DistanceOfIon":[[1.48,1.44,1.59],[1.44,1.41,1.56],[1.59,1.56,1.70]],"PsoRatio":[0.6],"PopSize":[5],"MaxStep":[3],"ICode":[13],"Split":"T","VSC":"T","MaxNumAtom":[31],
  "CtrlRange":[[1,10],[1,10],[1,10]],"PSTRESS":[0],"fmax":[0.01]}}
    

class TestCALYPSOScript(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_write_model_devi_out(self):
        #devi = write_model_devi_out(model_devi, 'calypso_test_path/model_devi.out')
        #ndevi = np.loadtxt('calypso_test_path/model_devi.out')
        devi = write_model_devi_out(model_devi, 'model_devi.out')
        ndevi = np.loadtxt('model_devi.out')
        self.assertEqual(ndevi[2,4],model_devi[2,4])
        os.remove('model_devi.out')

    def test_make_calypso_input(self):
        ret = make_calypso_input(["Mg","Al","Cu"],[1,1,1],[1,4],30,[
                              [1.48,1.44,1.59],[1.44,1.41,1.56],[1.59,1.56,1.70]
                              ],0.6,5,3,13,"T","T",31,[[1,10],[1,10],[1,10]],0,0.01)
        #with open('calypso_test_path/input.dat','w') as fin:
        with open('input.dat','w') as fin:
            fin.write(ret)
        f = open('input.dat')
        #f = open('calypso_test_path/input.dat')
        lines = f.readlines()
        f.close()
        for line in lines :
            if line[0] == '#':
                continue
            if 'PopSize' in line:
                temp_1 = line.split('=')[1].strip()
                self.assertEqual(int(temp_1),5)
            if 'MaxStep' in line:
                temp_2 = line.split('=')[1].strip()
                self.assertEqual(int(temp_2),3)
                os.remove('input.dat')
                break

    def test_parse_calypso_input(self):
        ret = make_calypso_input(["Mg","Al","Cu"],[1,1,1],[1,4],30,[
                              [1.48,1.44,1.59],[1.44,1.41,1.56],[1.59,1.56,1.70]
                              ],0.6,5,3,13,"T","T",31,[[1,10],[1,10],[1,10]],0,0.01)
        #with open('calypso_test_path/input.dat','w') as fin:
        with open('input.dat','w') as fin:
            fin.write(ret)
        formula = _parse_calypso_input('NumberOfFormula','input.dat').split()
        #formula = _parse_calypso_input('NumberOfFormula',calypso_data).split()
        formula = list(map(int,formula))
        self.assertEqual(formula,model_devi_jobs.get('model_devi_jobs').get('NumberOfFormula'))

        nameofatoms = _parse_calypso_input('NameOfAtoms','input.dat').split()
        #nameofatoms = _parse_calypso_input('NameOfAtoms',calypso_data).split()
        self.assertEqual(nameofatoms,model_devi_jobs.get('model_devi_jobs').get('NameOfAtoms'))
        
        min_dis = _parse_calypso_dis_mtx(len(nameofatoms),'input.dat')
        #min_dis = _parse_calypso_dis_mtx(len(nameofatoms),calypso_data)
        self.assertEqual(float(min_dis),np.nanmin(model_devi_jobs.get('model_devi_jobs').get('DistanceOfIon')))
        os.remove('input.dat')


if __name__ == '__main__':
    unittest.main(verbosity=2)
