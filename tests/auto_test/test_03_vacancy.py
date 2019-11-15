import numpy as np
import unittest
import json,re,os,shutil,glob
from .input_data import *
from .context import setUpModule
from dpgen.auto_test import gen_03_vacancy,cmpt_03_vacancy

class TestVacancy(unittest.TestCase):

    def test_gen_vacancy(self):
        conf_dir="confs/Cu/std-fcc"
        global_task_name='03.vacancy'
        task_path=os.path.abspath(re.sub('confs', global_task_name, conf_dir))
        with open (param_file, 'r') as fp :
            jdata = json.load (fp)

        supercell=jdata['supercell']

        gen_03_vacancy.make_vasp(jdata,conf_dir,supercell)
        kspacing = jdata['vasp_params']['kspacing']
        vasp_str='vasp-k%.2f' % kspacing
        vasp_path=os.path.join(task_path,vasp_str)
        vasp_check=[]
        struct_path=os.path.join(vasp_path,'struct-*')
        struct_task=glob.glob(struct_path)
        for ss in struct_task :
            vasp_check+=[os.path.join(ss,ii) for ii in vasp_input]
        for ii in vasp_check:
            if self.assertTrue(os.path.isfile(ii)):
                os.remove(ii)

        gen_03_vacancy.make_lammps(jdata,conf_dir,"deepmd",supercell)
        dp_path = os.path.join(task_path,'deepmd')
        dp_check=[]
        struct_path=os.path.join(dp_path,'struct-*')
        struct_task=glob.glob(struct_path)
        for ss in struct_task :
            dp_check+=[os.path.join(ss,ii) for ii in dp_input]
        for ii in dp_check:
            if self.assertTrue(os.path.isfile(ii)):
                os.remove(ii)

    def test_cmpt_vacancy(self):
        conf_dir="confs/Cu/std-fcc"
        with open (param_file, 'r') as fp :
            jdata = json.load (fp)
        supercell=jdata['supercell']
        cmpt_03_vacancy.cmpt_vasp(jdata, conf_dir,supercell)
        cmpt_03_vacancy.cmpt_deepmd_lammps(jdata, conf_dir,supercell,'deepmd')





if __name__== '__main__':
    unittest.main()
