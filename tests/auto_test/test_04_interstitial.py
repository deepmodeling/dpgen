import numpy as np
import unittest
import json,re,os,shutil,glob
from input_data import *
from dpgen.auto_test import gen_04_interstitial,cmpt_04_interstitial

class Testinterstitial(unittest,TestCase):

    def test_gen_interstitial(self):
        conf_dir="confs/Cu/std-fcc"
        global_task_name='04.interstitial'
        task_path=os.path.abspath(re.sub('confs', global_task_name, conf_dir))
        with open (param_file, 'r') as fp :
            jdata = json.load (fp)

        supercell=jdata['supercell']
        insert_ele=jdata['insert_ele']

        gen_04_interstitial.make_vasp(jdata,conf_dir,supercell,insert_ele)
        kspacing = jdata['vasp_params']['kspacing']
        vasp_str='vasp-k%.2f' % kspacing
        vasp_path=os.path.join(task_path,vasp_str)
        vasp_check=[]
        struct_path=os.path.join(vasp_path,'struct-*')
        struct_task=glob.glob(struct_path)
        for ss in struct_task :
            vasp_check+=[os.path.join(ss,ii) for ii in vasp_input]
        for ii in vasp_check:
            if os.path.isfile(ii):
                os.remove(ii)
            else:
                raise "error in gen_04_interstitial.make_vasp "

        gen_04_interstitial.make_lammps(jdata,conf_dir,supercell,insert_ele,'deepmd')
        dp_path = os.path.join(task_path,'deepmd')
        dp_check=[]
        struct_path=os.path.join(dp_path,'struct-*')
        struct_task=glob.glob(struct_path)
        for ss in struct_task :
            dp_check+=[os.path.join(ss,ii) for ii in dp_input]
        for ii in dp_check:
            if os.path.isfile(ii):
                os.remove(ii)
            else:
                raise "error in gen_04_interstitial.make_lammps "

    def test_cmpt_interstitial(self):
        conf_dir="confs/Cu/std-fcc"
        with open (param_file, 'r') as fp :
            jdata = json.load (fp)
        supercell=jdata['supercell']
        insert_ele=jdata['insert_ele']
        cmpt_04_interstitial.cmpt_vasp(jdata, conf_dir,supercell)
        cmpt_04_interstitial.cmpt_deepmd_lammps(jdata, conf_dir,supercell,insert_ele,'deepmd')





if __name__== '__main__':
    unittest.main()
