import numpy as np
import unittest
import json,re,os,filecmp,glob
from .input_data import *
from .context import setUpModule
from dpgen.auto_test import gen_05_surf,cmpt_05_surf

class TestSurf(unittest.TestCase):

    def test_gen_surf(self):
        conf_dir="confs/Cu/std-fcc"
        global_task_name='05.surf'
        task_path=os.path.abspath(re.sub('confs', global_task_name, conf_dir))
        with open (param_file, 'r') as fp :
            jdata = json.load (fp)
        max_miller=jdata['max_miller']

        gen_05_surf.make_vasp(jdata,conf_dir,max_miller,False,False)
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


        gen_05_surf.make_lammps(jdata,conf_dir,max_miller,False,False,'deepmd')
        dp_path = os.path.join(task_path,'deepmd')
        dp_check=[]
        struct_path=os.path.join(dp_path,'struct-*')
        struct_task=glob.glob(struct_path)
        for ss in struct_task :
            dp_check+=[os.path.join(ss,ii) for ii in dp_input]
        for ii in dp_check:
            if self.assertTrue(os.path.isfile(ii)):
                os.remove(ii)

    def test_cmpt_surf(self):
        conf_dir="confs/Cu/std-fcc"
        global_task_name='05.surf'
        task_path=os.path.abspath(re.sub('confs', global_task_name, conf_dir))

        with open (param_file, 'r') as fp :
            jdata = json.load (fp)
        cmpt_05_surf.cmpt_vasp(jdata, conf_dir,False)
        kspacing = jdata['vasp_params']['kspacing']
        vasp_str='vasp-k%.2f' % kspacing
        vasp_path=os.path.join(task_path,vasp_str)
        result =os.path.join(vasp_path,'result')
        ref = os.path.join(vasp_path,'ref')
        self.assertTrue(filecmp.cmp(result,ref))

        cmpt_05_surf.cmpt_deepmd_lammps(jdata, conf_dir,'deepmd',False)
        dp_path = os.path.join(task_path,'deepmd')
        result =os.path.join(dp_path,'result')
        ref = os.path.join(dp_path,'ref')
        self.assertTrue(filecmp.cmp(result,ref))





if __name__== '__main__':
    unittest.main()
