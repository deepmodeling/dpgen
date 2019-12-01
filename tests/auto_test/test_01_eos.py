import numpy as np
import unittest
import json,re,os,filecmp
from .input_data import *
from .context import setUpModule
from dpgen.auto_test import gen_01_eos,cmpt_01_eos

class TestEos(unittest.TestCase):

    def test_gen_eos(self):
        conf_dir="confs/Cu/std-fcc"
        global_task_name='01.eos'
        task_path=os.path.abspath(re.sub('confs', global_task_name, conf_dir))
        with open (param_file, 'r') as fp :
            jdata = json.load (fp)
        vol_start = jdata['vol_start']
        vol_end = jdata['vol_end']
        vol_step = jdata['vol_step']

        gen_01_eos.make_vasp(jdata,conf_dir)
        kspacing = jdata['vasp_params']['kspacing']
        vasp_str='vasp-k%.2f' % kspacing
        vasp_path=os.path.join(task_path,vasp_str)

        vasp_check=[]
        for vol in np.arange(vol_start, vol_end, vol_step) :
            vol_path = os.path.join(vasp_path, 'vol-%.2f' % vol)
            vasp_check+=[os.path.join(vol_path,ii) for ii in vasp_input]
        for ii in vasp_check:
            if self.assertTrue(os.path.isfile(ii)):
                os.remove(ii)

        gen_01_eos.make_lammps(jdata,conf_dir,"deepmd")
        dp_path = os.path.join(task_path,'deepmd')
        dp_check=[]
        for vol in np.arange(vol_start, vol_end, vol_step) :
            vol_path = os.path.join(dp_path, 'vol-%.2f' % vol)
            dp_check+=[os.path.join(vol_path,ii) for ii in dp_input]
        for ii in vasp_check:
            if self.assertTrue(os.path.isfile(ii)):
                os.remove(ii)

    def test_cmpt_eos(self):
        conf_dir="confs/Cu/std-fcc"
        global_task_name='01.eos'
        task_path=os.path.abspath(re.sub('confs', global_task_name, conf_dir))
        with open (param_file, 'r') as fp :
            jdata = json.load (fp)

        cmpt_01_eos.comput_vasp_eos(jdata, conf_dir)
        kspacing = jdata['vasp_params']['kspacing']
        vasp_str='vasp-k%.2f' % kspacing
        vasp_path=os.path.join(task_path,vasp_str)
        result =os.path.join(vasp_path,'result')
        ref = os.path.join(vasp_path,'ref')
        self.assertTrue(filecmp.cmp(result,ref))

        cmpt_01_eos.comput_lmp_eos(jdata, conf_dir,'deepmd')
        dp_path = os.path.join(task_path,'deepmd')
        result =os.path.join(dp_path,'result')
        ref = os.path.join(dp_path,'ref')
        self.assertTrue(filecmp.cmp(result,ref))





if __name__== '__main__':
    unittest.main()
