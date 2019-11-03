import numpy as np
import unittest
import json,re,os
from input_data import *
from dpgen.auto_test import gen_02_elastic,cmpt_02_elastic

class TestElastic(unittest,TestCase):

    def test_gen_elastic(self):
        conf_dir="confs/Cu/std-fcc"
        global_task_name='02.elastic'
        task_path=os.path.abspath(re.sub('confs', global_task_name, conf_dir))
        with open (param_file, 'r') as fp :
            jdata = json.load (fp)

        gen_02_elastic.make_vasp(jdata,conf_dir)
        kspacing = jdata['vasp_params']['kspacing']
        vasp_str='vasp-k%.2f' % kspacing
        vasp_path=os.path.join(task_path,vasp_str)

        vasp_check=[]
        for dfm in np.arange(24) :
            dfm_path = os.path.join(vasp_path, 'dfm-%03d' % dfm)
            vasp_check+=[os.path.join(dfm_path,ii) for ii in vasp_input]
        for ii in vasp_check:
            if os.path.isfile(ii):
                os.remove(ii)
            else:
                raise "error in gen_02_elastic.make_vasp "

        gen_02_elastic.make_lammps(jdata,conf_dir,"deepmd")
        dp_path = os.path.join(task_path,'deepmd')
        dp_check=[]
        for dfm in np.arange(24) :
            dfm_path = os.path.join(dp_path, 'dfm-%03d' % dfm)
            dp_check+=[os.path.join(dfm_path,ii) for ii in dp_input]
        for ii in dp_check:
            if os.path.isfile(ii):
                os.remove(ii)
            else:
                raise "error in gen_02_elastic.make_lammps "

    def test_cmpt_elastic(self):
        conf_dir="confs/Cu/std-fcc"
        with open (param_file, 'r') as fp :
            jdata = json.load (fp)
        cmpt_02_elastic.cmpt_vasp(jdata, conf_dir)
        cmpt_02_elastic.cmpt_deepmd_lammps(jdata, conf_dir,'deepmd')





if __name__== '__main__':
    unittest.main()
