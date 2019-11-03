import numpy as np
import unittest
import json,re,os,shutil,glob
from input_data import *
from dpgen.auto_test import gen_05_surf,cmpt_05_surf

class Testsurf(unittest,TestCase):

    def test_gen_surf(self):
        conf_dir="confs/Cu/std-fcc"
        global_task_name='05.surf'
        task_path=os.path.abspath(re.sub('confs', global_task_name, conf_dir))
        with open (param_file, 'r') as fp :
            jdata = json.load (fp)
        max_miller=jdata['max_miller']
        relax_box=jdata['relax_box']
        static=jdata['static-opt']

        gen_05_surf.make_vasp(jdata,conf_dir,max_miller,static,relax_box)
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
                raise "error in gen_05_surf.make_vasp "

        gen_05_surf.make_lammps(jdata,conf_dir,max_miller,static,relax_box,'deepmd')
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
                raise "error in gen_05_surf.make_lammps "

    def test_cmpt_surf(self):
        conf_dir="confs/Cu/std-fcc"
        with open (param_file, 'r') as fp :
            jdata = json.load (fp)
        static_opt=jdata['static-opt']
        cmpt_05_surf.cmpt_vasp(jdata, conf_dir,static_opt)
        cmpt_05_surf.cmpt_deepmd_lammps(jdata, conf_dir,'deepmd',static_opt)





if __name__== '__main__':
    unittest.main()
