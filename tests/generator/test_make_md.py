import os,json,glob,shutil
import dpdata
import numpy as np
import unittest

from .context import make_model_devi
from .context import parse_cur_job
from .context import param_file
from .context import machine_file
from .context import my_file_cmp
from .context import setUpModule
from .comp_sys import test_atom_names
from .comp_sys import test_atom_types
from .comp_sys import test_coord
from .comp_sys import test_cell

def _make_fake_models(idx, numb_models) :
    train_dir = os.path.join('iter.%06d' % idx, 
                             '00.train')
    os.makedirs(train_dir, exist_ok = True)
    pwd = os.getcwd()
    os.chdir(train_dir)
    for ii in range(numb_models) :
        os.makedirs('%03d' % ii, exist_ok = True)
        with open(os.path.join('%03d' % ii, 'forzen_model.pb'), 'w') as fp:
            fp.write(str(ii))
        if not os.path.isfile('graph.%03d.pb' % ii) :
            os.symlink(os.path.join('%03d' % ii, 'forzen_model.pb'), 
                       'graph.%03d.pb' % ii)
    os.chdir(pwd)


def _check_confs(testCase, idx, jdata) :
    md_dir = os.path.join('iter.%06d' % idx, 
                          '01.model_devi')
    tasks = glob.glob(os.path.join(md_dir, 'task.*'))
    tasks.sort()
    cur_job = jdata['model_devi_jobs'][idx]
    sys_idx = cur_job['sys_idx']
    sys_configs = jdata['sys_configs']
    poscars = []
    for ii in sys_idx :
        sys_poscars = []
        for ss in sys_configs[ii]:
            tmp_poscars = glob.glob(ss)
            sys_poscars += tmp_poscars
        sys_poscars.sort()
        poscars.append(sys_poscars)
    for ii in tasks :
        conf_file = os.path.join(ii, 'conf.lmp')
        l_conf_file = os.path.basename(os.readlink(conf_file))
        poscar_file = poscars[int(l_conf_file.split('.')[0])][int(l_conf_file.split('.')[1])]
        sys_0 = dpdata.System(conf_file, type_map = jdata['type_map'])
        sys_1 = dpdata.System(poscar_file)
        test_atom_names(testCase, sys_0, sys_1)
        test_atom_types(testCase, sys_0, sys_1)
        test_cell(testCase, sys_0, sys_1)
        test_coord(testCase, sys_0, sys_1)
        

def _check_pb(testCase, idx) :
    md_dir = os.path.join('iter.%06d' % idx, 
                          '01.model_devi')
    tr_dir = os.path.join('iter.%06d' % idx, 
                          '00.train')
    md_pb = glob.glob(os.path.join(md_dir, 'grapb*pb'))
    tr_pb = glob.glob(os.path.join(tr_dir, 'grapb*pb'))
    md_pb.sort()
    tr_pb.sort()
    for ii,jj in zip(md_pb, tr_pb) :
        my_file_cmp(testCase,ii,jj)

        
def _check_traj_dir(testCase, idx) :
    md_dir = os.path.join('iter.%06d' % idx, 
                          '01.model_devi')
    tasks = glob.glob(os.path.join(md_dir, 'task.*'))
    tasks.sort()
    for ii in tasks:
        testCase.assertTrue(os.path.isdir(os.path.join(ii, 'traj')))


def _get_lammps_pt(lmp_input) :
    with open(lmp_input) as fp: 
        for ii in fp:
            if 'variable' in ii and 'TEMP' in ii :
                lt = float(ii.split()[3])
            if 'variable' in ii and 'PRES' in ii :
                lp = float(ii.split()[3])
    return lt,lp

def _check_pt(testCase, idx, jdata) :
    md_dir = os.path.join('iter.%06d' % idx, 
                          '01.model_devi')
    tasks = glob.glob(os.path.join(md_dir, 'task.*'))
    tasks.sort()
    cur_job = jdata['model_devi_jobs'][idx]
    ensemble, nsteps, trj_freq, temps, press, pka_e, dt = parse_cur_job(cur_job)
    testCase.assertTrue(ensemble, 'npt')
    # get poscars
    sys_idx = cur_job['sys_idx']
    sys_configs = jdata['sys_configs']
    poscars = []
    for ii in sys_idx :
        sys_poscars = []
        for ss in sys_configs[ii]:
            tmp_poscars = glob.glob(ss)
            sys_poscars += tmp_poscars
        sys_poscars.sort()
        poscars.append(sys_poscars)
    for sidx,ii in enumerate(poscars) :
        count = 0
        for ss in ii:
            for tt in temps:
                for pp in press:
                    task_dir = os.path.join(md_dir, 'task.%03d.%06d' % (sidx, count))
                    lt, lp = _get_lammps_pt(os.path.join(task_dir, 'input.lammps'))
                    testCase.assertAlmostEqual(lt, tt)
                    testCase.assertAlmostEqual(lp, pp)
                    count += 1
                    

class TestMakeModelDevi(unittest.TestCase):
    def test_make_model_devi (self) :        
        if os.path.isdir('iter.000000') :
            shutil.rmtree('iter.000000')
        with open (param_file, 'r') as fp :
            jdata = json.load (fp)
        with open (machine_file, 'r') as fp:
            mdata = json.load (fp)
        _make_fake_models(0, jdata['numb_models'])
        make_model_devi(0, jdata, mdata)
        _check_pb(self, 0)
        _check_confs(self, 0, jdata)
        _check_traj_dir(self, 0)
        _check_pt(self, 0, jdata)
        shutil.rmtree('iter.000000')


if __name__ == '__main__':
    unittest.main()
