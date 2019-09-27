#!/usr/bin/env python3

import os,sys,json,glob,shutil
import numpy as np
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
__package__ = 'generator'
from .context import make_train
from .context import param_file
from .context import machine_file
from .context import setUpModule

def _comp_sys_files (sys0, sys1) :
    pwd = os.getcwd()
    os.chdir(sys0) 
    files = glob.glob('*.raw') 
    set_files = glob.glob('set.*/*npy') 
    # files += set_files
    os.chdir(pwd)
    for ii in files :
        with open(os.path.join(sys0, ii)) as fp0 :
            with open(os.path.join(sys1, ii)) as fp1:
                if fp0.read() != fp1.read() :
                    return False
    for ii in set_files:
        t0 = np.load(os.path.join(sys0, ii))
        t1 = np.load(os.path.join(sys1, ii))
        if np.linalg.norm(t0-t1) > 1e-12 :
            return False
    return True

def _comp_init_data(testCase, iter_idx, init_data_prefix, init_data_sys) :
    for ii in init_data_sys :
        sys0 = os.path.join(init_data_prefix, ii) 
        sys1 = os.path.join('iter.%06d' % iter_idx, 
                            '00.train',
                            'data.init', 
                            ii)
        testCase.assertTrue(_comp_sys_files(sys0, sys1),
                            'systems %s %s are not identical' % (sys0, sys1))

def _check_numb_models(testCase, iter_idx, numb_models) :
    models = glob.glob(os.path.join('iter.%06d' % iter_idx, 
                                    '00.train', 
                                    '[0-9][0-9][0-9]'))
    testCase.assertTrue(len(models), numb_models)


def _check_model_inputs(testCase, iter_idx, jdata) :
    train_param = jdata['train_param']
    numb_models = jdata['numb_models']
    default_training_param = jdata['default_training_param']
    init_data_sys = [os.path.join('..', 'data.init', ii) for ii in jdata['init_data_sys']]
    init_batch_size = jdata['init_batch_size']
    sys_batch_size = jdata['sys_batch_size']
    if iter_idx > 0 :
        systems = glob.glob(os.path.join('iter.*', '02.fp', 'data.*'))
        for ii in systems :
            init_data_sys.append(os.path.join('..', 'data.iters', ii))
            sys_idx = int(os.path.basename(ii).split('.')[1])
            init_batch_size.append(sys_batch_size[sys_idx])
    for kk in range(numb_models) :
        with open(os.path.join('iter.%06d' % iter_idx, 
                               '00.train', 
                               '%03d' % kk,
                               train_param)) as fp :
            jdata0 = json.load(fp)
        # keys except 'systems', 'batch_size', 'seed' should be identical
        for ii in jdata0.keys() :
            if ii == 'systems' :
                for jj,kk in zip(jdata0[ii], init_data_sys):
                    testCase.assertEqual(jj, kk)
            elif ii == 'batch_size' :
                for jj, kk in zip(jdata0[ii], init_batch_size) :
                    testCase.assertEqual(jj, kk)
            elif ii == 'seed':
                pass
            else :
                testCase.assertEqual(jdata0[ii], default_training_param[ii])

def _make_fake_fp(iter_idx, sys_idx, nframes):
    for ii in range(nframes) :
        dirname = os.path.join('iter.%06d' % iter_idx, 
                               '02.fp', 
                               'task.%03d.%06d' % (sys_idx, ii))
        os.makedirs(dirname, exist_ok = True)           
    dirname = os.path.join('iter.%06d' % iter_idx, 
                           '02.fp', 
                           'data.%03d' % sys_idx)
    os.makedirs(dirname, exist_ok = True)
    box_str = ['0' for ii in range(9)]
    box_str = ' '.join(box_str)
    with open(os.path.join(dirname, 'box.raw'), 'w') as fp :
        for ii in range(nframes) :
            fp.write(box_str + '\n')


def _check_pb_link(testCase, iter_idx, numb_models) :
    pwd = os.getcwd()
    os.chdir(os.path.join('iter.%06d' % iter_idx, 
                          '00.train'))
    for ii in range(numb_models) :
        lnk = os.readlink('graph.%03d.pb' % ii) 
        testCase.assertEqual(lnk, os.path.join('%03d' % ii, 'frozen_model.pb'))
    os.chdir(pwd)

class TestMakeTrain(unittest.TestCase):
    def test_0 (self) :        
        with open (param_file, 'r') as fp :
            jdata = json.load (fp)
        with open (machine_file, 'r') as fp:
            mdata = json.load (fp)
        fp_task_min = jdata['fp_task_min']
        make_train(0, jdata, mdata)
        # comp init data
        init_data_prefix = jdata['init_data_prefix']
        init_data_sys = jdata['init_data_sys']
        _comp_init_data(self, 0, init_data_prefix, init_data_sys)
        # check number of models
        _check_numb_models(self, 0, jdata['numb_models'])
        # check models inputs
        _check_model_inputs(self, 0, jdata)
        # remove iter
        shutil.rmtree('iter.000000')

    def test_1_data(self) :
        with open (param_file, 'r') as fp :
            jdata = json.load (fp)
        with open (machine_file, 'r') as fp:
            mdata = json.load (fp)
        make_train(0, jdata, mdata)
        # make fake fp results #data == fp_task_min
        _make_fake_fp(0, 0, jdata['fp_task_min'])
        # make iter1 train
        make_train(1, jdata, mdata)
        # check data is linked
        self.assertTrue(os.path.isdir(os.path.join('iter.000001', '00.train', 'data.iters', 'iter.000000', '02.fp')))
        # check models inputs
        _check_model_inputs(self, 1, jdata)
        # remove testing dirs
        shutil.rmtree('iter.000001')
        shutil.rmtree('iter.000000')


    def test_1_skip(self):
        with open (param_file, 'r') as fp :
            jdata = json.load (fp)
        with open (machine_file, 'r') as fp:
            mdata = json.load (fp)
        make_train(0, jdata, mdata)
        # make fake fp results #data == fp_task_min - 1
        _make_fake_fp(0, 0, jdata['fp_task_min'] - 1)
        # make iter1 train
        make_train(1, jdata, mdata)
        self.assertTrue(os.path.isfile(os.path.join('iter.000001', '00.train', 'copied')))
        # check pb file linked
        _check_pb_link(self, 1, jdata['numb_models'])
        # remove testing dirs
        shutil.rmtree('iter.000001')
        shutil.rmtree('iter.000000')
        

if __name__ == '__main__':
    unittest.main()
