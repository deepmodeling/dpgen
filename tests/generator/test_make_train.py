#!/usr/bin/env python3

import os,sys,json,glob,shutil,dpdata
import numpy as np
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
__package__ = 'generator'
import tempfile
from .context import make_train, run_train
from .context import param_file
from .context import param_file_v1
from .context import param_file_v1_et
from .context import machine_file
from .context import machine_file_v1
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
    train_param = jdata.get('train_param', 'input.json')
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

def _check_model_input_dict(testCase, input_dict, init_data_sys, init_batch_size, default_training_param):
    for ii in input_dict.keys() :
        if ii == 'systems' :
            for jj,kk in zip(input_dict[ii], init_data_sys):
                testCase.assertEqual(jj, kk)
        elif ii == 'batch_size' :
            for jj, kk in zip(input_dict[ii], init_batch_size) :
                testCase.assertEqual(jj, kk)
        elif ii == 'seed':
            # can be anything
            pass
        elif ii == 'numb_fparam':
            testCase.assertEqual(input_dict[ii], 1)
        elif ii == 'numb_aparam':
            testCase.assertEqual(input_dict[ii], 1)
        else :                        
            testCase.assertEqual(input_dict[ii], default_training_param[ii])


def _check_model_inputs_v1(testCase, iter_idx, jdata, reuse = False) :
    train_param = jdata.get('train_param', 'input.json')
    numb_models = jdata['numb_models']
    use_ele_temp = jdata.get('use_ele_temp', 0)
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
        if use_ele_temp == 1:
            testCase.assertTrue('numb_fparam' in jdata0['model']['fitting_net'])
            testCase.assertFalse('numb_aparam' in jdata0['model']['fitting_net'])
        if use_ele_temp == 2:
            testCase.assertTrue('numb_aparam' in jdata0['model']['fitting_net'])
            testCase.assertFalse('numb_fparam' in jdata0['model']['fitting_net'])
        _check_model_input_dict(testCase, jdata0['model']['descriptor'], init_data_sys, init_batch_size, default_training_param['model']['descriptor'])
        _check_model_input_dict(testCase, jdata0['model']['fitting_net'], init_data_sys, init_batch_size, default_training_param['model']['fitting_net'])
        _check_model_input_dict(testCase, jdata0['loss'], init_data_sys, init_batch_size, default_training_param['loss'])
        _check_model_input_dict(testCase, jdata0['learning_rate'], init_data_sys, init_batch_size, default_training_param['learning_rate'])
        _check_model_input_dict(testCase, jdata0['training'], init_data_sys, init_batch_size, default_training_param['training'])
        if reuse:
            testCase.assertEqual(jdata['training_reuse_stop_batch'], 
                                 jdata0['training']['stop_batch'])
            testCase.assertEqual(jdata['training_reuse_start_lr'], 
                                 jdata0['learning_rate']['start_lr'])
            testCase.assertEqual(jdata['training_reuse_start_pref_e'], 
                                 jdata0['loss']['start_pref_e'])
            testCase.assertEqual(jdata['training_reuse_start_pref_f'], 
                                 jdata0['loss']['start_pref_f'])
            old_ratio = jdata['training_reuse_old_ratio']
            testCase.assertEqual(jdata0['training']['auto_prob_style'], 
                                 "prob_sys_size; 0:1:%f; 1:2:%f" % (old_ratio, 1-old_ratio))


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
    tmp_sys = dpdata.LabeledSystem('out_data_post_fp_vasp/02.fp/task.000.000000/OUTCAR')
    tmp_sys1 = tmp_sys.sub_system([0])
    tmp_sys2 = tmp_sys1
    for ii in range(1, nframes):
        tmp_sys2.append(tmp_sys1)
    tmp_sys2.to('deepmd/npy', dirname)


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
        # No longer support for DeePMD-kit-0.x version. 
        return    
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
        # No longer support for DeePMD-kit-0.x version. 
        return
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
        # No longer support for DeePMD-kit-0.x version. 
        return
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


    def test_1_data_v1(self) :
        with open (param_file_v1, 'r') as fp :
            jdata = json.load (fp)
            jdata.pop('use_ele_temp', None)
        with open (machine_file_v1, 'r') as fp:
            mdata = json.load (fp)
        make_train(0, jdata, mdata)
        # make fake fp results #data == fp_task_min
        _make_fake_fp(0, 0, jdata['fp_task_min'])
        # make iter1 train
        make_train(1, jdata, mdata)
        # check data is linked
        self.assertTrue(os.path.isdir(os.path.join('iter.000001', '00.train', 'data.iters', 'iter.000000', '02.fp')))
        # check models inputs
        _check_model_inputs_v1(self, 1, jdata)
        # remove testing dirs
        shutil.rmtree('iter.000001')
        shutil.rmtree('iter.000000')


    def test_1_data_reuse_v1(self) :
        with open (param_file_v1, 'r') as fp :
            jdata = json.load (fp)
            jdata.pop('use_ele_temp', None)
            jdata['training_reuse_iter'] = 1
            jdata['training_reuse_old_ratio'] = 0.8
            jdata['training_reuse_stop_batch'] = 400000
            jdata['training_reuse_start_lr'] = 1e-4
            jdata['training_reuse_start_pref_e'] = 0.1
            jdata['training_reuse_start_pref_f'] = 100
        with open (machine_file_v1, 'r') as fp:
            mdata = json.load (fp)
        make_train(0, jdata, mdata)
        # make fake fp results #data == fp_task_min
        _make_fake_fp(0, 0, jdata['fp_task_min'])
        # make iter1 train
        make_train(1, jdata, mdata)
        # check data is linked
        self.assertTrue(os.path.isdir(os.path.join('iter.000001', '00.train', 'data.iters', 'iter.000000', '02.fp')))
        # check old models are linked
        self.assertTrue(os.path.isdir(os.path.join('iter.000001', '00.train', '000', 'old')))
        self.assertTrue(os.path.isdir(os.path.join('iter.000001', '00.train', '001', 'old')))
        self.assertTrue(os.path.isdir(os.path.join('iter.000001', '00.train', '002', 'old')))
        self.assertTrue(os.path.isdir(os.path.join('iter.000001', '00.train', '003', 'old')))
        # check models inputs
        _check_model_inputs_v1(self, 1, jdata, reuse = True)
        # remove testing dirs
        shutil.rmtree('iter.000001')
        shutil.rmtree('iter.000000')

        
    def test_1_data_v1_eletron_temp(self) :
        with open (param_file_v1_et, 'r') as fp :
            jdata = json.load (fp)
        with open (machine_file_v1, 'r') as fp:
            mdata = json.load (fp)
        make_train(0, jdata, mdata)
        # make fake fp results #data == fp_task_min
        _make_fake_fp(0, 0, jdata['fp_task_min'])
        # make iter1 train
        make_train(1, jdata, mdata)
        # check data is linked
        self.assertTrue(os.path.isdir(os.path.join('iter.000001', '00.train', 'data.iters', 'iter.000000', '02.fp')))
        # check models inputs
        _check_model_inputs_v1(self, 1, jdata)
        # remove testing dirs
        shutil.rmtree('iter.000001')
        shutil.rmtree('iter.000000')
        
    def test_1_data_v1_h5(self) :
        """Test HDF5 file as input data."""
        dpdata.LabeledSystem("data/deepmd", fmt='deepmd/npy').to_deepmd_hdf5('data/deepmd.hdf5')
        with open (param_file_v1, 'r') as fp :
            jdata = json.load (fp)
        jdata.pop('use_ele_temp', None)
        jdata['init_data_sys'].append('deepmd.hdf5')
        jdata['init_batch_size'].append('auto')
        with open (machine_file_v1, 'r') as fp:
            mdata = json.load (fp)
        make_train(0, jdata, mdata)
        # make fake fp results #data == fp_task_min
        _make_fake_fp(0, 0, jdata['fp_task_min'])
        # make iter1 train
        make_train(1, jdata, mdata)
        # check data is linked
        self.assertTrue(os.path.isdir(os.path.join('iter.000001', '00.train', 'data.iters', 'iter.000000', '02.fp')))
        # check models inputs
        with open(os.path.join('iter.%06d' % 1, 
                               '00.train', 
                               '%03d' % 0,
                               "input.json")) as fp:
            jdata0 = json.load(fp)
        self.assertEqual(jdata0['training']['systems'], [
            '../data.init/deepmd',
            '../data.init/deepmd.hdf5#',
            '../data.iters/iter.000000/02.fp/data.000',
        ])
        # test run_train -- confirm transferred files are correct
        with tempfile.TemporaryDirectory() as remote_root:
            run_train(1, jdata, {
                "api_version": "1.0",
                "train_command": (
                    "test -d ../data.init/deepmd"
                    "&& test -f ../data.init/deepmd.hdf5"
                    "&& test -d ../data.iters/iter.000000/02.fp/data.000"
                    "&& touch frozen_model.pb lcurve.out model.ckpt.meta model.ckpt.index model.ckpt.data-00000-of-00001 checkpoint"
                    "&& echo dp"
                ),
                "train_machine": {
                    "batch_type": "shell",
                    "local_root": "./",
                    "remote_root": remote_root,
                    "context_type": "local",
                },
                "train_resources": {
                    "group_size": 1,
                },
            })

        # remove testing dirs
        shutil.rmtree('iter.000001')
        shutil.rmtree('iter.000000')
        os.remove('data/deepmd.hdf5')


if __name__ == '__main__':
    unittest.main()
