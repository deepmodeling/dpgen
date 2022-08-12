import glob
import json
import os
import re

import numpy as np
from monty.serialization import loadfn, dumpfn

import dpgen.auto_test.lib.vasp as vasp
from dpgen import dlog
from dpgen.auto_test.Property import Property
from dpgen.auto_test.refine import make_refine
from dpgen.auto_test.reproduce import make_repro
from dpgen.auto_test.reproduce import post_repro

import dpgen.generator.lib.abacus_scf as abacus_scf
import dpgen.auto_test.lib.abacus as abacus

class EOS(Property):
    def __init__(self,
                 parameter,inter_param=None):
        parameter['reproduce'] = parameter.get('reproduce', False)
        self.reprod = parameter['reproduce']
        if not self.reprod:
            if not ('init_from_suffix' in parameter and 'output_suffix' in parameter):
                self.vol_start = parameter['vol_start']
                self.vol_end = parameter['vol_end']
                self.vol_step = parameter['vol_step']
                parameter['vol_abs'] = parameter.get('vol_abs', False)
                self.vol_abs = parameter['vol_abs']
            parameter['cal_type'] = parameter.get('cal_type', 'relaxation')
            self.cal_type = parameter['cal_type']
            default_cal_setting = {"relax_pos": True,
                                   "relax_shape": True,
                                   "relax_vol": False}
            if 'cal_setting' not in parameter:
                parameter['cal_setting'] = default_cal_setting
            else:
                if "relax_pos" not in parameter['cal_setting']:
                    parameter['cal_setting']['relax_pos'] = default_cal_setting['relax_pos']
                if "relax_shape" not in parameter['cal_setting']:
                    parameter['cal_setting']['relax_shape'] = default_cal_setting['relax_shape']
                if "relax_vol" not in parameter['cal_setting']:
                    parameter['cal_setting']['relax_vol'] = default_cal_setting['relax_vol']
            self.cal_setting = parameter['cal_setting']
        else:
            parameter['cal_type'] = 'static'
            self.cal_type = parameter['cal_type']
            default_cal_setting = {"relax_pos": False,
                                   "relax_shape": False,
                                   "relax_vol": False}
            if 'cal_setting' not in parameter:
                parameter['cal_setting'] = default_cal_setting
            else:
                if "relax_pos" not in parameter['cal_setting']:
                    parameter['cal_setting']['relax_pos'] = default_cal_setting['relax_pos']
                if "relax_shape" not in parameter['cal_setting']:
                    parameter['cal_setting']['relax_shape'] = default_cal_setting['relax_shape']
                if "relax_vol" not in parameter['cal_setting']:
                    parameter['cal_setting']['relax_vol'] = default_cal_setting['relax_vol']
            self.cal_setting = parameter['cal_setting']
            parameter['init_from_suffix'] = parameter.get('init_from_suffix', '00')
            self.init_from_suffix = parameter['init_from_suffix']
        self.parameter = parameter
        self.inter_param = inter_param if inter_param != None else {'type': 'vasp'}

    def make_confs(self,
                   path_to_work,
                   path_to_equi,
                   refine=False):
        path_to_work = os.path.abspath(path_to_work)
        if os.path.exists(path_to_work):
            dlog.warning('%s already exists' % path_to_work)
        else:
            os.makedirs(path_to_work)
        path_to_equi = os.path.abspath(path_to_equi)

        if 'start_confs_path' in self.parameter and os.path.exists(self.parameter['start_confs_path']):
            init_path_list = glob.glob(os.path.join(self.parameter['start_confs_path'], '*'))
            struct_init_name_list = []
            for ii in init_path_list:
                struct_init_name_list.append(ii.split('/')[-1])
            struct_output_name = path_to_work.split('/')[-2]
            assert struct_output_name in struct_init_name_list
            path_to_equi = os.path.abspath(os.path.join(self.parameter['start_confs_path'],
                                                        struct_output_name, 'relaxation', 'relax_task'))

        cwd = os.getcwd()
        task_list = []
        if self.reprod:
            print('eos reproduce starts')
            if 'init_data_path' not in self.parameter:
                raise RuntimeError("please provide the initial data path to reproduce")
            init_data_path = os.path.abspath(self.parameter['init_data_path'])
            task_list = make_repro(self.inter_param,init_data_path, self.init_from_suffix,
                                   path_to_work, self.parameter.get('reprod_last_frame', True))
            os.chdir(cwd)

        else:
            if refine:
                print('eos refine starts')
                task_list = make_refine(self.parameter['init_from_suffix'],
                                        self.parameter['output_suffix'],
                                        path_to_work)
                os.chdir(cwd)

                init_from_path = re.sub(self.parameter['output_suffix'][::-1],
                                        self.parameter['init_from_suffix'][::-1],
                                        path_to_work[::-1], count=1)[::-1]
                task_list_basename = list(map(os.path.basename, task_list))

                for ii in task_list_basename:
                    init_from_task = os.path.join(init_from_path, ii)
                    output_task = os.path.join(path_to_work, ii)
                    os.chdir(output_task)
                    if os.path.isfile('eos.json'):
                        os.remove('eos.json')
                    if os.path.islink('eos.json'):
                        os.remove('eos.json')
                    os.symlink(os.path.relpath(os.path.join(init_from_task, 'eos.json')), 'eos.json')
                os.chdir(cwd)

            else:
                print('gen eos from ' + str(self.vol_start) + ' to ' + str(self.vol_end) + ' by every ' + str(self.vol_step))
                if self.vol_abs : 
                    dlog.info('treat vol_start and vol_end as absolute volume')
                else : 
                    dlog.info('treat vol_start and vol_end as relative volume')

                if self.inter_param['type'] == 'abacus':
                    equi_contcar = os.path.join(path_to_equi,abacus.final_stru(path_to_equi))
                else:
                    equi_contcar = os.path.join(path_to_equi, 'CONTCAR')

                if not os.path.isfile(equi_contcar):
                    raise RuntimeError("Can not find %s, please do relaxation first" % equi_contcar)

                if self.inter_param['type'] == 'abacus':
                    stru_data = abacus_scf.get_abacus_STRU(equi_contcar)  
                    vol_to_poscar = abs(np.linalg.det(stru_data['cells'])) / np.array(stru_data['atom_numbs']).sum()
                else:
                    vol_to_poscar = vasp.poscar_vol(equi_contcar) / vasp.poscar_natoms(equi_contcar)
                self.parameter['scale2equi'] = []

                task_num = 0
                while self.vol_start + self.vol_step * task_num < self.vol_end:
                # for vol in np.arange(int(self.vol_start * 100), int(self.vol_end * 100), int(self.vol_step * 100)):
                    # vol = vol / 100.0
                    vol = self.vol_start + task_num * self.vol_step
                    #task_num = int((vol - self.vol_start) / self.vol_step)
                    output_task = os.path.join(path_to_work, 'task.%06d' % task_num)
                    os.makedirs(output_task, exist_ok=True)
                    os.chdir(output_task)
                    if self.inter_param['type'] == 'abacus':
                        POSCAR = 'STRU'
                        POSCAR_orig = 'STRU.orig'
                        scale_func = abacus.stru_scale
                    else:
                        POSCAR = 'POSCAR'
                        POSCAR_orig = 'POSCAR.orig'
                        scale_func = vasp.poscar_scale

                    for ii in ['INCAR', 'POTCAR', POSCAR_orig, POSCAR, 'conf.lmp', 'in.lammps']:
                        if os.path.exists(ii):
                            os.remove(ii)
                    task_list.append(output_task)
                    os.symlink(os.path.relpath(equi_contcar), POSCAR_orig)
                    # scale = (vol / vol_to_poscar) ** (1. / 3.)

                    if self.vol_abs :
                        scale = (vol / vol_to_poscar) ** (1. / 3.)
                        eos_params = {'volume': vol, 'scale': scale}
                    else :
                        scale = vol ** (1. / 3.)
                        eos_params = {'volume': vol * vol_to_poscar, 'scale': scale}
                    dumpfn(eos_params, 'eos.json', indent=4)
                    self.parameter['scale2equi'].append(scale)  # 06/22
                    scale_func(POSCAR_orig,POSCAR,scale)
                    task_num += 1
                os.chdir(cwd)
        return task_list

    def post_process(self, task_list):
        pass

    def task_type(self):
        return self.parameter['type']

    def task_param(self):
        return self.parameter

    def _compute_lower(self,
                       output_file,
                       all_tasks,
                       all_res):
        output_file = os.path.abspath(output_file)
        res_data = {}
        ptr_data = "conf_dir: " + os.path.dirname(output_file) + "\n"
        if not self.reprod:
            ptr_data += ' VpA(A^3)  EpA(eV)\n'
            for ii in range(len(all_tasks)):
                # vol = self.vol_start + ii * self.vol_step
                vol = loadfn(os.path.join(all_tasks[ii], 'eos.json'))['volume']
                task_result = loadfn(all_res[ii])
                res_data[vol] = task_result['energies'][-1] / sum(task_result['atom_numbs'])
                ptr_data += '%7.3f  %8.4f \n' % (vol, task_result['energies'][-1] / sum(task_result['atom_numbs']))
                # res_data[vol] = all_res[ii]['energy'] / len(all_res[ii]['force'])
                # ptr_data += '%7.3f  %8.4f \n' % (vol, all_res[ii]['energy'] / len(all_res[ii]['force']))

        else:
            if 'init_data_path' not in self.parameter:
                raise RuntimeError("please provide the initial data path to reproduce")
            init_data_path = os.path.abspath(self.parameter['init_data_path'])
            res_data, ptr_data = post_repro(init_data_path, self.parameter['init_from_suffix'],
                                            all_tasks, ptr_data, self.parameter.get('reprod_last_frame', True))

        with open(output_file, 'w') as fp:
            json.dump(res_data, fp, indent=4)

        return res_data, ptr_data
