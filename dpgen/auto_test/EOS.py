import glob
import json
import os

import numpy as np
from monty.serialization import loadfn, dumpfn

import dpgen.auto_test.lib.vasp as vasp
from dpgen import dlog
from dpgen.auto_test.Property import Property
from dpgen.auto_test.refine import make_refine
from dpgen.auto_test.reproduce import make_repro
from dpgen.auto_test.reproduce import post_repro


class EOS(Property):
    def __init__(self,
                 parameter):
        parameter['reprod-opt'] = parameter.get('reprod-opt', False)
        self.reprod = parameter['reprod-opt']
        if not self.reprod:
            self.vol_start = parameter['vol_start']
            self.vol_end = parameter['vol_end']
            self.vol_step = parameter['vol_step']
            parameter['cal_type'] = parameter.get('cal_type', 'relaxation')
            self.cal_type = parameter['cal_type']
            default_cal_setting = {"relax_pos": True,
                                   "relax_shape": True,
                                   "relax_vol": False}
            if 'cal_setting' not in parameter:
                parameter['cal_setting'] = default_cal_setting
            elif "relax_pos" not in parameter['cal_setting']:
                parameter['cal_setting']['relax_pos'] = default_cal_setting['relax_pos']
            elif "relax_shape" not in parameter['cal_setting']:
                parameter['cal_setting']['relax_shape'] = default_cal_setting['relax_shape']
            elif "relax_vol" not in parameter['cal_setting']:
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
            elif "relax_pos" not in parameter['cal_setting']:
                parameter['cal_setting']['relax_pos'] = default_cal_setting['relax_pos']
            elif "relax_shape" not in parameter['cal_setting']:
                parameter['cal_setting']['relax_shape'] = default_cal_setting['relax_shape']
            elif "relax_vol" not in parameter['cal_setting']:
                parameter['cal_setting']['relax_vol'] = default_cal_setting['relax_vol']
            self.cal_setting = parameter['cal_setting']
            parameter['init_from_suffix'] = parameter.get('init_from_suffix', '00')
            self.init_from_suffix = parameter['init_from_suffix']
        self.parameter = parameter

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
                                                        struct_output_name, 'relaxation'))

        cwd = os.getcwd()
        task_list = []
        if refine:
            print('EOS refine starts')
            task_list = make_refine(self.parameter['init_from_suffix'],
                                    self.parameter['output_suffix'],
                                    path_to_work)
            os.chdir(cwd)
        if self.reprod:
            print('eos reproduce starts')
            if 'vasp_lmp_path' not in self.parameter:
                raise RuntimeError("please provide the vasp_lmp_path for reproduction")
            vasp_lmp_path = os.path.abspath(self.parameter['vasp_lmp_path'])
            task_list = make_repro(vasp_lmp_path, self.init_from_suffix, path_to_work)
            os.chdir(cwd)
        else:
            print('gen eos from ' + str(self.vol_start) + ' to ' + str(self.vol_end) + ' by every ' + str(self.vol_step))
            equi_contcar = os.path.join(path_to_equi, 'CONTCAR')
            if not os.path.exists(equi_contcar):
                raise RuntimeError("please do relaxation first")
            vol_to_poscar = vasp.poscar_vol(equi_contcar) / vasp.poscar_natoms(equi_contcar)
            self.parameter['scale2equi'] = []
            for vol in np.arange(self.vol_start, self.vol_end, self.vol_step):
                task_num = (vol - self.vol_start) / self.vol_step
                output_task = os.path.join(path_to_work, 'task.%06d' % task_num)
                os.makedirs(output_task, exist_ok=True)
                os.chdir(output_task)
                for ii in ['INCAR', 'POTCAR', 'POSCAR.orig', 'POSCAR', 'conf.lmp', 'in.lammps']:
                    if os.path.exists(ii):
                        os.remove(ii)
                task_list.append(output_task)
                os.symlink(os.path.relpath(equi_contcar), 'POSCAR.orig')
                # scale = (vol / vol_to_poscar) ** (1. / 3.)
                scale = vol ** (1. / 3.)
                eos_params = {'volume': vol * vol_to_poscar, 'scale': scale}
                dumpfn(eos_params, 'eos.json', indent=4)
                self.parameter['scale2equi'].append(scale)  # 06/22
                vasp.poscar_scale('POSCAR.orig', 'POSCAR', scale)
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
                res_data[vol] = task_result['energies'][-1] / task_result['atom_numbs'][0]
                ptr_data += '%7.3f  %8.4f \n' % (vol, task_result['energies'][-1] / task_result['atom_numbs'][0])
                # res_data[vol] = all_res[ii]['energy'] / len(all_res[ii]['force'])
                # ptr_data += '%7.3f  %8.4f \n' % (vol, all_res[ii]['energy'] / len(all_res[ii]['force']))

        else:
            if 'vasp_lmp_path' not in self.parameter:
                raise RuntimeError("please provide the vasp_lmp_path for reproduction")
            vasp_lmp_path = os.path.abspath(self.parameter['vasp_lmp_path'])
            res_data, ptr_data = post_repro(vasp_lmp_path, self.parameter['init_from_suffix'], all_tasks, ptr_data)

        with open(output_file, 'w') as fp:
            json.dump(res_data, fp, indent=4)

        return res_data, ptr_data
