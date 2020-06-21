from dpgen.auto_test.Property import Property
from dpgen.auto_test.refine import make_refine
from dpgen.auto_test import reproduce
import dpgen.auto_test.lib.vasp as vasp
import numpy as np
import os, json


class EOS(Property):
    def __init__(self,
                 parameter):
        self.parameter = parameter
        self.vol_start = parameter['vol_start']
        self.vol_end = parameter['vol_end']
        self.vol_step = parameter['vol_step']
        self.change_box = parameter.get('change_box', True)
        self.reprod = parameter.get('reprod-opt', False)

    def make_confs(self,
                   path_to_work,
                   path_to_equi,
                   refine=False):
        path_to_work = os.path.abspath(path_to_work)
        path_to_equi = os.path.abspath(path_to_equi)
        cwd = os.getcwd()
        task_list = []
        if refine:
            task_list = make_refine(self.parameter['init_from_suffix'],
                                    self.parameter['output_suffix'],
                                    path_to_work,
                                    (self.vol_end - self.vol_start) / self.vol_step)
            os.chdir(cwd)
        if self.reprod:
            if 'vasp_lmp_path' not in self.parameter:
                raise RuntimeError("please provide the vasp_lmp_path for reproduction")
            vasp_lmp_path = os.path.abspath(self.parameter['vasp_lmp_path'])
            task_list = reproduce.make_repro(vasp_lmp_path, path_to_work)
            os.chdir(cwd)
        else:
            equi_contcar = os.path.join(path_to_equi, 'CONTCAR')
            if not os.path.exists(equi_contcar):
                raise RuntimeError("please do relaxation first")
            vol_to_poscar = vasp.poscar_vol(equi_contcar) / vasp.poscar_natoms(equi_contcar)
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
                scale = (vol / vol_to_poscar) ** (1. / 3.)
                self.parameter['scale2equi'] = scale  # 06/14
                vasp.poscar_scale('POSCAR.orig', 'POSCAR', scale)
            os.chdir(cwd)
        return task_list

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
                vol = self.vol_start + ii * self.vol_step
                res_data[vol] = all_res[ii]['energy'] / len(all_res[ii]['force']) * 3
                ptr_data += '%7.3f  %8.4f \n' % (vol, all_res[ii]['energy'] / len(all_res[ii]['force']) * 3)

        else:
            if 'vasp_lmp_path' not in self.parameter:
                raise RuntimeError("please provide the vasp_lmp_path for reproduction")
            vasp_lmp_path = os.path.abspath(self.parameter['vasp_lmp_path'])
            res_data, ptr_data = reproduce.post_repro(vasp_lmp_path, all_tasks, ptr_data)

        with open(output_file, 'w') as fp:
            json.dump(res_data, fp, indent=4)

        return res_data, ptr_data
