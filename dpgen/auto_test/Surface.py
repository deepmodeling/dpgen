from dpgen.auto_test.Property import Property
from dpgen.auto_test.refine import make_refine
from dpgen.auto_test import reproduce
import dpgen.auto_test.lib.vasp as vasp
import dpgen.auto_test.lib.lammps as lammps
from pymatgen.core.structure import Structure
from pymatgen.core.surface import generate_all_slabs
from monty.serialization import loadfn, dumpfn
from dpgen import dlog
import numpy as np
import os, json


class Surface(Property):
    def __init__(self,
                 parameter):
        parameter['reprod-opt'] = parameter.get('reprod-opt', False)
        self.reprod = parameter['reprod-opt']
        if not self.reprod:
            self.min_slab_size = parameter['min_slab_size']
            self.min_vacuum_size = parameter['min_vacuum_size']
            self.pert_xz = parameter['pert_xz']
            default_max_miller = 2
            self.miller = parameter.get('max_miller', default_max_miller)
            parameter['cal_type'] = parameter.get('cal_type', 'relaxation')
            self.cal_type = parameter['cal_type']
            default_cal_setting = {"relax_pos": True,
                                   "relax_shape": True,
                                   "relax_vol": False}
            parameter['cal_setting'] = parameter.get('cal_setting', default_cal_setting)
            self.cal_setting = parameter['cal_setting']
        else:
            parameter['cal_type'] = 'static'
            self.cal_type = parameter['cal_type']
            parameter['cal_setting'] = {"relax_pos": False,
                                        "relax_shape": False,
                                        "relax_vol": False}
            self.cal_setting = parameter['cal_setting']
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
        task_list = []
        cwd = os.getcwd()

        if self.reprod:
            print('surface reproduce starts')
            if 'vasp_lmp_path' not in self.parameter:
                raise RuntimeError("please provide the vasp_lmp_path for reproduction")
            vasp_lmp_path = os.path.abspath(self.parameter['vasp_lmp_path'])
            task_list = reproduce.make_repro(vasp_lmp_path, path_to_work)
            os.chdir(cwd)

        else:
            equi_contcar = os.path.join(path_to_equi, 'CONTCAR')
            if not os.path.exists(equi_contcar):
                raise RuntimeError("please do relaxation first")
            ptypes = vasp.get_poscar_types(equi_contcar)
            # gen structure
            ss = Structure.from_file(equi_contcar)
            # gen slabs
            all_slabs = generate_all_slabs(ss, self.miller, self.min_slab_size, self.min_vacuum_size)

            if refine:
                print('surface refine starts')
                task_list = make_refine(self.parameter['init_from_suffix'],
                                        self.parameter['output_suffix'],
                                        path_to_work)
                # record miller
                for ii in range(len(task_list)):
                    os.chdir(task_list[ii])
                    dumpfn(all_slabs[ii].miller_index, 'miller.json')
                os.chdir(cwd)

            else:
                os.chdir(path_to_work)
                if os.path.isfile('POSCAR'):
                    os.remove('POSCAR')
                os.symlink(os.path.relpath(equi_contcar), 'POSCAR')
                #           task_poscar = os.path.join(output, 'POSCAR')
                for ii in range(len(all_slabs)):
                    output_task = os.path.join(path_to_work, 'task.%06d' % ii)
                    os.makedirs(output_task, exist_ok=True)
                    os.chdir(output_task)
                    for jj in ['INCAR', 'POTCAR', 'POSCAR', 'conf.lmp', 'in.lammps']:
                        if os.path.exists(jj):
                            os.remove(jj)
                    task_list.append(output_task)
                    print("# %03d generate " % ii, output_task, " \t %d atoms" % len(all_slabs[ii].sites))
                    # make confs
                    all_slabs[ii].to('POSCAR', 'POSCAR.tmp')
                    vasp.regulate_poscar('POSCAR.tmp', 'POSCAR')
                    vasp.sort_poscar('POSCAR', 'POSCAR', ptypes)
                    vasp.perturb_xz('POSCAR', 'POSCAR', self.pert_xz)
                    # record miller
                    dumpfn(all_slabs[ii].miller_index, 'miller.json')
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
        ptr_data = os.path.dirname(output_file) + '\n'

        if not self.reprod:
            ptr_data += "Miller_Indices: \tSurf_E(J/m^2) EpA(eV) equi_EpA(eV)\n"
            for ii in all_tasks:
                task_result = loadfn(os.path.join(ii, 'result_task.json'))
                natoms = task_result['atom_numbs'][0]
                epa = task_result['energies'][-1] / task_result['atom_numbs'][0]
                AA = np.linalg.norm(np.cross(task_result['cells'][0][0], task_result['cells'][0][1]))

                equi_path = os.path.abspath(os.path.join(os.path.dirname(output_file), '../relaxation'))
                equi_result = loadfn(os.path.join(equi_path, 'result.json'))
                equi_epa = equi_result['energies'][-1] / equi_result['atom_numbs'][0]
                structure_dir = os.path.basename(ii)

                Cf = 1.60217657e-16 / (1e-20 * 2) * 0.001
                evac = (task_result['energies'][-1] - equi_epa * natoms) / AA * Cf

                miller_index = loadfn(os.path.join(ii, 'miller.json'))
                ptr_data += "%-25s     %7.3f    %8.3f %8.3f\n" % (
                    str(miller_index) + '-' + structure_dir + ':', evac, epa, equi_epa)
                res_data[str(miller_index) + '-' + structure_dir] = [evac, epa, equi_epa]

        else:
            if 'vasp_lmp_path' not in self.parameter:
                raise RuntimeError("please provide the vasp_lmp_path for reproduction")
            vasp_lmp_path = os.path.abspath(self.parameter['vasp_lmp_path'])
            res_data, ptr_data = reproduce.post_repro(vasp_lmp_path, all_tasks, ptr_data)

        with open(output_file, 'w') as fp:
            json.dump(res_data, fp, indent=4)

        return res_data, ptr_data
