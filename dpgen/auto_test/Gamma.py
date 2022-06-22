import glob
import json
import os
import re

import numpy as np
from monty.serialization import loadfn, dumpfn
from pymatgen.core.structure import Structure
from pymatgen.core.surface import SlabGenerator

import dpgen.auto_test.lib.vasp as vasp
from dpgen import dlog
from dpgen.auto_test.Property import Property
from dpgen.auto_test.refine import make_refine
from dpgen.auto_test.reproduce import make_repro
from dpgen.auto_test.reproduce import post_repro


class Gamma(Property):
    """
    Calculation of gamma line for bcc and fcc (v1.1 add half z judgement)
    """

    def __init__(self,
                 parameter):
        parameter['reproduce'] = parameter.get('reproduce', False)
        self.reprod = parameter['reproduce']
        if not self.reprod:
            if not ('init_from_suffix' in parameter and 'output_suffix' in parameter):
                self.miller_index = parameter['miller_index']
                parameter['min_slab_size'] = parameter.get('min_slab_size', 10)
                self.min_slab_size = parameter['min_slab_size']
                parameter['min_supercell_size'] = parameter.get('min_supercell_size', (5,5,10))
                self.min_supercell_size = parameter['min_supercell_size']
                self.min_vacuum_size = parameter['min_vacuum_size']
                parameter['n_steps'] = parameter.get('n_steps', 10)
                self.n_steps = parameter['n_steps']
                self.atom_num = None
            #                parameter['pert_xz'] = parameter.get('pert_xz', 0.01)
            #                self.pert_xz = parameter['pert_xz']
            #                default_max_miller = 2
            #                parameter['max_miller'] = parameter.get('max_miller', default_max_miller)
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

        task_list = []
        cwd = os.getcwd()

        if self.reprod:
            print('surface reproduce starts')
            if 'init_data_path' not in self.parameter:
                raise RuntimeError("please provide the initial data path to reproduce")
            init_data_path = os.path.abspath(self.parameter['init_data_path'])
            task_list = make_repro(init_data_path, self.init_from_suffix,
                                   path_to_work, self.parameter.get('reprod_last_frame', True))
            os.chdir(cwd)

        else:
            if refine:
                print('surface refine starts')
                task_list = make_refine(self.parameter['init_from_suffix'],
                                        self.parameter['output_suffix'],
                                        path_to_work)
                os.chdir(cwd)
                # record miller
                init_from_path = re.sub(self.parameter['output_suffix'][::-1],
                                        self.parameter['init_from_suffix'][::-1],
                                        path_to_work[::-1], count=1)[::-1]
                task_list_basename = list(map(os.path.basename, task_list))

                for ii in task_list_basename:
                    init_from_task = os.path.join(init_from_path, ii)
                    output_task = os.path.join(path_to_work, ii)
                    os.chdir(output_task)
                    if os.path.isfile('miller.json'):
                        os.remove('miller.json')
                    if os.path.islink('miller.json'):
                        os.remove('miller.json')
                    os.symlink(os.path.relpath(os.path.join(init_from_task, 'miller.json')), 'miller.json')
                os.chdir(cwd)

            else:
                equi_contcar = os.path.join(path_to_equi, 'CONTCAR')
                if not os.path.exists(equi_contcar):
                    raise RuntimeError("please do relaxation first")
                ptypes = vasp.get_poscar_types(equi_contcar)
                # read structure from relaxed CONTCAR
                ss = Structure.from_file(equi_contcar)
                # rewrite new CONTCAR with direct coords
                os.chdir(path_to_equi)
                ss.to('POSCAR', 'CONTCAR.direct')
                # re-read new CONTCAR
                ss = Structure.from_file('CONTCAR.direct')
                # gen initial slab for fp calculation
                slabGen_fp = SlabGenerator(ss, miller_index=self.miller_index,
                                           min_slab_size=self.min_slab_size,
                                           min_vacuum_size=self.min_vacuum_size,
                                           center_slab=True, in_unit_planes=True,
                                           lll_reduce=False, primitive=True,
                                           max_normal_search=5)
                slab_fp = slabGen_fp.get_slab()
                # gen initial slab for md calculation
                slabGen_md = SlabGenerator(ss, miller_index=self.miller_index,
                                           min_slab_size=self.min_supercell_size[2],
                                           min_vacuum_size=0, center_slab=True,
                                           in_unit_planes=True, lll_reduce=True,
                                           max_normal_search=5, primitive=False,
                                           reorient_lattice=False)
                slab_md = slabGen_md.get_slab()
                # make supercell for md calculation
                slab_md.make_supercell(scaling_matrix=[self.min_supercell_size[0],
                                                       self.min_supercell_size[1], 1])
                # define displace vectors
                disp_vector_fp = (1, 0, 0)
                disp_vector_md = (0.5/self.min_supercell_size[0], 0.5/self.min_supercell_size[1], 0)
                # displace structure
                all_slabs_fp = self.__displace_slab(slab_fp, disp_vector=disp_vector_fp)
                all_slabs_md = self.__displace_slab(slab_md, disp_vector=disp_vector_md)
                self.atom_num = len(all_slabs_fp[0].sites)

                os.chdir(path_to_work)
                if os.path.isfile('POSCAR'):
                    os.remove('POSCAR')
                if os.path.islink('POSCAR'):
                    os.remove('POSCAR')
                os.symlink(os.path.relpath(equi_contcar), 'POSCAR')
                #           task_poscar = os.path.join(output, 'POSCAR')
                for ii in range(len(all_slabs_fp)):
                    output_task = os.path.join(path_to_work, 'task.%06d' % ii)
                    os.makedirs(output_task, exist_ok=True)
                    os.chdir(output_task)
                    for jj in ['INCAR', 'POTCAR', 'POSCAR', 'conf.lmp', 'in.lammps']:
                        if os.path.exists(jj):
                            os.remove(jj)
                    task_list.append(output_task)
                    print("# %03d generate " % ii, output_task)
                    # print("# %03d generate " % ii, output_task, " \t %d atoms" % self.atom_num)
                    # make confs
                    all_slabs_fp[ii].to('POSCAR', 'POSCAR.fp')
                    all_slabs_md[ii].to('POSCAR', 'POSCAR.tmp')
                    vasp.regulate_poscar('POSCAR.tmp', 'POSCAR')
                    vasp.sort_poscar('POSCAR', 'POSCAR', ptypes)
                    # vasp.perturb_xz('POSCAR', 'POSCAR', self.pert_xz)
                    # record miller
                    dumpfn(all_slabs_fp[ii].miller_index, 'miller.json')
                os.chdir(cwd)

        return task_list

    def __direction_dict(self):
        pass

    def __displace_slab(self,
                        slab, disp_vector):
        """
        return a list of displaced slab objects
        """
        all_slabs = [slab.copy()]
        for ii in list(range(self.n_steps)):
            frac_disp = 1 / self.n_steps
            unit_vector = frac_disp * np.array(disp_vector)
            # return list of atoms number to be displaced which above 0.5 z
            disp_atoms_list = np.where(slab.frac_coords[:,2]>0.5)[0]
            slab.translate_sites(indices=disp_atoms_list, vector=unit_vector,
                                 frac_coords=True, to_unit_cell=True)
            all_slabs.append(slab.copy())
        return all_slabs

    def __pos_fix(self,
                  poscar):
        """
        add position fix condition of x and y in POSCAR
        """
        insert_pos = -self.atom_num
        with open(poscar, 'r') as fin1:
            contents = fin1.readlines()
            contents.insert(insert_pos - 1, 'Selective dynamics\n')
            for ii in range(insert_pos, 0, 1):
                contents[ii] = contents[ii].replace('\n', '')
                contents[ii] += ' ' + 'F F T' + '\n'
        with open(poscar, 'w') as fin2:
            for ii in range(len(contents)):
                fin2.write(contents[ii])

    def post_process(self,
                     task_list):
        if True:
            for ii in task_list:
                inter = os.path.join(ii, 'inter.json')
                poscar_fp = os.path.join(ii, 'POSCAR.fp')
                poscar = os.path.join(ii, 'POSCAR')
                calc_type = loadfn(inter)['type']
                if calc_type == 'vasp':
                    self.__pos_fix(poscar_fp)
                    os.renames(poscar_fp, poscar)
                else:
                    os.remove(poscar_fp)

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
            ptr_data += "No_steps: \tStacking_Fault_E(J/m^2) EpA(eV) equi_EpA(eV)\n"
            for ii in all_tasks:
                task_result = loadfn(os.path.join(ii, 'result_task.json'))
                natoms = np.sum(task_result['atom_numbs'])
                epa = task_result['energies'][-1] / natoms
                AA = np.linalg.norm(np.cross(task_result['cells'][0][0], task_result['cells'][0][1]))

                equi_path = os.path.abspath(os.path.join(os.path.dirname(output_file), '../relaxation/relax_task'))
                equi_result = loadfn(os.path.join(equi_path, 'result.json'))
                equi_epa = equi_result['energies'][-1] / np.sum(equi_result['atom_numbs'])
                structure_dir = os.path.basename(ii)

                Cf = 1.60217657e-16 / (1e-20 * 2) * 0.001
                evac = (task_result['energies'][-1] - equi_epa * natoms) / AA * Cf

                miller_index = loadfn(os.path.join(ii, 'miller.json'))
                ptr_data += "%-25s     %7.3f    %8.3f %8.3f\n" % (
                    str(miller_index) + '-' + structure_dir + ':', evac, epa, equi_epa)
                res_data[str(miller_index) + '-' + structure_dir] = [evac, epa, equi_epa]

        else:
            if 'init_data_path' not in self.parameter:
                raise RuntimeError("please provide the initial data path to reproduce")
            init_data_path = os.path.abspath(self.parameter['init_data_path'])
            res_data, ptr_data = post_repro(init_data_path, self.parameter['init_from_suffix'],
                                            all_tasks, ptr_data, self.parameter.get('reprod_last_frame', True))

        with open(output_file, 'w') as fp:
            json.dump(res_data, fp, indent=4)

        return res_data, ptr_data
