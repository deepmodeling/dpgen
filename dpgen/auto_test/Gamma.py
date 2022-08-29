import glob
import json
import os
import re

import dpdata
import numpy as np
from monty.serialization import loadfn, dumpfn
from pymatgen.core.structure import Structure
from pymatgen.core.surface import SlabGenerator
from pymatgen.io.ase import AseAtomsAdaptor
from ase.lattice.cubic import BodyCenteredCubic as bcc
from ase.lattice.cubic import FaceCenteredCubic as fcc
from ase.lattice.hexagonal import HexagonalClosedPacked as hcp

import dpgen.auto_test.lib.vasp as vasp
from dpgen import dlog
from dpgen.auto_test.Property import Property
from dpgen.auto_test.refine import make_refine
from dpgen.auto_test.reproduce import make_repro
from dpgen.auto_test.reproduce import post_repro

import dpgen.auto_test.lib.abacus as abacus


class Gamma(Property):
    """
    Calculation of common gamma lines for bcc and fcc
    """
    def __init__(self,
                 parameter,inter_param=None):
        parameter['reproduce'] = parameter.get('reproduce', False)
        self.reprod = parameter['reproduce']
        if not self.reprod:
            if not ('init_from_suffix' in parameter and 'output_suffix' in parameter):
                self.miller_index = parameter['miller_index']
                self.displace_direction = parameter['displace_direction']
                self.lattice_type = parameter['lattice_type']
                parameter['supercell_size'] = parameter.get('supercell_size', (1,1,5))
                self.supercell_size = parameter['supercell_size']
                parameter['min_vacuum_size'] = parameter.get('min_vacuum_size', 20)
                self.min_vacuum_size = parameter['min_vacuum_size']
                parameter['add_fix'] = parameter.get('add_fix', ['true','true','false']) # standard method
                self.add_fix = parameter['add_fix']
                parameter['n_steps'] = parameter.get('n_steps', 10)
                self.n_steps = parameter['n_steps']
                self.atom_num = None
            parameter['cal_type'] = parameter.get('cal_type', 'relaxation')
            self.cal_type = parameter['cal_type']
            default_cal_setting = {"relax_pos": True,
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

        task_list = []
        cwd = os.getcwd()

        if self.reprod:
            print('gamma line reproduce starts')
            if 'init_data_path' not in self.parameter:
                raise RuntimeError("please provide the initial data path to reproduce")
            init_data_path = os.path.abspath(self.parameter['init_data_path'])
            task_list = make_repro(init_data_path, self.init_from_suffix,
                                   path_to_work, self.parameter.get('reprod_last_frame', True))
            os.chdir(cwd)

        else:
            if refine:
                print('gamma line refine starts')
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
                if self.inter_param['type'] == 'abacus':
                    CONTCAR = abacus.final_stru(path_to_equi)
                    POSCAR = 'STRU'
                else:
                    CONTCAR = 'CONTCAR'
                    POSCAR = 'POSCAR'

                equi_contcar = os.path.join(path_to_equi, CONTCAR)
                if not os.path.exists(equi_contcar):
                    raise RuntimeError("please do relaxation first")
                print('we now only support gamma line calculation for BCC and FCC metals')
                print('supported slip systems are planes/direction: 100/010, 110/111, 111/110, 111/112, 112/111, and 123/111')

                if self.inter_param['type'] == 'abacus':
                    stru = dpdata.System(equi_contcar, fmt="stru")
                    stru.to('contcar','CONTCAR.tmp')
                    ptypes = vasp.get_poscar_types('CONTCAR.tmp')
                    ss = Structure.from_file('CONTCAR.tmp')
                    os.remove('CONTCAR.tmp')
                else:
                    ptypes = vasp.get_poscar_types(equi_contcar)
                    # read structure from relaxed CONTCAR
                    ss = Structure.from_file(equi_contcar)

                # rewrite new CONTCAR with direct coords
                os.chdir(path_to_equi)
                ss.to('POSCAR', 'CONTCAR.direct')
                # re-read new CONTCAR
                ss = Structure.from_file('CONTCAR.direct')
                relax_a = ss.lattice.a
                relax_b = ss.lattice.b
                relax_c = ss.lattice.c
                # gen initial slab
                slab = self.__gen_slab_ase(symbol=ptypes[0],
                                           lat_param=[relax_a,relax_b,relax_c])
                # define displace vectors
                disp_vector = (1/self.supercell_size[0], 0, 0)
                # displace structure
                all_slabs = self.__displace_slab(slab, disp_vector=disp_vector)
                self.atom_num = len(all_slabs[0].sites)

                os.chdir(path_to_work)
                if os.path.isfile(POSCAR):
                    os.remove(POSCAR)
                if os.path.islink(POSCAR):
                    os.remove(POSCAR)
                os.symlink(os.path.relpath(equi_contcar), POSCAR)
                #           task_poscar = os.path.join(output, 'POSCAR')
                for ii in range(len(all_slabs)):
                    output_task = os.path.join(path_to_work, 'task.%06d' % ii)
                    os.makedirs(output_task, exist_ok=True)
                    os.chdir(output_task)
                    for jj in ['INCAR', 'POTCAR', POSCAR, 'conf.lmp', 'in.lammps']:
                        if os.path.exists(jj):
                            os.remove(jj)
                    task_list.append(output_task)
                    #print("# %03d generate " % ii, output_task)
                    print("# %03d generate " % ii, output_task, " \t %d atoms" % self.atom_num)
                    # make confs
                    all_slabs[ii].to('POSCAR', 'POSCAR.tmp')
                    vasp.regulate_poscar('POSCAR.tmp', 'POSCAR')
                    vasp.sort_poscar('POSCAR', 'POSCAR', ptypes)
                    if self.inter_param['type'] == 'abacus':
                        abacus.poscar2stru("POSCAR",self.inter_param,"STRU")
                        os.remove('POSCAR')
                    # vasp.perturb_xz('POSCAR', 'POSCAR', self.pert_xz)
                    # record miller
                    dumpfn(self.miller_index, 'miller.json')
                os.chdir(cwd)

        return task_list

    @staticmethod
    def centralize_slab(slab) -> None:
        z_pos_list = list(set([site.position[2] for site in slab]))
        z_pos_list.sort()
        central_atoms = (z_pos_list[-1] - z_pos_list[0])/2
        central_cell = slab.cell[2][2]/2
        disp_length = central_cell - central_atoms
        for site in slab:
            site.position[2] += disp_length

    def return_direction(self):
        miller_str = ''
        direct_str = ''
        for ii in range(len(self.miller_index)):
            miller_str += str(self.miller_index[ii])
        for ii in range(len(self.displace_direction)):
            direct_str += str(self.displace_direction[ii])
        search_key = miller_str + '/' + direct_str
        # define specific cell vectors
        dict_directions = {
            '100/010': [(0,1,0), (0,0,1), (1,0,0)],
            '110/111': [(-1,1,1), (1,-1,1), (1,1,0)],
            '111/110': [(-1,1,0), (-1,-1,2), (1,1,1)],
            '111/112': [(1,1,-2), (-1,1,0), (1,1,1)],
            '112/111': [(-1,-1,1), (1,-1,0), (1,1,2)],
            '123/111': [(-1,-1,1), (2,-1,0), (1,2,3)]
        }
        try:
            directions = dict_directions[search_key]
        except KeyError:
            raise RuntimeError(f'Unsupported input combination of miller index and displacement direction: '
                               f'{miller_str}:{direct_str}')
        return directions

    def __gen_slab_ase(self,
                       symbol, lat_param):
        if not self.lattice_type:
            raise RuntimeError('Error! Please provide the input lattice type!')
        elif self.lattice_type == 'bcc':
            slab_ase = bcc(symbol=symbol, size=self.supercell_size, latticeconstant=lat_param[0],
                           directions=self.return_direction())
        elif self.lattice_type == 'fcc':
            slab_ase = fcc(symbol=symbol, size=self.supercell_size, latticeconstant=lat_param[0],
                           directions=self.return_direction())
        elif self.lattice_type == 'hcp':
            pass
        else:
            raise RuntimeError(f'unsupported lattice type: {self.lattice_type}')
        self.centralize_slab(slab_ase)
        if self.min_vacuum_size > 0:
            slab_ase.center(vacuum=self.min_vacuum_size/2, axis=2)
        slab_pymatgen = AseAtomsAdaptor.get_structure(slab_ase)
        return slab_pymatgen

    # leave this function to later use
    #def __gen_slab_pmg(self,
    #                   pmg_struc):
    #    slabGen = SlabGenerator(pmg_struc, miller_index=self.miller_index,
    #                            min_slab_size=self.supercell_size[2],
    #                            min_vacuum_size=self.min_vacuum_size,
    #                            center_slab=True, in_unit_planes=True, lll_reduce=False,
    #                            primitive=True, max_normal_search=5)
    #    slab_pmg = slabGen.get_slab()
    #    slab_pmg.make_supercell(scaling_matrix=[self.supercell_size[0],self.supercell_size[1],1])
    #    return slab_pmg

    def __displace_slab(self,
                        slab, disp_vector):
        # return a list of displaced slab objects
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

    def __poscar_fix(self, poscar) -> None:
        # add position fix condition of x and y in POSCAR
        insert_pos = -self.atom_num
        fix_dict = {
            'true': 'F',
            'false': 'T'
        }
        add_fix_str = ' ' + fix_dict[self.add_fix[0]] + \
                      ' ' + fix_dict[self.add_fix[1]] + \
                      ' ' + fix_dict[self.add_fix[2]] + '\n'
        with open(poscar, 'r') as fin1:
            contents = fin1.readlines()
            contents.insert(insert_pos-1, 'Selective dynamics\n')
            for ii in range(insert_pos, 0, 1):
                contents[ii] = contents[ii].replace('\n', '')
                contents[ii] += add_fix_str
        with open(poscar, 'w') as fin2:
            for ii in range(len(contents)):
                fin2.write(contents[ii])

    def __stru_fix(self,stru) -> None:
        fix_dict = {
            'true': True, 
            'false': False
        }
        fix_xyz = [fix_dict[i] for i in self.addfix]
        abacus.stru_fix_atom(stru,fix_atom=fix_xyz)

    def __inLammpes_fix(self, inLammps) -> None:
        # add position fix condition of x and y of in.lammps
        fix_dict = {
            'true': '0',
            'false': 'NULL'
        }
        add_fix_str = 'fix             1 all setforce' + \
                      ' ' + fix_dict[self.add_fix[0]] + \
                      ' ' + fix_dict[self.add_fix[1]] + \
                      ' ' + fix_dict[self.add_fix[2]] + '\n'
        with open(inLammps, 'r') as fin1:
            contents = fin1.readlines()
            for ii in range(len(contents)):
                upper = re.search("variable        N equal count\(all\)", contents[ii])
                lower = re.search("min_style       cg", contents[ii])
                if lower:
                    lower_id = ii
                    #print(lower_id)
                elif upper:
                    upper_id = ii
                    #print(upper_id)
            del contents[lower_id+1:upper_id-1]
            contents.insert(lower_id+1, add_fix_str)
        with open(inLammps, 'w') as fin2:
            for ii in range(len(contents)):
                fin2.write(contents[ii])

    def post_process(self,
                     task_list):
        if self.add_fix:
            count = 0
            for ii in task_list:
                count += 1
                inter = os.path.join(ii, 'inter.json')
                poscar = os.path.join(ii, 'POSCAR')
                calc_type = loadfn(inter)['type']
                if calc_type == 'vasp':
                    self.__poscar_fix(poscar)
                elif calc_type == 'abacus':
                    self.__stru_fix(os.path.join(ii, 'STRU'))
                else:
                    inLammps = os.path.join(ii, 'in.lammps')
                    if count == 1:
                        self.__inLammpes_fix(inLammps)


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
            ptr_data += str(tuple(self.miller_index)) + ' plane along ' + str(self.displace_direction)
            ptr_data += "No_task: \tDisplacement \tStacking_Fault_E(J/m^2) EpA(eV) slab_equi_EpA(eV)\n"
            all_tasks.sort()
            task_result_slab_equi = loadfn(os.path.join(all_tasks[0], 'result_task.json'))
            for ii in all_tasks:
                task_result = loadfn(os.path.join(ii, 'result_task.json'))
                natoms = np.sum(task_result['atom_numbs'])
                epa = task_result['energies'][-1] / natoms
                equi_epa_slab = task_result_slab_equi['energies'][-1] / natoms
                AA = np.linalg.norm(np.cross(task_result['cells'][0][0], task_result['cells'][0][1]))

                equi_path = os.path.abspath(os.path.join(os.path.dirname(output_file), '../relaxation/relax_task'))
                equi_result = loadfn(os.path.join(equi_path, 'result.json'))
                equi_epa = equi_result['energies'][-1] / np.sum(equi_result['atom_numbs'])
                structure_dir = os.path.basename(ii)

                Cf = 1.60217657e-16 / 1e-20 * 0.001
                sfe = (task_result['energies'][-1] - task_result_slab_equi['energies'][-1]) / AA * Cf

                miller_index = loadfn(os.path.join(ii, 'miller.json'))
                ptr_data += "%-25s     %7.2f   %7.3f    %8.3f %8.3f\n" % (
                    str(miller_index) + '-' + structure_dir + ':', int(ii[-4:])/self.n_steps, sfe, epa, equi_epa_slab)
                res_data[int(ii[-4:])/self.n_steps] = [sfe, epa, equi_epa]


        else:
            if 'init_data_path' not in self.parameter:
                raise RuntimeError("please provide the initial data path to reproduce")
            init_data_path = os.path.abspath(self.parameter['init_data_path'])
            res_data, ptr_data = post_repro(init_data_path, self.parameter['init_from_suffix'],
                                            all_tasks, ptr_data, self.parameter.get('reprod_last_frame', True))

        with open(output_file, 'w') as fp:
            json.dump(res_data, fp, indent=4)

        return res_data, ptr_data
