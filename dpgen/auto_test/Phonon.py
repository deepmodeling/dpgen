import os , argparse , json, glob, shutil
import subprocess as sp
import numpy as np
import dpgen.auto_test.lib.vasp as vasp
import dpgen.auto_test.lib.lammps as lammps

from dpgen import ROOT_PATH
from dpgen.auto_test.Property import Property
from dpgen import dlog
from pymatgen.io.vasp import Incar
from dpgen.generator.lib.vasp import incar_upper

class Phonon(Property):
    def __init__(self,parameter,inter_param = None):
        parameter['reproduce'] = parameter.get('reproduce', False)
        parameter['primitive'] = parameter.get('primitive', False)
        self.reprod = parameter['reproduce']
        if not self.reprod:
            if not ('init_from_suffix' in parameter and 'output_suffix' in parameter):
                self.band_path = parameter['band_path']
                self.supercell_matrix = parameter['supercell_matrix']
                self.primitive = parameter['primitive']
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

    def make_confs(self,path_to_work,path_to_equi,refine=False):
        supercell_matrix = self.supercell_matrix
        band_path = self.band_path 
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
        equi_contcar = os.path.join(path_to_equi, 'CONTCAR')

        os.chdir(path_to_work)
        if os.path.isfile('POSCAR'):
            os.remove('POSCAR')
        if os.path.islink('POSCAR'):
            os.remove('POSCAR')
        os.symlink(os.path.relpath(equi_contcar), 'POSCAR')
        if self.primitive:
            os.system('phonopy --symmetry')
            os.system('cp PPOSCAR POSCAR')
        os.chdir(cwd)

        if refine:
            raise RuntimeError('Refine is not supported for phonon calculation')
        else:
            if not os.path.exists(equi_contcar):
                raise RuntimeError("please do relaxation first")
            n_task = 1

            for ii in range(n_task):
                output_task = os.path.join(path_to_work,'task.%06d' % ii)
                os.makedirs(output_task,exist_ok=True)
                os.chdir(output_task)
                for jj in ['INCAR', 'POTCAR', 'POSCAR', 'conf.lmp', 'in.lammps','POSCAR-unitcell','SPOSCAR']:
                    if os.path.exists(jj):
                        os.remove(jj)
                task_list.append(output_task)
                os.symlink(os.path.relpath(equi_contcar), 'POSCAR-unitcell')
                
                #gen band.conf
                with open("POSCAR-unitcell",'r') as fp :
                    lines = fp.read().split('\n')
                    ele_list = lines[5].split()
                with open('band.conf','w') as fp:
                    fp.write('ATOM_NAME = ')
                    for ii in ele_list:
                        fp.write(ii)
                        fp.write(' ')
                    fp.write('\n')
                    fp.write('DIM = %s %s %s\n'%(supercell_matrix[0],supercell_matrix[1],supercell_matrix[2]))
                    fp.write('BAND = %s\n'%band_path)
                    fp.write('FORCE_CONSTANTS=READ\n')
            
                #gen POSCAR
                if (self.inter_param['type'] == "vasp"):
                    os.system('phonopy -d --dim="%d %d %d" -c POSCAR-unitcell'%(int(supercell_matrix[0]),int(supercell_matrix[1]),int(supercell_matrix[2])))
                    os.symlink('SPOSCAR','POSCAR')
                else:
                    os.system('cp POSCAR-unitcell POSCAR')
                
            os.chdir(cwd)
        return task_list

    def post_process(self,task_list):
        pass

    def task_type(self):
        return self.parameter['type']

    def task_param(self):
        return self.parameter

    def _compute_lower(self,
                       output_file,
                       all_tasks,
                       all_res):
        cwd = os.getcwd()
        supercell_matrix = self.supercell_matrix
        output_file = os.path.abspath(output_file)
        path_to_work = os.path.dirname(output_file)

        res_data = {}
        res_data['supercell_matrix'] = supercell_matrix
        ptr_data = "conf_dir: " + os.path.dirname(output_file) + "\n"
        if not self.reprod:
            for ii in range(len(all_tasks)):
                os.chdir(all_tasks[ii])
            
            if (self.inter_param['type'] == 'vasp'):
                if os.path.isfile('vasprun.xml'):
                    os.system('phonopy --fc vasprun.xml')
                    if os.path.isfile('FORCE_CONSTANTS'):
                        os.system('phonopy --dim="%s %s %s" -c POSCAR-unitcell band.conf'%(supercell_matrix[0],supercell_matrix[1],supercell_matrix[2]))
                        os.system('phonopy-bandplot --gnuplot band.yaml > band.dat')
                        print('band.dat is created')
                    else:
                        print('FORCE_CONSTANTS No such file')
                else:
                    print('vasprun.xml No such file')
            else:
                if os.path.isfile('FORCE_CONSTANTS'):
                    os.system('phonopy --dim="%s %s %s" -c POSCAR band.conf'%(supercell_matrix[0],supercell_matrix[1],supercell_matrix[2]))
                    os.system('phonopy-bandplot --gnuplot band.yaml > band.dat')
                else:
                    print('FORCE_CONSTANTS No such file')
            with open ("band.dat","r") as f:
                ptr_data += '%s'%(f.read())
        else:
            if 'init_data_path' not in self.parameter:
                raise RuntimeError("please provide the initial data path to reproduce")
            init_data_path = os.path.abspath(self.parameter['init_data_path'])
            res_data, ptr_data = post_repro(init_data_path, self.parameter['init_from_suffix'],
                                            all_tasks, ptr_data, self.parameter.get('reprod_last_frame', True))

        with open(output_file, 'w') as fp:
            json.dump(res_data, fp, indent=4)
        
        os.chdir(cwd)
        return res_data, ptr_data