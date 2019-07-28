#!/usr/bin/env python3

import os, re, argparse, filecmp, json, glob
import subprocess as sp
import numpy as np
import dpgen.auto_test.lib.vasp as vasp
import dpgen.auto_test.lib.lammps as lammps
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms
import yaml
import phonopy


global_equi_name = '00.equi'
global_task_name = '06.phonon'

'''
link poscar
link potcar
make incar
''' 

def get_structure_from_poscar(file_name, number_of_dimensions=3):
    """
    Read crystal structure from a VASP POSCAR type file
    :param file_name: POSCAR filename
    :param number_of_dimensions: number of dimensions of the crystal structure
    :return: Atoms (phonopy) type object containing the crystal structure
    """
    # Check file exists
    if not os.path.isfile(file_name):
        print('Structure file does not exist!')
        exit()

    # Read from VASP POSCAR file
    poscar_file = open(file_name, 'r')
    data_lines = poscar_file.read().split('\n')
    poscar_file.close()

    multiply = float(data_lines[1])
    direct_cell = np.array([data_lines[i].split()
                            for i in range(2, 2+number_of_dimensions)], dtype=float)
    direct_cell *= multiply
    scaled_positions = None
    positions = None

    try:
        number_of_types = np.array(data_lines[3+number_of_dimensions].split(),dtype=int)

        coordinates_type = data_lines[4+number_of_dimensions][0]
        if coordinates_type == 'D' or coordinates_type == 'd' :

            scaled_positions = np.array([data_lines[8+k].split()[0:3]
                                         for k in range(np.sum(number_of_types))],dtype=float)
        else:
            positions = np.array([data_lines[8+k].split()[0:3]
                                  for k in range(np.sum(number_of_types))],dtype=float)

        atomic_types = []
        for i,j in enumerate(data_lines[5].split()):
            atomic_types.append([j]*number_of_types[i])
        atomic_types = [item for sublist in atomic_types for item in sublist]

    # Old style POSCAR format
    except ValueError:
        number_of_types = np.array(data_lines[5].split(), dtype=int)
        coordinates_type = data_lines[6][0]
        if coordinates_type == 'D' or coordinates_type == 'd':
            scaled_positions = np.array([data_lines[7+k].split()[0:3]
                                         for k in range(np.sum(number_of_types))], dtype=float)
        else:
            positions = np.array([data_lines[7+k].split()[0:3]
                                  for k in range(np.sum(number_of_types))], dtype=float)

        atomic_types = []
        for i,j in enumerate(data_lines[0].split()):
            atomic_types.append([j]*number_of_types[i])
        atomic_types = [item for sublist in atomic_types for item in sublist]

    return PhonopyAtoms(symbols=atomic_types,
                        scaled_positions=scaled_positions,
                        cell=direct_cell)

def make_vasp(jdata, conf_dir) :
    fp_params = jdata['vasp_params']
    ecut = fp_params['ecut']
    ediff = fp_params['ediff']
    npar = fp_params['npar']
    kpar = fp_params['kpar']
    kspacing = fp_params['kspacing']
    kgamma = fp_params['kgamma']
    supercell_matrix=jdata['supercell_matrix']
    band_path=jdata['band']

    conf_path = os.path.abspath(conf_dir)
    conf_poscar = os.path.join(conf_path, 'POSCAR')
    equi_path = re.sub('confs', global_equi_name, conf_path)
    equi_path = os.path.join(equi_path, 'vasp-k%.2f' % kspacing)
    equi_contcar = os.path.join(equi_path, 'CONTCAR')
    task_path = re.sub('confs', global_task_name, conf_path)
    task_path = os.path.join(task_path, 'vasp-k%.2f' % kspacing)
    os.makedirs(task_path, exist_ok=True)
    cwd = os.getcwd()
    os.chdir(task_path)
    print(task_path)
    if os.path.isfile('POSCAR-unitcell') :
        os.remove('POSCAR-unitcell')
    if os.path.isfile('POSCAR') :
        os.remove('POSCAR')
    os.symlink(os.path.relpath(equi_contcar), 'POSCAR-unitcell')
    os.chdir(cwd)
    task_poscar = os.path.join(task_path, 'POSCAR-unitcell')   
    # gen incar
    fc = vasp.make_vasp_phonon_incar(ecut, ediff, npar, kpar, kspacing = None, kgamma = None)
    with open(os.path.join(task_path, 'INCAR'), 'w') as fp :
        fp.write(fc)
    # gen potcar
    with open(task_poscar,'r') as fp :
        lines = fp.read().split('\n')
        ele_list = lines[5].split()
    potcar_map = jdata['potcar_map']
    potcar_list = []
    for ii in ele_list :
        assert(os.path.exists(potcar_map[ii]))
        potcar_list.append(potcar_map[ii])
    with open(os.path.join(task_path,'POTCAR'), 'w') as outfile:
        for fname in potcar_list:
            with open(fname) as infile:
                outfile.write(infile.read())
    # gen kpoints
    fc = vasp.make_kspacing_kpoints(task_poscar, kspacing, kgamma)
    with open(os.path.join(task_path,'KPOINTS'), 'w') as fp:
        fp.write(fc)
    # gen band.conf
    os.chdir(task_path)
    with open('band.conf','w') as fp:
        fp.write('ATOM_NAME = ')
        for ii in ele_list:
            fp.write(ii)
            fp.write(' ')
        fp.write('\n')
        fp.write('DIM = %d %d %d\n'%(supercell_matrix[0],supercell_matrix[1],supercell_matrix[2]))
        fp.write('BAND = %s\n'%band_path)
        fp.write('FORCE_CONSTANTS=READ')
    # gen POSCAR
    os.system('phonopy -d --dim="%d %d %d" -c POSCAR-unitcell'%(supercell_matrix[0],supercell_matrix[1],supercell_matrix[2]))
    os.symlink('SPOSCAR', 'POSCAR')

def make_lammps(jdata, conf_dir,task_type) :
    fp_params = jdata['lammps_params']
    model_dir = fp_params['model_dir']
    type_map = fp_params['type_map'] 
    model_dir = os.path.abspath(model_dir)
    model_name =fp_params['model_name']
    if not model_name :
        models = glob.glob(os.path.join(model_dir, '*pb'))
        model_name = [os.path.basename(ii) for ii in models]
    else:
        models = [os.path.join(model_dir,ii) for ii in model_name]

    model_param = {'model_name' :      fp_params['model_name'],
                  'param_type':          fp_params['model_param_type']}

    supercell_matrix=jdata['supercell_matrix']
    band_path=jdata['band']

    conf_path = os.path.abspath(conf_dir)
    conf_poscar = os.path.join(conf_path, 'POSCAR')
    # get equi poscar
    equi_path = re.sub('confs', global_equi_name, conf_path)
    equi_path = os.path.join(equi_path, task_type)
    task_path = re.sub('confs', global_task_name, conf_path)
    task_path = os.path.join(task_path, task_type)
    os.makedirs(task_path, exist_ok=True)
    
    task_poscar = os.path.join(task_path, 'POSCAR')
    cwd = os.getcwd()
    os.chdir(task_path)
    if os.path.isfile('POSCAR') :
        os.remove('POSCAR')
    os.symlink(os.path.relpath(conf_poscar), 'POSCAR')
    os.chdir(cwd)    
    with open(task_poscar,'r') as fp :
        lines = fp.read().split('\n')
        ele_list = lines[5].split()
    
    print(task_path)
    # make conf.lmp  
    conf_file = os.path.join(task_path, 'conf.lmp')
    lammps.cvt_lammps_conf(task_poscar, os.path.relpath(conf_file))
    ptypes = vasp.get_poscar_types(task_poscar)
    lammps.apply_type_map(conf_file, type_map, ptypes) 
    # make lammps.in
    ntypes=len(ele_list)
    unitcell=get_structure_from_poscar(task_poscar)
    if task_type=='deepmd':
        fc = lammps.make_lammps_phonon('conf.lmp', 
                                    unitcell.masses, 
                                    lammps.inter_deepmd,
                                    model_name)
    if task_type=='meam':
        fc = lammps.make_lammps_phonon('conf.lmp', 
                                    unitcell.masses, 
                                    lammps.inter_meam,
                                    model_param)      
    f_lammps_in = os.path.join(task_path, 'lammps.in')
    with open(f_lammps_in, 'w') as fp :
        fp.write(fc)
    cwd = os.getcwd()
    # link models
    os.chdir(task_path)
    for ii in model_name :
        if os.path.exists(ii) :
            os.remove(ii)
    for (ii,jj) in zip(models, model_name) :
        os.symlink(os.path.relpath(ii), jj)
    # gen band.conf
    os.chdir(task_path)
    with open('band.conf','w') as fp:
        fp.write('ATOM_NAME = ')
        for ii in ele_list:
            fp.write(ii)
            fp.write(' ')
        fp.write('\n')
        fp.write('DIM = %d %d %d\n'%(supercell_matrix[0],supercell_matrix[1],supercell_matrix[2]))
        fp.write('BAND = %s\n'%band_path)
        fp.write('FORCE_CONSTANTS=READ\n')
    os.system('phonolammps lammps.in --dim %d %d %d -c POSCAR'%(supercell_matrix[0],supercell_matrix[1],supercell_matrix[2]))


def _main() :
    parser = argparse.ArgumentParser(
        description="gen 06.phonon")
    parser.add_argument('TASK', type=str,
                        help='the task of generation, vasp or lammps')
    parser.add_argument('PARAM', type=str,
                        help='json parameter file')
    parser.add_argument('CONF', type=str,
                        help='the path to conf')
    args = parser.parse_args()

    with open (args.PARAM, 'r') as fp :
        jdata = json.load (fp)

#    print('generate %s task with conf %s' % (args.TASK, args.CONF))
    if args.TASK == 'vasp':
        make_vasp(jdata, args.CONF)               
    elif args.TASK == 'deepmd' or args.TASK =='meam' :
        make_lammps(jdata, args.CONF,args.TASK)
    else :
        raise RuntimeError("unknow task ", args.TASK)
    
if __name__ == '__main__' :
    _main()

