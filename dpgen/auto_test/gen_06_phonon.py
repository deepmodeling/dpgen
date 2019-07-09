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
#from phonolammps import Phonolammps


global_equi_name = '00.equi'
global_task_name = '06.phonon'

'''
link poscar
link potcar
make incar
'''
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
    fc = vasp.make_vasp_phonon_incar(ecut, ediff, 1, 1, kspacing = None, kgamma = None)
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
    '''
    phonon = phonopy.load(supercell_matrix=supercell_matrix,
                      primitive_matrix='auto',
                      unitcell_filename='POSCAR-unitcell')
    phonon.save(filename='phonopy_disp.yaml')       
    with open('phonopy_disp.yaml', 'r') as f:
        temp = yaml.load(f.read())
    with open('POSCAR', 'w') as fp:
        for ii in ele_list:
            fp.write(ii)
            fp.write(' ')
        fp.write('\n')
        data=open('POSCAR-unitcell', 'r')
        next(data)
        fp.write(data.readline())
        for ii in temp['supercell']['lattice']:
            fp.write(str(ii).replace(',', '').replace('[', '').replace(']','\n'))
        for ii in ele_list:
            fp.write(str(str(temp['supercell']['points']).count(ii)))
            fp.write(' ')
        fp.write('\n')
        fp.write('Direct\n')
        for ii in temp['supercell']['points']:
            fp.write(str(ii['coordinates']).replace(',', '').replace('[', '').replace(']', '\n'))
    '''
    os.system('phonopy -d --dim="%d %d %d" -c POSCAR-unitcell'%(supercell_matrix[0],supercell_matrix[1],supercell_matrix[2]))
    os.system('cp SPOSCAR POSCAR')

def make_deepmd_lammps(jdata, conf_dir) :
    deepmd_model_dir = jdata['deepmd_model_dir']
    deepmd_type_map = jdata['deepmd_type_map']
    ntypes = len(deepmd_type_map)    
    deepmd_model_dir = os.path.abspath(deepmd_model_dir)
    deepmd_models = glob.glob(os.path.join(deepmd_model_dir, '*pb'))
    deepmd_models_name = [os.path.basename(ii) for ii in deepmd_models]
    supercell_matrix=jdata['supercell_matrix']
    band_path=jdata['band']

    conf_path = os.path.abspath(conf_dir)
    conf_poscar = os.path.join(conf_path, 'POSCAR')
    # get equi poscar
    equi_path = re.sub('confs', global_equi_name, conf_path)
    equi_path = os.path.join(equi_path, 'deepmd')
    task_path = re.sub('confs', global_task_name, conf_path)
    task_path = os.path.join(task_path, 'deepmd')
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
    lammps.apply_type_map(conf_file, deepmd_type_map, ptypes) 
    # make lammps.in
    ntypes=len(ele_list)
    unitcell=PhonopyAtoms(symbols=ele_list,cell=(np.eye(3)),scaled_positions=np.zeros((ntypes,3)))
    fc = lammps.make_lammps_phonon('conf.lmp', 
                                    unitcell.masses, 
                                    lammps.inter_deepmd,
                                    deepmd_models_name)        
    f_lammps_in = os.path.join(task_path, 'lammps.in')
    with open(f_lammps_in, 'w') as fp :
        fp.write(fc)
    cwd = os.getcwd()
    # link models
    os.chdir(task_path)
    for ii in deepmd_models_name :
        if os.path.exists(ii) :
            os.remove(ii)
    for (ii,jj) in zip(deepmd_models, deepmd_models_name) :
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
        fp.write('FORCE_CONSTANTS=READ')
    # gen task
    ''' 
    phlammps = Phonolammps('lammps.in',supercell_matrix=supercell_matrix)
    unitcell = phlammps.get_unitcell()
    phonon = Phonopy(unitcell,supercell_matrix)
    phonon.save(filename='phonopy_disp.yaml')       
    with open('phonopy_disp.yaml', 'r') as f:
        temp = yaml.load(f.read())
    with open('POSCAR-unitcell', 'w') as fp:
        for ii in ele_list:
            fp.write(ii)
            fp.write(' ')
        fp.write('\n')
        data=open('POSCAR', 'r')
        next(data)
        fp.write(data.readline())
        for ii in temp['unit_cell']['lattice']:
            fp.write(str(ii).replace(',', '').replace('[', '').replace(']','\n'))
        for ii in ele_list:
            fp.write(str(str(temp['unit_cell']['points']).count(ii)))
            fp.write(' ')
        fp.write('\n')
        fp.write('Direct\n')
        for ii in temp['unit_cell']['points']:
            fp.write(str(ii['coordinates']).replace(',', '').replace('[', '').replace(']', '\n'))
    '''
    os.system('phonolammps lammps.in --dim %d %d %d -c POSCAR-unitcell'%(supercell_matrix[0],supercell_matrix[1],supercell_matrix[2]))


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
    elif args.TASK == 'deepmd' :
        make_deepmd_lammps(jdata, args.CONF)
    elif args.TASK == 'meam' :
        make_meam_lammps(jdata, args.CONF)
    else :
        raise RuntimeError("unknow task ", args.TASK)
    
if __name__ == '__main__' :
    _main()

