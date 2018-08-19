#!/usr/bin/env python3

import os, re, argparse, filecmp, json, glob
import subprocess as sp
import numpy as np
import lib.vasp as vasp
import lib.lammps as lammps
from pymatgen.core.surface import generate_all_slabs, Structure

global_equi_name = '00.equi'
global_task_name = '05.surf'

def make_vasp(jdata, conf_dir, max_miller = 2, relax_box = False) :
    fp_params = jdata['vasp_params']
    ecut = fp_params['ecut']
    ediff = fp_params['ediff']
    npar = fp_params['npar']
    kpar = fp_params['kpar']
    kspacing = fp_params['kspacing']
    kgamma = fp_params['kgamma']
    min_slab_size = jdata['min_slab_size']
    min_vacuum_size = jdata['min_vacuum_size']
    pert_xz = jdata['pert_xz']

    # get conf poscar
    # conf_path = os.path.abspath(conf_dir)
    # conf_poscar = os.path.join(conf_path, 'POSCAR')
    equi_path = re.sub('confs', global_equi_name, conf_dir)
    equi_path = os.path.join(equi_path, 'vasp-k%.2f' % kspacing)
    equi_path = os.path.abspath(equi_path)    
    equi_contcar = os.path.join(equi_path, 'CONTCAR')
    task_path = re.sub('confs', global_task_name, conf_dir)
    task_path = os.path.abspath(task_path)
    task_path = os.path.join(task_path, 'vasp-k%.2f' % kspacing)
    os.makedirs(task_path, exist_ok=True)
    cwd = os.getcwd()
    os.chdir(task_path)
    if os.path.isfile('POSCAR') :
        os.remove('POSCAR')
    os.symlink(os.path.relpath(equi_contcar), 'POSCAR')
    os.chdir(cwd)
    task_poscar = os.path.join(task_path, 'POSCAR')
    # gen strcture
    ss = Structure.from_file(task_poscar)
    # gen slabs
    all_slabs = generate_all_slabs(ss, max_miller, min_slab_size, min_vacuum_size)
    # gen incar
    fc = vasp.make_vasp_relax_incar(ecut, ediff, True, relax_box, False, 1, 1, kspacing = kspacing, kgamma = kgamma)
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
    # gen tasks    
    cwd = os.getcwd()
    for ii in range(len(all_slabs)) :
        slab = all_slabs[ii]
        miller_str = "m%d.%d.%dm" % (slab.miller_index[0], slab.miller_index[1], slab.miller_index[2])
        # make dir
        struct_path = os.path.join(task_path, 'struct-%03d-%s' % (ii, miller_str))
        os.makedirs(struct_path, exist_ok=True)
        os.chdir(struct_path)
        for jj in ['POSCAR', 'POTCAR', 'INCAR'] :
            if os.path.isfile(jj):
                os.remove(jj)
        print("# %03d generate " % ii, struct_path, " \t %d atoms" % len(slab.sites))
        # make conf
        slab.to('POSCAR', 'POSCAR')
        vasp.regulate_poscar('POSCAR', 'POSCAR')
        vasp.perturb_xz('POSCAR', 'POSCAR', pert_xz)
        # record miller
        np.savetxt('miller.out', slab.miller_index, fmt='%d')
        # link incar, potcar, kpoints
        os.symlink(os.path.relpath(os.path.join(task_path, 'INCAR')), 'INCAR')
        os.symlink(os.path.relpath(os.path.join(task_path, 'POTCAR')), 'POTCAR')
    cwd = os.getcwd()

def make_deepmd_lammps(jdata, conf_dir, max_miller = 2, relax_box = False) :
    fp_params = jdata['vasp_params']
    kspacing = fp_params['kspacing']
    deepmd_model_dir = jdata['deepmd_model_dir']
    deepmd_type_map = jdata['deepmd_type_map']
    ntypes = len(deepmd_type_map)    
    deepmd_model_dir = os.path.abspath(deepmd_model_dir)
    deepmd_models = glob.glob(os.path.join(deepmd_model_dir, '*pb'))
    deepmd_models_name = [os.path.basename(ii) for ii in deepmd_models]
    min_slab_size = jdata['min_slab_size']
    min_vacuum_size = jdata['min_vacuum_size']

    # get equi poscar
    # conf_path = os.path.abspath(conf_dir)
    # conf_poscar = os.path.join(conf_path, 'POSCAR')
    equi_path = re.sub('confs', global_equi_name, conf_dir)
    equi_path = os.path.join(equi_path, 'vasp-k%.2f' % kspacing)
    equi_path = os.path.abspath(equi_path)    
    equi_contcar = os.path.join(equi_path, 'CONTCAR')    
    task_path = re.sub('confs', global_task_name, conf_dir)
    task_path = os.path.abspath(task_path)
    task_path = os.path.join(task_path, 'lmp')
    os.makedirs(task_path, exist_ok=True)
    cwd = os.getcwd()
    os.chdir(task_path)
    if os.path.isfile('POSCAR') :
        os.remove('POSCAR')
    os.symlink(os.path.relpath(equi_contcar), 'POSCAR')
    os.chdir(cwd)
    task_poscar = os.path.join(task_path, 'POSCAR')
    # gen strcture
    ss = Structure.from_file(task_poscar)
    # gen slabs
    all_slabs = generate_all_slabs(ss, max_miller, min_slab_size, min_vacuum_size)
    # make lammps.in
    fc = lammps.make_lammps_equi('conf.lmp', ntypes, deepmd_models_name, change_box = relax_box)
    f_lammps_in = os.path.join(task_path, 'lammps.in')
    with open(f_lammps_in, 'w') as fp :
        fp.write(fc)
    cwd = os.getcwd()
    for ii in range(len(all_slabs)) :
        slab = all_slabs[ii]
        miller_str = "m%d.%d.%dm" % (slab.miller_index[0], slab.miller_index[1], slab.miller_index[2])
        # make dir
        struct_path = os.path.join(task_path, 'struct-%03d-%s' % (ii, miller_str))
        os.makedirs(struct_path, exist_ok=True)
        os.chdir(struct_path)
        for jj in ['conf.lmp', 'lammps.in'] + deepmd_models_name :
            if os.path.isfile(jj):
                os.remove(jj)
        print("# %03d generate " % ii, struct_path, " \t %d atoms" % len(slab.sites))
        # make conf
        slab.to('POSCAR', 'POSCAR')
        vasp.regulate_poscar('POSCAR', 'POSCAR')
        lammps.cvt_lammps_conf('POSCAR', 'conf.lmp')
        ptypes = vasp.get_poscar_types('POSCAR')
        lammps.apply_type_map('conf.lmp', deepmd_type_map, ptypes)    
        # record miller
        np.savetxt('miller.out', slab.miller_index, fmt='%d')
        # link lammps.in
        os.symlink(os.path.relpath(f_lammps_in), 'lammps.in')
        # link models
        for (ii,jj) in zip(deepmd_models, deepmd_models_name) :
            os.symlink(os.path.relpath(ii), jj)
    cwd = os.getcwd()
    
def _main() :
    parser = argparse.ArgumentParser(
        description="gen 01.eos")
    parser.add_argument('TASK', type=str,
                        help='the task of generation, vasp or lammps')
    parser.add_argument('PARAM', type=str,
                        help='json parameter file')
    parser.add_argument('CONF', type=str,
                        help='the path to conf')
    parser.add_argument('MAX_MILLER', type=int,
                        help='the maximum miller index')
    parser.add_argument('-r', '--relax-box', action = 'store_true',
                        help='set if the box is relaxed, otherwise only relax atom positions')
    args = parser.parse_args()

    with open (args.PARAM, 'r') as fp :
        jdata = json.load (fp)

    print('# generate %s task with conf %s' % (args.TASK, args.CONF))
    if args.TASK == 'vasp':
        make_vasp(jdata, args.CONF, args.MAX_MILLER, args.relax_box)
    elif args.TASK == 'lammps' :
        make_deepmd_lammps(jdata, args.CONF, args.MAX_MILLER, args.relax_box)
    else :
        raise RuntimeError("unknow task ", args.TASK)
    
if __name__ == '__main__' :
    _main()

    
