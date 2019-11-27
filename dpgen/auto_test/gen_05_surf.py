#!/usr/bin/env python3

import os, re, argparse, filecmp, json, glob
import subprocess as sp
import numpy as np
import dpgen.auto_test.lib.vasp as vasp
import dpgen.auto_test.lib.lammps as lammps
from pymatgen.core.surface import generate_all_slabs, Structure

global_equi_name = '00.equi'
global_task_name = '05.surf'

def make_vasp(jdata, conf_dir, max_miller = 2, relax_box = False, static = False) :

    min_slab_size = jdata['min_slab_size']
    min_vacuum_size = jdata['min_vacuum_size']
    pert_xz = jdata['pert_xz']

    if 'relax_incar' in jdata.keys():
        vasp_str='vasp-relax_incar'
    else:
        kspacing = jdata['vasp_params']['kspacing']
        vasp_str='vasp-k%.2f' % (kspacing)

    # get conf poscar
    # conf_path = os.path.abspath(conf_dir)
    # conf_poscar = os.path.join(conf_path, 'POSCAR')
    equi_path = re.sub('confs', global_equi_name, conf_dir)
    equi_path = os.path.join(equi_path, vasp_str)
    equi_path = os.path.abspath(equi_path)
    equi_contcar = os.path.join(equi_path, 'CONTCAR')
    assert os.path.exists(equi_contcar),"Please compute the equilibrium state using vasp first"
    task_path = re.sub('confs', global_task_name, conf_dir)
    task_path = os.path.abspath(task_path)
    if static:
        if 'scf_incar' in jdata.keys():
            vasp_static_str='vasp-static-scf_incar'
        else:
            vasp_static_str='vasp-static-k%.2f' % (kspacing)
        task_path = os.path.join(task_path, vasp_static_str)
    else :
        task_path = os.path.join(task_path, vasp_str)
    os.makedirs(task_path, exist_ok=True)
    cwd = os.getcwd()
    os.chdir(task_path)
    if os.path.isfile('POSCAR') :
        os.remove('POSCAR')
    os.symlink(os.path.relpath(equi_contcar), 'POSCAR')
    os.chdir(cwd)
    task_poscar = os.path.join(task_path, 'POSCAR')
    ptypes = vasp.get_poscar_types(task_poscar)
    # gen strcture
    ss = Structure.from_file(task_poscar)
    # gen slabs
    all_slabs = generate_all_slabs(ss, max_miller, min_slab_size, min_vacuum_size)
    # gen incar
    if static :
        if  'scf_incar' in jdata.keys():
            scf_incar_path = jdata['scf_incar']
            assert(os.path.exists(scf_incar_path))
            scf_incar_path = os.path.abspath(scf_incar_path)
            fc = open(scf_incar_path).read()
        else :
            fp_params = jdata['vasp_params']
            ecut = fp_params['ecut']
            ediff = fp_params['ediff']
            npar = fp_params['npar']
            kpar = fp_params['kpar']
            kspacing = fp_params['kspacing']
            kgamma = fp_params['kgamma']
            fc = vasp.make_vasp_static_incar(ecut, ediff, npar=npar,kpar=kpar, kspacing = kspacing, kgamma = kgamma)
    else :
        if  'relax_incar' in jdata.keys():
            relax_incar_path = jdata['relax_incar']
            assert(os.path.exists(relax_incar_path))
            relax_incar_path = os.path.abspath(relax_incar_path)
            fc = open(relax_incar_path).read()
        else :
            fp_params = jdata['vasp_params']
            ecut = fp_params['ecut']
            ediff = fp_params['ediff']
            npar = fp_params['npar']
            kpar = fp_params['kpar']
            kspacing = fp_params['kspacing']
            kgamma = fp_params['kgamma']
            fc = vasp.make_vasp_relax_incar(ecut, ediff, True, relax_box, False, npar=npar,kpar=kpar, kspacing = kspacing, kgamma = kgamma)
    with open(os.path.join(task_path, 'INCAR'), 'w') as fp :
        fp.write(fc)
    # gen potcar
    with open(task_poscar,'r') as fp :
        lines = fp.read().split('\n')
        ele_list = lines[5].split()
    potcar_map = jdata['potcar_map']
    potcar_list = []
    for ii in ele_list :
        assert os.path.exists(os.path.abspath(potcar_map[ii])),"No POTCAR in the potcar_map of %s"%(ii)
        potcar_list.append(os.path.abspath(potcar_map[ii]))
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
        slab.to('POSCAR', 'POSCAR.tmp')
        vasp.regulate_poscar('POSCAR.tmp', 'POSCAR')
        vasp.sort_poscar('POSCAR', 'POSCAR', ptypes)
        vasp.perturb_xz('POSCAR', 'POSCAR', pert_xz)
        # record miller
        np.savetxt('miller.out', slab.miller_index, fmt='%d')
        # link incar, potcar, kpoints
        os.symlink(os.path.relpath(os.path.join(task_path, 'INCAR')), 'INCAR')
        os.symlink(os.path.relpath(os.path.join(task_path, 'POTCAR')), 'POTCAR')
    cwd = os.getcwd()

def make_lammps(jdata, conf_dir, max_miller = 2, static = False, relax_box = False, task_type = 'wrong-task') :
    #kspacing = jdata['vasp_params']['kspacing']
    fp_params = jdata['lammps_params']
    model_dir = fp_params['model_dir']
    type_map = fp_params['type_map']
    model_dir = os.path.abspath(model_dir)
    model_name =fp_params['model_name']
    if not model_name and task_type=='deepmd':
        models = glob.glob(os.path.join(model_dir, '*pb'))
        model_name = [os.path.basename(ii) for ii in models]
        assert len(model_name)>0,"No deepmd model in the model_dir"
    else:
        models = [os.path.join(model_dir,ii) for ii in model_name]

    model_param = {'model_name' :      fp_params['model_name'],
                  'param_type':          fp_params['model_param_type']}

    ntypes = len(type_map)

    min_slab_size = jdata['min_slab_size']
    min_vacuum_size = jdata['min_vacuum_size']

    # get equi poscar
    # conf_path = os.path.abspath(conf_dir)
    # conf_poscar = os.path.join(conf_path, 'POSCAR')
    if 'relax_incar' in jdata.keys():
        vasp_str='vasp-relax_incar'
    else:
        kspacing = jdata['vasp_params']['kspacing']
        vasp_str='vasp-k%.2f' % (kspacing)

    equi_path = re.sub('confs', global_equi_name, conf_dir)
    equi_path = os.path.join(equi_path, vasp_str)
    equi_path = os.path.abspath(equi_path)
    equi_contcar = os.path.join(equi_path, 'CONTCAR')
    assert os.path.exists(equi_contcar),"Please compute the equilibrium state using vasp first"
    task_path = re.sub('confs', global_task_name, conf_dir)
    task_path = os.path.abspath(task_path)
    if static:
        task_path = os.path.join(task_path, task_type+'-static')
    else:
        task_path = os.path.join(task_path, task_type)
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
    if task_type =='deepmd':
        if static :
            fc = lammps.make_lammps_eval('conf.lmp',
                                     ntypes,
                                     lammps.inter_deepmd,
                                     model_name)
        else :
            fc = lammps.make_lammps_equi('conf.lmp',
                                     ntypes,
                                     lammps.inter_deepmd,
                                     model_name,
                                     change_box = relax_box)
    elif task_type =='meam':
        if static :
            fc = lammps.make_lammps_eval('conf.lmp',
                                     ntypes,
                                     lammps.inter_meam,
                                     model_param)
        else :
            fc = lammps.make_lammps_equi('conf.lmp',
                                     ntypes,
                                     lammps.inter_meam,
                                     model_param,
                                     change_box = relax_box)
    f_lammps_in = os.path.join(task_path, 'lammps.in')
    with open(f_lammps_in, 'w') as fp :
        fp.write(fc)
    cwd = os.getcwd()

    os.chdir(task_path)
    for ii in model_name :
        if os.path.exists(ii) :
            os.remove(ii)
    for (ii,jj) in zip(models, model_name) :
        os.symlink(os.path.relpath(ii), jj)
    share_models = [os.path.join(task_path,ii) for ii in model_name]

    for ii in range(len(all_slabs)) :
        slab = all_slabs[ii]
        miller_str = "m%d.%d.%dm" % (slab.miller_index[0], slab.miller_index[1], slab.miller_index[2])
        # make dir
        struct_path = os.path.join(task_path, 'struct-%03d-%s' % (ii, miller_str))
        os.makedirs(struct_path, exist_ok=True)
        os.chdir(struct_path)
        for jj in ['conf.lmp', 'lammps.in'] + model_name :
            if os.path.isfile(jj):
                os.remove(jj)
        print("# %03d generate " % ii, struct_path, " \t %d atoms" % len(slab.sites))
        # make conf
        slab.to('POSCAR', 'POSCAR')
        vasp.regulate_poscar('POSCAR', 'POSCAR')
        lammps.cvt_lammps_conf('POSCAR', 'conf.lmp')
        ptypes = vasp.get_poscar_types('POSCAR')
        lammps.apply_type_map('conf.lmp', type_map, ptypes)
        # record miller
        np.savetxt('miller.out', slab.miller_index, fmt='%d')
        # link lammps.in
        os.symlink(os.path.relpath(f_lammps_in), 'lammps.in')
        # link models
        for (ii,jj) in zip(share_models, model_name) :
            os.symlink(os.path.relpath(ii), jj)
    cwd = os.getcwd()

def _main() :
    parser = argparse.ArgumentParser(
        description="gen 05.surf")
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
        make_vasp(jdata, args.CONF, args.MAX_MILLER, static = False, relax_box = args.relax_box)
    elif args.TASK == 'vasp-static':
        make_vasp(jdata, args.CONF, args.MAX_MILLER, static = True, relax_box = args.relax_box)
    elif args.TASK == 'deepmd' or args.TASK =='meam':
        make_lammps(jdata, args.CONF, args.MAX_MILLER, static = False, relax_box = args.relax_box, task_type = args.TASK)
    elif args.TASK == 'deepmd-static' or args.TASK == 'meam-static':
        args.TASK=args.TASK.replace('-static','')
        make_lammps(jdata, args.CONF, args.MAX_MILLER, static = True, relax_box = args.relax_box, task_type = args.TASK)
    else :
        raise RuntimeError("unknow task ", args.TASK)

if __name__ == '__main__' :
    _main()
