#!/usr/bin/env python3

import os, re, argparse, filecmp, json, glob
import subprocess as sp
import numpy as np
import dpgen.auto_test.lib.vasp as vasp
import dpgen.auto_test.lib.lammps as lammps
from pymatgen.core.structure import Structure
from pymatgen.analysis.defects.core import Vacancy
from pymatgen.analysis.defects.generators import VacancyGenerator

global_equi_name = '00.equi'
global_task_name = '03.vacancy'

def make_vasp(jdata, conf_dir, supercell = [1,1,1]) :
    fp_params = jdata['vasp_params']
    ecut = fp_params['ecut']
    ediff = fp_params['ediff']
    npar = fp_params['npar']
    kpar = fp_params['kpar']
    kspacing = fp_params['kspacing']
    kgamma = fp_params['kgamma']
    conf_path = os.path.abspath(conf_dir)
    conf_poscar = os.path.join(conf_path, 'POSCAR')
    # get equi poscar
    equi_path = re.sub('confs', global_equi_name, conf_path)
    equi_path = os.path.join(equi_path, 'vasp-k%.2f' % kspacing)
    equi_contcar = os.path.join(equi_path, 'CONTCAR')
    task_path = re.sub('confs', global_task_name, conf_path)
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
    # gen defects
    vds = VacancyGenerator(ss)
    dss = []
    for jj in vds :
        dss.append(jj.generate_defect_structure(supercell))
    # gen incar
    fc = vasp.make_vasp_relax_incar(ecut, ediff, True, True, True, 1, 1, kspacing = kspacing, kgamma = kgamma)
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
    copy_str = "%sx%sx%s" % (supercell[0], supercell[1], supercell[2])
    cwd = os.getcwd()
    for ii in range(len(dss)) :
        struct_path = os.path.join(task_path, 'struct-%s-%03d' % (copy_str,ii))
        print('# generate %s' % (struct_path))
        os.makedirs(struct_path, exist_ok=True)
        os.chdir(struct_path)
        for jj in ['POSCAR', 'POTCAR', 'INCAR'] :
            if os.path.isfile(jj):
                os.remove(jj)
        # make conf
        dss[ii].to('POSCAR', 'POSCAR')
        # link incar, potcar, kpoints
        os.symlink(os.path.relpath(os.path.join(task_path, 'INCAR')), 'INCAR')
        os.symlink(os.path.relpath(os.path.join(task_path, 'POTCAR')), 'POTCAR')
        # save supercell
        np.savetxt('supercell.out', supercell, fmt='%d')
    os.chdir(cwd)

def make_deepmd_lammps(jdata, conf_dir, supercell) :
    fp_params = jdata['vasp_params']
    kspacing = fp_params['kspacing']
    deepmd_model_dir = jdata['deepmd_model_dir']
    deepmd_type_map = jdata['deepmd_type_map']
    ntypes = len(deepmd_type_map)    
    deepmd_model_dir = os.path.abspath(deepmd_model_dir)
    deepmd_models = glob.glob(os.path.join(deepmd_model_dir, '*pb'))
    deepmd_models_name = [os.path.basename(ii) for ii in deepmd_models]

    conf_path = os.path.abspath(conf_dir)
    conf_poscar = os.path.join(conf_path, 'POSCAR')
    # get equi poscar
    equi_path = re.sub('confs', global_equi_name, conf_path)
    equi_path = os.path.join(equi_path, 'vasp-k%.2f' % kspacing)
    equi_contcar = os.path.join(equi_path, 'CONTCAR')
    # equi_path = re.sub('confs', global_equi_name, conf_path)
    # equi_path = os.path.join(equi_path, 'lmp')
    # equi_dump = os.path.join(equi_path, 'dump.relax')
    task_path = re.sub('confs', global_task_name, conf_path)
    task_path = os.path.join(task_path, 'deepmd')
    os.makedirs(task_path, exist_ok=True)
    # gen task poscar
    task_poscar = os.path.join(task_path, 'POSCAR')
    # lammps.poscar_from_last_dump(equi_dump, task_poscar, deepmd_type_map)
    cwd = os.getcwd()
    os.chdir(task_path)
    if os.path.isfile('POSCAR') :
        os.remove('POSCAR')
    os.symlink(os.path.relpath(equi_contcar), 'POSCAR')
    os.chdir(cwd)
    # gen structure from equi poscar
    ss = Structure.from_file(task_poscar)
    # gen defects
    vds = VacancyGenerator(ss)
    dss = []
    for jj in vds :
        dss.append(jj.generate_defect_structure(supercell))
    # gen tasks    
    cwd = os.getcwd()
    # make lammps.in, relax at 0 bar (scale = 1)
    fc = lammps.make_lammps_press_relax('conf.lmp', 
                                        ntypes, 
                                        1, 
                                        lammps.inter_deepmd,
                                        deepmd_models_name)
    f_lammps_in = os.path.join(task_path, 'lammps.in')
    with open(f_lammps_in, 'w') as fp :
        fp.write(fc)
    # gen tasks    
    copy_str = "%sx%sx%s" % (supercell[0], supercell[1], supercell[2])
    cwd = os.getcwd()
    for ii in range(len(dss)) :
        struct_path = os.path.join(task_path, 'struct-%s-%03d' % (copy_str,ii))
        print('# generate %s' % (struct_path))
        os.makedirs(struct_path, exist_ok=True)
        os.chdir(struct_path)
        for jj in ['conf.lmp', 'lammps.in'] + deepmd_models_name :
            if os.path.isfile(jj):
                os.remove(jj)
        # make conf
        dss[ii].to('POSCAR', 'POSCAR')
        lammps.cvt_lammps_conf('POSCAR', 'conf.lmp')
        ptypes = vasp.get_poscar_types('POSCAR')
        lammps.apply_type_map('conf.lmp', deepmd_type_map, ptypes)    
        # link lammps.in
        os.symlink(os.path.relpath(f_lammps_in), 'lammps.in')
        # link models
        for (ii,jj) in zip(deepmd_models, deepmd_models_name) :
            os.symlink(os.path.relpath(ii), jj)
        # save supercell
        np.savetxt('supercell.out', supercell, fmt='%d')
    os.chdir(cwd)

def make_meam_lammps(jdata, conf_dir, supercell) :
    fp_params = jdata['vasp_params']
    kspacing = fp_params['kspacing']

    meam_potfile_dir = jdata['meam_potfile_dir']
    meam_potfile_dir = os.path.abspath(meam_potfile_dir)
    meam_potfile = jdata['meam_potfile']
    meam_potfile = [os.path.join(meam_potfile_dir,ii) for ii in meam_potfile]
    meam_potfile_name = jdata['meam_potfile']
    type_map = jdata['meam_type_map']
    ntypes = len(type_map)
    meam_param = {'meam_potfile' :      jdata['meam_potfile'],
                  'meam_type':          jdata['meam_param_type']}

    conf_path = os.path.abspath(conf_dir)
    conf_poscar = os.path.join(conf_path, 'POSCAR')
    # get equi poscar
    equi_path = re.sub('confs', global_equi_name, conf_path)
    equi_path = os.path.join(equi_path, 'vasp-k%.2f' % kspacing)
    equi_contcar = os.path.join(equi_path, 'CONTCAR')
    # equi_path = re.sub('confs', global_equi_name, conf_path)
    # equi_path = os.path.join(equi_path, 'lmp')
    # equi_dump = os.path.join(equi_path, 'dump.relax')
    task_path = re.sub('confs', global_task_name, conf_path)
    task_path = os.path.join(task_path, 'meam')
    os.makedirs(task_path, exist_ok=True)
    # gen task poscar
    task_poscar = os.path.join(task_path, 'POSCAR')
    # lammps.poscar_from_last_dump(equi_dump, task_poscar, deepmd_type_map)
    cwd = os.getcwd()
    os.chdir(task_path)
    if os.path.isfile('POSCAR') :
        os.remove('POSCAR')
    os.symlink(os.path.relpath(equi_contcar), 'POSCAR')
    os.chdir(cwd)
    # gen structure from equi poscar
    ss = Structure.from_file(task_poscar)
    # gen defects
    vds = VacancyGenerator(ss)
    dss = []
    for jj in vds :
        dss.append(jj.generate_defect_structure(supercell))
    # gen tasks    
    cwd = os.getcwd()
    # make lammps.in, relax at 0 bar (scale = 1)
    fc = lammps.make_lammps_press_relax('conf.lmp', 
                                        ntypes, 
                                        1, 
                                        lammps.inter_meam,
                                        meam_param)
    f_lammps_in = os.path.join(task_path, 'lammps.in')
    with open(f_lammps_in, 'w') as fp :
        fp.write(fc)
    # gen tasks    
    copy_str = "%sx%sx%s" % (supercell[0], supercell[1], supercell[2])
    cwd = os.getcwd()
    for ii in range(len(dss)) :
        struct_path = os.path.join(task_path, 'struct-%s-%03d' % (copy_str,ii))
        print('# generate %s' % (struct_path))
        os.makedirs(struct_path, exist_ok=True)
        os.chdir(struct_path)
        for jj in ['conf.lmp', 'lammps.in'] + meam_potfile_name :
            if os.path.isfile(jj):
                os.remove(jj)
        # make conf
        dss[ii].to('POSCAR', 'POSCAR')
        lammps.cvt_lammps_conf('POSCAR', 'conf.lmp')
        ptypes = vasp.get_poscar_types('POSCAR')
        lammps.apply_type_map('conf.lmp', type_map, ptypes)    
        # link lammps.in
        os.symlink(os.path.relpath(f_lammps_in), 'lammps.in')
        # link models
        for (ii,jj) in zip(meam_potfile, meam_potfile_name) :
            os.symlink(os.path.relpath(ii), jj)
        # save supercell
        np.savetxt('supercell.out', supercell, fmt='%d')
    os.chdir(cwd)
    
def _main() :
    parser = argparse.ArgumentParser(
        description="gen 03.vacancy")
    parser.add_argument('TASK', type=str,
                        help='the task of generation, vasp or lammps')
    parser.add_argument('PARAM', type=str,
                        help='json parameter file')
    parser.add_argument('CONF', type=str,
                        help='the path to conf')
    parser.add_argument('COPY', type=int, nargs = 3,
                        help='the path to conf')
    args = parser.parse_args()

    with open (args.PARAM, 'r') as fp :
        jdata = json.load (fp)

#    print('# generate %s task with conf %s' % (args.TASK, args.CONF))
    if args.TASK == 'vasp':
        make_vasp(jdata, args.CONF, args.COPY)
    elif args.TASK == 'deepmd' :
        make_deepmd_lammps(jdata, args.CONF, args.COPY)
    elif args.TASK == 'meam' :
        make_meam_lammps(jdata, args.CONF, args.COPY)
    else :
        raise RuntimeError("unknow task ", args.TASK)
    
if __name__ == '__main__' :
    _main()

    
