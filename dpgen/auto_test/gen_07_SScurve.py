#!/usr/bin/env python3

import os, re, argparse, filecmp, json, glob
import subprocess as sp
import numpy as np
import dpgen.auto_test.lib.vasp as vasp
import dpgen.auto_test.lib.lammps as lammps
from pymatgen.core.structure import Structure
from pymatgen.analysis.elasticity.strain import Deformation, DeformedStructureSet, Strain

global_equi_name = '00.equi'
global_task_name = '07.SScurve'

def make_vasp(jdata, conf_dir, norm_def = 2e-3, shear_def = 5e-3) :
    fp_params = jdata['vasp_params']
    ecut = fp_params['ecut']
    ediff = fp_params['ediff']
    npar = fp_params['npar']
    kpar = fp_params['kpar']
    kspacing = fp_params['kspacing']
    kgamma = fp_params['kgamma']
    strain_start=jdata['strain_start']
    strain_end=jdata['strain_end']
    strain_step=jdata['strain_step']
    strain_direct=jdata['strain_direct']

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
    # stress
    equi_outcar = os.path.join(equi_path, 'OUTCAR')
    stress = vasp.get_stress(equi_outcar)
    np.savetxt(os.path.join(task_path, 'equi.stress.out'), stress)
    # gen strcture
    ss = Structure.from_file(task_poscar)
    # gen defomations
    norm_strains=np.arange(strain_start,strain_end,strain_step)
    print('gen with norm '+str(norm_strains))
    deformations=[]
    for ii in norm_strains:
        strain = Strain.from_index_amount(strain_direct, ii)
        deformations.append(strain.get_deformation_matrix())
    deformed_structures = [defo.apply_to_structure(ss)
                                    for defo in deformations]
    n_dfm = len(deformed_structures)
    # gen incar
    fc = vasp.make_vasp_relax_incar(ecut, ediff, True, False, False, npar=npar, kpar=kpar, kspacing = None, kgamma = None)
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
    # gen tasks    
    cwd = os.getcwd()
    for ii in range(n_dfm) :
        # make dir
        dfm_path = os.path.join(task_path, 'dfm-%03d' % ii)
        os.makedirs(dfm_path, exist_ok=True)
        os.chdir(dfm_path)
        for jj in ['POSCAR', 'POTCAR', 'INCAR', 'KPOINTS'] :
            if os.path.isfile(jj):
                os.remove(jj)
        # make conf
        deformed_structures[ii].to('POSCAR', 'POSCAR')
        # record strain
        strain = Strain.from_deformation(deformations[ii])
        np.savetxt('strain.out', strain)
        # link incar, potcar, kpoints
        os.symlink(os.path.relpath(os.path.join(task_path, 'INCAR')), 'INCAR')
        os.symlink(os.path.relpath(os.path.join(task_path, 'POTCAR')), 'POTCAR')
        os.symlink(os.path.relpath(os.path.join(task_path, 'KPOINTS')), 'KPOINTS')
    cwd = os.getcwd()

def make_lammps(jdata, conf_dir,task_type) :
    fp_params = jdata['lammps_params']
    model_dir = fp_params['model_dir']
    type_map = fp_params['type_map'] 
    model_dir = os.path.abspath(model_dir)
    model_name =fp_params['model_name']
    deepmd_version = fp_params.get("deepmd_version", "0.12")
    if not model_name :
        models = glob.glob(os.path.join(model_dir, '*pb'))
        model_name = [os.path.basename(ii) for ii in models]
    else:
        models = [os.path.join(model_dir,ii) for ii in model_name]

    model_param = {'model_name' :      model_name,
                  'param_type':          fp_params['model_param_type'],
                  'deepmd_version' : deepmd_version}
    ntypes = len(type_map)
    strain_start=jdata['strain_start']
    strain_end=jdata['strain_end']
    strain_step=jdata['strain_step']
    strain_direct=jdata['strain_direct']

    conf_path = os.path.abspath(conf_dir)
    conf_poscar = os.path.join(conf_path, 'POSCAR')
    # get equi poscar
    equi_path = re.sub('confs', global_equi_name, conf_path)
    equi_path = os.path.join(equi_path, task_type)
    equi_dump = os.path.join(equi_path, 'dump.relax')
    task_path = re.sub('confs', global_task_name, conf_path)
    task_path = os.path.join(task_path, task_type)
    os.makedirs(task_path, exist_ok=True)
    task_poscar = os.path.join(task_path, 'POSCAR')
    lammps.poscar_from_last_dump(equi_dump, task_poscar, type_map)
    # get equi stress
    equi_log = os.path.join(equi_path, 'log.lammps')
    stress = lammps.get_stress(equi_log)
    np.savetxt(os.path.join(task_path, 'equi.stress.out'), stress)
    # gen strcture
    ss = Structure.from_file(task_poscar)
    # gen defomations
    norm_strains=np.arange(strain_start,strain_end,strain_step)
    print('gen with norm '+str(norm_strains))
    deformations=[]
    for ii in norm_strains:
        strain = Strain.from_index_amount(strain_direct, ii)
        deformations.append(strain.get_deformation_matrix())
    deformed_structures = [defo.apply_to_structure(ss)
                                    for defo in deformations]
    n_dfm = len(deformed_structures)
    # gen tasks    
    cwd = os.getcwd()
    # make lammps.in
    if task_type=='deepmd':
        fc = lammps.make_lammps_elastic('conf.lmp', 
                                    ntypes, 
                                    lammps.inter_deepmd,
                                    model_param)
    elif task_type =='meam':
        fc = lammps.make_lammps_elastic('conf.lmp', 
                                    ntypes, 
                                    lammps.inter_meam,
                                    model_param) 


    f_lammps_in = os.path.join(task_path, 'lammps.in')
    with open(f_lammps_in, 'w') as fp :
        fp.write(fc)
    cwd = os.getcwd()
    if task_type =='deepmd':
        os.chdir(task_path)
        for ii in model_name :
            if os.path.exists(ii) :
                os.remove(ii)
        for (ii,jj) in zip(models, model_name) :
                os.symlink(os.path.relpath(ii), jj)
        share_models = glob.glob(os.path.join(task_path, '*pb'))
    else:
        share_models = models

    for ii in range(n_dfm) :
        # make dir
        dfm_path = os.path.join(task_path, 'dfm-%03d' % ii)
        os.makedirs(dfm_path, exist_ok=True)
        os.chdir(dfm_path)
        for jj in ['conf.lmp', 'lammps.in'] + model_name :
            if os.path.isfile(jj):
                os.remove(jj)
        # make conf
        deformed_structures[ii].to('POSCAR', 'POSCAR')
        lammps.cvt_lammps_conf('POSCAR', 'conf.lmp')
        ptypes = vasp.get_poscar_types('POSCAR')
        lammps.apply_type_map('conf.lmp', type_map, ptypes)    
        # record strain
        strain = Strain.from_deformation(deformations[ii])
        np.savetxt('strain.out', strain)
        # link lammps.in
        os.symlink(os.path.relpath(f_lammps_in), 'lammps.in')
        # link models
        for (ii,jj) in zip(share_models, model_name) :
            os.symlink(os.path.relpath(ii), jj)
    cwd = os.getcwd()

def make_meam_lammps(jdata, conf_dir) :
    meam_potfile_dir = jdata['meam_potfile_dir']
    meam_potfile_dir = os.path.abspath(meam_potfile_dir)
    meam_potfile = jdata['meam_potfile']
    meam_potfile = [os.path.join(meam_potfile_dir,ii) for ii in meam_potfile]
    meam_potfile_name = jdata['meam_potfile']
    type_map = jdata['meam_type_map']
    ntypes = len(type_map)
    meam_param = {'meam_potfile' :      jdata['meam_potfile'],
                  'meam_type':          jdata['meam_param_type']}

    norm_def = jdata['norm_deform']
    shear_def = jdata['shear_deform']

    conf_path = os.path.abspath(conf_dir)
    conf_poscar = os.path.join(conf_path, 'POSCAR')
    # get equi poscar
    equi_path = re.sub('confs', global_equi_name, conf_path)
    equi_path = os.path.join(equi_path, 'meam')
    equi_dump = os.path.join(equi_path, 'dump.relax')
    task_path = re.sub('confs', global_task_name, conf_path)
    task_path = os.path.join(task_path, 'meam')
    os.makedirs(task_path, exist_ok=True)
    task_poscar = os.path.join(task_path, 'POSCAR')
    lammps.poscar_from_last_dump(equi_dump, task_poscar, type_map)
    # get equi stress
    equi_log = os.path.join(equi_path, 'log.lammps')
    stress = lammps.get_stress(equi_log)
    np.savetxt(os.path.join(task_path, 'equi.stress.out'), stress)
    # gen strcture
    # ss = Structure.from_file(conf_poscar)
    # print(ss)
    # ss = ss.from_file(task_poscar)
    # print(ss)
    ss = Structure.from_file(task_poscar)
    # gen defomations
    norm_strains = [-norm_def, -0.5*norm_def, 0.5*norm_def, norm_def]
    shear_strains = [-shear_def, -0.5*shear_def, 0.5*shear_def, shear_def]
    print('gen with norm '+str(norm_strains))
    print('gen with shear '+str(shear_strains))
    dfm_ss = DeformedStructureSet(ss, 
                                  symmetry = False, 
                                  norm_strains = norm_strains,
                                  shear_strains = shear_strains)
    n_dfm = len(dfm_ss)
    # gen tasks    
    cwd = os.getcwd()
    # make lammps.in
    fc = lammps.make_lammps_elastic('conf.lmp', 
                                    ntypes, 
                                    lammps.inter_meam,
                                    meam_param)        
    f_lammps_in = os.path.join(task_path, 'lammps.in')
    with open(f_lammps_in, 'w') as fp :
        fp.write(fc)
    cwd = os.getcwd()
    for ii in range(n_dfm) :
        # make dir
        dfm_path = os.path.join(task_path, 'dfm-%03d' % ii)
        os.makedirs(dfm_path, exist_ok=True)
        os.chdir(dfm_path)
        for jj in ['conf.lmp', 'lammps.in'] + meam_potfile_name :
            if os.path.isfile(jj):
                os.remove(jj)
        # make conf
        dfm_ss.deformed_structures[ii].to('POSCAR', 'POSCAR')
        lammps.cvt_lammps_conf('POSCAR', 'conf.lmp')
        ptypes = vasp.get_poscar_types('POSCAR')
        lammps.apply_type_map('conf.lmp', type_map, ptypes)    
        # record strain
        strain = Strain.from_deformation(dfm_ss.deformations[ii])
        np.savetxt('strain.out', strain)
        # link lammps.in
        os.symlink(os.path.relpath(f_lammps_in), 'lammps.in')
        # link models
        for (ii,jj) in zip(meam_potfile, meam_potfile_name) :
            os.symlink(os.path.relpath(ii), jj)
    cwd = os.getcwd()

    
def _main() :
    parser = argparse.ArgumentParser(
        description="gen 07.SScurve")
    parser.add_argument('TASK', type=str,
                        help='the task of generation, vasp or lammps')
    parser.add_argument('PARAM', type=str,
                        help='json parameter file')
    parser.add_argument('CONF', type=str,
                        help='the path to conf')
    args = parser.parse_args()

    with open (args.PARAM, 'r') as fp :
        jdata = json.load (fp)

    print('generate %s task with conf %s' % (args.TASK, args.CONF))
    if args.TASK == 'vasp':
        make_vasp(jdata, args.CONF)               
    elif args.TASK == 'deepmd' or args.TASK == 'meam' :
        make_lammps(jdata, args.CONF,args.TASK)
    #elif args.TASK == 'meam' :
    #    make_meam_lammps(jdata, args.CONF)
    else :
        raise RuntimeError("unknow task ", args.TASK)
    
if __name__ == '__main__' :
    _main()

