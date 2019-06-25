#!/usr/bin/env python3

import os, re, argparse, filecmp, json, glob, warnings
import subprocess as sp
import numpy as np
import lib.vasp as vasp
import lib.lammps as lammps
from pymatgen.core.structure import Structure
from pymatgen.analysis.defects.core import Interstitial
from pymatgen.analysis.defects.generators import InterstitialGenerator

global_equi_name = '00.equi'
global_task_name = '04.interstitial'

def make_vasp(json, conf_dir, supercell, insert_ele) :
    for ii in insert_ele :
        _make_vasp(json, conf_dir, supercell, ii)

def _make_vasp(jdata, conf_dir, supercell, insert_ele) :
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
    print("task poscar: ", task_poscar)
    ss = Structure.from_file(task_poscar)
    # gen defects
    vds = InterstitialGenerator(ss, insert_ele)
    dss = []
    for jj in vds :
        dss.append(jj.generate_defect_structure(supercell))
    # gen incar
    fc = vasp.make_vasp_relax_incar(ecut, ediff, True, True, True, 1, 1, kspacing = kspacing, kgamma = kgamma)
    with open(os.path.join(task_path, 'INCAR'), 'w') as fp :
        fp.write(fc)
    # gen tasks    
    copy_str = "%sx%sx%s" % (supercell[0], supercell[1], supercell[2])
    cwd = os.getcwd()
    for ii in range(len(dss)) :
        struct_path = os.path.join(task_path, 'struct-%s-%s-%03d' % (insert_ele,copy_str,ii))
        print('# generate %s' % (struct_path))
        os.makedirs(struct_path, exist_ok=True)
        os.chdir(struct_path)
        for jj in ['POSCAR', 'POTCAR', 'INCAR'] :
            if os.path.isfile(jj):
                os.remove(jj)
        # make conf
        dss[ii].to('POSCAR', 'POSCAR')
        # gen potcar
        _gen_potcar(jdata, 'POSCAR', 'POTCAR')
        # link incar
        os.symlink(os.path.relpath(os.path.join(task_path, 'INCAR')), 'INCAR')
        # save supercell
        np.savetxt('supercell.out', supercell, fmt='%d')
    os.chdir(cwd)

def _gen_potcar (jdata, task_poscar, filename) :
    # gen potcar
    with open(task_poscar,'r') as fp :
        lines = fp.read().split('\n')
        ele_list = lines[5].split()
    potcar_map = jdata['potcar_map']
    potcar_list = []
    for ii in ele_list :
        assert(os.path.exists(potcar_map[ii]))
        potcar_list.append(potcar_map[ii])
    with open(filename, 'w') as outfile:
        for fname in potcar_list:
            with open(fname) as infile:
                outfile.write(infile.read())

def make_deepmd_reprod_traj(jdata, conf_dir, supercell, insert_ele, task_name) : 
    for ii in insert_ele :
        _make_deepmd_reprod_traj(jdata, conf_dir, supercell, ii, task_name)

def _make_deepmd_reprod_traj(jdata, conf_dir, supercell, insert_ele, task_name) : 
    fp_params = jdata['vasp_params']
    kspacing = fp_params['kspacing']
    deepmd_model_dir = jdata['deepmd_model_dir']
    deepmd_type_map = jdata['deepmd_type_map']
    ntypes = len(deepmd_type_map)    
    deepmd_model_dir = os.path.abspath(deepmd_model_dir)
    deepmd_models = glob.glob(os.path.join(deepmd_model_dir, '*pb'))
    deepmd_models_name = [os.path.basename(ii) for ii in deepmd_models]

    conf_path = os.path.abspath(conf_dir)
    task_path = re.sub('confs', global_task_name, conf_path)
    vasp_path = os.path.join(task_path, 'vasp-k%.2f' % kspacing)
    lmps_path = os.path.join(task_path, task_name + '-k%.2f' % kspacing)    
    os.makedirs(lmps_path, exist_ok = True)
    copy_str = "%sx%sx%s" % (supercell[0], supercell[1], supercell[2])
    struct_widecard = os.path.join(vasp_path, 'struct-%s-%s-*' % (insert_ele,copy_str))
    vasp_struct = glob.glob(struct_widecard)
    vasp_struct.sort()
    cwd=os.getcwd()
    
    # make lammps.in
    fc = lammps.make_lammps_eval('conf.lmp', 
                                 ntypes, 
                                 lammps.inter_deepmd,
                                 deepmd_models_name)
    f_lammps_in = os.path.join(lmps_path, 'lammps.in')
    with open(f_lammps_in, 'w') as fp :
        fp.write(fc)

    for vs in vasp_struct :
        # get vasp energy
        outcar = os.path.join(vs, 'OUTCAR')
        energies = vasp.get_energies(outcar)
        # get xdat
        xdatcar = os.path.join(vs, 'XDATCAR')
        struct_basename  = os.path.basename(vs)
        ls = os.path.join(lmps_path, struct_basename)
        print(ls)
        os.makedirs(ls, exist_ok = True)
        os.chdir(ls)
        if os.path.exists('XDATCAR') :
            os.remove('XDATCAR')
        os.symlink(os.path.relpath(xdatcar), 'XDATCAR')
        xdat_lines = open('XDATCAR', 'r').read().split('\n')
        natoms = vasp.poscar_natoms('XDATCAR')
        xdat_secsize = natoms + 8
        xdat_nframes = len(xdat_lines) // xdat_secsize
        if xdat_nframes > len(energies) :
            warnings.warn('nframes %d in xdat is larger than energy %d, use the last %d frames' % (xdat_nframes, len(energies), len(energies)))
            xdat_nlines = len(energies) * xdat_secsize
            xdat_lines = xdat_lines[xdat_nlines:]
        xdat_nframes = len(xdat_lines) // xdat_secsize
        print(xdat_nframes, len(energies))
        # loop over frames
        for ii in range(xdat_nframes) :
            frame_path = 'frame.%06d' % ii
            os.makedirs(frame_path, exist_ok=True)
            os.chdir(frame_path)
            # clear dir
            for jj in ['conf.lmp'] :
                if os.path.isfile(jj):
                    os.remove(jj)            
            for jj in ['lammps.in'] + deepmd_models_name :
                if os.path.islink(jj):
                    os.unlink(jj)            
            # link lammps in
            os.symlink(os.path.relpath(f_lammps_in), 'lammps.in')
            # make conf
            with open('POSCAR', 'w') as fp :
                fp.write('\n'.join(xdat_lines[ii*xdat_secsize:(ii+1)*xdat_secsize]))
            lammps.cvt_lammps_conf('POSCAR', 'conf.lmp')
            ptypes = vasp.get_poscar_types('POSCAR')
            lammps.apply_type_map('conf.lmp', deepmd_type_map, ptypes)
            # link models
            for (kk,ll) in zip(deepmd_models, deepmd_models_name) :
                os.symlink(os.path.relpath(kk), ll)
            os.chdir(ls)
        os.chdir(cwd)

def make_meam_reprod_traj(jdata, conf_dir, supercell, insert_ele, task_name) : 
    for ii in insert_ele :
        _make_meam_reprod_traj(jdata, conf_dir, supercell, ii, task_name)

def _make_meam_reprod_traj(jdata, conf_dir, supercell, insert_ele, task_name) : 
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
    task_path = re.sub('confs', global_task_name, conf_path)
    vasp_path = os.path.join(task_path, 'vasp-k%.2f' % kspacing)
    lmps_path = os.path.join(task_path, task_name + '-k%.2f' % kspacing)    
    os.makedirs(lmps_path, exist_ok = True)
    copy_str = "%sx%sx%s" % (supercell[0], supercell[1], supercell[2])
    struct_widecard = os.path.join(vasp_path, 'struct-%s-%s-*' % (insert_ele,copy_str))
    vasp_struct = glob.glob(struct_widecard)
    vasp_struct.sort()
    cwd=os.getcwd()
    
    # make lammps.in
    fc = lammps.make_lammps_eval('conf.lmp', 
                                 ntypes, 
                                 lammps.inter_meam,
                                 meam_param)
    f_lammps_in = os.path.join(lmps_path, 'lammps.in')
    with open(f_lammps_in, 'w') as fp :
        fp.write(fc)

    for vs in vasp_struct :
        # get vasp energy
        outcar = os.path.join(vs, 'OUTCAR')
        energies = vasp.get_energies(outcar)
        # get xdat
        xdatcar = os.path.join(vs, 'XDATCAR')
        struct_basename  = os.path.basename(vs)
        ls = os.path.join(lmps_path, struct_basename)
        print(ls)
        os.makedirs(ls, exist_ok = True)
        os.chdir(ls)
        if os.path.exists('XDATCAR') :
            os.remove('XDATCAR')
        os.symlink(os.path.relpath(xdatcar), 'XDATCAR')
        xdat_lines = open('XDATCAR', 'r').read().split('\n')
        natoms = vasp.poscar_natoms('XDATCAR')
        xdat_secsize = natoms + 8
        xdat_nframes = len(xdat_lines) // xdat_secsize
        if xdat_nframes > len(energies) :
            warnings.warn('nframes %d in xdat is larger than energy %d, use the last %d frames' % (xdat_nframes, len(energies), len(energies)))
            xdat_nlines = len(energies) * xdat_secsize
            xdat_lines = xdat_lines[xdat_nlines:]
        xdat_nframes = len(xdat_lines) // xdat_secsize
        print(xdat_nframes, len(energies))
        # loop over frames
        for ii in range(xdat_nframes) :
            frame_path = 'frame.%06d' % ii
            os.makedirs(frame_path, exist_ok=True)
            os.chdir(frame_path)
            # clear dir
            for jj in ['conf.lmp'] :
                if os.path.isfile(jj):
                    os.remove(jj)            
            for jj in ['lammps.in'] + meam_potfile_name :
                if os.path.islink(jj):
                    os.unlink(jj)            
            # link lammps in
            os.symlink(os.path.relpath(f_lammps_in), 'lammps.in')
            # make conf
            with open('POSCAR', 'w') as fp :
                fp.write('\n'.join(xdat_lines[ii*xdat_secsize:(ii+1)*xdat_secsize]))
            lammps.cvt_lammps_conf('POSCAR', 'conf.lmp')
            ptypes = vasp.get_poscar_types('POSCAR')
            lammps.apply_type_map('conf.lmp', type_map, ptypes)
            # link models
            for (kk,ll) in zip(meam_potfile, meam_potfile_name) :
                os.symlink(os.path.relpath(kk), ll)
            os.chdir(ls)
        os.chdir(cwd)


def make_deepmd_lammps(jdata, conf_dir, supercell, insert_ele, task_name) :
    for ii in insert_ele:
        _make_deepmd_lammps(jdata, conf_dir, supercell, ii, task_name)

def _make_deepmd_lammps(jdata, conf_dir, supercell, insert_ele, task_name) :
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
    equi_path = os.path.join(equi_path, task_name)
    equi_dump = os.path.join(equi_path, 'dump.relax')
    task_path = re.sub('confs', global_task_name, conf_path)
    task_path = os.path.join(task_path, task_name)
    os.makedirs(task_path, exist_ok=True)
    task_poscar = os.path.join(task_path, 'POSCAR')
    cwd = os.getcwd()
    os.chdir(task_path)
    lammps.poscar_from_last_dump(equi_dump, task_poscar, deepmd_type_map)
    os.chdir(cwd)
    # gen structure from equi poscar
    print("task poscar: ", task_poscar)
    ss = Structure.from_file(task_poscar)
    # gen defects
    vds = InterstitialGenerator(ss, insert_ele)
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
        struct_path = os.path.join(task_path, 'struct-%s-%s-%03d' % (insert_ele,copy_str,ii))
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

def make_meam_lammps(jdata, conf_dir, supercell, insert_ele, task_name) :
    for ii in insert_ele:
        _make_meam_lammps(jdata, conf_dir, supercell, ii, task_name)

def _make_meam_lammps(jdata, conf_dir, supercell, insert_ele, task_name) :
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
    equi_path = os.path.join(equi_path, task_name)
    equi_dump = os.path.join(equi_path, 'dump.relax')
    task_path = re.sub('confs', global_task_name, conf_path)
    task_path = os.path.join(task_path, task_name)
    os.makedirs(task_path, exist_ok=True)
    task_poscar = os.path.join(task_path, 'POSCAR')
    cwd = os.getcwd()
    os.chdir(task_path)
    lammps.poscar_from_last_dump(equi_dump, task_poscar, type_map)
    os.chdir(cwd)
    # gen structure from equi poscar
    ss = Structure.from_file(task_poscar)
    # gen defects
    vds = InterstitialGenerator(ss, insert_ele)
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
        struct_path = os.path.join(task_path, 'struct-%s-%s-%03d' % (insert_ele,copy_str,ii))
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
        description="gen 04.interstitial")
    parser.add_argument('TASK', type=str,
                        help='the task of generation, vasp or lammps')
    parser.add_argument('PARAM', type=str,
                        help='json parameter file')
    parser.add_argument('CONF', type=str,
                        help='the path to conf')
    parser.add_argument('COPY', type=int, nargs = 3,
                        help='define the supercell')
    parser.add_argument('ELEMENT', type=str, nargs = '+',
                        help='the inserted element')
    args = parser.parse_args()

    with open (args.PARAM, 'r') as fp :
        jdata = json.load (fp)

#    print('# generate %s task with conf %s' % (args.TASK, args.CONF))
    if args.TASK == 'vasp':
        make_vasp(jdata, args.CONF, args.COPY, args.ELEMENT)
    elif args.TASK == 'deepmd' :
        make_deepmd_lammps(jdata, args.CONF, args.COPY, args.ELEMENT, args.TASK)
    elif args.TASK == 'deepmd-reprod' :
        make_deepmd_reprod_traj(jdata, args.CONF, args.COPY, args.ELEMENT, args.TASK)
    elif args.TASK == 'meam' :
        make_meam_lammps(jdata, args.CONF, args.COPY, args.ELEMENT, args.TASK)
    elif args.TASK == 'meam-reprod' :
        make_meam_reprod_traj(jdata, args.CONF, args.COPY, args.ELEMENT, args.TASK)
    else :
        raise RuntimeError("unknow task ", args.TASK)
    
if __name__ == '__main__' :
    _main()

    
