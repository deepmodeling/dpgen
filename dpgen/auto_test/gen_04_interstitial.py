#!/usr/bin/env python3

import os, re, argparse, filecmp, json, glob, warnings, shutil
import subprocess as sp
import numpy as np
import dpgen.auto_test.lib.vasp as vasp
import dpgen.auto_test.lib.lammps as lammps
from pymatgen.core.structure import Structure
from pymatgen.analysis.defects.core import Interstitial
from pymatgen.analysis.defects.generators import InterstitialGenerator

from dpgen import ROOT_PATH
from pymatgen.io.vasp import Incar
from dpgen.generator.lib.vasp import incar_upper
cvasp_file=os.path.join(ROOT_PATH,'generator/lib/cvasp.py')


global_equi_name = '00.equi'
global_task_name = '04.interstitial'

def make_vasp(json, conf_dir, supercell, insert_ele) :
    for ii in insert_ele :
        _make_vasp(json, conf_dir, supercell, ii)

def _make_vasp(jdata, conf_dir, supercell, insert_ele) :

    conf_path = os.path.abspath(conf_dir)
    conf_poscar = os.path.join(conf_path, 'POSCAR')
    # get equi poscar
    if 'relax_incar' in jdata.keys():
        vasp_str='vasp-relax_incar'
    else:
        kspacing = jdata['vasp_params']['kspacing']
        vasp_str='vasp-k%.2f' % (kspacing)
    equi_path = re.sub('confs', global_equi_name, conf_path)
    equi_path = os.path.join(equi_path, vasp_str)
    equi_contcar = os.path.join(equi_path, 'CONTCAR')
    assert os.path.exists(equi_contcar),"Please compute the equilibrium state using vasp first"
    task_path = re.sub('confs', global_task_name, conf_path)
    task_path = os.path.join(task_path, vasp_str)
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
    if  'relax_incar' in jdata.keys():
        relax_incar_path = jdata['relax_incar']
        assert(os.path.exists(relax_incar_path))
        relax_incar_path = os.path.abspath(relax_incar_path)
        incar = incar_upper(Incar.from_file(relax_incar_path))
        fc = incar.get_string()
        kspacing = incar['KSPACING']
        kgamma = incar['KGAMMA']

    else :
        fp_params = jdata['vasp_params']
        ecut = fp_params['ecut']
        ediff = fp_params['ediff']
        npar = fp_params['npar']
        kpar = fp_params['kpar']
        kspacing = fp_params['kspacing']
        kgamma = fp_params['kgamma']
        fc = vasp.make_vasp_relax_incar(ecut, ediff, True, True, True, npar=npar,kpar=kpar, kspacing = kspacing, kgamma = kgamma)
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
        with open('POSCAR','r') as fp :
            lines = fp.read().split('\n')
            ele_list = lines[5].split()

        os.chdir(cwd)
        potcar_map = jdata['potcar_map']
        potcar_list = []
        for ii in ele_list :
            assert os.path.exists(os.path.abspath(potcar_map[ii])),"No POTCAR in the potcar_map of %s"%(ii)
            potcar_list.append(os.path.abspath(potcar_map[ii]))
        os.chdir(struct_path)

        with open('POTCAR', 'w') as outfile:
            for fname in potcar_list:
                with open(fname) as infile:
                    outfile.write(infile.read())

        # link incar
        os.symlink(os.path.relpath(os.path.join(task_path, 'INCAR')), 'INCAR')
        # save supercell
        np.savetxt('supercell.out', supercell, fmt='%d')

        # write kp
        fc = vasp.make_kspacing_kpoints('POSCAR', kspacing, kgamma)
        with open('KPOINTS', 'w') as fp: fp.write(fc)

        #copy cvasp
        if ('cvasp' in jdata) and (jdata['cvasp'] == True):
           shutil.copyfile(cvasp_file, os.path.join(struct_path,'cvasp.py'))
    os.chdir(cwd)


def make_reprod_traj(jdata, conf_dir, supercell, insert_ele, task_type) :
    for ii in insert_ele :
        _make_reprod_traj(jdata, conf_dir, supercell, ii, task_type)

def _make_reprod_traj(jdata, conf_dir, supercell, insert_ele, task_type) :

    fp_params = jdata['lammps_params']
    model_dir = fp_params['model_dir']
    type_map = fp_params['type_map']
    model_dir = os.path.abspath(model_dir)
    model_name =fp_params['model_name']
    deepmd_version = fp_params.get("deepmd_version", "0.12")
    if not model_name and task_type=='deepmd':
        models = glob.glob(os.path.join(model_dir, '*pb'))
        model_name = [os.path.basename(ii) for ii in models]
        assert len(model_name)>0,"No deepmd model in the model_dir"
    else:
        models = [os.path.join(model_dir,ii) for ii in model_name]

    model_param = {'model_name' :      model_name,
                  'param_type':          fp_params['model_param_type'],
                  'deepmd_version' : deepmd_version}

    ntypes = len(type_map)

    conf_path = os.path.abspath(conf_dir)
    task_path = re.sub('confs', global_task_name, conf_path)
    if 'relax_incar' in jdata.keys():
        vasp_str='vasp-relax_incar'
        lmps_str= task_type + '-reprod-relax_incar'
    else:
        kspacing = jdata['vasp_params']['kspacing']
        vasp_str = 'vasp-k%.2f' % (kspacing)
        lmps_str = task_type + '-reprod-k%.2f' % (kspacing)

    vasp_path = os.path.join(task_path, vasp_str)
    lmps_path = os.path.join(task_path, lmps_str)

    os.makedirs(lmps_path, exist_ok = True)
    copy_str = "%sx%sx%s" % (supercell[0], supercell[1], supercell[2])
    struct_widecard = os.path.join(vasp_path, 'struct-%s-%s-*' % (insert_ele,copy_str))
    vasp_struct = glob.glob(struct_widecard)
    assert len(vasp_struct)>0 ,"Please compute the interstitial defect using vasp first"
    vasp_struct.sort()
    cwd=os.getcwd()

    # make lammps.in
    if task_type =='deepmd':
        fc = lammps.make_lammps_eval('conf.lmp',
                                 ntypes,
                                 lammps.inter_deepmd,
                                 model_param)
    elif task_type =='meam':
        fc = lammps.make_lammps_eval('conf.lmp',
                                 ntypes,
                                 lammps.inter_meam,
                                 model_param)
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
        #link lammps.in and model
        for jj in ['lammps.in'] + model_name :
            if os.path.islink(jj):
                os.unlink(jj)
        os.symlink(os.path.relpath(f_lammps_in), 'lammps.in')

        for ii in model_name :
            if os.path.exists(ii) :
                os.remove(ii)
        for (ii,jj) in zip(models, model_name) :
            os.symlink(os.path.relpath(ii), jj)
        share_models = [os.path.join(ls,ii) for ii in model_name]

        # loop over frames
        for ii in range(xdat_nframes) :
            frame_path = 'frame.%06d' % ii
            os.makedirs(frame_path, exist_ok=True)
            os.chdir(frame_path)
            # clear dir
            for jj in ['conf.lmp'] :
                if os.path.isfile(jj):
                    os.remove(jj)
            for jj in ['lammps.in'] + model_name :
                if os.path.islink(jj):
                    os.unlink(jj)
            # link lammps in
            os.symlink(os.path.relpath('../lammps.in'), 'lammps.in')
            # make conf
            with open('POSCAR', 'w') as fp :
                fp.write('\n'.join(xdat_lines[ii*xdat_secsize:(ii+1)*xdat_secsize]))
            lammps.cvt_lammps_conf('POSCAR', 'conf.lmp')
            ptypes = vasp.get_poscar_types('POSCAR')
            lammps.apply_type_map('conf.lmp', type_map, ptypes)
            # link models
            for (kk,ll) in zip(share_models, model_name) :
                os.symlink(os.path.relpath(kk), ll)
            os.chdir(ls)
        os.chdir(cwd)



def make_lammps(jdata, conf_dir, supercell, insert_ele, task_type) :
    for ii in insert_ele:
        _make_lammps(jdata, conf_dir, supercell, ii, task_type)

def _make_lammps(jdata, conf_dir, supercell, insert_ele, task_type) :

    fp_params = jdata['lammps_params']
    model_dir = fp_params['model_dir']
    type_map = fp_params['type_map']
    model_dir = os.path.abspath(model_dir)
    model_name =fp_params['model_name']
    deepmd_version = fp_params.get("deepmd_version", "0.12")
    if not model_name and task_type=='deepmd':
        models = glob.glob(os.path.join(model_dir, '*pb'))
        model_name = [os.path.basename(ii) for ii in models]
        assert len(model_name)>0,"No deepmd model in the model_dir"
    else:
        models = [os.path.join(model_dir,ii) for ii in model_name]

    model_param = {'model_name' :      model_name,
                  'param_type':          fp_params['model_param_type'],
                  'deepmd_version' : deepmd_version}

    ntypes = len(type_map)

    conf_path = os.path.abspath(conf_dir)
    conf_poscar = os.path.join(conf_path, 'POSCAR')
    # get equi poscar
    equi_path = re.sub('confs', global_equi_name, conf_path)
    if 'relax_incar' in jdata.keys():
        vasp_str='vasp-relax_incar'
    else:
        kspacing = jdata['vasp_params']['kspacing']
        vasp_str='vasp-k%.2f' % (kspacing)
    equi_path = os.path.join(equi_path, vasp_str)
    equi_contcar = os.path.join(equi_path, 'CONTCAR')
    #equi_path = os.path.join(equi_path, task_type)
    #equi_dump = os.path.join(equi_path, 'dump.relax')
    assert os.path.exists(equi_contcar),"Please compute the equilibrium state using vasp first"
    task_path = re.sub('confs', global_task_name, conf_path)
    task_path = os.path.join(task_path, task_type)
    os.makedirs(task_path, exist_ok=True)
    task_poscar = os.path.join(task_path, 'POSCAR')
    cwd = os.getcwd()
    os.chdir(task_path)
    if os.path.isfile('POSCAR') :
        os.remove('POSCAR')
    os.symlink(os.path.relpath(equi_contcar), 'POSCAR')
    #lammps.poscar_from_last_dump(equi_dump, task_poscar, type_map)
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
    if task_type=='deepmd':
        fc = lammps.make_lammps_press_relax('conf.lmp',
                                        ntypes,
                                        1,
                                        lammps.inter_deepmd,
                                        model_param)
    elif task_type =='meam':
        fc = lammps.make_lammps_press_relax('conf.lmp',
                                        ntypes,
                                        1,
                                        lammps.inter_meam,
                                        model_param)
    f_lammps_in = os.path.join(task_path, 'lammps.in')
    with open(f_lammps_in, 'w') as fp :
        fp.write(fc)
    # gen tasks
    copy_str = "%sx%sx%s" % (supercell[0], supercell[1], supercell[2])
    cwd = os.getcwd()

    os.chdir(task_path)
    for ii in model_name :
        if os.path.exists(ii) :
            os.remove(ii)
    for (ii,jj) in zip(models, model_name) :
        os.symlink(os.path.relpath(ii), jj)
    share_models = [os.path.join(task_path,ii) for ii in model_name]

    for ii in range(len(dss)) :
        struct_path = os.path.join(task_path, 'struct-%s-%s-%03d' % (insert_ele,copy_str,ii))
        print('# generate %s' % (struct_path))
        os.makedirs(struct_path, exist_ok=True)
        os.chdir(struct_path)
        for jj in ['conf.lmp', 'lammps.in'] + model_name :
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
        for (ii,jj) in zip(share_models, model_name) :
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
    elif args.TASK == 'deepmd' or args.TASK=='meam':
        make_lammps(jdata, args.CONF, args.COPY, args.ELEMENT, args.TASK)
    elif args.TASK == 'deepmd-reprod' or args.TASK == 'meam-reprod':
        args.TASK=args.TASK.replace('-reprod','')
        make_reprod_traj(jdata, args.CONF, args.COPY, args.ELEMENT, args.TASK)
    else :
        raise RuntimeError("unknow task ", args.TASK)

if __name__ == '__main__' :
    _main()
