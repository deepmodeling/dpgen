#!/usr/bin/env python3

import random, os, sys
import dpdata
import subprocess as sp
import dpgen.auto_test.lib.util as util
from distutils.version import LooseVersion
from dpdata.periodic_table import Element


def cvt_lammps_conf(fin, fout, type_map, ofmt='lammps/data'):
    """
    Format convert from fin to fout, specify the output format by ofmt
    Imcomplete situation
    """
    supp_ofmt = ['lammps/dump', 'lammps/data', 'vasp/poscar']
    supp_exts = ['dump', 'lmp', 'poscar/POSCAR']

    if 'dump' in fout:
        ofmt = 'lammps/dump'
    elif 'lmp' in fout:
        ofmt = 'lammps/data'
    elif 'poscar' in fout or 'POSCAR' in fout:
        ofmt = 'vasp/poscar'
    if not ofmt in supp_ofmt:
        raise RuntimeError("output format " + ofmt + " is not supported. use one of " + str(supp_ofmt))

    if 'lmp' in fout:
        d_poscar = dpdata.System(fin, fmt='vasp/poscar', type_map=type_map)
        d_poscar.to_lammps_lmp(fout, frame_idx=0)
    elif 'poscar' in fout or 'POSCAR' in fout:
        d_dump = dpdata.System(fin, fmt='lammps/dump', type_map=type_map)
        d_dump.to_vasp_poscar(fout, frame_idx=-1)


def apply_type_map(conf_file, deepmd_type_map, ptypes):
    """
    apply type map.
    conf_file:          conf file converted from POSCAR
    deepmd_type_map:    deepmd atom type map
    ptypes:             atom types defined in POSCAR
    """
    natoms = _get_conf_natom(conf_file)
    ntypes = len(deepmd_type_map)
    with open(conf_file, 'r') as fp:
        lines = fp.read().split('\n')
    # with open(conf_file+'.bk', 'w') as fp:
    #     fp.write("\n".join(lines))
    new_lines = lines
    # revise ntypes
    idx_ntypes = -1
    for idx, ii in enumerate(lines):
        if 'atom types' in ii:
            idx_ntypes = idx
    if idx_ntypes == -1:
        raise RuntimeError("cannot find the entry 'atom types' in ", conf_file)
    words = lines[idx_ntypes].split()
    words[0] = str(ntypes)
    new_lines[idx_ntypes] = " ".join(words)
    # find number of atoms
    idx_atom_entry = -1
    for idx, ii in enumerate(lines):
        if 'Atoms' in ii:
            idx_atom_entry = idx
    if idx_atom_entry == -1:
        raise RuntimeError("cannot find the entry 'Atoms' in ", conf_file)
    # revise atom type
    for idx in range(idx_atom_entry + 2, idx_atom_entry + 2 + natoms):
        ii = lines[idx]
        words = ii.split()
        assert (len(words) >= 5)
        old_id = int(words[1])
        new_id = deepmd_type_map.index(ptypes[old_id - 1]) + 1
        words[1] = str(new_id)
        ii = " ".join(words)
        new_lines[idx] = ii
    with open(conf_file, 'w') as fp:
        fp.write("\n".join(new_lines))


def _get_ntype(conf):
    with open(conf, 'r') as fp:
        lines = fp.read().split('\n')
    for ii in lines:
        if "atom types" in ii:
            return int(ii.split()[0])
    raise RuntimeError("cannot find line indicate atom types in ", conf)


def _get_conf_natom(conf):
    with open(conf, 'r') as fp:
        lines = fp.read().split('\n')
    for ii in lines:
        if "atoms" in ii:
            return int(ii.split()[0])
    raise RuntimeError("cannot find line indicate atom types in ", conf)


def inter_deepmd(param):
    models = param["model_name"]
    deepmd_version = param["deepmd_version"]
    ret = "pair_style deepmd "
    model_list = ""
    for ii in models:
        model_list += ii + " "
    if LooseVersion(deepmd_version) < LooseVersion('1'):
        ## DeePMD-kit version == 0.x
        if len(models) > 1:
            ret += '%s 10 model_devi.out\n' % model_list
        else:
            ret += models[0] + '\n'
    else:
        ## DeePMD-kit version >= 1
        if len(models) > 1:
            ret += "%s out_freq 10 out_file model_devi.out\n" % model_list
        else:
            ret += models[0] + '\n'
    ret += "pair_coeff * *\n"
    return ret


def inter_meam(param):
    ret = ""
    line = "pair_style      meam \n"
    line += "pair_coeff      * * %s " % param['model_name'][0]
    for ii in param['param_type']:
        line += ii + ' '
    line += "%s " % param['model_name'][1]
    for ii in param['param_type']:
        line += ii + ' '
    line += '\n'
    ret += line
    return ret


def inter_eam_fs(param):  # 06/08 eam.fs interaction
    ret = ""
    line = "pair_style      eam/fs \n"
    line += "pair_coeff      * * %s " % param['model_name'][0]
    for ii in param['param_type']:
        line += ii + ' '
    line += '\n'
    ret += line
    return ret


def inter_eam_alloy(param):  # 06/08 eam.alloy interaction
    ret = ""
    line = "pair_style      eam/alloy \n"
    line += "pair_coeff      * * %s " % param['model_name']
    for ii in param['param_type']:
        line += ii + ' '
    line += '\n'
    ret += line
    return ret


def element_list(type_map):
    type_map_reverse = {k: v for v, k in type_map.items()}
    type_map_list = []
    tmp_list = list(type_map_reverse.keys())
    tmp_list.sort()
    for ii in tmp_list:
        type_map_list.append(type_map_reverse[ii])
    return type_map_list


def make_lammps_eval(conf, type_map, interaction, param):
    type_map_list = element_list(type_map)

    """
    make lammps input for static calcualtion
    """
    ret = ""
    ret += "clear\n"
    ret += "units 	metal\n"
    ret += "dimension	3\n"
    ret += "boundary	p p p\n"
    ret += "atom_style	atomic\n"
    ret += "box         tilt large\n"
    ret += "read_data   %s\n" % conf
    for ii in range(len(type_map)):
        ret += "mass            %d %.3f\n" % (ii + 1, Element(type_map_list[ii]).mass)
    ret += "neigh_modify    every 1 delay 0 check no\n"
    ret += interaction(param)
    ret += "compute         mype all pe\n"
    ret += "thermo          100\n"
    ret += "thermo_style    custom step pe pxx pyy pzz pxy pxz pyz lx ly lz vol c_mype\n"
    ret += "dump            1 all custom 100 dump.relax id type xs ys zs fx fy fz\n"  # 06/09 give dump.relax
    ret += "run    0\n"
    ret += "variable        N equal count(all)\n"
    ret += "variable        V equal vol\n"
    ret += "variable        E equal \"c_mype\"\n"
    ret += "variable        tmplx equal lx\n"
    ret += "variable        tmply equal ly\n"
    ret += "variable        Pxx equal pxx\n"
    ret += "variable        Pyy equal pyy\n"
    ret += "variable        Pzz equal pzz\n"
    ret += "variable        Pxy equal pxy\n"
    ret += "variable        Pxz equal pxz\n"
    ret += "variable        Pyz equal pyz\n"
    ret += "variable        Epa equal ${E}/${N}\n"
    ret += "variable        Vpa equal ${V}/${N}\n"
    ret += "variable        AA equal (${tmplx}*${tmply})\n"
    ret += "print \"All done\"\n"
    ret += "print \"Total number of atoms = ${N}\"\n"
    ret += "print \"Final energy per atoms = ${Epa}\"\n"
    ret += "print \"Final volume per atoms = ${Vpa}\"\n"
    ret += "print \"Final Base area = ${AA}\"\n"
    ret += "print \"Final Stress (xx yy zz xy xz yz) = ${Pxx} ${Pyy} ${Pzz} ${Pxy} ${Pxz} ${Pyz}\"\n"
    return ret


def make_lammps_equi(conf, type_map, interaction, param,
                     etol=0, ftol=1e-10,
                     maxiter=5000, maxeval=500000,
                     change_box=True):
    type_map_list = element_list(type_map)

    """
    make lammps input for equilibritation
    """
    ret = ""
    ret += "clear\n"
    ret += "units 	metal\n"
    ret += "dimension	3\n"
    ret += "boundary	p p p\n"
    ret += "atom_style	atomic\n"
    ret += "box         tilt large\n"
    ret += "read_data   %s\n" % conf
    for ii in range(len(type_map)):
        ret += "mass            %d %.3f\n" % (ii + 1, Element(type_map_list[ii]).mass)
    ret += "neigh_modify    every 1 delay 0 check no\n"
    ret += interaction(param)
    ret += "compute         mype all pe\n"
    ret += "thermo          100\n"
    ret += "thermo_style    custom step pe pxx pyy pzz pxy pxz pyz lx ly lz vol c_mype\n"
    ret += "dump            1 all custom 100 dump.relax id type xs ys zs fx fy fz\n"
    ret += "min_style       cg\n"
    if change_box:
        ret += "fix             1 all box/relax iso 0.0 \n"
        ret += "minimize        %e %e %d %d\n" % (etol, ftol, maxiter, maxeval)
        ret += "fix             1 all box/relax aniso 0.0 \n"
        ret += "minimize        %e %e %d %d\n" % (etol, ftol, maxiter, maxeval)
        ret += "fix             1 all box/relax tri 0.0 \n"
    ret += "minimize        %e %e %d %d\n" % (etol, ftol, maxiter, maxeval)
    ret += "variable        N equal count(all)\n"
    ret += "variable        V equal vol\n"
    ret += "variable        E equal \"c_mype\"\n"
    ret += "variable        tmplx equal lx\n"
    ret += "variable        tmply equal ly\n"
    ret += "variable        Pxx equal pxx\n"
    ret += "variable        Pyy equal pyy\n"
    ret += "variable        Pzz equal pzz\n"
    ret += "variable        Pxy equal pxy\n"
    ret += "variable        Pxz equal pxz\n"
    ret += "variable        Pyz equal pyz\n"
    ret += "variable        Epa equal ${E}/${N}\n"
    ret += "variable        Vpa equal ${V}/${N}\n"
    ret += "variable        AA equal (${tmplx}*${tmply})\n"
    ret += "print \"All done\"\n"
    ret += "print \"Total number of atoms = ${N}\"\n"
    ret += "print \"Final energy per atoms = ${Epa}\"\n"
    ret += "print \"Final volume per atoms = ${Vpa}\"\n"
    ret += "print \"Final Base area = ${AA}\"\n"
    ret += "print \"Final Stress (xx yy zz xy xz yz) = ${Pxx} ${Pyy} ${Pzz} ${Pxy} ${Pxz} ${Pyz}\"\n"
    return ret


def make_lammps_elastic(conf, type_map, interaction, param,
                        etol=0, ftol=1e-10,
                        maxiter=5000, maxeval=500000):
    type_map_list = element_list(type_map)

    """
    make lammps input for elastic calculation
    """
    ret = ""
    ret += "clear\n"
    ret += "units 	metal\n"
    ret += "dimension	3\n"
    ret += "boundary	p p p\n"
    ret += "atom_style	atomic\n"
    ret += "box         tilt large\n"
    ret += "read_data   %s\n" % conf
    for ii in range(len(type_map)):
        ret += "mass            %d %.3f\n" % (ii + 1, Element(type_map_list[ii]).mass)
    ret += "neigh_modify    every 1 delay 0 check no\n"
    ret += interaction(param)
    ret += "compute         mype all pe\n"
    ret += "thermo          100\n"
    ret += "thermo_style    custom step pe pxx pyy pzz pxy pxz pyz lx ly lz vol c_mype\n"
    ret += "dump            1 all custom 100 dump.relax id type xs ys zs fx fy fz\n"
    ret += "min_style       cg\n"
    ret += "minimize        %e %e %d %d\n" % (etol, ftol, maxiter, maxeval)
    ret += "variable        N equal count(all)\n"
    ret += "variable        V equal vol\n"
    ret += "variable        E equal \"c_mype\"\n"
    ret += "variable        Pxx equal pxx\n"
    ret += "variable        Pyy equal pyy\n"
    ret += "variable        Pzz equal pzz\n"
    ret += "variable        Pxy equal pxy\n"
    ret += "variable        Pxz equal pxz\n"
    ret += "variable        Pyz equal pyz\n"
    ret += "variable        Epa equal ${E}/${N}\n"
    ret += "variable        Vpa equal ${V}/${N}\n"
    ret += "print \"All done\"\n"
    ret += "print \"Total number of atoms = ${N}\"\n"
    ret += "print \"Final energy per atoms = ${Epa}\"\n"
    ret += "print \"Final volume per atoms = ${Vpa}\"\n"
    ret += "print \"Final Stress (xx yy zz xy xz yz) = ${Pxx} ${Pyy} ${Pzz} ${Pxy} ${Pxz} ${Pyz}\"\n"
    return ret


def make_lammps_press_relax(conf, type_map, scale2equi, interaction, param,
                            B0=70, bp=0,
                            etol=0, ftol=1e-10,
                            maxiter=5000, maxeval=500000):
    type_map_list = element_list(type_map)

    """
    make lammps input for relaxation at a certain volume
    scale2equi: the volume scale with respect to equilibrium volume
    """
    ret = ""
    ret += "clear\n"
    ret += "variable        GPa2bar	equal 1e4\n"
    ret += "variable        B0		equal %f\n" % B0
    ret += "variable        bp		equal %f\n" % bp
    ret += "variable	    xx		equal %f\n" % scale2equi
    ret += "variable        yeta	equal 1.5*(${bp}-1)\n"
    ret += "variable        Px0		equal 3*${B0}*(1-${xx})/${xx}^2*exp(${yeta}*(1-${xx}))\n"
    ret += "variable        Px		equal ${Px0}*${GPa2bar}\n"
    ret += "units       metal\n"
    ret += "dimension   3\n"
    ret += "boundary	p p p\n"
    ret += "atom_style	atomic\n"
    ret += "box         tilt large\n"
    ret += "read_data   %s\n" % conf
    for ii in range(len(type_map)):
        ret += "mass            %d %.3f\n" % (ii + 1, Element(type_map_list[ii]).mass)
    ret += "neigh_modify    every 1 delay 0 check no\n"
    ret += interaction(param)
    ret += "compute         mype all pe\n"
    ret += "thermo          100\n"
    ret += "thermo_style    custom step pe pxx pyy pzz pxy pxz pyz lx ly lz vol c_mype\n"
    ret += "dump            1 all custom 100 dump.relax id type xs ys zs fx fy fz\n"
    ret += "min_style       cg\n"
    ret += "fix             1 all box/relax iso ${Px} \n"
    ret += "minimize        %e %e %d %d\n" % (etol, ftol, maxiter, maxeval)
    ret += "fix             1 all box/relax aniso ${Px} \n"
    ret += "minimize        %e %e %d %d\n" % (etol, ftol, maxiter, maxeval)
    ret += "variable        N equal count(all)\n"
    ret += "variable        V equal vol\n"
    ret += "variable        E equal \"c_mype\"\n"
    ret += "variable        Pxx equal pxx\n"
    ret += "variable        Pyy equal pyy\n"
    ret += "variable        Pzz equal pzz\n"
    ret += "variable        Pxy equal pxy\n"
    ret += "variable        Pxz equal pxz\n"
    ret += "variable        Pyz equal pyz\n"
    ret += "variable        Epa equal ${E}/${N}\n"
    ret += "variable        Vpa equal ${V}/${N}\n"
    ret += "print \"All done\"\n"
    ret += "print \"Total number of atoms  = ${N}\"\n"
    ret += "print \"Relax at Press         = ${Px} Bar\"\n"
    ret += "print \"Final energy per atoms = ${Epa} eV\"\n"
    ret += "print \"Final volume per atoms = ${Vpa} A^3\"\n"
    ret += "print \"Final Stress (xx yy zz xy xz yz) = ${Pxx} ${Pyy} ${Pzz} ${Pxy} ${Pxz} ${Pyz}\"\n"
    return ret


def make_lammps_phonon(conf, masses, interaction, param,
                       etol=0, ftol=1e-10,
                       maxiter=5000, maxeval=500000):
    """
    make lammps input for elastic calculation
    """
    ret = ""
    ret += "clear\n"
    ret += "units 	metal\n"
    ret += "dimension	3\n"
    ret += "boundary	p p p\n"
    ret += "atom_style	atomic\n"
    ret += "box         tilt large\n"
    ret += "read_data   %s\n" % conf
    ntypes = len(masses)
    for ii in range(ntypes):
        ret += "mass            %d %f\n" % (ii + 1, masses[ii])
    ret += "neigh_modify    every 1 delay 0 check no\n"
    ret += interaction(param)
    return ret


def _get_epa(lines):
    for ii in lines:
        if ("Final energy per atoms" in ii) and (not 'print' in ii):
            return float(ii.split('=')[1].split()[0])
    raise RuntimeError("cannot find key \"Final energy per atoms\" in lines, something wrong")


def _get_vpa(lines):
    for ii in lines:
        if ("Final volume per atoms" in ii) and (not 'print' in ii):
            return float(ii.split('=')[1].split()[0])
    raise RuntimeError("cannot find key \"Final volume per atoms\" in lines, something wrong")


def _get_natoms(lines):
    for ii in lines:
        if ("Total number of atoms" in ii) and (not 'print' in ii):
            return int(ii.split('=')[1].split()[0])
    raise RuntimeError("cannot find key \"Total number of atoms\" in lines, something wrong")


def get_nev(log):
    """
    get natoms, energy_per_atom and volume_per_atom from lammps log
    """
    with open(log, 'r') as fp:
        lines = fp.read().split('\n')
    epa = _get_epa(lines)
    vpa = _get_vpa(lines)
    natoms = _get_natoms(lines)
    return natoms, epa, vpa


def get_base_area(log):
    """
    get base area
    """
    with open(log, 'r') as fp:
        lines = fp.read().split('\n')
    for ii in lines:
        if ("Final Base area" in ii) and (not 'print' in ii):
            return float(ii.split('=')[1].split()[0])


def get_stress(log):
    """
    get stress from lammps log
    """
    with open(log, 'r') as fp:
        lines = fp.read().split('\n')
    for ii in lines:
        if ('Final Stress' in ii) and (not 'print' in ii):
            vstress = [float(jj) for jj in ii.split('=')[1].split()]
    stress = util.voigt_to_stress(vstress)
    return stress


def poscar_from_last_dump(dump, poscar_out, deepmd_type_map):
    """
    get poscar from the last frame of a lammps MD traj (dump format)
    """
    with open(dump, 'r') as fp:
        lines = fp.read().split('\n')
    step_idx = -1
    for idx, ii in enumerate(lines):
        if 'ITEM: TIMESTEP' in ii:
            step_idx = idx
    if step_idx == -1:
        raise RuntimeError("cannot find timestep in lammps dump, something wrong")
    with open('tmp_dump', 'w') as fp:
        fp.write("\n".join(lines[step_idx:]))
    cvt_lammps_conf('tmp_dump', poscar_out, ofmt='vasp')
    os.remove('tmp_dump')
    with open(poscar_out, 'r') as fp:
        lines = fp.read().split('\n')
    types = [deepmd_type_map[int(ii.split('_')[1])] for ii in lines[5].split()]
    lines[5] = " ".join(types)
    with open(poscar_out, 'w') as fp:
        lines = fp.write("\n".join(lines))


def check_finished_new(fname, keyword):
    with open(fname, 'r') as fp:
        lines = fp.read().split('\n')
    flag = False
    for jj in lines:
        if (keyword in jj) and (not 'print' in jj):
            flag = True
    return flag


def check_finished(fname):
    with open(fname, 'r') as fp:
        return 'Total wall time:' in fp.read()
