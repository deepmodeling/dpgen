#!/usr/bin/env python3

import random, os, sys
import subprocess as sp
import lib.util as util

def cvt_lammps_conf (fin, 
                     fout, 
                     ofmt = 'lammps_data') :
    """
    Format convert from fin to fout, specify the output format by ofmt
    """
    thisfile = os.path.abspath(__file__)
    thisdir = os.path.dirname(thisfile)
    cmd = os.path.join(thisdir, 'ovito_file_convert.py')
    cmd_opt = '-m '+ofmt
    cmd_line = cmd + ' ' + cmd_opt + ' ' + fin + ' ' + fout
    sp.check_call(cmd_line, shell = True)    
    # sp.check_call([cmd, cmd_opt, fin, fout])

def apply_type_map(conf_file, deepmd_type_map, ptypes) :
    """
    apply type map.
    conf_file:          conf file converted from POSCAR
    deepmd_type_map:    deepmd atom type map
    ptypes:             atom types defined in POSCAR
    """
    natoms = _get_natom(conf_file)
    ntypes = len(deepmd_type_map)
    with open(conf_file, 'r') as fp:
        lines = fp.read().split('\n')
    # with open(conf_file+'.bk', 'w') as fp:
    #     fp.write("\n".join(lines))
    new_lines = lines
    # revise ntypes
    idx_ntypes = -1
    for idx, ii in enumerate(lines) :
        if 'atom types' in ii :
            idx_ntypes = idx
    if idx_ntypes == -1 :
        raise RuntimeError("cannot find the entry 'atom types' in ", conf_file)
    words = lines[idx_ntypes].split()
    words[0] = str(ntypes)
    new_lines[idx_ntypes] = " ".join(words)    
    # find number of atoms
    idx_atom_entry = -1
    for idx, ii in enumerate(lines) :
        if 'Atoms' in ii :
            idx_atom_entry = idx
    if idx_atom_entry == -1 :
        raise RuntimeError("cannot find the entry 'Atoms' in ", conf_file)
    # revise atom type
    for idx in range(idx_atom_entry+2, idx_atom_entry+2+natoms) :
        ii = lines[idx]
        words = ii.split()
        assert(len(words) >= 5)
        old_id = int(words[1])
        new_id = deepmd_type_map.index(ptypes[old_id-1])+1
        words[1] = str(new_id)
        ii = " ".join(words)
        new_lines[idx] = ii
    with open(conf_file, 'w') as fp:
        fp.write("\n".join(new_lines))        

def _get_ntype(conf) :
    with open(conf, 'r') as fp:
        lines = fp.read().split('\n')
    for ii in lines :
        if "atom types" in ii :
            return int(ii.split()[0])
    raise RuntimeError("cannot find line indicate atom types in ", conf)

def _get_natom(conf) :
    with open(conf, 'r') as fp:
        lines = fp.read().split('\n')
    for ii in lines :
        if "atoms" in ii :
            return int(ii.split()[0])
    raise RuntimeError("cannot find line indicate atom types in ", conf)

def _make_lammps_deepmd_model(models) :
    ret = ""
    line = "pair_style deepmd "
    for ii in models :
        line += ii + ' '
    line += ' 10 model_devi.out\npair_coeff\n'
    ret += line
    return ret

def make_lammps_equi(conf, ntypes, models, 
                     etol=1e-12, ftol=1e-6, 
                     maxiter=5000, maxeval=500000, 
                     change_box = True) :
    """
    make lammps input for equilibritation
    """
    ret = ""
    ret += "clear\n"
    ret += "units 	metal\n"
    ret += "dimension	3\n"
    ret += "boundary	p	p    p\n"
    ret += "atom_style	atomic\n"
    ret += "box         tilt large\n"
    ret += "read_data   %s\n" % conf
    for ii in range(ntypes) :
        ret += "mass            %d 1\n" % (ii+1)            
    ret += "neigh_modify    every 1 delay 0 check no\n"
    ret += _make_lammps_deepmd_model(models)
    ret += "compute         mype all pe\n"
    ret += "thermo          100\n"
    ret += "thermo_style    custom step pe pxx pyy pzz pxy pxz pyz lx ly lz vol c_mype\n"
    ret += "dump            1 all custom 100 dump.relax id type xs ys zs fx fy fz\n"
    ret += "min_style       cg\n"
    if change_box :
        ret += "fix             1 all box/relax iso 0.0 \n"
        ret += "minimize        %e %e %d %d\n" % (etol, ftol, maxiter, maxeval)
        ret += "fix             1 all box/relax aniso 0.0 \n"
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

def make_lammps_elastic(conf, ntypes, models, 
                        etol=1e-12, ftol=1e-6, 
                        maxiter=5000, maxeval=500000) :
    """
    make lammps input for elastic calculation
    """
    ret = ""
    ret += "clear\n"
    ret += "units 	metal\n"
    ret += "dimension	3\n"
    ret += "boundary	p	p    p\n"
    ret += "atom_style	atomic\n"
    ret += "box         tilt large\n"
    ret += "read_data   %s\n" % conf
    for ii in range(ntypes) :
        ret += "mass            %d 1\n" % (ii+1)            
    ret += "neigh_modify    every 1 delay 0 check no\n"
    ret += _make_lammps_deepmd_model(models)
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

def make_lammps_press_relax(conf, ntypes, scale2equi, models,
                            B0 = 70, bp = 0, 
                            etol=1e-12, ftol=1e-6, 
                            maxiter=5000, maxeval=500000) :
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
    for ii in range(ntypes) :
        ret += "mass            %d 1\n" % (ii+1)            
    ret += "neigh_modify    every 1 delay 0 check no\n"
    ret += _make_lammps_deepmd_model(models)
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

def _get_epa (lines) :
    for ii in lines:
        if ("Final energy per atoms" in ii) and (not 'print' in ii):
            return float(ii.split('=')[1].split()[0])
    raise RuntimeError("cannot find key \"Final energy per atoms\" in lines, something wrong")

def _get_vpa (lines) :
    for ii in lines :
        if ("Final volume per atoms" in ii) and (not 'print' in ii):
            return float(ii.split('=')[1].split()[0])
    raise RuntimeError("cannot find key \"Final volume per atoms\" in lines, something wrong")

def _get_natoms (lines) :
    for ii in lines:
        if ("Total number of atoms" in ii) and (not 'print' in ii):
            return int(ii.split('=')[1].split()[0])
    raise RuntimeError("cannot find key \"Total number of atoms\" in lines, something wrong")

def get_nev (log) :
    """
    get natoms, energy_per_atom and volume_per_atom from lammps log
    """
    with open(log, 'r') as fp:
        lines = fp.read().split('\n') 
    epa = _get_epa(lines)
    vpa = _get_vpa(lines)
    natoms = _get_natoms(lines)
    return natoms, epa, vpa

def get_base_area (log) :
    """
    get base area
    """
    with open(log, 'r') as fp:
        lines = fp.read().split('\n') 
    for ii in lines:
        if ("Final Base area" in ii) and (not 'print' in ii):
            return float(ii.split('=')[1].split()[0])

def get_stress(log) :
    """
    get stress from lammps log
    """
    with open(log, 'r') as fp :
        lines = fp.read().split('\n')
    for ii in lines :
        if ('Final Stress' in ii) and (not 'print' in ii):
            vstress = [float(jj) for jj in ii.split('=')[1].split()]
    stress = util.voigt_to_stress(vstress)
    return stress

def poscar_from_last_dump(dump, poscar_out, deepmd_type_map) :
    """
    get poscar from the last frame of a lammps MD traj (dump format)
    """
    with open(dump, 'r') as fp :
        lines = fp.read().split('\n')
    step_idx = -1
    for idx,ii in enumerate(lines) :
        if 'ITEM: TIMESTEP' in ii :
            step_idx = idx
    if step_idx == -1 :
        raise RuntimeError("cannot find timestep in lammps dump, something wrong")    
    with open('tmp_dump', 'w') as fp:
        fp.write("\n".join(lines[step_idx:]))
    cvt_lammps_conf('tmp_dump', poscar_out, ofmt='vasp')
    os.remove('tmp_dump')
    with open(poscar_out, 'r') as fp:
        lines = fp.read().split('\n')
    types = [ deepmd_type_map[int(ii.split('_')[1])-1] for ii in lines[5].split()]
    lines[5] = " ".join(types)
    with open(poscar_out, 'w') as fp:
        lines = fp.write("\n".join(lines))



    
 
