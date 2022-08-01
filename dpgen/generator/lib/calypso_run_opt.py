#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os,sys,glob,time 
import numpy as np
from ase.io import read 
from ase.optimize import BFGS,QuasiNewton,LBFGS
from ase.constraints import UnitCellFilter, ExpCellFilter
from deepmd.calculator import DP 
'''
structure optimization with DP model and ASE
PSTRESS and fmax should exist in input.dat
'''

def Get_Element_Num(elements):
    '''Using the Atoms.symples to Know Element&Num'''
    element = []
    ele = {}
    element.append(elements[0])
    for x in elements: 
        if x not in element :
            element.append(x)
    for x in element: 
        ele[x] = elements.count(x)
    return element, ele 
        
def Write_Contcar(element, ele, lat, pos):
    '''Write CONTCAR''' 
    f = open('CONTCAR','w')
    f.write('ASE-DPKit-Optimization\n')
    f.write('1.0\n') 
    for i in range(3):
        f.write('%15.10f %15.10f %15.10f\n' % tuple(lat[i]))
    for x in element: 
        f.write(x + '  ')
    f.write('\n') 
    for x in element:
        f.write(str(ele[x]) + '  ') 
    f.write('\n') 
    f.write('Direct\n')
    na = sum(ele.values())
    dpos = np.dot(pos,np.linalg.inv(lat))
    for i in range(na): 
        f.write('%15.10f %15.10f %15.10f\n' % tuple(dpos[i]))
        
def Write_Outcar(element, ele, volume, lat, pos, ene, force, stress,pstress):
    '''Write OUTCAR'''
    f = open('OUTCAR','w')
    for x in element: 
        f.write('VRHFIN =' + str(x) + '\n')
    f.write('ions per type =')
    for x in element:
        f.write('%5d' % ele[x])
    f.write('\nDirection     XX             YY             ZZ             XY             YZ             ZX\n') 
    f.write('in kB') 
    f.write('%15.6f' % stress[0])
    f.write('%15.6f' % stress[1])
    f.write('%15.6f' % stress[2])
    f.write('%15.6f' % stress[3])
    f.write('%15.6f' % stress[4])
    f.write('%15.6f' % stress[5])
    f.write('\n') 
    ext_pressure = np.sum(stress[0] + stress[1] + stress[2])/3.0 - pstress
    f.write('external pressure = %20.6f kB    Pullay stress = %20.6f  kB\n'% (ext_pressure, pstress))
    f.write('volume of cell : %20.6f\n' % volume)
    f.write('direct lattice vectors\n')
    for i in range(3):
        f.write('%10.6f %10.6f %10.6f\n' % tuple(lat[i]))
    f.write('POSITION                                       TOTAL-FORCE(eV/Angst)\n')  
    f.write('-------------------------------------------------------------------\n')
    na = sum(ele.values()) 
    for i in range(na): 
        f.write('%15.6f %15.6f %15.6f' % tuple(pos[i])) 
        f.write('%15.6f %15.6f %15.6f\n' % tuple(force[i]))
    f.write('-------------------------------------------------------------------\n')
    f.write('energy  without entropy= %20.6f %20.6f\n' % (ene, ene/na))
    enthalpy = ene + pstress * volume / 1602.17733      
    f.write('enthalpy is  TOTEN    = %20.6f %20.6f\n' % (enthalpy, enthalpy/na)) 
        
def read_stress_fmax():
    pstress = 0
    fmax = 0.01
    #assert os.path.exists('./input.dat'), 'input.dat does not exist!'
    try:
        f = open('input.dat','r')
    except:
        assert os.path.exists('../input.dat'),' now we are in %s, do not find ../input.dat'%(os.getcwd())
        f = open('../input.dat','r')
    lines = f.readlines()
    f.close()
    for line in lines:
        if line[0] == '#':
            continue
        if 'PSTRESS' in line or 'pstress' in line:
            pstress = float(line.split('=')[1])
        if 'fmax' in line:
            fmax = float(line.split('=')[1])
    return fmax,pstress
        
def run_opt(fmax,stress): 
    '''Using the ASE&DP to Optimize Configures''' 
    
    calc = DP(model='../graph.000.pb')    # init the model before iteration
    os.system('mv OUTCAR OUTCAR-last')

    print('Start to Optimize Structures by DP----------')
        
    Opt_Step = 1000 
    start = time.time() 
    # pstress kbar
    pstress = stress
    # kBar to eV/A^3
    # 1 eV/A^3 = 160.21766028 GPa
    # 1 / 160.21766028 ~ 0.006242
    aim_stress = 1.0 * pstress* 0.01 * 0.6242 / 10.0 
    to_be_opti = read('POSCAR') 
    to_be_opti.calc = calc 
    ucf = UnitCellFilter(to_be_opti, scalar_pressure=aim_stress)
    atoms_vol_2 = to_be_opti.get_volume() 
    # opt
    opt = LBFGS(ucf,trajectory='traj.traj') 
    opt.run(fmax=fmax,steps=Opt_Step) 
    
    atoms_lat = to_be_opti.cell 
    atoms_pos = to_be_opti.positions
    atoms_force = to_be_opti.get_forces() 
    atoms_stress = to_be_opti.get_stress() 
    # eV/A^3 to GPa
    atoms_stress = atoms_stress/(0.01*0.6242)
    atoms_symbols = to_be_opti.get_chemical_symbols() 
    atoms_ene = to_be_opti.get_potential_energy() 
    atoms_vol = to_be_opti.get_volume() 
    element, ele = Get_Element_Num(atoms_symbols) 

    Write_Contcar(element, ele, atoms_lat, atoms_pos)
    Write_Outcar(element, ele, atoms_vol, atoms_lat, atoms_pos,atoms_ene, atoms_force, atoms_stress * -10.0, pstress)

    stop = time.time()
    _cwd = os.getcwd()
    _cwd = os.path.basename(_cwd)
    print('%s is done, time: %s' % (_cwd,stop-start))

def run():
    fmax, stress = read_stress_fmax()
    run_opt(fmax, stress)

if __name__=='__main__':
    run()
