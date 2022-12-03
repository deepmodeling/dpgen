#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import os,sys,glob,time
from deepmd.calculator import DP
from ase.io import read 

'''
check if structure optimization worked well
if not, this script will generate a fake outcar
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
    f.write('ASE-DPKit-FAILED-nan\n')
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
    f.write('energy  without entropy= %20.6f %20.6f\n' % (ene, ene)) 
    enthalpy = ene + pstress * volume / 1602.17733 
    f.write('enthalpy is  TOTEN    = %20.6f %20.6f\n' % (enthalpy, enthalpy))

def check():

    from deepmd.calculator import DP
    from ase.io import read
    calc = DP(model='../graph.000.pb')    # init the model before iteration

    to_be_opti = read('POSCAR')
    to_be_opti.calc = calc     
    # ---------------------------------
    # for failed outcar 
    atoms_symbols_f = to_be_opti.get_chemical_symbols()
    element_f, ele_f = Get_Element_Num(atoms_symbols_f)
    atoms_vol_f = to_be_opti.get_volume() 
    atoms_stress_f = to_be_opti.get_stress()
    atoms_stress_f = atoms_stress_f/(0.01*0.6242)
    atoms_lat_f = to_be_opti.cell 
    atoms_pos_f = to_be_opti.positions
    atoms_force_f = to_be_opti.get_forces()
    atoms_ene_f =  610612509 
    # --------------------------------- 
    Write_Contcar(element_f, ele_f, atoms_lat_f, atoms_pos_f)
    Write_Outcar(element_f, ele_f, atoms_vol_f, atoms_lat_f, atoms_pos_f,atoms_ene_f, atoms_force_f, atoms_stress_f * -10.0, 0) 
 
cwd = os.getcwd()
if not os.path.exists(os.path.join(cwd,'OUTCAR')):
    check()
