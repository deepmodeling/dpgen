#!/usr/bin/env python3

import os, re, argparse
import dpgen.auto_test.lib.crys as crys
from pymatgen.ext.matproj import MPRester, Composition
from pymatgen.analysis.structure_matcher import StructureMatcher

global_std_crystal = {
    'fcc' : crys.fcc,
    'hcp' : crys.hcp,
    'dhcp' : crys.dhcp,
    'bcc' : crys.bcc,
    'diamond' : crys.diamond,
    'sc' : crys.sc
}

def test_fit(struct, data) :
    m = StructureMatcher() 
    for ii in data :
        if m.fit(ii['structure'], struct) :
            return True
    return False

def make_path_mp(ii) :
    pf = ii['pretty_formula']
    pf = re.sub('\d+', '', pf)
    task_id = ii['task_id']
    work_path = 'confs'
    work_path = os.path.join(work_path, pf)
    work_path = os.path.join(work_path, task_id)
    return work_path

def gen_ele_std(ele_name, ctype):
    struct = global_std_crystal[ctype](ele_name)
    work_path = 'confs'
    work_path = os.path.join(work_path, ele_name)
    work_path = os.path.join(work_path, 'std-'+ctype)
    os.makedirs(work_path, exist_ok = True)
    fposcar = os.path.join(work_path, 'POSCAR')
    fjson = os.path.join(work_path, 'data.json')
    struct.to('poscar', fposcar)    
    return struct

def gen_element(ele_name,key) :
    assert(type(ele_name) == str)
    mpr = MPRester(key)
    data = mpr.query({'elements':[ele_name], 'nelements':1}, 
                     properties=["task_id", 
                                 "pretty_formula", 
                                 'formula', 
                                 "anonymous_formula",                             
                                 'formation_energy_per_atom',
                                 'energy_per_atom',
                                 'structure'])
    for ii in data :
        work_path = make_path_mp(ii)
        os.makedirs(work_path, exist_ok = True)
        fposcar = os.path.join(work_path, 'POSCAR')
        fjson = os.path.join(work_path, 'data.json')
        ii['structure'].to('poscar', fposcar)
        ii['structure'].to('json', fjson)

    m = StructureMatcher() 
    for ii in global_std_crystal.keys() :
        ss = gen_ele_std(ele_name, ii)
        find = False
        for jj in data:
            if m.fit(ss,jj['structure']) :
                find = True
                break
        if find :            
            work_path = make_path_mp(jj)
            with open(os.path.join(work_path,'std-crys'), 'w') as fp :
                fp.write(ii+'\n')

def gen_element_std(ele_name) :
    assert(type(ele_name) == str)
    for ii in global_std_crystal.keys() :
        ss = gen_ele_std(ele_name, ii)

def gen_alloy(eles,key) :
    
    mpr = MPRester(key)

    data = mpr.query({'elements':{'$all': eles}, 'nelements':len(eles)}, 
                     properties=["task_id", 
                                 "pretty_formula", 
                                 'formula', 
                                 "anonymous_formula",                             
                                 'formation_energy_per_atom',
                                 'energy_per_atom',
                                 'structure'])
    if len(data) == 0 :
        return
    
    alloy_file = make_path_mp(data[0])
    os.makedirs(alloy_file, exist_ok = True)
    alloy_file = os.path.join(alloy_file, '..')
    alloy_file = os.path.join(alloy_file, 'alloy')
    with open(alloy_file, 'w') as fp :
        None    
    
    for ii in data :
        work_path = make_path_mp(ii)
        os.makedirs(work_path, exist_ok = True)
        fposcar = os.path.join(work_path, 'POSCAR')
        fjson = os.path.join(work_path, 'data.json')
        ii['structure'].to('poscar', fposcar)
        ii['structure'].to('json', fjson)

def _main() :
    parser = argparse.ArgumentParser(
        description="gen structures")
    parser.add_argument('key', type=str,
                        help='key id of material project')
    parser.add_argument('elements', 
                        type=str,
                        nargs = '+',
                        help="the list of appeared elements")
    args = parser.parse_args()

    print('generate %s' % (args.elements))
    if len(args.elements) == 1 :
        gen_element(args.elements[0],key)
        # gen_element_std(args.elements[0])
    else :
        gen_alloy(args.elements,key)

if __name__ == '__main__' :
    _main()

