#!/usr/bin/python3
import os,sys
from unicodedata import numeric
import dpdata
import dpgen.generator.lib.abacus_scf as abacus_scf
import numpy as np
from pymatgen.core.structure import Structure

A2BOHR = 1.8897261254578281
MASS_DICT = {"H":1.0079,"He":4.0026,"Li":6.941,"Be":9.0122,"B":10.811,"C":12.0107,"N":14.0067,\
             "O":15.9994,"F":18.9984,"Ne":20.1797,"Na":22.9897,"Mg":24.305,"Al":26.9815,"Si":28.0855,\
             "P":30.9738,"S":32.065,"Cl":35.453,"K":39.0983,"Ar":39.948,"Ca":40.078,"Sc":44.9559,\
             "Ti":47.867,"V":50.9415,"Cr":51.9961,"Mn":54.938,"Fe":55.845,"Ni":58.6934,"Co":58.9332,\
             "Cu":63.546,"Zn":65.39,"Ga":69.723,"Ge":72.64,"As":74.9216,"Se":78.96,"Br":79.904,"Kr":83.8,\
             "Rb":85.4678,"Sr":87.62,"Y":88.9059,"Zr":91.224,"Nb":92.9064,"Mo":95.94,"Tc":98,"Ru":101.07,\
             "Rh":102.9055,"Pd":106.42,"Ag":107.8682,"Cd":112.411,"In":114.818,"Sn":118.71,"Sb":121.76,"I":126.9045,\
             "Te":127.6,"Xe":131.293,"Cs":132.9055,"Ba":137.327,"La":138.9055,"Ce":140.116,"Pr":140.9077,"Nd":144.24,\
             "Pm":145,"Sm":150.36,"Eu":151.964,"Gd":157.25,"Tb":158.9253,"Dy":162.5,"Ho":164.9303,"Er":167.259,\
             "Tm":168.9342,"Yb":173.04,"Lu":174.967,"Hf":178.49,"Ta":180.9479,"W":183.84,"Re":186.207,"Os":190.23,\
             "Ir":192.217,"Pt":195.078,"Au":196.9665,"Hg":200.59,"Tl":204.3833,"Pb":207.2,"Bi":208.9804,"Po":209,\
             "At":210,"Rn":222,"Fr":223,"Ra":226,"Ac":227,"Pa":231.0359,"Th":232.0381,"Np":237,"U":238.0289,"Am":243,\
             "Pu":244,"Cm":247,"Bk":247,"Cf":251,"Es":252,"Fm":257,"Md":258,"No":259,"Rf":261,"Lr":262,"Db":262,"Bh":264,\
             "Sg":266,"Mt":268,"Rg":272,"Hs":277,"H":1.0079,"He":4.0026,"Li":6.941,"Be":9.0122,"B":10.811,"C":12.0107,\
             "N":14.0067,"O":15.9994,"F":18.9984,"Ne":20.1797,"Na":22.9897,"Mg":24.305,"Al":26.9815,"Si":28.0855,"P":30.9738,\
             "S":32.065,"Cl":35.453,"K":39.0983,"Ar":39.948,"Ca":40.078,"Sc":44.9559,"Ti":47.867,"V":50.9415,"Cr":51.9961,\
             "Mn":54.938,"Fe":55.845,"Ni":58.6934,"Co":58.9332,"Cu":63.546,"Zn":65.39,"Ga":69.723,"Ge":72.64,"As":74.9216,\
             "Se":78.96,"Br":79.904,"Kr":83.8,"Rb":85.4678,"Sr":87.62,"Y":88.9059,"Zr":91.224,"Nb":92.9064,"Mo":95.94,"Tc":98,\
             "Ru":101.07,"Rh":102.9055,"Pd":106.42,"Ag":107.8682,"Cd":112.411,"In":114.818,"Sn":118.71,"Sb":121.76,\
             "I":126.9045,"Te":127.6,"Xe":131.293,"Cs":132.9055,"Ba":137.327,"La":138.9055,"Ce":140.116,"Pr":140.9077,\
             "Nd":144.24,"Pm":145,"Sm":150.36,"Eu":151.964,"Gd":157.25,"Tb":158.9253,"Dy":162.5,"Ho":164.9303,"Er":167.259,\
             "Tm":168.9342,"Yb":173.04,"Lu":174.967,"Hf":178.49,"Ta":180.9479,"W":183.84,"Re":186.207,"Os":190.23,"Ir":192.217,\
             "Pt":195.078,"Au":196.9665,"Hg":200.59,"Tl":204.3833,"Pb":207.2,"Bi":208.9804,"Po":209,"At":210,"Rn":222,"Fr":223,\
             "Ra":226,"Ac":227,"Pa":231.0359,"Th":232.0381,"Np":237,"U":238.0289,"Am":243,"Pu":244,"Cm":247,"Bk":247,"Cf":251,\
             "Es":252,"Fm":257,"Md":258,"No":259,"Rf":261,"Lr":262,"Db":262,"Bh":264,"Sg":266,"Mt":268,"Rg":272,"Hs":277}
key_words_list = ["ATOMIC_SPECIES", "NUMERICAL_ORBITAL", "LATTICE_CONSTANT", "LATTICE_VECTORS", "ATOMIC_POSITIONS", "NUMERICAL_DESCRIPTOR"]

def poscar2stru(poscar,inter_param,stru):
    '''
        - poscar:           POSCAR for input
        - inter_param:      dictionary of 'interaction' from param.json
                            some key words for ABACUS are:
                                - atom_masses:  a dictionary of atoms' masses
                                - orb_files:     a dictionary of orbital files
                                - deepks_desc:  a string of deepks descriptor file
        - stru:            output filename, usally is 'STRU'
    '''
    stru = dpdata.System(poscar, fmt = 'vasp/poscar')
    stru_data = stru.data
    atom_mass = []
    pseudo = None
    orb = None
    deepks_desc = None 

    if 'atom_masses' not in inter_param:
        atom_mass_dict = {i:1.0 if i not in MASS_DICT else MASS_DICT[i] for i in stru_data['atom_names']}
    else:
        atom_mass_dict = inter_param['atom_masses']
    for atom in stru_data['atom_names']:
        assert(atom in atom_mass_dict), "the mass of %s is not defined in interaction:atom_masses" % atom
        atom_mass.append(atom_mass_dict[atom]) 

    if 'potcars' in inter_param:
        pseudo = []
        for atom in stru_data['atom_names']:
            assert(atom in inter_param['potcars']), "the pseudopotential of %s is not defined in interaction:potcars" % atom
            pseudo.append("./pp_orb/" + inter_param['potcars'][atom].split('/')[-1])

    if 'orb_files' in inter_param:
        orb = []
        for atom in stru_data['atom_names']:
            assert(atom in inter_param['orb_files']), "orbital file of %s is not defined in interaction:orb_files" % atom
            orb.append("./pp_orb/" + inter_param['orb_files'][atom].split('/')[-1])

    if 'deepks_desc' in inter_param:
        deepks_desc ="./pp_orb/%s\n" % inter_param['deepks_desc']

    stru.to("stru", "STRU", mass = atom_mass, pp_file = pseudo, numerical_orbital = orb, numerical_descriptor = deepks_desc)


def stru_fix_atom(struf,fix_atom = [True,True,True]):
    '''
...
ATOMIC_POSITIONS
Cartesian               #Cartesian(Unit is LATTICE_CONSTANT)
Si                      #Name of element
0.0                     #Magnetic for this element.
2                       #Number of atoms
0.00 0.00 0.00 0 0 0    #x,y,z, move_x, move_y, move_z
0.25 0.25 0.25 0 0 0
    '''
    fix_xyz = ['0' if i else '1' for i in fix_atom ]
    if os.path.isfile(struf):
        with open(struf) as f1: lines = f1.readlines()
        for i in range(len(lines)):
            if "ATOMIC_POSITIONS" in lines[i]: break
        i += 1
        flag_read_coord_type = False
        flag_read_atom_number = 2
        flag_atom_number = 0
        while i < len(lines):
            if lines[i].strip() == '':pass 
            elif lines[i].split()[0] in key_words_list: break
            elif not flag_read_coord_type: 
                flag_read_coord_type = True
            elif flag_atom_number:
                flag_atom_number -= 1
                x,y,z = lines[i].split()[:3]
                lines[i] = "%s %s %s %s %s %s\n" % tuple([x,y,z] + fix_xyz)
            elif flag_read_coord_type and flag_read_atom_number:
                flag_read_atom_number -= 1
            elif not flag_read_atom_number:
                flag_read_atom_number = 2
                flag_atom_number = int(lines[i].split()[0])
            i += 1
        with open(struf,'w') as f1: f1.writelines(lines)
    else:
        raise RuntimeError("Error: Try to modify struc file %s, but can not find it" % struf)

def stru_scale (stru_in, stru_out, scale) :
    with open(stru_in, 'r') as fin : lines = fin.readlines()
    for i in range(len(lines)):
        if "LATTICE_CONSTANT" in lines[i]:
            lines[i+1] = str(float(lines[i+1].strip())*scale) + '\n'
            break
    with open(stru_out,'w') as f1: f1.writelines(lines)


def write_kpt(kptf,kptlist):
    context = "K_POINTS\n0\nGamma\n"
    for i in kptlist: context += str(i) + " "
    with open(kptf,'w') as f1:  f1.write(context)     

def write_input(inputf,inputdict):
    context = "INPUT_PARAMETERS\n"
    for key in inputdict.keys():
        if key[0] in ['_','#']:continue
        context += key + " " + str(inputdict[key]) + "\n"
    with open(inputf,'w') as f1:  f1.write(context)       

def make_kspacing_kpt(struf,kspacing):
    stru_data = abacus_scf.get_abacus_STRU(struf)
    cell = stru_data['cells'] / abacus_scf.bohr2ang
    volume = abs(cell[0].dot(np.cross(cell[1],cell[2])))
    coef = 2 * np.pi / volume / kspacing
    kpt = [max(1,int(np.linalg.norm(np.cross(cell[x],cell[y]))*coef+1)) for x,y in [[1,2],[2,0],[0,1]] ]
    return kpt

def check_finished(fname):
    with open(fname, 'r') as fp:
        return 'Total  Time  :' in fp.read()

def final_stru(abacus_path):
    with open(os.path.join(abacus_path, 'INPUT')) as f1: lines = f1.readlines()
    suffix = 'ABACUS'
    calculation = 'scf'
    out_stru = False
    for line in lines:
        if 'suffix' in line and line.split()[0] == 'suffix': 
            suffix = line.split()[1]
        elif 'calculation' in line and line.split()[0] == 'calculation':
            calculation = line.split()[1]
        elif 'out_stru' in line and line.split()[0] == 'out_stru':
            out_stru = bool(line.split()[1])
    logf = os.path.join(abacus_path, 'OUT.%s/running_%s.log'%(suffix,calculation))
    if calculation in ['relax','cell-relax']:
        if not out_stru:
            return 'OUT.%s/STRU_ION_D' % suffix
        else:
            with open(logf) as f1: lines = f1.readlines()
            for i in range(1,len(lines)):
                if lines[-i][36:41] == 'istep':
                    max_step = int(lines[-i].split()[-1])
                    break
            return 'OUT.%s/STRU_ION%d_D' % (suffix,max_step)
    elif calculation == 'md':
        with open(logf) as f1: lines = f1.readlines()
        for i in range(1,len(lines)):
            if lines[-i][1:27] == 'STEP OF MOLECULAR DYNAMICS':
                max_step = int(lines[-i].split()[-1])
                break
        return 'OUT.%s/STRU_MD_%d' % (suffix,max_step)
    elif calculation == 'scf':
        return 'STRU'
    else:
        print("Unrecognized calculation type in %s/INPUT" % abacus_path)
        return 'STRU'

def stru2Structure(struf):
    stru = dpdata.System(struf, fmt="stru")
    stru.to('poscar','POSCAR.tmp')
    ss = Structure.from_file('POSCAR.tmp')
    os.remove('POSCAR.tmp')
    return ss

def check_stru_fixed(struf,fixed):
    block = {}
    with open(struf) as f1: lines = f1.readlines()
    for line in lines:
        if line.strip() == '':continue
        elif line.split()[0] in key_words_list:
            key = line.split()[0]
            block[key] = []
        else:
            block[key].append(line)
    i = 3
    while i < len(block['ATOMIC_POSITIONS']):
        natom = int(block['ATOMIC_POSITIONS'][i])
        for j in range(natom):
            i += 1
            for k in block['ATOMIC_POSITIONS'][i].split()[3:6]:
                if fixed and bool(int(k)):return False
                elif not fixed and not bool(int(k)): return False
        i += 1
    return True

def modify_stru_path(strucf,tpath):
    if tpath[-1] != '/':tpath += '/'
    with open(strucf) as f1: lines = f1.readlines()
    for i,line in enumerate(lines):
        if "ATOMIC_SPECIES" in line and line.split()[0] == "ATOMIC_SPECIES":
            file_numb = 2
        elif ("NUMERICAL_ORBITAL" in line and line.split()[0] == "NUMERICAL_ORBITAL") or \
             ("NUMERICAL_DESCRIPTOR" in line and line.split()[0] == "NUMERICAL_DESCRIPTOR"): 
             file_numb = 0 
        else:continue

        for j in range(i+1,len(lines)):
            if lines[j].strip() in key_words_list: break
            elif lines[j].strip() == '':continue
            ppfile = tpath + os.path.split(lines[j].split()[file_numb])[1]
            tmp_line = ''
            for k in range(file_numb): tmp_line += lines[j].split()[k] + ' '
            lines[j] = tmp_line + ppfile + '\n'

    with open(strucf,'w') as f1: f1.writelines(lines)

