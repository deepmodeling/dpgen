#!/usr/bin/env python3
import os
import numpy as np

def make_run_opt_script(fmax ) :

    ret ="import os,sys,glob,time                                                                 \n"
    ret+="import numpy as np                                                                      \n"
    ret+="from ase.io import read                                                                 \n"
    ret+="from ase.optimize import BFGS,QuasiNewton,LBFGS                                         \n"
    ret+="from ase.constraints import UnitCellFilter, ExpCellFilter\n"
    ret+="from deepmd.calculator import DP                                                        \n"
    ret+="                                                                                        \n"
    ret+="def Get_Element_Num(elements):                                                          \n"
    ret+="    '''Using the Atoms.symples to Know Element&Num'''                                   \n"
    ret+="    element = []                                                                        \n"
    ret+="    ele = {}                                                                            \n"
    ret+="    element.append(elements[0])                                                         \n"
    ret+="    for x  in elements:                                                                 \n"
    ret+="        if x not in element :                                                           \n"
    ret+="            element.append(x)                                                           \n"
    ret+="    for x in element:                                                                   \n"
    ret+="        ele[x] = elements.count(x)                                                      \n"
    ret+="    return element, ele                                                                 \n"
    ret+="                                                                                        \n"
    ret+="def Write_Contcar(element, ele, lat, pos):                                              \n"
    ret+="    '''Write CONTCAR'''                                                                 \n"
    ret+="    f = open('CONTCAR','w')                                                             \n"
    ret+="    f.write('ASE-DPKit-Optimization\\n')                                                \n"
    ret+="    f.write('1.0\\n')                                                                   \n"
    ret+="    for i in range(3):                                                                  \n"
    ret+="        f.write('%15.10f %15.10f %15.10f\\n' % tuple(lat[i]))                           \n"
    ret+="    for x in element:                                                                   \n"
    ret+="        f.write(x + '  ')                                                               \n"
    ret+="    f.write('\\n')                                                                      \n"
    ret+="    for x in element:                                                                   \n"
    ret+="        f.write(str(ele[x]) + '  ')                                                     \n"
    ret+="    f.write('\\n')                                                                      \n"
    ret+="    f.write('Direct\\n')                                                                \n"
    ret+="    na = sum(ele.values())                                                              \n"
    ret+="    dpos = np.dot(pos,np.linalg.inv(lat))                                               \n"
    ret+="    for i in range(na):                                                                 \n"
    ret+="        f.write('%15.10f %15.10f %15.10f\\n' % tuple(dpos[i]))                          \n"
    ret+="                                                                                        \n"
    ret+="def Write_Outcar(element, ele, volume, lat, pos, ene, force, stress,pstress):           \n"
    ret+="    '''Write OUTCAR'''                                                                  \n"
    ret+="    f = open('OUTCAR','w')                                                              \n"
    ret+="    for x in element:                                                                   \n"
    ret+="        f.write('VRHFIN =' + str(x) + '\\n')                                            \n"
    ret+="    f.write('ions per type =')                                                          \n"
    ret+="    for x in element:                                                                   \n"
    ret+="        f.write('%5d' % ele[x])                                                         \n"
    ret+="    #f.write('\\nvolume of cell :\\n')                                                  \n"
    ret+="    f.write('\\nDirection     XX             YY             ZZ             XY             YZ             ZX\\n') \n"
    ret+="    f.write('in kB')                                                                    \n"
    ret+="    f.write('%15.6f' % stress[0])                                                       \n"
    ret+="    f.write('%15.6f' % stress[1])                                                       \n"
    ret+="    f.write('%15.6f' % stress[2])                                                       \n"
    ret+="    f.write('%15.6f' % stress[3])                                                       \n"
    ret+="    f.write('%15.6f' % stress[4])                                                       \n"
    ret+="    f.write('%15.6f' % stress[5])                                                       \n"
    ret+="    f.write('\\n')                                                                      \n"
    ret+="    ext_pressure = np.sum(stress[0] + stress[1] + stress[2])/3.0 - pstress              \n"
    ret+="    f.write('external pressure = %20.6f kB    Pullay stress = %20.6f  kB\\n'% (ext_pressure, pstress))\n"
    ret+="    f.write('volume of cell : %20.6f\\n' % volume)                                      \n"
    ret+="    f.write('direct lattice vectors\\n')                                                \n"
    ret+="    for i in range(3):                                                                  \n"
    ret+="        f.write('%10.6f %10.6f %10.6f\\n' % tuple(lat[i]))                              \n"
    ret+="    f.write('POSITION                                       TOTAL-FORCE(eV/Angst)\\n')  \n"
    ret+="    f.write('-------------------------------------------------------------------\\n')   \n"
    ret+="    na = sum(ele.values())                                                              \n"
    ret+="    for i in range(na):                                                                 \n"
    ret+="        f.write('%15.6f %15.6f %15.6f' % tuple(pos[i]))                                 \n"
    ret+="        f.write('%15.6f %15.6f %15.6f\\n' % tuple(force[i]))                            \n"
    ret+="    f.write('-------------------------------------------------------------------\\n')   \n"
    ret+="    f.write('energy  without entropy= %20.6f %20.6f\\n' % (ene, ene/na))                \n"
    ret+="    enthalpy = ene + pstress * volume / 1602.17733                                      \n"
    ret+="    f.write('enthalpy is  TOTEN    = %20.6f %20.6f\\n' % (enthalpy, enthalpy/na))       \n"
    ret+="                                                                                        \n"
    ret+="def Write_calylog(templog):                                                             \n"
    ret+="    '''For Write Evolve Structures Log into caly.log'''                                 \n"
    ret+="    f = open('opt.log','a+')                                                           \n"
    ret+="    f.write(templog+'\\n')                                                              \n"
    ret+="    f.close()                                                                           \n"
    ret+="def read_stress():\n"
    ret+="    pstress = 0\n"
    ret+="    #assert os.path.exists('./input.dat'), 'input.dat does not exist!'\n"
    ret+="    try:\n"
    ret+="        f = open('input.dat','r')\n"
    ret+="    except:\n"
    ret+="        assert os.path.exists('../input.dat'),' now we are in %s, do not find ../input.dat'%(os.getcwd())\n"
    ret+="        f = open('../input.dat','r')\n"
    ret+="    lines = f.readlines()\n"
    ret+="    f.close()\n"
    ret+="    for line in lines:\n"
    ret+="        if line[0] == '#':\n"
    ret+="            continue\n"
    ret+="        if 'PSTRESS' in line:\n"
    ret+="            pstress = float(line.split('=')[1])\n"
    ret+="    return pstress\n"
    ret+="                                                                                        \n"
    ret+="                                                                                        \n"
    ret+="def run_opt(stress):                                                                          \n"
    ret+="    '''Using the ASE&DP to Optimize Configures'''                                       \n"
    ret+="    ''' > 600 Steps Called Failure Optimization'''                                      \n"
    ret+="    \n"
    ret+="    os.system('mv OUTCAR OUTCAR-last')\n"
    ret+="    is_dpgen = True\n"
    ret+="    nn = 1\n"
    ret+="    try:\n"
    ret+="        model_path = sys.argv[1]                                                            \n"
    ret+="        Model_List = glob.glob('%s/graph*pb'%model_path)                                    \n"
    ret+="        calc = DP(model=Model_List[1])    # init the model before iteration          \n"
    ret+="    except:\n"
    ret+="        assert os.path.exists('graph.pb'), 'did not found graph.pb in this directory %s, or did you forget to add args?'%(os.getcwd())\n"
    ret+="        calc = DP(model='graph.pb')    # init the model before iteration          \n"
    ret+="        is_dpgen = False\n"
    ret+="        nn = 3\n"
    ret+="                                                                                        \n"
    ret+="    print('Start to Optimize Structures by DP----------')                               \n"
    ret+="                                                                                        \n"
    ret+="    Opt_Step = 200                                                                      \n"
    ret+="    start = time.time()                                                                 \n"
    ret+="    # pstress kbar\n"
    ret+="    pstress = stress\n"
    ret+="    # kBar to eV/A^3\n"
    ret+="    # 1 eV/A^3 = 160.21766028 GPa\n"
    ret+="    # 1 / 160.21766028 ~ 0.006242\n"
    ret+="    aim_stress = 1.0 * pstress* 0.01 * 0.6242  / 10.0 \n"
    ret+="    to_be_opti = read('POSCAR')                                                         \n"
    ret+="    to_be_opti.calc = calc                                                              \n"
    ret+="    ucf = UnitCellFilter(to_be_opti, scalar_pressure=aim_stress)\n"
    ret+="    atoms_vol_2 = to_be_opti.get_volume()                                               \n"
    ret+="    for i in range(nn):                                                                  \n"
    ret+="        if is_dpgen:\n"
    ret+="            opt = LBFGS(ucf,trajectory='traj.traj')                                                         \n"
    ret+="        else:\n"
    ret+="            opt = LBFGS(ucf)                                                         \n"
    ret+="        #opt = QuasiNewton(to_be_opti)                                                  \n"
    ret+="        #opt = BFGS(to_be_opti)                                                         \n"
    ret+="        #opt = BFGS(to_be_opti,trajectory='traj.traj',logfile='opt.log')                \n"
    ret+="        opt.run(fmax=%s,steps=200)                       \n"%str(fmax)
    ret+="                                                                                        \n"
    ret+="        atoms_lat = to_be_opti.cell                                                     \n"
    ret+="        atoms_pos = to_be_opti.positions                                                \n"
    ret+="        atoms_force = to_be_opti.get_forces()                                           \n"
    ret+="        atoms_stress = to_be_opti.get_stress()                                          \n"
    ret+="        # eV/A^3 to GPa\n"
    ret+="        atoms_stress = atoms_stress/(0.01*0.6242)\n"
    ret+="        #atoms_num = to_be_opti.get_atomic_numbers()                                    \n"
    ret+="        atoms_symbols = to_be_opti.get_chemical_symbols()                               \n"
    ret+="        #atoms_formula = to_be_opti.get_chemical_formula()                              \n"
    ret+="        atoms_ene = to_be_opti.get_potential_energy()                                   \n"
    ret+="        atoms_vol = to_be_opti.get_volume()                                             \n"
    ret+="                                                                                        \n"
    ret+="        element, ele = Get_Element_Num(atoms_symbols)                                   \n"
    ret+="                                                                                        \n"
    ret+="        Write_Contcar(element, ele, atoms_lat, atoms_pos)                               \n"
    ret+="        Write_Outcar(element, ele, atoms_vol, atoms_lat, atoms_pos,atoms_ene, atoms_force, atoms_stress * -10.0, pstress)                   \n"
    ret+="                                                                                        \n"
    ret+="                                                                                        \n"
    ret+="    stop = time.time()                                                                  \n"
    ret+="    _cwd = os.getcwd()                                                                  \n"
    ret+="    _cwd = os.path.basename(_cwd)                                                       \n"
    ret+="    print('%s is done, time: %s' % (_cwd,stop-start))                                   \n"
    ret+="stress = read_stress()\n"
    ret+="run_opt(stress)                                                                               \n"
    return ret                                                                                     
                                                                                                   

def make_check_outcar_script( ) :

    ret = "import numpy as np                \n" 
    ret+= "import os,sys,glob,time           \n" 
    ret+= "from deepmd.calculator import DP  \n"
    ret+= "from ase.io import read           \n" 
    ret+= "                                                      \n" 
    ret+= "def Get_Element_Num(elements):                        \n" 
    ret+= "    '''Using the Atoms.symples to Know Element&Num''' \n"     
    ret+= "    element = []                                      \n"     
    ret+= "    ele = {}                                          \n"
    ret+= "    element.append(elements[0])                       \n"
    ret+= "    for x  in elements:                               \n"
    ret+= "        if x not in element :                         \n"
    ret+= "            element.append(x)                         \n"
    ret+= "    for x in element:                                 \n"
    ret+= "        ele[x] = elements.count(x)                    \n"
    ret+= "    return element, ele                               \n"
    ret+= "                                                      \n" 
    ret+= "def Write_Contcar(element, ele, lat, pos):                                              \n"
    ret+= "    '''Write CONTCAR'''                                                                 \n"
    ret+= "    f = open('CONTCAR','w')                                                             \n"
    ret+= "    f.write('ASE-DPKit-FAILED-nan\\n')                                                \n" 
    ret+= "    f.write('1.0\\n')                                                                   \n"
    ret+= "    for i in range(3):                                                                  \n"
    ret+= "        f.write('%15.10f %15.10f %15.10f\\n' % tuple(lat[i]))                           \n"
    ret+= "    for x in element:                                                                   \n"
    ret+= "        f.write(x + '  ')                                                               \n"
    ret+= "    f.write('\\n')                                                                      \n"
    ret+= "    for x in element:                                                                   \n"
    ret+= "        f.write(str(ele[x]) + '  ')                                                     \n"
    ret+= "    f.write('\\n')                                                                      \n"
    ret+= "    f.write('Direct\\n')                                                                \n"
    ret+= "    na = sum(ele.values())                                                              \n"
    ret+= "    dpos = np.dot(pos,np.linalg.inv(lat))                                               \n"
    ret+= "    for i in range(na):                                                                 \n"
    ret+= "        f.write('%15.10f %15.10f %15.10f\\n' % tuple(dpos[i]))                          \n"
    ret+= "                                                                                        \n"
    ret+= "def Write_Outcar(element, ele, volume, lat, pos, ene, force, stress,pstress):           \n"
    ret+= "    '''Write OUTCAR'''                                                                  \n"
    ret+= "    f = open('OUTCAR','w')                                                              \n"
    ret+= "    for x in element:                                                                   \n" 
    ret+= "        f.write('VRHFIN =' + str(x) + '\\n')                                            \n"
    ret+= "    f.write('ions per type =')                                                          \n"
    ret+= "    for x in element:                                                                   \n"
    ret+= "        f.write('%5d' % ele[x])                                                         \n"
    ret+= "    #f.write('\\nvolume of cell :\\n')                                                  \n"
    ret+= "    f.write('\\nDirection     XX             YY             ZZ             XY             YZ             ZX\\n') \n"
    ret+= "    f.write('in kB')                                                                    \n"
    ret+= "    f.write('%15.6f' % stress[0])                                                       \n"
    ret+= "    f.write('%15.6f' % stress[1])                                                       \n"
    ret+= "    f.write('%15.6f' % stress[2])                                                       \n"
    ret+= "    f.write('%15.6f' % stress[3])                                                       \n"
    ret+= "    f.write('%15.6f' % stress[4])                                                       \n"
    ret+= "    f.write('%15.6f' % stress[5])                                                       \n"
    ret+= "    f.write('\\n')                                                                      \n"
    ret+= "    ext_pressure = np.sum(stress[0] + stress[1] + stress[2])/3.0 - pstress              \n"
    ret+= "    f.write('external pressure = %20.6f kB    Pullay stress = %20.6f  kB\\n'% (ext_pressure, pstress))\n"
    ret+= "    f.write('volume of cell : %20.6f\\n' % volume)                                      \n"
    ret+= "    f.write('direct lattice vectors\\n')                                                \n"
    ret+= "    for i in range(3):                                                                  \n"
    ret+= "        f.write('%10.6f %10.6f %10.6f\\n' % tuple(lat[i]))                              \n"
    ret+= "    f.write('POSITION                                       TOTAL-FORCE(eV/Angst)\\n')  \n"
    ret+= "    f.write('-------------------------------------------------------------------\\n')   \n"
    ret+= "    na = sum(ele.values())                                                              \n"
    ret+= "    for i in range(na):                                                                 \n"
    ret+= "        f.write('%15.6f %15.6f %15.6f' % tuple(pos[i]))                                 \n"
    ret+= "        f.write('%15.6f %15.6f %15.6f\\n' % tuple(force[i]))                            \n"
    ret+= "    f.write('-------------------------------------------------------------------\\n')   \n"
    ret+= "    f.write('energy  without entropy= %20.6f %20.6f\\n' % (ene, ene))                \n"
    ret+= "    enthalpy = ene + pstress * volume / 1602.17733                                      \n"
    ret+= "    f.write('enthalpy is  TOTEN    = %20.6f %20.6f\\n' % (enthalpy, enthalpy))       \n"
    ret+= "                                                                                                                                \n"
    ret+= "def check():                                                                                                                    \n"
    ret+= "                                                                                                                                \n"
    ret+= "    from deepmd.calculator import DP                                                                                            \n"
    ret+= "    from ase.io import read                                                                                                     \n"
    ret+= "    model_path = sys.argv[1]                                                                                                    \n" 
    ret+= "    Model_List = glob.glob('%s/graph*pb'%model_path)                                                                            \n"
    ret+= "    calc = DP(model='%s'%(Model_List[0]))    # init the model before iteration                                                  \n"
    ret+= "                                                                                                                                \n"
    ret+= "    to_be_opti = read('POSCAR')                                                                                                 \n"
    ret+= "    to_be_opti.calc = calc                                                                                                      \n"
    ret+= "    # ---------------------------------                                                                                         \n"
    ret+= "    # for failed outcar                                                                                                         \n"
    ret+= "    atoms_symbols_f = to_be_opti.get_chemical_symbols()                                                                         \n"
    ret+= "    element_f, ele_f = Get_Element_Num(atoms_symbols_f)                                                                         \n"
    ret+= "    atoms_vol_f = to_be_opti.get_volume()                                                                                       \n"
    ret+= "    atoms_stress_f = to_be_opti.get_stress()                                                                                    \n"
    ret+= "    atoms_stress_f = atoms_stress_f/(0.01*0.6242)                                                                               \n"
    ret+= "    atoms_lat_f = to_be_opti.cell                                                                                               \n"
    ret+= "    atoms_pos_f = to_be_opti.positions                                                                                          \n"
    ret+= "    atoms_force_f = to_be_opti.get_forces()                                                                                     \n"
    ret+= "    atoms_ene_f =  610612509                                                                                                    \n"
    ret+= "    # ---------------------------------                                                                                         \n"
    ret+= "    Write_Contcar(element_f, ele_f, atoms_lat_f, atoms_pos_f)                                                                   \n"
    ret+= "    Write_Outcar(element_f, ele_f, atoms_vol_f, atoms_lat_f, atoms_pos_f,atoms_ene_f, atoms_force_f, atoms_stress_f * -10.0, 0) \n" 
    ret+= " \n"
    ret+= "cwd = os.getcwd()                                   \n"
    ret+= "if not os.path.exists(os.path.join(cwd,'OUTCAR')):  \n"
    ret+= "    check()  \n"
    return ret                                                                                     
                                                                                                   

def make_calypso_input(nameofatoms,numberofatoms,
                           numberofformula,volume,
                           distanceofion,psoratio,popsize,
                           maxstep,icode,split,vsc,
                           maxnumatom,ctrlrange,pstress,fmax):
    assert len(distanceofion) == len(nameofatoms) #"check distance of ions and the number of atoms"
    assert len(distanceofion[0]) == len(nameofatoms)
    ret = "################################ The Basic Parameters of CALYPSO ################################\n"
    ret+= "# A string of one or several words contain a descriptive name of the system (max. 40 characters).\n"
    ret+= "SystemName = %s\n"%(''.join(nameofatoms))
    ret+= "# Number of different atomic species in the simulation.\n"
    ret+= "NumberOfSpecies = %d\n"%(len(nameofatoms))
    ret+= "# Element symbols of the different chemical species.\n"
    ret+= "NameOfAtoms = %s\n"%(' '.join(nameofatoms))
    ret+= "# Number of atoms for each chemical species in one formula unit. \n"
    ret+= "NumberOfAtoms = %s\n"%(' '.join(list(map(str,numberofatoms))))
    ret+= "# The range of formula unit per cell in your simulation. \n"
    ret+= "NumberOfFormula = %s\n"%(' '.join(list(map(str,numberofformula))))
    ret+= "# The volume per formula unit. Unit is in angstrom^3.\n"
    ret+= "Volume = %s\n"%(volume[0])
    ret+= "# Minimal distance between atoms of each chemical species. Unit is in angstrom.\n"
    ret+= "@DistanceOfIon \n"
    for temp in distanceofion:
        ret+="%4s \n"%(' '.join(list(map(str,temp))))
    ret+= "@End\n"
    ret+= "# It determines which algorithm should be adopted in the simulation.\n"
    ret+= "Ialgo = 2\n"
    ret+= "# Ialgo = 1 for Global PSO\n"
    ret+= "# Ialgo = 2 for Local PSO (default value)\n"
    ret+= "# The proportion of the structures generated by PSO.\n"
    ret+= "PsoRatio = %s\n"%(psoratio[0])
    ret+= "# The population size. Normally, it has a larger number for larger systems.\n"
    ret+= "PopSize = %d\n"%(popsize[0])
    ret+= "# The Max step for iteration\n"
    ret+= "MaxStep = %d\n"%(maxstep[0])
    ret+= "#It determines which method should be adopted in generation the random structure. \n"                             
    ret+= "GenType= 1 \n"
    ret+= "# 1 under symmetric constraints\n"
    ret+= "# 2 grid method for large system\n"
    ret+= "# 3 and 4 core grow method \n"
    ret+= "# 0 combination of all method\n"
    ret+= "# If GenType=3 or 4, it determined the small unit to grow the whole structure\n"
    ret+= "# It determines which local optimization method should be interfaced in the simulation.\n"
    ret+= "ICode= %d\n"%(icode[0])
    ret+= "# ICode= 1 interfaced with VASP\n"
    ret+= "# ICode= 2 interfaced with SIESTA\n"
    ret+= "# ICode= 3 interfaced with GULP\n"
    ret+= "# The number of lbest for local PSO\n"
    ret+= "NumberOfLbest=4\n"
    ret+= "# The Number of local optimization for each structure.\n"
    ret+= "NumberOfLocalOptim= 3\n"
    ret+= "# The command to perform local optimiztion calculation (e.g., VASP, SIESTA) on your computer.\n"
    ret+= "Command = sh submit.sh\n"
    ret+= "MaxTime = 9000 \n"
    ret+= "# If True, a previous calculation will be continued.\n"
    ret+= "PickUp = F\n"
    ret+= "# At which step will the previous calculation be picked up.\n"
    ret+= "PickStep = 1\n"
    ret+= "# If True, the local optimizations performed by parallel\n"
    ret+= "Parallel = F\n"
    ret+= "# The number node for parallel \n"
    ret+= "NumberOfParallel = 4\n"
    ret+= "Split = %s\n"%(split)
    ret+= "PSTRESS = %s\n"%(str(pstress[0]))
    ret+= "fmax = %s\n"%(str(fmax[0]))
    ret+= "################################ End of The Basic Parameters of CALYPSO #######################\n"
    if vsc == 'T':
        assert len(ctrlrange) == len(nameofatoms) #'check distance of ions and the number of atoms'
        ret+= "##### The Parameters For Variational Stoichiometry  ##############\n"
        ret+= "## If True, Variational Stoichiometry structure prediction is performed\n"
        ret+= "VSC = %s\n"%(vsc)
        ret+= "# The Max Number of Atoms in unit cell\n"
        ret+= "MaxNumAtom = %s\n"%(maxnumatom[0])
        ret+= "# The Variation Range for each type atom \n"
        ret+= "@CtrlRange\n"                                                                                                      
        for ttemp in ctrlrange:
            ret+="%4s \n"%(' '.join(list(map(str,ttemp))))
        ret+= "@end\n"
        ret+= "###################End Parameters for VSC ##########################\n"
    return ret


def _make_model_devi_native_calypso(jdata, cur_job, calypso_run_opt_path):

    # Crystal Parameters
    nameofatoms = cur_job.get('NameOfAtoms')
    numberofatoms = cur_job.get('NumberOfAtoms')
    numberofformula = cur_job.get('NumberOfFormula',[1,1])
    volume = cur_job.get('Volume')
    distanceofion = cur_job.get('DistanceOfIon')
    psoratio = cur_job.get('PsoRatio')
    popsize = cur_job.get('PopSize')
    maxstep = cur_job.get('MaxStep')
    icode = cur_job.get('ICode',[1])
    split = cur_job.get('Split','T')
    # VSC Control
    maxnumatom = None
    ctrlrange = None
    vsc = cur_job.get('VSC','F')
    if vsc == 'T':
        maxnumatom = cur_job.get('MaxNumAtom')
        ctrlrange = cur_job.get('CtrlRange')

    # Optimization
    pstress = cur_job.get('PSTRESS',[0.001])
    fmax = cur_job.get('fmax',[0.01])

    # Cluster

    # 2D

    file_c = make_calypso_input(nameofatoms,numberofatoms,
                               numberofformula,volume,
                               distanceofion,psoratio,popsize,
                               maxstep,icode,split,vsc,
                               maxnumatom,ctrlrange,pstress,fmax)
    with open(os.path.join(calypso_run_opt_path, 'input.dat'), 'w') as cin :
        cin.write(file_c)

def write_model_devi_out(devi, fname):
    assert devi.shape[1] == 8
    #assert devi.shape[1] == 7
    header = "%5s" % "step"
    for item in 'vf':
        header += "%16s%16s%16s" % (f"max_devi_{item}", f"min_devi_{item}",f"avg_devi_{item}")
    header += "%16s"%str('min_dis')
    np.savetxt(fname,
               devi,
               fmt=['%5d'] + ['%17.6e' for _ in range(7)],
               delimiter='',
               header=header)
    return devi

