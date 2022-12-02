import os,sys,json,glob,shutil
import dpdata
import numpy as np
import unittest
import importlib

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
__package__ = 'generator'
from .context import make_fp
from .context import detect_multiplicity
from .context import parse_cur_job
from .context import param_file
from .context import param_file_merge_traj
from .context import param_old_file
from .context import param_pwscf_file
from .context import param_pwscf_old_file
from .context import param_abacus_post_file
from .context import param_diy_abacus_post_file
from .context import param_siesta_file
from .context import param_gaussian_file
from .context import param_cp2k_file
from .context import param_cp2k_file_exinput
from .context import param_amber_file
from .context import ref_cp2k_file_input
from .context import ref_cp2k_file_exinput
from .context import machine_file
from .context import param_diy_file
from .context import param_multiple_trust_file
from .context import make_kspacing_kpoints
from .context import my_file_cmp
from .context import setUpModule
from .comp_sys import test_atom_names
from .comp_sys import test_atom_types
from .comp_sys import test_coord
from .comp_sys import test_cell
from pymatgen.io.vasp import Kpoints,Incar
from .context import param_pwmat_file
import scipy.constants as pc

vasp_incar_ref = "PREC=A\n\
ENCUT=600\n\
ISYM=0\n\
ALGO=fast\n\
EDIFF=1e-05\n\
LREAL=A\n\
NPAR=1\n\
KPAR=1\n\
NELMIN=4\n\
ISIF=2\n\
ISMEAR=1\n\
SIGMA=0.25\n\
IBRION=-1\n\
NSW=0\n\
LWAVE=F\n\
LCHARG=F\n\
PSTRESS=0\n\
KSPACING=0.16\n\
KGAMMA=F\n";

vasp_incar_ele_temp_ref = "PREC=A\n\
ENCUT=600\n\
ISYM=0\n\
ALGO=fast\n\
EDIFF=1e-05\n\
LREAL=A\n\
NPAR=1\n\
KPAR=1\n\
NELMIN=4\n\
ISIF=2\n\
ISMEAR=-1\n\
SIGMA=%.10f\n\
IBRION=-1\n\
NSW=0\n\
LWAVE=F\n\
LCHARG=F\n\
PSTRESS=0\n\
KSPACING=0.16\n\
KGAMMA=F\n";

pwscf_input_ref="&control\n\
calculation='scf',\n\
restart_mode='from_scratch',\n\
outdir='./OUT',\n\
tprnfor=.TRUE.,\n\
tstress=.TRUE.,\n\
disk_io='none',\n\
pseudo_dir='./',\n\
/\n\
&system\n\
vdw_corr='TS',\n\
ecutwfc=110,\n\
ts_vdw_econv_thr=1e-08,\n\
nosym=.TRUE.,\n\
ibrav=0,\n\
nat=6,\n\
ntyp=3,\n\
/\n\
&electrons\n\
conv_thr=1e-08,\n\
/\n"

siesta_input_ref="\
SystemName        system\n\
SystemLabel       system\n\
NumberOfAtoms     6\n\
NumberOfSpecies   3\n\
\n\
WriteForces       T\n\
WriteCoorStep     T\n\
WriteCoorXmol     T\n\
WriteMDXmol       T\n\
WriteMDHistory    T\n\
\n\
MeshCutoff            300 Ry\n\
DM.Tolerance          1.000000e-04\n\
DM.NumberPulay         5\n\
DM.UseSaveDM          true\n\
XC.functional          GGA\n\
XC.authors             PBE\n\
MD.UseSaveXV           T\n\
\n\
DM.UseSaveDM           F\n\
WriteDM                F\n\
WriteDM.NetCDF         F\n\
WriteDMHS.NetCDF       F\n\
\n\
%block Chemical_Species_label\n\
1	6	C\n\
2	1	H\n\
3	7	N\n"


gaussian_input_ref="""%nproc=14
#force b3lyp/6-31g*

DPGEN
"""


pwmat_input_ref = "4 1\n\
in.atom=atom.config\n\
ecut=50\n\
e_error=0.0001\n\
rho_error=0.0001\n\
scf_iter0_1=6 4 3 0.0000 0.025 2\n\
scf_iter0_2=94 4 3 1.0000 0.025 2\n\
mp_n123=147 57 39 0 0 0 2\n\
out.wg=F\n\
out.rho=F\n\
out.mlmd=T\n\
job=scf\n\
IN.PSP1 = C.SG15.PBE.UPF\n\
IN.PSP2 = H.SG15.PBE.UPF\n\
IN.PSP3 = N.SG15.PBE.UPF\n";

abacus_input_ref = "INPUT_PARAMETERS\n\
calculation scf\n\
ntype 2\n\
ecutwfc 80.000000\n\
scf_thr 1.000000e-07\n\
scf_nmax 50\n\
basis_type pw\n\
gamma_only 1\n\
dft_functional pbe\n\
mixing_type pulay\n\
mixing_beta 0.400000\n\
symmetry 1\n\
nbands 5\n\
nspin 1\n\
ks_solver cg\n\
smearing_method fixed\n\
smearing_sigma 0.001000\n\
cal_force 1\n\
cal_stress 1\n\
deepks_out_labels 0\n\
deepks_descriptor_lmax 0\n\
deepks_scf 0\n\
deepks_model model.ptg\n"

abacus_kpt_ref = "K_POINTS\n\
0\n\
Gamma\n\
1 1 1 0 0 0\n"


def _box2lmpbox(orig, box) :
    lohi = np.zeros([3,2])
    for dd in range(3) :
        lohi[dd][0] = orig[dd]
    tilt = np.zeros(3)
    tilt[0] = box[1][0]
    tilt[1] = box[2][0]
    tilt[2] = box[2][1]
    lens = np.zeros(3)
    lens[0] = box[0][0]
    lens[1] = box[1][1]
    lens[2] = box[2][2]
    for dd in range(3) :
        lohi[dd][1] = lohi[dd][0] + lens[dd]
    return lohi, tilt

def _box2dumpbox(orig, box) :
    lohi, tilt = _box2lmpbox(orig, box)
    xy = tilt[0]
    xz = tilt[1]
    yz = tilt[2]
    bounds = np.zeros([3,2])
    bounds[0][0] = lohi[0][0] + min(0.0,xy,xz,xy+xz)
    bounds[0][1] = lohi[0][1] + max(0.0,xy,xz,xy+xz)
    bounds[1][0] = lohi[1][0] + min(0.0,yz)
    bounds[1][1] = lohi[1][1] + max(0.0,yz)
    bounds[2][0] = lohi[2][0]
    bounds[2][1] = lohi[2][1]
    return bounds, tilt


def _write_lammps_dump(sys, dump_file, f_idx = 0) :
    cell = sys['cells'][f_idx].reshape([3,3])
    coord = sys['coords'][f_idx].reshape([-1,3])
    bd, tilt = _box2dumpbox(np.zeros(3), cell)
    atype = sys['atom_types']
    natoms = len(sys['atom_types'])
    with open(dump_file, 'a') as fp:
        fp.write('ITEM: TIMESTEP\n')
        fp.write('0\n')
        fp.write('ITEM: NUMBER OF ATOMS\n')
        fp.write(str(natoms) + '\n')
        fp.write('ITEM: BOX BOUNDS xy xz yz pp pp pp\n')
        for ii in range(3):
            fp.write('%f %f %f\n' % (bd[ii][0], bd[ii][1], tilt[ii]))
        fp.write('ITEM: ATOMS id type x y z\n')
        for ii in range(natoms) :
            fp.write('%d %d %f %f %f\n' % (ii+1, atype[ii]+1, coord[ii][0], coord[ii][1], coord[ii][2]))


def _make_fake_md(idx, md_descript, atom_types, type_map, ele_temp = None) :
    """
    md_descript: list of dimension
                 [n_sys][n_MD][n_frame]
    ele_temp: list of dimension
                 [n_sys][n_MD]
    """
    natoms = len(atom_types)
    ntypes = len(type_map)
    atom_types = np.array(atom_types, dtype = int)
    atom_numbs = [np.sum(atom_types == ii) for ii in range(ntypes)]
    sys = dpdata.System()
    sys.data['atom_names'] = type_map
    sys.data['atom_numbs'] = atom_numbs
    sys.data['atom_types'] = atom_types
    for sidx,ss in enumerate(md_descript) :
        for midx,mm in enumerate(ss) :
            nframes = len(mm)
            cells  = np.random.random([nframes,3,3])
            coords = np.random.random([nframes,natoms,3])
            sys.data['coords'] = coords
            sys.data['cells'] = cells
            task_dir = os.path.join('iter.%06d' % idx,
                                    '01.model_devi',
                                    'task.%03d.%06d' % (sidx, midx))
            os.makedirs(os.path.join(task_dir, 'traj'), exist_ok = True)
            for ii in range(nframes) :
                _write_lammps_dump(sys,
                                   os.path.join(task_dir,
                                                'traj',
                                                '%d.lammpstrj' % ii))
            md_out = np.zeros([nframes, 7])
            md_out[:,0] = np.arange(nframes)
            md_out[:,4] = mm
            np.savetxt(os.path.join(task_dir, 'model_devi.out'), md_out)
            if ele_temp is not None:
                with open(os.path.join(task_dir, 'job.json'), 'w') as fp:
                    json.dump({"ele_temp": ele_temp[sidx][midx]}, fp)

def _make_fake_md_merge_traj(idx, md_descript, atom_types, type_map, ele_temp = None) :
    """
    md_descript: list of dimension
                 [n_sys][n_MD][n_frame]
    ele_temp: list of dimension
                 [n_sys][n_MD]
    """
    natoms = len(atom_types)
    ntypes = len(type_map)
    atom_types = np.array(atom_types, dtype = int)
    atom_numbs = [np.sum(atom_types == ii) for ii in range(ntypes)]
    sys = dpdata.System()
    sys.data['atom_names'] = type_map
    sys.data['atom_numbs'] = atom_numbs
    sys.data['atom_types'] = atom_types
    for sidx,ss in enumerate(md_descript) :
        for midx,mm in enumerate(ss) :
            nframes = len(mm)
            cells  = np.random.random([nframes,3,3])
            coords = np.random.random([nframes,natoms,3])
            sys.data['coords'] = coords
            sys.data['cells'] = cells
            task_dir = os.path.join('iter.%06d' % idx,
                                    '01.model_devi',
                                    'task.%03d.%06d' % (sidx, midx))
            cwd = os.getcwd()
            os.makedirs(task_dir,exist_ok = True)
            for ii in range(nframes):
                _write_lammps_dump(sys,os.path.join(task_dir,'all.lammpstrj'),ii)
            file_content = """\
0.000000000000000000e+01 0.000000000000000000e+00 0.000000000000000000e+00 0.000000000000000000e+00 0.000000000000000000e+00 0.000000000000000000e+00 0.000000000000000000e+00
1.000000000000000000e+01 0.000000000000000000e+00 0.000000000000000000e+00 0.000000000000000000e+00 2.899999999999999800e-02 0.000000000000000000e+00 0.000000000000000000e+00
2.000000000000000000e+01 0.000000000000000000e+00 0.000000000000000000e+00 0.000000000000000000e+00 5.799999999999999600e-02 0.000000000000000000e+00 0.000000000000000000e+00
3.000000000000000000e+01 0.000000000000000000e+00 0.000000000000000000e+00 0.000000000000000000e+00 8.699999999999999400e-02 0.000000000000000000e+00 0.000000000000000000e+00
4.000000000000000000e+01 0.000000000000000000e+00 0.000000000000000000e+00 0.000000000000000000e+00 1.159999999999999920e-01 0.000000000000000000e+00 0.000000000000000000e+00
5.000000000000000000e+01 0.000000000000000000e+00 0.000000000000000000e+00 0.000000000000000000e+00 1.449999999999999900e-01 0.000000000000000000e+00 0.000000000000000000e+00
6.000000000000000000e+01 0.000000000000000000e+00 0.000000000000000000e+00 0.000000000000000000e+00 1.739999999999999880e-01 0.000000000000000000e+00 0.000000000000000000e+00
7.000000000000000000e+01 0.000000000000000000e+00 0.000000000000000000e+00 0.000000000000000000e+00 2.029999999999999860e-01 0.000000000000000000e+00 0.000000000000000000e+00
8.000000000000000000e+01 0.000000000000000000e+00 0.000000000000000000e+00 0.000000000000000000e+00 2.319999999999999840e-01 0.000000000000000000e+00 0.000000000000000000e+00
9.000000000000000000e+01 0.000000000000000000e+00 0.000000000000000000e+00 0.000000000000000000e+00 2.610000000000000098e-01 0.000000000000000000e+00 0.000000000000000000e+00
"""
            with open(os.path.join(task_dir, 'model_devi.out') , 'w') as fp:
                fp.write(file_content)
            if ele_temp is not None:
                with open(os.path.join(task_dir, 'job.json'), 'w') as fp:
                    json.dump({"ele_temp": ele_temp[sidx][midx]}, fp)

def _check_poscars(testCase, idx, fp_task_max, type_map) :
    fp_path = os.path.join('iter.%06d' % idx, '02.fp')
    candi_files = glob.glob(os.path.join(fp_path, 'candidate.shuffled.*.out'))
    candi_files.sort()
    sys_idx = [str(os.path.basename(ii).split('.')[2]) for ii in candi_files]
    for sidx,ii in zip(sys_idx, candi_files) :
        md_task = []
        f_idx = []
        with open(ii) as fp:
            for ii in fp :
                md_task.append(ii.split()[0])
                f_idx.append(ii.split()[1])
        md_task = md_task[:fp_task_max]
        f_idx = f_idx[:fp_task_max]
        cc = 0
        for tt,ff in zip(md_task, f_idx) :
            traj_file = os.path.join(tt, 'traj', '%d.lammpstrj' % int(ff))
            poscar_file = os.path.join(fp_path,
                                       'task.%03d.%06d' % (int(sidx), cc),
                                       'POSCAR')
            cc += 1
            sys0 = dpdata.System(traj_file, fmt = 'lammps/dump', type_map = type_map)
            sys1 = dpdata.System(poscar_file, fmt = 'vasp/poscar')
            test_atom_names(testCase, sys0, sys1)

def _check_poscars_merge_traj(testCase, idx, fp_task_max, type_map ) :
    fp_path = os.path.join('iter.%06d' % idx, '02.fp')
    candi_files = glob.glob(os.path.join(fp_path, 'candidate.shuffled.*.out'))
    candi_files.sort()
    sys_idx = [str(os.path.basename(ii).split('.')[2]) for ii in candi_files]
    for sidx,ii in zip(sys_idx, candi_files) :
        md_task = []
        f_idx = []
        with open(ii) as fp:
            for ii in fp :
                md_task.append(ii.split()[0])
                f_idx.append(ii.split()[1])
        md_task = md_task[:fp_task_max]
        f_idx = f_idx[:fp_task_max]
        cc = 0
        label_0 = 0
        label_1 = 0

        for tt,ff in zip(md_task, f_idx) :
            traj_file = os.path.join(tt, 'all.lammpstrj')
            poscar_file = os.path.join(fp_path,
                                       'task.%03d.%06d' % (int(sidx), cc),
                                       'POSCAR')
            cc += 1
            sys0 = dpdata.System(traj_file, fmt = 'lammps/dump', type_map = type_map)
            sys1 = dpdata.System(poscar_file, fmt = 'vasp/poscar')
            new_coords_0 = float(sys1["coords"][0][0][0])
            new_coords_1 = float(sys1["coords"][0][1][0])
            if (label_0 == new_coords_0 and label_1 == new_coords_1):
                raise RuntimeError("The exact same POSCAR is generated under different first-principles calculation catalogs")
            label_0 = new_coords_0
            label_1 = new_coords_1
            test_atom_names(testCase, sys0[int(int(ff)/10)], sys1)

def _check_kpoints_exists(testCase, idx) :
    fp_path = os.path.join('iter.%06d' % idx, '02.fp')
    tasks = glob.glob(os.path.join(fp_path, 'task.*'))
    for ii in tasks :
        testCase.assertTrue(os.path.isfile(os.path.join(ii, 'KPOINTS')))

def _check_kpoints(testCase, idx) :
    fp_path = os.path.join('iter.%06d' % idx, '02.fp')
    tasks = glob.glob(os.path.join(fp_path, 'task.*'))
    for ii in tasks :
        kpoints=Kpoints.from_file(os.path.join(os.path.join(ii, 'KPOINTS')))
        incar=Incar.from_file(os.path.join(os.path.join(ii, 'INCAR')))
        kspacing = incar['KSPACING']
        gamma = incar['KGAMMA']
        if isinstance(gamma,bool):
           pass
        else:
           if gamma[0].upper()=="T":
              gamma=True
           else:
              gamma=False
        ret=make_kspacing_kpoints(os.path.join(os.path.join(ii, 'POSCAR')), kspacing, gamma)
        kpoints_ref=Kpoints.from_string(ret)
        testCase.assertEqual(repr(kpoints), repr(kpoints_ref))


def _check_incar_exists(testCase, idx) :
    fp_path = os.path.join('iter.%06d' % idx, '02.fp')
    # testCase.assertTrue(os.path.isfile(os.path.join(fp_path, 'INCAR')))
    tasks = glob.glob(os.path.join(fp_path, 'task.*'))
    for ii in tasks :
        my_file_cmp(testCase,
                    os.path.join(fp_path, 'INCAR'),
                    os.path.join(ii, 'INCAR'))


def _check_potcar(testCase, idx, fp_pp_path, fp_pp_files) :
    nfile = len(fp_pp_files)
    fp_path = os.path.join('iter.%06d' % idx, '02.fp')
    for ii in range(nfile):
        testCase.assertTrue(os.path.isfile(os.path.join(fp_pp_path, fp_pp_files[ii])))
    tasks = glob.glob(os.path.join(fp_path, 'task.*'))
    for ii in tasks :
        for jj in range(nfile):
            my_file_cmp(testCase,
                        os.path.join(fp_pp_path, fp_pp_files[jj]),
                        os.path.join(ii, fp_pp_files[jj]))


def _check_sel(testCase, idx, fp_task_max, flo, fhi):
    
    def _trust_limitation_check(sys_idx, lim):
        if isinstance(lim, list):
            sys_lim = lim[sys_idx]
        elif isinstance(lim, dict):
            sys_lim = lim[str(sys_idx)]
        else:
            sys_lim = lim
        return sys_lim
    
    fp_path = os.path.join('iter.%06d' % idx, '02.fp')
    candi_files = glob.glob(os.path.join(fp_path, 'candidate.shuffled.*.out'))
    candi_files.sort()
    sys_idx = [str(os.path.basename(ii).split('.')[2]) for ii in candi_files]
    for sidx,ii in zip(sys_idx, candi_files) :
        md_task = []
        f_idx = []
        with open(ii) as fp:
            for ii in fp :
                md_task.append(ii.split()[0])
                f_idx.append(ii.split()[1])
        md_task = md_task[:fp_task_max]
        f_idx = f_idx[:fp_task_max]
        flo = _trust_limitation_check(int(sidx), flo)
        fhi = _trust_limitation_check(int(sidx), fhi)
        for tt,ff in zip(md_task, f_idx):
            md_value = np.loadtxt(os.path.join(tt, 'model_devi.out'))
            fvalue = md_value[int(ff)][4]
            testCase.assertTrue(fvalue >= flo)
            testCase.assertTrue(fvalue <  fhi)


def _check_incar(testCase, idx):
    fp_path = os.path.join('iter.%06d' % idx, '02.fp')
    tasks = glob.glob(os.path.join(fp_path, 'task.*'))
    cwd = os.getcwd()
    for ii in tasks :
        os.chdir(ii)
        with open('INCAR') as fp:
            incar = fp.read()
            testCase.assertEqual(incar.strip(), vasp_incar_ref.strip())
        os.chdir(cwd)

def _check_incar_ele_temp(testCase, idx, ele_temp):
    fp_path = os.path.join('iter.%06d' % idx, '02.fp')
    tasks = glob.glob(os.path.join(fp_path, 'task.*'))
    cwd = os.getcwd()
    for ii in tasks :
        os.chdir(ii)
        bname = os.path.basename(ii)
        sidx = int(bname.split('.')[1])
        tidx = int(bname.split('.')[2])
        with open('INCAR') as fp:
            incar = fp.read()
            incar0 = Incar.from_string(incar)
            # make_fake_md: the frames in a system shares the same ele_temp
            incar1 = Incar.from_string(vasp_incar_ele_temp_ref%(ele_temp[sidx][0] * pc.Boltzmann / pc.electron_volt))
            for ii in incar0.keys():
                # skip checking nbands...
                if ii == 'NBANDS':
                    continue
                testCase.assertAlmostEqual(incar0[ii], incar1[ii], msg = 'key %s differ' % (ii), places = 5)
        os.chdir(cwd)

def _check_pwscf_input_head(testCase, idx) :
    fp_path = os.path.join('iter.%06d' % idx, '02.fp')
    tasks = glob.glob(os.path.join(fp_path, 'task.*'))
    for ii in tasks :
        ifile = os.path.join(ii, 'input')
        testCase.assertTrue(os.path.isfile(ifile))
        with open(ifile) as fp:
            lines = fp.read().split('\n')
        for idx, jj in enumerate(lines) :
            if 'ATOMIC_SPECIES' in jj :
                break
        lines = lines[:idx]
        testCase.assertEqual(('\n'.join(lines)).strip(), pwscf_input_ref.strip())

def _check_abacus_input(testCase, idx) :
    fp_path = os.path.join('iter.%06d' % idx, '02.fp')
    tasks = glob.glob(os.path.join(fp_path, 'task.*'))
    for ii in tasks :
        ifile = os.path.join(ii, 'INPUT')
        testCase.assertTrue(os.path.isfile(ifile))
        with open(ifile) as fp:
            lines = fp.read().split('\n')
        testCase.assertEqual(('\n'.join(lines)).strip(), abacus_input_ref.strip())

def _check_abacus_kpt(testCase, idx) :
    fp_path = os.path.join('iter.%06d' % idx, '02.fp')
    tasks = glob.glob(os.path.join(fp_path, 'task.*'))
    for ii in tasks :
        ifile = os.path.join(ii, 'KPT')
        testCase.assertTrue(os.path.isfile(ifile))
        with open(ifile) as fp:
            lines = fp.read().split('\n')
        testCase.assertEqual(('\n'.join(lines)).strip(), abacus_kpt_ref.strip())

def _check_siesta_input_head(testCase, idx) :
    fp_path = os.path.join('iter.%06d' % idx, '02.fp')
    tasks = glob.glob(os.path.join(fp_path, 'task.*'))
    for ii in tasks :
        ifile = os.path.join(ii, 'input')
        testCase.assertTrue(os.path.isfile(ifile))
        with open(ifile) as fp:
            lines = fp.read().split('\n')
        for idx, jj in enumerate(lines) :
            if '%endblock Chemical_Species_label' in jj :
                break
        lines = lines[:idx]
        testCase.assertEqual(('\n'.join(lines)).strip(), siesta_input_ref.strip())


def _check_gaussian_input_head(testCase, idx) :
    fp_path = os.path.join('iter.%06d' % idx, '02.fp')
    tasks = glob.glob(os.path.join(fp_path, 'task.*'))
    for ii in tasks :
        ifile = os.path.join(ii, 'input')
        testCase.assertTrue(os.path.isfile(ifile))
        with open(ifile) as fp:
            lines = fp.read().split('\n')
        for idx, jj in enumerate(lines) :
            if '0 1' in jj :
                break
        lines = lines[:idx]
        testCase.assertEqual(('\n'.join(lines)).strip(), gaussian_input_ref.strip())


def _check_cp2k_input_head(testCase, idx, ref_out) :
    fp_path = os.path.join('iter.%06d' % idx, '02.fp')
    tasks = glob.glob(os.path.join(fp_path, 'task.*'))
    for ii in tasks :
        ifile = os.path.join(ii, 'input.inp')
        testCase.assertTrue(os.path.isfile(ifile))
        with open(ifile) as fp:
            lines = fp.read().split('\n')
        for idx, jj in enumerate(lines) :
            if '&CELL' in jj :
                cell_start_idx = idx
            if '&END CELL' in jj :
                cell_end_idx = idx
        lines_check = lines[:cell_start_idx+1] + lines[cell_end_idx:]
        testCase.assertEqual(('\n'.join(lines_check)).strip(), ref_out.strip())



def _check_pwmat_input(testCase, idx):
    fp_path = os.path.join('iter.%06d' % idx, '02.fp')
    tasks = glob.glob(os.path.join(fp_path, 'task.*'))
    cwd = os.getcwd()
    for ii in tasks :
        os.chdir(ii)
        os.system("sed -i '8c mp_n123=147 57 39 0 0 0 2' etot.input")
        with open('etot.input') as fp:
            lines = fp.read()
            testCase.assertEqual(lines.strip(), pwmat_input_ref.strip())
        os.chdir(cwd)

def _check_symlink_user_forward_files(testCase, idx, file):
    fp_path = os.path.join('iter.%06d' % idx, '02.fp')
    tasks = glob.glob(os.path.join(fp_path, 'task.*'))
    cwd = os.getcwd()
    for ii in tasks:
        os.chdir(ii)
        testCase.assertEqual(os.path.isfile("vdw_kernel.bindat"), True)
        os.chdir(cwd)

class TestMakeFPPwscf(unittest.TestCase):
    def test_make_fp_pwscf(self):
        setUpModule()
        if os.path.isdir('iter.000000') :
            shutil.rmtree('iter.000000')
        with open (param_pwscf_file, 'r') as fp :
            jdata = json.load (fp)
        with open (machine_file, 'r') as fp:
            mdata = json.load (fp)
        md_descript = []
        nsys = 2
        nmd = 3
        n_frame = 10
        for ii in range(nsys) :
            tmp = []
            for jj in range(nmd) :
                tmp.append(np.arange(0, 0.29, 0.29/10))
            md_descript.append(tmp)
        atom_types = [0, 1, 2, 2, 0, 1]
        type_map = jdata['type_map']
        _make_fake_md(0, md_descript, atom_types, type_map)
        make_fp(0, jdata, {})
        _check_sel(self, 0, jdata['fp_task_max'], jdata['model_devi_f_trust_lo'], jdata['model_devi_f_trust_hi'])
        _check_poscars(self, 0, jdata['fp_task_max'], jdata['type_map'])
        _check_pwscf_input_head(self, 0)
        _check_potcar(self, 0, jdata['fp_pp_path'], jdata['fp_pp_files'])
        shutil.rmtree('iter.000000')

    def test_make_fp_pwscf_old(self):
        setUpModule()
        if os.path.isdir('iter.000000') :
            shutil.rmtree('iter.000000')
        with open (param_pwscf_old_file, 'r') as fp :
            jdata = json.load (fp)
        with open (machine_file, 'r') as fp:
            mdata = json.load (fp)
        md_descript = []
        nsys = 2
        nmd = 3
        n_frame = 10
        for ii in range(nsys) :
            tmp = []
            for jj in range(nmd) :
                tmp.append(np.arange(0, 0.29, 0.29/10))
            md_descript.append(tmp)
        atom_types = [0, 1, 2, 2, 0, 1]
        type_map = jdata['type_map']
        _make_fake_md(0, md_descript, atom_types, type_map)
        make_fp(0, jdata, {})
        _check_sel(self, 0, jdata['fp_task_max'], jdata['model_devi_f_trust_lo'], jdata['model_devi_f_trust_hi'])
        _check_poscars(self, 0, jdata['fp_task_max'], jdata['type_map'])
        _check_pwscf_input_head(self, 0)
        _check_potcar(self, 0, jdata['fp_pp_path'], jdata['fp_pp_files'])
        shutil.rmtree('iter.000000')

class TestMakeFPABACUS(unittest.TestCase):
    def test_make_fp_abacus(self):
        setUpModule()
        if os.path.isdir('iter.000000') :
            shutil.rmtree('iter.000000')
        with open (param_abacus_post_file, 'r') as fp :
            jdata = json.load (fp)
        fp.close()
        with open (machine_file, 'r') as fp:
            mdata = json.load (fp)
        fp.close()
        md_descript = []
        nsys = 2
        nmd = 3
        n_frame = 10
        for ii in range(nsys) :
            tmp = []
            for jj in range(nmd) :
                tmp.append(np.arange(0, 0.29, 0.29/10))
            md_descript.append(tmp)
        atom_types = [0, 0, 0, 0, 1]
        type_map = jdata['type_map']
        _make_fake_md(0, md_descript, atom_types, type_map)
        make_fp(0, jdata, {})
        _check_sel(self, 0, jdata['fp_task_max'], jdata['model_devi_f_trust_lo'], jdata['model_devi_f_trust_hi'])
        _check_poscars(self, 0, jdata['fp_task_max'], jdata['type_map'])
        _check_abacus_input(self, 0)
        _check_abacus_kpt(self, 0)
        _check_potcar(self, 0, jdata['fp_pp_path'], jdata['fp_pp_files'])
        shutil.rmtree('iter.000000')

    def test_make_fp_abacus_from_input(self):
        ## Verify if user chooses to diy ABACUS INPUT totally.
        setUpModule()
        if os.path.isdir('iter.000000') :
            shutil.rmtree('iter.000000')
        with open (param_diy_abacus_post_file, 'r') as fp :
            jdata = json.load (fp)
        fp.close()
        with open (machine_file, 'r') as fp:
            mdata = json.load (fp)
        fp.close()
        md_descript = []
        nsys = 2
        nmd = 3
        n_frame = 10
        for ii in range(nsys) :
            tmp = []
            for jj in range(nmd) :
                tmp.append(np.arange(0, 0.29, 0.29/10))
            md_descript.append(tmp)
        atom_types = [0, 1, 0, 1]
        type_map = jdata['type_map']
        _make_fake_md(0, md_descript, atom_types, type_map)
        make_fp(0, jdata, {})
        _check_sel(self, 0, jdata['fp_task_max'], jdata['model_devi_f_trust_lo'], jdata['model_devi_f_trust_hi'])
        _check_poscars(self, 0, jdata['fp_task_max'], jdata['type_map'])
        _check_abacus_input(self, 0)
        _check_abacus_kpt(self, 0)
        _check_potcar(self, 0, jdata['fp_pp_path'], jdata['fp_pp_files'])
        shutil.rmtree('iter.000000')

class TestMakeFPAMBERDiff(unittest.TestCase):
    def test_make_fp_amber_diff(self):
        setUpModule()
        if os.path.isdir('iter.000000') :
            shutil.rmtree('iter.000000')
        with open(param_amber_file, 'r') as fp:
            jdata = json.load(fp)
        jdata['mdin_prefix'] = os.path.abspath(jdata['mdin_prefix'])
        task_dir = os.path.join('iter.%06d' % 0,
                        '01.model_devi',
                        'task.%03d.%06d' % (0, 0))
        os.makedirs(task_dir, exist_ok = True)
        with open(os.path.join(task_dir, "rc.mdout"), 'w') as f:
            f.write("Active learning frame written with max. frc. std.:     3.29037 kcal/mol/A")
        import ase
        from ase.io.netcdftrajectory import write_netcdftrajectory
        write_netcdftrajectory(os.path.join(task_dir, 'rc.nc'), ase.Atoms("C", positions=np.zeros((1, 3))))
        make_fp(0, jdata, {})


class TestMakeFPSIESTA(unittest.TestCase):
    def test_make_fp_siesta(self):
        setUpModule()
        if os.path.isdir('iter.000000') :
            shutil.rmtree('iter.000000')
        with open (param_siesta_file, 'r') as fp :
            jdata = json.load (fp)
        with open (machine_file, 'r') as fp:
            mdata = json.load (fp)
        md_descript = []
        nsys = 2
        nmd = 3
        n_frame = 10
        for ii in range(nsys) :
            tmp = []
            for jj in range(nmd) :
                tmp.append(np.arange(0, 0.29, 0.29/10))
            md_descript.append(tmp)
        atom_types = [0, 1, 2, 2, 0, 1]
        type_map = jdata['type_map']
        _make_fake_md(0, md_descript, atom_types, type_map)
        make_fp(0, jdata, {})
        _check_sel(self, 0, jdata['fp_task_max'], jdata['model_devi_f_trust_lo'], jdata['model_devi_f_trust_hi'])
        _check_poscars(self, 0, jdata['fp_task_max'], jdata['type_map'])
        _check_siesta_input_head(self, 0)
        _check_potcar(self, 0, jdata['fp_pp_path'], jdata['fp_pp_files'])
        shutil.rmtree('iter.000000')

class TestMakeFPVasp(unittest.TestCase):
    def test_make_fp_vasp(self):
        setUpModule()
        if os.path.isdir('iter.000000') :
            shutil.rmtree('iter.000000')
        with open (param_file, 'r') as fp :
            jdata = json.load (fp)
        with open (machine_file, 'r') as fp:
            mdata = json.load (fp)
        md_descript = []
        nsys = 2
        nmd = 3
        n_frame = 10
        for ii in range(nsys) :
            tmp = []
            for jj in range(nmd) :
                tmp.append(np.arange(0, 0.29, 0.29/10))
            md_descript.append(tmp)
        atom_types = [0, 1, 0, 1]
        type_map = jdata['type_map']
        _make_fake_md(0, md_descript, atom_types, type_map)
        make_fp(0, jdata, {"fp_user_forward_files" : ["vdw_kernel.bindat"] })
        _check_sel(self, 0, jdata['fp_task_max'], jdata['model_devi_f_trust_lo'], jdata['model_devi_f_trust_hi'])
        _check_poscars(self, 0, jdata['fp_task_max'], jdata['type_map'])
        # _check_incar_exists(self, 0)
        _check_incar(self, 0)
        _check_kpoints_exists(self, 0)
        _check_kpoints(self,0)
        # checked elsewhere
        # _check_potcar(self, 0, jdata['fp_pp_path'], jdata['fp_pp_files'])
        shutil.rmtree('iter.000000')
    
    def test_make_fp_vasp_merge_traj(self):
        setUpModule()
        if os.path.isdir('iter.000000') :
            shutil.rmtree('iter.000000')
        with open (param_file_merge_traj, 'r') as fp :
            jdata = json.load (fp)
        with open (machine_file, 'r') as fp:
            mdata = json.load (fp)
        md_descript = []
        nsys = 2
        nmd = 3
        n_frame = 10
        for ii in range(nsys) :
            tmp = []
            for jj in range(nmd) :
                tmp.append(np.arange(0, 0.29, 0.29/10))
            md_descript.append(tmp)
        atom_types = [0, 1, 0, 1]
        type_map = jdata['type_map']

        _make_fake_md_merge_traj(0, md_descript, atom_types, type_map)
        make_fp(0, jdata, {"fp_user_forward_files" : ["vdw_kernel.bindat"] })
        _check_poscars_merge_traj(self, 0, jdata['fp_task_max'], jdata['type_map'])
        #_check_incar_exists(self, 0)
        _check_incar(self, 0)
        _check_kpoints_exists(self, 0)
        _check_kpoints(self,0)
        # checked elsewhere
        # _check_potcar(self, 0, jdata['fp_pp_path'], jdata['fp_pp_files'])
        shutil.rmtree('iter.000000')

    def test_make_fp_vasp_old(self):
        setUpModule()
        if os.path.isdir('iter.000000') :
            shutil.rmtree('iter.000000')
        with open (param_old_file, 'r') as fp :
            jdata = json.load (fp)
        with open (machine_file, 'r') as fp:
            mdata = json.load (fp)
        md_descript = []
        nsys = 2
        nmd = 3
        n_frame = 10
        for ii in range(nsys) :
            tmp = []
            for jj in range(nmd) :
                tmp.append(np.arange(0, 0.29, 0.29/10))
            md_descript.append(tmp)
        atom_types = [0, 1, 0, 1]
        type_map = jdata['type_map']
        _make_fake_md(0, md_descript, atom_types, type_map)
        make_fp(0, jdata, {})
        _check_sel(self, 0, jdata['fp_task_max'], jdata['model_devi_f_trust_lo'], jdata['model_devi_f_trust_hi'])
        _check_poscars(self, 0, jdata['fp_task_max'], jdata['type_map'])
        # _check_incar_exists(self, 0)
        _check_incar(self, 0)
        _check_kpoints_exists(self, 0)
        _check_kpoints(self,0)
        # checked elsewhere
        # _check_potcar(self, 0, jdata['fp_pp_path'], jdata['fp_pp_files'])
        shutil.rmtree('iter.000000')

    def test_make_fp_vasp_less_sel(self):
        if os.path.isdir('iter.000000') :
            shutil.rmtree('iter.000000')
        with open (param_file, 'r') as fp :
            jdata = json.load (fp)
        with open (machine_file, 'r') as fp:
            mdata = json.load (fp)
        md_descript = []
        nsys = 1
        nmd = 1
        n_frame = 8
        for ii in range(nsys) :
            tmp = []
            for jj in range(nmd) :
                tmp.append(np.arange(0, 0.29, 0.29/10))
            md_descript.append(tmp)
        atom_types = [0, 1, 0, 1]
        type_map = jdata['type_map']
        _make_fake_md(0, md_descript, atom_types, type_map)
        make_fp(0, jdata, {})
        _check_sel(self, 0, jdata['fp_task_max'], jdata['model_devi_f_trust_lo'], jdata['model_devi_f_trust_hi'])
        _check_poscars(self, 0, jdata['fp_task_max'], jdata['type_map'])
        # _check_incar_exists(self, 0)
        _check_incar(self, 0)
        _check_kpoints_exists(self, 0)
        _check_kpoints(self,0)
        # checked elsewhere
        # _check_potcar(self, 0, jdata['fp_pp_path'], jdata['fp_pp_files'])
        shutil.rmtree('iter.000000')


    def test_make_fp_vasp_from_incar(self):
        ## Verify if user chooses to diy VASP INCAR totally.
        setUpModule()
        if os.path.isdir('iter.000000') :
            shutil.rmtree('iter.000000')
        with open (param_diy_file, 'r') as fp :
            jdata = json.load (fp)
        fp.close()
        with open (machine_file, 'r') as fp:
            mdata = json.load (fp)
        fp.close()
        md_descript = []
        nsys = 2
        nmd = 3
        n_frame = 10
        for ii in range(nsys) :
            tmp = []
            for jj in range(nmd) :
                tmp.append(np.arange(0, 0.29, 0.29/10))
            md_descript.append(tmp)
        atom_types = [0, 1, 0, 1]
        type_map = jdata['type_map']
        _make_fake_md(0, md_descript, atom_types, type_map)
        make_fp(0, jdata, {})
        _check_sel(self, 0, jdata['fp_task_max'], jdata['model_devi_f_trust_lo'], jdata['model_devi_f_trust_hi'])
        _check_poscars(self, 0, jdata['fp_task_max'], jdata['type_map'])
        # _check_incar_exists(self, 0)
        _check_incar(self, 0)
        _check_kpoints_exists(self, 0)
        _check_kpoints(self,0)
        # checked elsewhere
        # _check_potcar(self, 0, jdata['fp_pp_path'], jdata['fp_pp_files'])
        shutil.rmtree('iter.000000')

    def test_make_fp_vasp_ele_temp(self):
        ## Verify if user chooses to diy VASP INCAR totally.
        setUpModule()
        if os.path.isdir('iter.000000') :
            shutil.rmtree('iter.000000')
        with open (param_diy_file, 'r') as fp :
            jdata = json.load (fp)
        fp.close()
        with open (machine_file, 'r') as fp:
            mdata = json.load (fp)
        fp.close()
        md_descript = []
        ele_temp = []
        nsys = 2
        nmd = 3
        n_frame = 10
        for ii in range(nsys) :
            tmp = []
            for jj in range(nmd) :
                tmp.append(np.arange(0, 0.29, 0.29/10))
            md_descript.append(tmp)
            ele_temp.append([np.random.random() * 100000] * nmd)
        atom_types = [0, 1, 0, 1]
        type_map = jdata['type_map']
        _make_fake_md(0, md_descript, atom_types, type_map, ele_temp = ele_temp)
        make_fp(0, jdata, {})
        _check_sel(self, 0, jdata['fp_task_max'], jdata['model_devi_f_trust_lo'], jdata['model_devi_f_trust_hi'])
        _check_poscars(self, 0, jdata['fp_task_max'], jdata['type_map'])
        _check_incar_ele_temp(self, 0, ele_temp)
        _check_kpoints_exists(self, 0)
        _check_kpoints(self,0)
        # checked elsewhere
        # _check_potcar(self, 0, jdata['fp_pp_path'], jdata['fp_pp_files'])
        shutil.rmtree('iter.000000')

    def test_make_fp_vasp_multiple_trust_level(self):
        # Verify if sys_idx dependent trust level could be read.
        setUpModule()
        if os.path.isdir('iter.000000') :
            shutil.rmtree('iter.000000')
        with open (param_multiple_trust_file, 'r') as fp :
            jdata = json.load (fp)
        fp.close()
        with open (machine_file, 'r') as fp:
            mdata = json.load (fp)
        fp.close()
        md_descript = []
        ele_temp = []
        nsys = 2
        nmd = 3
        n_frame = 10
        for ii in range(nsys) :
            tmp = []
            for jj in range(nmd) :
                tmp.append(np.arange(0, 0.29, 0.29/10))
            md_descript.append(tmp)
            ele_temp.append([np.random.random() * 100000] * nmd)
        atom_types = [0, 1, 0, 1]
        type_map = jdata['type_map']
        _make_fake_md(0, md_descript, atom_types, type_map, ele_temp = ele_temp)
        make_fp(0, jdata, {})
        _check_sel(self, 0, jdata['fp_task_max'], jdata['model_devi_f_trust_lo'], jdata['model_devi_f_trust_hi'])
        _check_poscars(self, 0, jdata['fp_task_max'], jdata['type_map'])
        _check_incar_ele_temp(self, 0, ele_temp)
        _check_kpoints_exists(self, 0)
        _check_kpoints(self,0)
        # checked elsewhere
        # _check_potcar(self, 0, jdata['fp_pp_path'], jdata['fp_pp_files'])
        shutil.rmtree('iter.000000')
    

class TestMakeFPGaussian(unittest.TestCase):
    def make_fp_gaussian(self, multiplicity="auto"):
        setUpModule()
        if os.path.isdir('iter.000000') :
            shutil.rmtree('iter.000000')
        with open (param_gaussian_file, 'r') as fp :
            jdata = json.load (fp)
        jdata['user_fp_params']['multiplicity'] = multiplicity
        with open (machine_file, 'r') as fp:
            mdata = json.load (fp)
        md_descript = []
        nsys = 2
        nmd = 3
        n_frame = 10
        for ii in range(nsys) :
            tmp = []
            for jj in range(nmd) :
                tmp.append(np.arange(0, 0.29, 0.29/10))
            md_descript.append(tmp)
        atom_types = [0, 1, 2, 2, 0, 1]
        type_map = jdata['type_map']
        _make_fake_md(0, md_descript, atom_types, type_map)
        make_fp(0, jdata, {})
        _check_sel(self, 0, jdata['fp_task_max'], jdata['model_devi_f_trust_lo'], jdata['model_devi_f_trust_hi'])
        _check_poscars(self, 0, jdata['fp_task_max'], jdata['type_map'])
        _check_gaussian_input_head(self, 0)
        _check_potcar(self, 0, jdata['fp_pp_path'], jdata['fp_pp_files'])
        shutil.rmtree('iter.000000')

    @unittest.skipIf(importlib.util.find_spec("openbabel") is None, "requires openbabel")
    def test_make_fp_gaussian(self):
        self.make_fp_gaussian()

    def test_make_fp_gaussian_multiplicity_one(self):
        self.make_fp_gaussian(multiplicity=1)

    def test_detect_multiplicity(self):
        # oxygen O2 3
        self._check_multiplicity(['O', 'O'], 3)
        # methane CH4 1
        self._check_multiplicity(['C', 'H', 'H', 'H', 'H'], 1)
        # CH3 2
        self._check_multiplicity(['C', 'H', 'H', 'H'], 2)
        # CH2 1
        self._check_multiplicity(['C', 'H', 'H'], 1)
        # CH 2
        self._check_multiplicity(['C', 'H'], 2)

    def _check_multiplicity(self, symbols, multiplicity):
        self.assertEqual(detect_multiplicity(np.array(symbols)), multiplicity)

class TestMakeFPCP2K(unittest.TestCase):
    def test_make_fp_cp2k(self):
        setUpModule()
        if os.path.isdir('iter.000000') :
            shutil.rmtree('iter.000000')
        with open (param_cp2k_file, 'r') as fp :
            jdata = json.load (fp)
        with open (machine_file, 'r') as fp:
            mdata = json.load (fp)
        md_descript = []
        nsys = 2
        nmd = 3
        n_frame = 10
        for ii in range(nsys) :
            tmp = []
            for jj in range(nmd) :
                tmp.append(np.arange(0, 0.29, 0.29/10))
            md_descript.append(tmp)
        atom_types = [0, 1, 2, 2, 0, 1]
        type_map = jdata['type_map']
        _make_fake_md(0, md_descript, atom_types, type_map)
        make_fp(0, jdata, {})
        _check_sel(self, 0, jdata['fp_task_max'], jdata['model_devi_f_trust_lo'], jdata['model_devi_f_trust_hi'])
        _check_poscars(self, 0, jdata['fp_task_max'], jdata['type_map'])
        with open(ref_cp2k_file_input, 'r') as f:
            cp2k_input_ref = ''.join(f.readlines())
        _check_cp2k_input_head(self, 0, cp2k_input_ref)
        _check_potcar(self, 0, jdata['fp_pp_path'], jdata['fp_pp_files'])
        shutil.rmtree('iter.000000')
    def test_make_fp_cp2k_exinput(self):
        setUpModule()
        if os.path.isdir('iter.000000') :
            shutil.rmtree('iter.000000')
        with open (param_cp2k_file_exinput, 'r') as fp :
            jdata = json.load (fp)
        with open (machine_file, 'r') as fp:
            mdata = json.load (fp)
        md_descript = []
        nsys = 2
        nmd = 3
        n_frame = 10
        for ii in range(nsys) :
            tmp = []
            for jj in range(nmd) :
                tmp.append(np.arange(0, 0.29, 0.29/10))
            md_descript.append(tmp)
        atom_types = [0, 1, 2, 2, 0, 1]
        type_map = jdata['type_map']
        _make_fake_md(0, md_descript, atom_types, type_map)
        make_fp(0, jdata, {})
        _check_sel(self, 0, jdata['fp_task_max'], jdata['model_devi_f_trust_lo'], jdata['model_devi_f_trust_hi'])
        _check_poscars(self, 0, jdata['fp_task_max'], jdata['type_map'])
        with open(ref_cp2k_file_exinput, 'r') as f:
            cp2k_exinput_ref = ''.join(f.readlines())
        _check_cp2k_input_head(self, 0, cp2k_exinput_ref)
        _check_potcar(self, 0, jdata['fp_pp_path'], jdata['fp_pp_files'])
        shutil.rmtree('iter.000000')



class TestMakeFPPWmat(unittest.TestCase):
    def test_make_fp_pwmat(self):
        setUpModule()
        if os.path.isdir('iter.000000') :
            shutil.rmtree('iter.000000')
        with open (param_pwmat_file, 'r') as fp :
            jdata = json.load (fp)
        with open (machine_file, 'r') as fp:
            mdata = json.load (fp)
        md_descript = []
        nsys = 2
        nmd = 3
        n_frame = 10
        for ii in range(nsys) :
            tmp = []
            for jj in range(nmd) :
                tmp.append(np.arange(0, 0.29, 0.29/10))
            md_descript.append(tmp)
        atom_types = [0, 1, 2, 2, 0, 1]
        type_map = jdata['type_map']
        _make_fake_md(0, md_descript, atom_types, type_map)
        make_fp(0, jdata, {})
        _check_sel(self, 0, jdata['fp_task_max'], jdata['model_devi_f_trust_lo'], jdata['model_devi_f_trust_hi'])
        _check_poscars(self, 0, jdata['fp_task_max'], jdata['type_map'])
        _check_pwmat_input(self, 0)
        _check_potcar(self, 0, jdata['fp_pp_path'], jdata['fp_pp_files'])
        os.system('rm -r iter.000000')
        #shutil.rmtree('iter.000000')

if __name__ == '__main__':
    unittest.main()


