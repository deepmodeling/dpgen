import os,sys,json,glob,shutil
import dpdata
import numpy as np
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
__package__ = 'generator'
from .context import make_fp
from .context import detect_multiplicity
from .context import parse_cur_job
from .context import param_file
from .context import param_old_file
from .context import param_pwscf_file
from .context import param_pwscf_old_file
from .context import param_siesta_file
from .context import param_gaussian_file
from .context import param_cp2k_file
from .context import machine_file
from .context import param_diy_file
from .context import make_kspacing_kpoints
from .context import my_file_cmp
from .context import setUpModule
from .comp_sys import test_atom_names
from .comp_sys import test_atom_types
from .comp_sys import test_coord
from .comp_sys import test_cell
from pymatgen.io.vasp import Kpoints,Incar
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

cp2k_input_ref="\
&GLOBAL\n\
PROJECT DPGEN\n\
&END GLOBAL\n\
&FORCE_EVAL\n\
METHOD QS\n\
STRESS_TENSOR ANALYTICAL\n\
&DFT\n\
BASIS_SET_FILE_NAME ./cp2k_basis_pp_file/BASIS_MOLOPT\n\
POTENTIAL_FILE_NAME ./cp2k_basis_pp_file/GTH_POTENTIALS\n\
CHARGE 0\n\
UKS F\n\
MULTIPLICITY 1\n\
&MGRID\n\
CUTOFF 400\n\
REL_CUTOFF 50\n\
NGRIDS 4\n\
&END MGRID\n\
&QS\n\
EPS_DEFAULT 1.0E-12\n\
&END QS\n\
&SCF\n\
SCF_GUESS ATOMIC\n\
EPS_SCF 1.0E-6\n\
MAX_SCF 50\n\
&OT\n\
MINIMIZER DIIS\n\
PRECONDITIONER FULL_SINGLE_INVERSE\n\
&END OT\n\
&END SCF\n\
&XC\n\
&XC_FUNCTIONAL PBE\n\
&END XC_FUNCTIONAL\n\
&VDW_POTENTIAL\n\
DISPERSION_FUNCTIONAL PAIR_POTENTIAL\n\
&PAIR_POTENTIAL\n\
TYPE DFTD3\n\
PARAMETER_FILE_NAME ./cp2k_basis_pp_file/dftd3.dat\n\
REFERENCE_FUNCTIONAL PBE\n\
&END PAIR_POTENTIAL\n\
&END VDW_POTENTIAL\n\
&END XC\n\
&END DFT\n\
&SUBSYS\n\
&CELL\n\
&END CELL\n\
&COORD\n\
@include coord.xyz\n\
&END COORD\n\
&KIND H\n\
BASIS_SET DZVP-MOLOPT-GTH\n\
POTENTIAL GTH-PBE-q1\n\
&END KIND\n\
&KIND C\n\
BASIS_SET DZVP-MOLOPT-GTH\n\
POTENTIAL GTH-PBE-q4\n\
&END KIND\n\
&KIND N\n\
BASIS_SET DZVP-MOLOPT-GTH\n\
POTENTIAL GTH-PBE-q5\n\
&END KIND\n\
&END SUBSYS\n\
&PRINT\n\
&FORCES ON\n\
&END FORCES\n\
&END PRINT\n\
&END FORCE_EVAL\n"


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
    with open(dump_file, 'w') as fp:
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


def _check_cp2k_input_head(testCase, idx) :
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
        testCase.assertEqual(('\n'.join(lines_check)).strip(), cp2k_input_ref.strip())


class TestMakeFPPwscf(unittest.TestCase):
    def test_make_fp_pwscf(self):
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

class TestMakeFPSIESTA(unittest.TestCase):
    def test_make_fp_siesta(self):
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

    def test_make_fp_vasp_old(self):
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


class TestMakeFPGaussian(unittest.TestCase):
    def test_make_fp_gaussian(self):
        if os.path.isdir('iter.000000') :
            shutil.rmtree('iter.000000')
        with open (param_gaussian_file, 'r') as fp :
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
        _check_gaussian_input_head(self, 0)
        _check_potcar(self, 0, jdata['fp_pp_path'], jdata['fp_pp_files'])
        shutil.rmtree('iter.000000')

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
        _check_cp2k_input_head(self, 0)
        _check_potcar(self, 0, jdata['fp_pp_path'], jdata['fp_pp_files'])
        shutil.rmtree('iter.000000')


if __name__ == '__main__':
    unittest.main()


