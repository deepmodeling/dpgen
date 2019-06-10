import os,json,glob,shutil,filecmp
import dpdata
import numpy as np
import unittest

from context import make_fp_vasp
from context import make_fp_pwscf
from context import parse_cur_job
from context import param_file
from context import machine_file
from comp_sys import test_atom_names
from comp_sys import test_atom_types
from comp_sys import test_coord
from comp_sys import test_cell

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
    

def _make_fake_md(idx, md_descript, atom_types, type_map) :
    """
    md_descript: list of dimension
                 [n_sys][n_MD][n_frame]
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
            

def _check_incar(testCase, idx) :
    fp_path = os.path.join('iter.%06d' % idx, '02.fp')
    testCase.assertTrue(os.path.isfile(os.path.join(fp_path, 'INCAR')))
    tasks = glob.glob(os.path.join(fp_path, 'task.*'))
    for ii in tasks :
        testCase.assertTrue(filecmp.cmp(
            os.path.join(fp_path, 'INCAR'), 
            os.path.join(ii, 'INCAR')))
    

def _check_potcar(testCase, idx, fp_pp_path, fp_pp_files) :
    testCase.assertEqual(len(fp_pp_files), 1)
    fp_path = os.path.join('iter.%06d' % idx, '02.fp')
    testCase.assertTrue(os.path.isfile(os.path.join(fp_pp_path, fp_pp_files[0])))
    tasks = glob.glob(os.path.join(fp_path, 'task.*'))
    for ii in tasks :
        testCase.assertTrue(filecmp.cmp(
            os.path.join(fp_pp_path, fp_pp_files[0]), 
            os.path.join(ii, fp_pp_files[0])))
    

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
        make_fp_vasp(0, jdata)
        _check_sel(self, 0, jdata['fp_task_max'], jdata['model_devi_f_trust_lo'], jdata['model_devi_f_trust_hi'])
        _check_poscars(self, 0, jdata['fp_task_max'], jdata['type_map'])
        _check_incar(self, 0)
        _check_potcar(self, 0, jdata['fp_pp_path'], jdata['fp_pp_files'])
        shutil.rmtree('iter.000000')

if __name__ == '__main__':
    unittest.main()

            
