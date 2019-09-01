import os,json,glob,shutil,filecmp
import unittest
from pymatgen import Structure

from context_surf import *


class TestCollVasp(unittest.TestCase):
    def setUp(self):
        self.surfs=["surf-100"] #,"surf-110","surf-111"]
        self.elongs=["elong-0.500", "elong-1.000", "elong-1.500", "elong-2.000", "elong-2.500",\
             "elong-3.000", "elong-3.500", "elong-4.000", "elong-5.000", "elong-6.000",\
             "elong-7.000", "elong-8.000" ]
        with open (param_file, 'r') as fp :
            jdata = json.load (fp)
        out_dir = out_dir_name(jdata)
        jdata['out_dir'] = out_dir
        self.root_dir= out_dir
        create_path(out_dir)
        make_super_cell_pymatgen(jdata)
        place_element(jdata)
        make_vasp_relax(jdata)
        make_scale(jdata)
        pert_scaled(jdata)

    def tearDown(self):
        shutil.rmtree(self.root_dir)

    def test(self):
        surfs=glob.glob("surf.al.fcc.01x01x01/00.place_ele/surf*")
        surfs=[ii.split('/')[-1] for ii in surfs]
        surfs.sort()
        self.assertEqual(surfs,self.surfs)
        poscars=glob.glob("surf.al.fcc.01x01x01/00.place_ele/surf*/sys*/POSCAR")
        for poscar in poscars:
            surf=poscar.split('/')[-3]
            st1=Structure.from_file(surf+'.POSCAR')
            st2=Structure.from_file(poscar)
            self.assertEqual(st1,st2)
        
        for surf in self.surfs:
            elongs=glob.glob("surf.al.fcc.01x01x01/01.scale_pert/"+surf+"/sys-*/scale-1.000/el*")
            elongs=[ii.split('/')[-1] for ii in elongs]
            elongs.sort()
            self.assertEqual(elongs,self.elongs)
             
