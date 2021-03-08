import os,sys,json,glob,shutil
import unittest
from pymatgen.core import Structure, Element

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
__package__ = 'data'
from .context import setUpModule
from .context_surf import *

class TestGenSurf(unittest.TestCase):
    def setUp(self):
        self.surfs=["surf-100"] 
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
        self.jdata=jdata

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
            vacuum_size=float(Element(self.jdata['elements'][0]).atomic_radius*2)
            self.assertTrue(st1.lattice.c+vacuum_size-st2.lattice.c<0.01)
        
        for surf in self.surfs:
            elongs=glob.glob("surf.al.fcc.01x01x01/01.scale_pert/"+surf+"/sys-*/scale-1.000/el*")
            elongs=[ii.split('/')[-1] for ii in elongs]
            elongs.sort()
            self.assertEqual(elongs,self.elongs)
             
