import os,sys,json,glob,shutil
import unittest
from pymatgen.core import Structure, Composition

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
__package__ = 'data'
from .context import setUpModule
from .context_bulk import *

class TestGenBulk(unittest.TestCase):
    def setUp(self):
        self.alloy=[] 
        with open (param_file, 'r') as fp :
            jdata = json.load (fp)
        if "init_fp_style" not in jdata:
            jdata["init_fp_style"] = "VASP"
        out_dir = out_dir_name(jdata)
        self.out_dir= out_dir
        jdata['out_dir'] = out_dir
        self.elements=jdata["elements"]
        self.scale_numb=len(jdata["scale"])
        self.pert_numb=jdata["pert_numb"]
        self.root_dir= out_dir
        create_path(out_dir)
        make_unit_cell(jdata)
        make_super_cell(jdata)
        place_element(jdata)
        make_vasp_relax(jdata,{"fp_resources":{}})
        make_scale(jdata)
        pert_scaled(jdata)

    def tearDown(self):
        shutil.rmtree(self.root_dir)

    def test(self):
        path=self.out_dir+"/00.place_ele"
        struct0=Structure.from_file(os.path.join(path,"POSCAR"))
        alloys=glob.glob(os.path.join(path,"sys-*"))
        self.assertEqual(len(alloys),len(struct0)-1)
        for ii in alloys:
            elem_numb=[int(i) for i in ii.split('/')[-1].split('-')[1:]]
            comp=''
            for num, el in zip(elem_numb,self.elements):
                comp+=el+str(num)
            comp=Composition(comp)
            struct=Structure.from_file(os.path.join(ii,"POSCAR"))
            self.assertEqual(struct.composition,comp) 
        path=self.out_dir+"/01.scale_pert"
        alloys=glob.glob(os.path.join(path,"sys-*"))
        self.assertEqual(len(alloys),len(struct0)-1)
        for ii in alloys:
            scales=glob.glob(os.path.join(ii,"scale-*"))
            self.assertEqual(len(scales),self.scale_numb)
            for scale in scales:
                perts=glob.glob(os.path.join(scale,"[0-9]*")) 
                self.assertEqual(len(perts),self.pert_numb+1)


if __name__ == '__main__':
    unittest.main()