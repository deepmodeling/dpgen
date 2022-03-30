import os,sys,json,glob,shutil
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
__package__ = 'data'
from .context import setUpModule
from .context_bulk import *

class TestGenBulkABACUS(unittest.TestCase):
    def setUp(self):
        self.alloy=[] 
        with open (abacus_param_file, 'r') as fp :
            jdata = json.load (fp)
        out_dir = out_dir_name(jdata)
        self.out_dir= out_dir
        jdata['out_dir'] = out_dir
        self.elements=jdata["elements"]
        self.scale_numb=len(jdata["scale"])
        self.pert_numb=jdata["pert_numb"]
        self.root_dir= out_dir
        create_path(out_dir)
        stru_data = make_unit_cell_ABACUS(jdata)
        supercell_stru = make_super_cell_ABACUS(jdata, stru_data)            
        place_element_ABACUS(jdata, supercell_stru)
        make_abacus_relax(jdata, {"fp_resources":{}})
        make_scale_ABACUS(jdata)
        pert_scaled(jdata)

    def tearDown(self):
        shutil.rmtree(self.root_dir)

    def test(self):
        path=self.out_dir+"/00.place_ele"
        #struct0=Structure.from_file(os.path.join(path,"STRU"))
        alloys=glob.glob(os.path.join(path,"sys-*"))
        stru0 = get_abacus_STRU(os.path.join(alloys[0], "STRU"))
        self.assertEqual(len(alloys),stru0['coords'].shape[0]-1)
        for ii in alloys:
            elem_numb=[int(i) for i in ii.split('/')[-1].split('-')[1:]]
            struct=get_abacus_STRU(os.path.join(ii,"STRU"))
            self.assertEqual(struct["atom_names"],self.elements)
            self.assertEqual(struct["atom_numbs"], elem_numb)
        path=self.out_dir+"/01.scale_pert"
        alloys=glob.glob(os.path.join(path,"sys-*"))
        self.assertEqual(len(alloys), stru0['coords'].shape[0]-1)
        for ii in alloys:
            scales=glob.glob(os.path.join(ii,"scale-*"))
            self.assertEqual(len(scales),self.scale_numb)
            for scale in scales:
                perts=glob.glob(os.path.join(scale,"[0-9]*")) 
                self.assertEqual(len(perts),self.pert_numb+1)


if __name__ == '__main__':
    unittest.main()