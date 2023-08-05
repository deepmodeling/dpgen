import glob
import json
import os
import shutil
import sys
import unittest

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
__package__ = "data"
from .context import setUpModule  # noqa: F401
from .context_bulk import (
    abacus_param_file,
    abacus_ref_Cu_coord,
    create_path,
    get_abacus_STRU,
    make_abacus_relax,
    make_scale_ABACUS,
    make_super_cell_ABACUS,
    make_super_cell_STRU,
    make_unit_cell_ABACUS,
    out_dir_name,
    pert_scaled,
    place_element_ABACUS,
)


class TestGenBulkABACUS(unittest.TestCase):
    def setUp(self):
        self.alloy = []
        with open(abacus_param_file) as fp:
            jdata = json.load(fp)
        out_dir = out_dir_name(jdata)
        self.out_dir = out_dir
        jdata["out_dir"] = out_dir
        self.elements = jdata["elements"]
        self.scale_numb = len(jdata["scale"])
        self.pert_numb = jdata["pert_numb"]
        self.root_dir = out_dir
        self.jdata = jdata
        create_path(out_dir)

    def tearDown(self):
        # pass
        shutil.rmtree(self.root_dir)

    def test(self):
        jdata = self.jdata
        stru_data = make_unit_cell_ABACUS(jdata)
        supercell_stru = make_super_cell_ABACUS(jdata, stru_data)
        place_element_ABACUS(jdata, supercell_stru)
        make_abacus_relax(jdata, {"fp_resources": {}})
        make_scale_ABACUS(jdata)
        pert_scaled(jdata)
        path = self.out_dir + "/00.place_ele"
        # struct0=Structure.from_file(os.path.join(path,"STRU"))
        alloys = glob.glob(os.path.join(path, "sys-*"))
        stru0 = get_abacus_STRU(os.path.join(alloys[0], "STRU"))
        self.assertEqual(len(alloys), stru0["coords"].shape[0] - 1)
        for ii in alloys:
            elem_numb = [int(i) for i in ii.split("/")[-1].split("-")[1:]]
            struct = get_abacus_STRU(os.path.join(ii, "STRU"))
            self.assertEqual(struct["atom_names"], self.elements)
            self.assertEqual(struct["atom_numbs"], elem_numb)
        path = self.out_dir + "/01.scale_pert"
        alloys = glob.glob(os.path.join(path, "sys-*"))
        self.assertEqual(len(alloys), stru0["coords"].shape[0] - 1)
        for ii in alloys:
            scales = glob.glob(os.path.join(ii, "scale-*"))
            self.assertEqual(len(scales), self.scale_numb)
            for scale in scales:
                perts = glob.glob(os.path.join(scale, "[0-9]*"))
                self.assertEqual(len(perts), self.pert_numb + 1)

    def testSTRU(self):
        jdata = self.jdata
        jdata["from_poscar_path"] = "./Cu.STRU"
        make_super_cell_STRU(jdata)
        make_abacus_relax(jdata, {"fp_resources": {}})
        make_scale_ABACUS(jdata)
        pert_scaled(jdata)
        path = self.out_dir + "/00.place_ele"
        # struct0=Structure.from_file(os.path.join(path,"STRU"))
        alloys = glob.glob(os.path.join(path, "sys-*"))
        stru0 = get_abacus_STRU(os.path.join(alloys[0], "STRU"))
        self.assertEqual(
            jdata["super_cell"][0] * jdata["super_cell"][1] * jdata["super_cell"][2],
            stru0["coords"].shape[0],
        )
        for ii in alloys:
            elem_numb = [int(i) for i in ii.split("/")[-1].split("-")[1:]]
            struct = get_abacus_STRU(os.path.join(ii, "STRU"))
            self.assertEqual(struct["atom_numbs"], elem_numb)
            if os.path.basename(ii) == "sys-0004":
                np.testing.assert_almost_equal(
                    struct["coords"], abacus_ref_Cu_coord, decimal=3
                )
        path = self.out_dir + "/01.scale_pert"
        alloys = glob.glob(os.path.join(path, "sys-*"))
        for ii in alloys:
            scales = glob.glob(os.path.join(ii, "scale-*"))
            self.assertEqual(len(scales), self.scale_numb)
            for scale in scales:
                perts = glob.glob(os.path.join(scale, "[0-9]*"))
                self.assertEqual(len(perts), self.pert_numb + 1)


if __name__ == "__main__":
    unittest.main()
