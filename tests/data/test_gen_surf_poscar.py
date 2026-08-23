import glob
import json
import os
import shutil
import sys
import tempfile
import unittest

from pymatgen.core import Element, Structure

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
__package__ = "data"
from .context import setUpModule  # noqa: F401
from .context_surf_poscar import (
    create_path,
    make_scale,
    make_super_cell_pymatgen,
    make_vasp_relax,
    out_dir_name,
    param_file,
    pert_scaled,
    place_element,
)


class TestGenSurfPOSCAR(unittest.TestCase):
    def setUp(self):
        self.surfs = ["surf-100"]
        self.elongs = [
            "elong-0.500",
            "elong-1.000",
            "elong-1.500",
            "elong-2.000",
            "elong-4.000",
        ]
        with open(param_file) as fp:
            jdata = json.load(fp)
        self.input_dir = tempfile.mkdtemp()
        from_poscar_path = os.path.join(self.input_dir, "nested.POSCAR")
        shutil.copy2("POSCAR", from_poscar_path)
        jdata["from_poscar_path"] = from_poscar_path
        out_dir = out_dir_name(jdata)
        jdata["out_dir"] = out_dir
        self.root_dir = out_dir
        create_path(out_dir)
        make_super_cell_pymatgen(jdata)
        place_element(jdata)
        make_vasp_relax(jdata)
        make_scale(jdata)
        pert_scaled(jdata)
        self.jdata = jdata

    def tearDown(self):
        shutil.rmtree(self.root_dir)
        shutil.rmtree(self.input_dir)

    def test(self):
        surfs = glob.glob(os.path.join(self.root_dir, "01.scale_pert", "surf*"))
        surfs = [ii.split("/")[-1] for ii in surfs]
        surfs.sort()
        self.assertEqual(surfs, self.surfs)
        poscars = glob.glob(
            os.path.join(self.root_dir, "00.place_ele", "surf*", "sys*", "POSCAR")
        )
        for poscar in poscars:
            surf = poscar.split("/")[-3]
            st1 = Structure.from_file(surf + ".POSCAR")
            st2 = Structure.from_file(poscar)
            vacuum_size = float(Element(self.jdata["elements"][0]).atomic_radius * 2)
            self.assertTrue(st1.lattice.c + vacuum_size - st2.lattice.c < 0.01)

        for surf in self.surfs:
            elongs = glob.glob(
                os.path.join(
                    self.root_dir,
                    "01.scale_pert",
                    surf,
                    "sys-*",
                    "scale-1.000",
                    "el*",
                )
            )
            elongs = [ii.split("/")[-1] for ii in elongs]
            elongs.sort()
            self.assertEqual(elongs, self.elongs)
