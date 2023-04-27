import os

import numpy as np
from dpdata import LabeledSystem
from monty.serialization import dumpfn

import dpgen.auto_test.lib.abacus as abacus
import dpgen.generator.lib.abacus_scf as abacus_scf
from dpgen import dlog
from dpgen.auto_test.Task import Task
from dpgen.util import sepline


class ABACUS(Task):
    def __init__(self, inter_parameter, path_to_poscar):
        self.inter = inter_parameter
        self.inter_type = inter_parameter["type"]
        self.incar = inter_parameter.get("incar", {})
        self.potcar_prefix = inter_parameter.get("potcar_prefix", "")
        self.potcars = inter_parameter.get("potcars", None)
        self.orbfile = inter_parameter.get("orb_files", None)
        self.deepks = inter_parameter.get("deepks_desc", None)
        self.deepks_model = inter_parameter.get("deepks_model", None)
        self.path_to_poscar = path_to_poscar
        self.if_define_orb_file = False if self.orbfile == None else True

    def make_potential_files(self, output_dir):
        stru = os.path.abspath(os.path.join(output_dir, "STRU"))
        if not os.path.isfile(stru):
            raise FileNotFoundError("No file %s" % stru)
        stru_data = abacus_scf.get_abacus_STRU(stru)
        atom_names = stru_data["atom_names"]
        orb_files = stru_data["orb_files"]
        pp_files = stru_data["pp_files"]
        dpks_descriptor = stru_data["dpks_descriptor"]

        if os.path.islink(os.path.join(output_dir, "STRU")):
            stru_path, tmpf = os.path.split(
                os.readlink(os.path.join(output_dir, "STRU"))
            )
        else:
            stru_path = output_dir

        if pp_files == None:
            raise RuntimeError("No pseudopotential information in STRU file")

        pp_dir = os.path.abspath(self.potcar_prefix)
        cwd = os.getcwd()
        os.chdir(output_dir)
        if not os.path.isdir("./pp_orb"):
            os.mkdir("./pp_orb")

        pp_orb_file = []
        for iatom, atomname in enumerate(atom_names):
            # pseudopotential file
            if not self.potcars:
                raise RuntimeError(
                    "please specify the pseudopotential file for each atom type in 'potcars'"
                )
            if atomname not in self.potcars:
                raise RuntimeError(
                    "please specify the pseudopotential file of '%s'" % atomname
                )
            pp_orb_file.append([pp_files[iatom], self.potcars[atomname]])

            # orbital file
            if orb_files:
                if not self.orbfile:
                    raise RuntimeError(
                        "Orbital file is defined in STRU, so please specify the orbital file for each atom type in parameter setting file by 'orb_files'"
                    )
                if atomname not in self.orbfile:
                    raise RuntimeError(
                        "please specify the orbital file of '%s'" % atomname
                    )
                pp_orb_file.append([orb_files[iatom], self.orbfile[atomname]])
            elif self.orbfile:
                dlog.warning(
                    "Orbital is not needed by STRU, so ignore the setting of 'orb_files' in parameter setting file"
                )

        # dpks_descriptor
        if dpks_descriptor:
            if not self.deepks:
                raise RuntimeError(
                    "Deepks descriptor file is defined in STRU, so please specify in parameter setting file by 'deepks_desc'"
                )
            pp_orb_file.append([dpks_descriptor, self.deepks])
        elif self.deepks:
            dlog.warning(
                "Deepks descriptor is not needed by STRU, so ignore the setting of 'deepks_desc' in parameter setting file"
            )

        # dpks model
        if self.deepks_model:
            pp_orb_file.append([self.deepks_model, self.deepks_model])

        # link the files
        for file_stru, file_param in pp_orb_file:
            filename_in_stru = os.path.split(file_stru)[1]
            filename_in_para = os.path.split(file_param)[1]
            if filename_in_stru != filename_in_para:
                dlog.warning(
                    "file name in STRU is not match that defined in parameter setting file: '%s', '%s'."
                    % (filename_in_stru, filename_in_para)
                )

            src_file = os.path.join(pp_dir, file_param)
            if not os.path.isfile(src_file):
                raise RuntimeError("Can not find file %s" % src_file)
            tar_file = os.path.join("pp_orb", filename_in_stru)
            if os.path.isfile(tar_file):
                os.remove(tar_file)
            os.symlink(src_file, tar_file)

        os.chdir(cwd)

        dumpfn(self.inter, os.path.join(output_dir, "inter.json"), indent=4)

    def modify_input(self, incar, x, y):
        if x in incar and incar[x] != y:
            dlog.info("setting %s to %s" % (x, y))
        incar[x] = y

    def make_input_file(self, output_dir, task_type, task_param):
        sepline(ch=output_dir)
        dumpfn(task_param, os.path.join(output_dir, "task.json"), indent=4)

        assert os.path.exists(self.incar), "no INPUT file for relaxation"
        relax_incar_path = os.path.abspath(self.incar)
        incar_relax = abacus_scf.get_abacus_input_parameters(relax_incar_path)

        # deal with relaxation

        cal_type = task_param["cal_type"]
        cal_setting = task_param["cal_setting"]

        # user input INCAR for property calculation
        if "input_prop" in cal_setting and os.path.isfile(cal_setting["input_prop"]):
            incar_prop = os.path.abspath(cal_setting["input_prop"])
            incar = abacus_scf.get_abacus_input_parameters(incar_prop)
            dlog.info(
                "Detected 'input_prop' in 'relaxation', use %s as INPUT, and ignore 'cal_setting'"
                % incar_prop
            )

        # revise INCAR based on the INCAR provided in the "interaction"
        else:
            incar = incar_relax
            for key in cal_setting:
                if key in ["relax_pos", "relax_shape", "relax_vol", "K_POINTS", ""]:
                    continue
                if key[0] == "_":
                    continue
                if "interaction" in key.lower():
                    continue
                incar[key.lower()] = cal_setting[key]

            fix_atom = [False, False, False]
            if cal_type == "relaxation":
                relax_pos = cal_setting["relax_pos"]
                relax_shape = cal_setting["relax_shape"]
                relax_vol = cal_setting["relax_vol"]
                if [relax_pos, relax_shape, relax_vol] == [True, False, False]:
                    self.modify_input(incar, "calculation", "relax")
                elif [relax_pos, relax_shape, relax_vol] == [True, True, True]:
                    self.modify_input(incar, "calculation", "cell-relax")
                elif [relax_pos, relax_shape, relax_vol] == [True, True, False]:
                    self.modify_input(incar, "calculation", "cell-relax")
                    self.modify_input(incar, "fixed_axes", "volume")
                elif [relax_pos, relax_shape, relax_vol] == [False, True, False]:
                    self.modify_input(incar, "calculation", "cell-relax")
                    self.modify_input(incar, "fixed_axes", "volume")
                    fix_atom = [True, True, True]
                elif [relax_pos, relax_shape, relax_vol] == [False, True, True]:
                    self.modify_input(incar, "calculation", "cell-relax")
                    fix_atom = [True, True, True]
                elif [relax_pos, relax_shape, relax_vol] == [False, False, True]:
                    raise RuntimeError(
                        "relax volume but fix shape is not supported for ABACUS"
                    )
                elif [relax_pos, relax_shape, relax_vol] == [False, False, False]:
                    self.modify_input(incar, "calculation", "scf")
                else:
                    raise RuntimeError("not supported calculation setting for ABACUS")

            elif cal_type == "static":
                self.modify_input(incar, "calculation", "scf")

            else:
                raise RuntimeError("not supported calculation type for ABACUS")

            # modify STRU file base on the value of fix_atom
            abacus.stru_fix_atom(os.path.join(output_dir, "STRU"), fix_atom)

        if "basis_type" not in incar:
            dlog.info("'basis_type' is not defined, set to be 'pw'!")
            self.modify_input(incar, "basis_type", "pw")
        if "lcao" in incar["basis_type"].lower() and not self.if_define_orb_file:
            mess = (
                "The basis_type is %s, but not define orbital file!!!"
                % incar["basis_type"]
            )
            raise RuntimeError(mess)
        if "deepks_model" in incar:
            model_file = os.path.split(incar["deepks_model"])[1]
            self.modify_input(incar, "deepks_model", os.path.join("pp_orb", model_file))
        abacus.write_input(os.path.join(output_dir, "../INPUT"), incar)
        cwd = os.getcwd()
        os.chdir(output_dir)
        if not os.path.islink("INPUT"):
            os.symlink("../INPUT", "INPUT")
        elif not "../INPUT" == os.readlink("INPUT"):
            os.remove("INPUT")
            os.symlink("../INPUT", "INPUT")
        os.chdir(cwd)

        if "kspacing" in incar:
            if isinstance(incar["kspacing"],str):
                kspacing = [float(i) for i in incar["kspacing"].split()]
            elif isinstance(incar["kspacing"],(int,float)):
                kspacing = [incar["kspacing"]]
            else:
                kspacing = incar["kspacing"]
            if len(kspacing) == 1:
                kspacing = 3 * kspacing
                
            if os.path.isfile(os.path.join(output_dir, "STRU")):
                kpt = abacus.make_kspacing_kpt(
                    os.path.join(output_dir, "STRU"), kspacing
                )
                kpt += [0, 0, 0]
            else:
                kpt = [1, 1, 1, 0, 0, 0]
        elif "K_POINTS" in cal_setting:
            kpt = cal_setting["K_POINTS"]
        else:
            mess = "K point information is not defined\n"
            mess += "You can set key word 'kspacing' (unit in 1/bohr) as one or three float value in INPUT\n"
            mess += "or set key word 'K_POINTS' as a list in 'cal_setting', e.g. [1,2,3,0,0,0]\n"
            raise RuntimeError(mess)
        abacus.write_kpt(os.path.join(output_dir, "KPT"), kpt)

    def compute(self, output_dir):
        if not os.path.isfile(os.path.join(output_dir, "INPUT")):
            dlog.warning("cannot find INPUT in " + output_dir + " skip")
            return None
        ls = LabeledSystem(output_dir, fmt="abacus/relax")
        outcar_dict = ls.as_dict()
        return outcar_dict

    def forward_files(self, property_type="relaxation"):
        return ["INPUT", "STRU", "KPT", "pp_orb"]

    def forward_common_files(self, property_type="relaxation"):
        return []

    def backward_files(self, property_type="relaxation"):
        return []
