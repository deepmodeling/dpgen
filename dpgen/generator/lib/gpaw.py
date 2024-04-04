#!/usr/bin/env python3

### This file may not need

import random

import dpdata
import numpy as np
import scipy.constants as pc
from packaging.version import Version


def make_gpaw_input_aimd(ensemble,
                         conf_file,
                         graphs,
                         nsteps,
                         dt,
                         neidelay,
                         trj_freq,
                         mass_map,
                         temp,
                         jdata,
                         tau_t=0.1,
                         pres=None,
                         tau_p=0.5,
                         pka_e=None,
                         ele_temp_f=None,
                         ele_temp_a=None,
                         max_seed=1000000,
                         nopbc=False,
                         deepmd_version="0.1",
                         nbeads=None,):
    if (ele_temp_f is not None or ele_temp_a is not None) and Version(
        deepmd_version
    ) < Version("1"):
        raise RuntimeError(
            "the electron temperature is only supported by deepmd-kit >= 1.0.0, please upgrade your deepmd-kit"
        )
    if ele_temp_f is not None and ele_temp_a is not None:
        raise RuntimeError(
            "the frame style ele_temp and atom style ele_temp should not be set at the same time"
        )
    ret = "variable        NSTEPS          equal %d\n" % nsteps
    if nbeads is not None:
        if nbeads <= 0:
            raise ValueError(
                "The number of beads should be positive. Check your nbeads setting."
            )
        power = 1
        while power < nbeads:
            power *= 10
        ret += "variable        ibead           uloop %d pad\n" % (power - 1)
    if nbeads is not None:
        ret += "atom_modify        map yes\n"
    ret += "variable        THERMO_FREQ     equal %d\n" % trj_freq
    ret += "variable        DUMP_FREQ       equal %d\n" % trj_freq
    ret += "variable        TEMP            equal %f\n" % temp
    if nbeads is not None:
        ret += "variable        TEMP_NBEADS            equal %f\n" % (temp * nbeads)
    if ele_temp_f is not None:
        ret += "variable        ELE_TEMP        equal %f\n" % ele_temp_f
    if ele_temp_a is not None:
        ret += "variable        ELE_TEMP        equal %f\n" % ele_temp_a
    ret += "variable        PRES            equal %f\n" % pres
    ret += "variable        TAU_T           equal %f\n" % tau_t
    ret += "variable        TAU_P           equal %f\n" % tau_p
    ret += "\n"
    ret += "units           metal\n"
    if nopbc:
        ret += "boundary        f f f\n"
    else:
        ret += "boundary        p p p\n"
    ret += "atom_style      atomic\n"
    ret += "\n"
    ret += "neighbor        1.0 bin\n"
    if neidelay is not None:
        ret += "neigh_modify    delay %d\n" % neidelay
    ret += "\n"
    ret += "box          tilt large\n"
    if nbeads is None:
        ret += (
            'if "${restart} > 0" then "read_restart dpgen.restart.*" else "read_data %s"\n'
            % conf_file
        )
    else:
        ret += (
            'if "${restart} > 0" then "read_restart dpgen.restart${ibead}.*" else "read_data %s"\n'
            % conf_file
        )
    ret += "change_box   all triclinic\n"
    for jj in range(len(mass_map)):
        ret += "mass            %d %f\n" % (jj + 1, mass_map[jj])
    graph_list = ""
    for ii in graphs:
        graph_list += ii + " "
    if Version(deepmd_version) < Version("1"):
        # 0.x
        ret += "pair_style      deepmd %s ${THERMO_FREQ} model_devi.out\n" % graph_list
    else:
        # 1.x
        keywords = ""
        if jdata.get("use_clusters", False):
            keywords += "atomic "
        if jdata.get("use_relative", False):
            keywords += "relative %s " % jdata["epsilon"]
        if jdata.get("use_relative_v", False):
            keywords += "relative_v %s " % jdata["epsilon_v"]
        if ele_temp_f is not None:
            keywords += "fparam ${ELE_TEMP}"
        if ele_temp_a is not None:
            keywords += "aparam ${ELE_TEMP}"
        if nbeads is None:
            ret += f"pair_style      deepmd {graph_list} out_freq ${{THERMO_FREQ}} out_file model_devi.out {keywords}\n"
        else:
            ret += f"pair_style      deepmd {graph_list} out_freq ${{THERMO_FREQ}} out_file model_devi${{ibead}}.out {keywords}\n"
    ret += "pair_coeff      * *\n"
    ret += "\n"
    ret += "thermo_style    custom step temp pe ke etotal press vol lx ly lz xy xz yz\n"
    ret += "thermo          ${THERMO_FREQ}\n"
    model_devi_merge_traj = jdata.get("model_devi_merge_traj", False)
    if nbeads is None:
        if model_devi_merge_traj is True:
            ret += "dump            1 all custom ${DUMP_FREQ} all.lammpstrj id type x y z fx fy fz\n"
            ret += 'if "${restart} > 0" then "dump_modify     1 append yes"\n'
        else:
            ret += "dump            1 all custom ${DUMP_FREQ} traj/*.lammpstrj id type x y z fx fy fz\n"
    else:
        if model_devi_merge_traj is True:
            ret += "dump            1 all custom ${DUMP_FREQ} all.lammpstrj${ibead} id type x y z fx fy fz\n"
            ret += 'if "${restart} > 0" then "dump_modify     1 append yes"\n'
        else:
            ret += "dump            1 all custom ${DUMP_FREQ} traj/*.lammpstrj${ibead} id type x y z fx fy fz\n"
    if nbeads is None:
        ret += "restart         10000 dpgen.restart\n"
    else:
        ret += "restart         10000 dpgen.restart${ibead}\n"
    ret += "\n"
    if ensemble.split("-")[0] == "npt":
        assert pres is not None
        if nopbc:
            raise RuntimeError("ensemble %s is conflicting with nopbc" % ensemble)
    if nbeads is None:
        if ensemble == "npt" or ensemble == "npt-i" or ensemble == "npt-iso":
            ret += "fix             1 all npt temp ${TEMP} ${TEMP} ${TAU_T} iso ${PRES} ${PRES} ${TAU_P}\n"
        elif ensemble == "npt-a" or ensemble == "npt-aniso":
            ret += "fix             1 all npt temp ${TEMP} ${TEMP} ${TAU_T} aniso ${PRES} ${PRES} ${TAU_P}\n"
        elif ensemble == "npt-t" or ensemble == "npt-tri":
            ret += "fix             1 all npt temp ${TEMP} ${TEMP} ${TAU_T} tri ${PRES} ${PRES} ${TAU_P}\n"
        elif ensemble == "nvt":
            ret += "fix             1 all nvt temp ${TEMP} ${TEMP} ${TAU_T}\n"
        elif ensemble == "nve":
            ret += "fix             1 all nve\n"
        else:
            raise RuntimeError("unknown emsemble " + ensemble)
    else:
        if ensemble == "npt" or ensemble == "npt-i" or ensemble == "npt-iso":
            ret += "fix 1 all pimd/langevin fmmode physical ensemble npt integrator obabo thermostat PILE_L ${ibead} temp ${TEMP} tau ${TAU_T} scale 1.0 barostat BZP iso ${PRES} taup ${TAU_P}\n"
        elif ensemble == "npt-a" or ensemble == "npt-aniso":
            ret += "fix 1 all pimd/langevin fmmode physical ensemble npt integrator obabo thermostat PILE_L ${ibead} temp ${TEMP} tau ${TAU_T} scale 1.0 barostat BZP aniso ${PRES} taup ${TAU_P}\n"
        elif ensemble == "nvt":
            ret += "fix 1 all pimd/langevin fmmode physical ensemble nvt integrator obabo thermostat PILE_L ${ibead} temp ${TEMP} tau ${TAU_T} scale 1.0\n"
        elif ensemble == "nve":
            ret += "fix 1 all pimd/langevin fmmode physical ensemble nve integrator obabo temp ${TEMP}\n"
        else:
            raise RuntimeError(
                "unknown emsemble "
                + ensemble
                + " for fix pimd/langevin\nrefer to https://docs.lammps.org/fix_pimd.html for more information"
            )
    if nopbc:
        ret += "velocity        all zero linear\n"
        ret += "fix             fm all momentum 1 linear 1 1 1\n"
    ret += "\n"
    ret += "timestep        %f\n" % dt
    ret += "run             ${NSTEPS} upto\n"
    return ret
