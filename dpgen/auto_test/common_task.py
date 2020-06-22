import os
import dpgen.auto_test.lib.crys as crys
import glob, warnings, json
import dpgen.auto_test.lib.util as util

from dpgen.auto_test.VASP import VASP
from dpgen.auto_test.Lammps import  Lammps


def make_task(inter_parameter,
              path_to_poscar):
    """
    Make an instance of Task
    """
    inter_type = inter_parameter['type']
    if inter_type == 'vasp':
        return VASP(inter_parameter, path_to_poscar)
    elif inter_type in ['deepmd', 'meam', 'eam_fs', 'eam_alloy']:
        return Lammps(inter_parameter, path_to_poscar)
    else:
        raise RuntimeError(f'unsupported interaction {inter_type}')


