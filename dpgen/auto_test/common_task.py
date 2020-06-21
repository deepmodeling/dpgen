import os
import dpgen.auto_test.lib.crys as crys
import glob, warnings, json
import dpgen.auto_test.lib.util as util

from dpgen.auto_test.VASP import VASP
from dpgen.auto_test.DEEPMD_LMP import DEEPMD_LMP
from dpgen.auto_test.MEAM_LMP import MEAM_LMP
from dpgen.auto_test.EAM_FS_LMP import EAM_FS_LMP
from dpgen.auto_test.EAM_ALLOY_LMP import EAM_ALLOY_LMP


def make_task(inter_parameter,
              path_to_poscar):
    """
    Make an instance of Task
    """
    inter_type = inter_parameter['type']
    if inter_type == 'vasp':
        return VASP(inter_parameter, path_to_poscar)
    elif inter_type == 'deepmd':
        return DEEPMD_LMP(inter_parameter, path_to_poscar)
    elif inter_type == 'meam':
        return MEAM_LMP(inter_parameter, path_to_poscar)
    elif inter_type == 'eam_fs':
        return EAM_FS_LMP(inter_parameter, path_to_poscar)
    elif inter_type == 'eam_alloy':
        return EAM_ALLOY_LMP(inter_parameter, path_to_poscar)
    else:
        raise RuntimeError(f'unknown interaction {inter_type}')


def make_task_trans_files(inter_parameter):
    """
    Make the forward and backward file of an Task
    """
    inter_type = inter_parameter['type']
    if inter_type == 'vasp':
        return VASP.forward_files, VASP.forward_common_files, VASP.backward_files
    elif inter_type == 'deepmd':
        return DEEPMD_LMP.forward_files, DEEPMD_LMP.forward_common_files, DEEPMD_LMP.backward_files
    elif inter_type == 'meam':
        return MEAM_LMP.forward_files, MEAM_LMP.forward_common_files, MEAM_LMP.backward_files
    elif inter_type == 'eam_fs':
        return EAM_FS_LMP.forward_files, EAM_FS_LMP.forward_common_files, EAM_FS_LMP.backward_files
    elif inter_type == 'eam_alloy':
        return EAM_ALLOY_LMP.forward_files, EAM_ALLOY_LMP.forward_common_files, EAM_ALLOY_LMP.backward_files
    else:
        raise RuntimeError(f'unknown interaction {inter_type}')


