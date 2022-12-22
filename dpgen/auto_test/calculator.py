from dpgen.auto_test.VASP import VASP
from dpgen.auto_test.ABACUS import ABACUS
from dpgen.auto_test.Lammps import Lammps


def make_calculator(inter_parameter,
                    path_to_poscar):
    """
    Make an instance of Task
    """
    inter_type = inter_parameter['type']
    if inter_type == 'vasp':
        return VASP(inter_parameter, path_to_poscar)
    elif inter_type == 'abacus':
        return ABACUS(inter_parameter, path_to_poscar)
    elif inter_type in ['deepmd', 'meam', 'eam_fs', 'eam_alloy']:
        return Lammps(inter_parameter, path_to_poscar)
    #    if inter_type == 'siesta':
    #        return Siesta(inter_parameter, path_to_poscar)
    #        pass
    else:
        raise RuntimeError(f'unsupported interaction {inter_type}')
