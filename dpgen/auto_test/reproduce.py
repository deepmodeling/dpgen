import glob
import os

import numpy as np
from monty.serialization import loadfn


def make_repro(vasp_lmp_path, init_from_suffix, path_to_work):
    path_to_work = os.path.abspath(path_to_work)
    property_type = path_to_work.split('/')[-1].split('_')[0]
    vasp_lmp_path = os.path.join(vasp_lmp_path, '*', property_type + '_' + init_from_suffix)
    vasp_lmp_path_list = glob.glob(vasp_lmp_path)
    vasp_lmp_path_list.sort()
    cwd = os.getcwd()
    struct_init_name_list = []
    for ii in vasp_lmp_path_list:
        struct_init_name_list.append(ii.split('/')[-2])
    struct_output_name = path_to_work.split('/')[-2]

    assert struct_output_name in struct_init_name_list

    for idx, ii in enumerate(struct_init_name_list):
        if ii == struct_output_name:
            label = idx

    vasp_lmp_path_todo = vasp_lmp_path_list[label]

    if not os.path.exists(vasp_lmp_path_todo):
        raise RuntimeError("please do VASP or LAMMPS calcualtions first")

    vasp_lmp_task_todo = glob.glob(os.path.join(vasp_lmp_path_todo, 'task.[0-9]*[0-9]'))
    assert len(vasp_lmp_task_todo) > 0, "Please do VASP or LAMMPS calculations first"
    vasp_lmp_task_todo.sort()

    task_list = []
    task_num = 0
    for ii in vasp_lmp_task_todo:
        # get frame number
        task_result = loadfn(os.path.join(ii, 'result_task.json'))
        nframe = len(task_result['energies'])
        for jj in range(nframe):
            output_task = os.path.join(path_to_work, 'task.%06d' % task_num)
            task_num += 1
            task_list.append(output_task)
            os.makedirs(output_task, exist_ok=True)
            os.chdir(output_task)
            # clear dir
            for kk in ['INCAR', 'POTCAR', 'POSCAR.orig', 'POSCAR', 'conf.lmp', 'in.lammps']:
                if os.path.exists(kk):
                    os.remove(kk)
            # make conf
            task_result.to('vasp/poscar', 'POSCAR', frame_idx=jj)
    os.chdir(cwd)

    return task_list


def post_repro(vasp_lmp_path, init_from_suffix, all_tasks, ptr_data):
    ptr_data += "Reproduce: Initial_path Init_E(eV/atom)  Reprod_E(eV/atom)  Difference(eV/atom)\n"
    struct_output_name = all_tasks[0].split('/')[-3]
    property_type = all_tasks[0].split('/')[-2].split('_')[0]
    vasp_lmp_path = os.path.join(vasp_lmp_path, '*', property_type + '_' + init_from_suffix)
    vasp_lmp_path_list = glob.glob(vasp_lmp_path)
    vasp_lmp_path_list.sort()
    # cwd = os.getcwd()
    struct_init_name_list = []
    for ii in vasp_lmp_path_list:
        struct_init_name_list.append(ii.split('/')[-2])

    assert struct_output_name in struct_init_name_list

    for idx, ii in enumerate(struct_init_name_list):
        if ii == struct_output_name:
            label = idx

    vasp_lmp_path_todo = vasp_lmp_path_list[label]

    vasp_lmp_task_todo = glob.glob(os.path.join(vasp_lmp_path_todo, 'task.[0-9]*[0-9]'))
    assert len(vasp_lmp_task_todo) > 0, "Please do VASP or LAMMPS calcualtions first"
    vasp_lmp_task_todo.sort()

    idid = 0
    init_ener_tot = []
    output_ener_tot = []
    res_data = {}

    for ii in vasp_lmp_task_todo:
        init_task_result = loadfn(os.path.join(ii, 'result_task.json'))
        nframe = len(init_task_result['energies'])
        # idid += nframe
        natoms = init_task_result['atom_numbs'][0]
        init_ener = init_task_result['energies']
        init_ener_tot.extend(list(init_ener))
        output_ener = []
        for jj in range(idid, idid + nframe):
            output_task_result = loadfn(os.path.join(all_tasks[jj], 'result_task.json'))
            output_epa = output_task_result['energies'] / natoms
            output_ener.append(output_epa)
            output_ener_tot.extend(output_task_result['energies'])

            init_epa = init_ener[jj - idid] / natoms
            ptr_data += '%s %7.3f  %7.3f  %7.3f\n' % (ii, init_epa,
                                                      output_epa, output_epa - init_epa)
        idid += nframe
        output_ener = np.array(output_ener)
        output_ener = np.reshape(output_ener, [-1, 1])
        init_ener = np.reshape(init_ener, [-1, 1]) / natoms
        error_start = 1
        output_ener -= output_ener[-1] - init_ener[-1]
        diff = output_ener - init_ener
        diff = diff[error_start:]
        error = np.linalg.norm(diff) / np.sqrt(np.size(output_ener) - 1)
        res_data[ii] = {'nframes': len(init_ener), 'error': error}

    if not len(init_ener_tot) == len(output_ener_tot):
        raise RuntimeError("reproduce tasks not equal to init")
    #    for ii in range(len(lmp_ener_tot)):
    #        ptr_data += '%7.3f  %7.3f  %7.3f\n' % (vasp_ener_tot[ii], lmp_ener_tot[ii],
    #                                               lmp_ener_tot[ii] - vasp_ener_tot[ii])
    return res_data, ptr_data
