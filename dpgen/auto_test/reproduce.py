import glob
import os

import numpy as np
from monty.serialization import loadfn
import dpgen.auto_test.lib.abacus as abacus

def make_repro(inter_param,init_data_path, init_from_suffix, path_to_work, reprod_last_frame=True):
    path_to_work = os.path.abspath(path_to_work)
    property_type = path_to_work.split('/')[-1].split('_')[0]
    init_data_path = os.path.join(init_data_path, '*', property_type + '_' + init_from_suffix)
    init_data_path_list = glob.glob(init_data_path)
    init_data_path_list.sort()
    cwd = os.getcwd()
    struct_init_name_list = []
    for ii in init_data_path_list:
        struct_init_name_list.append(ii.split('/')[-2])
    struct_output_name = path_to_work.split('/')[-2]

    assert struct_output_name in struct_init_name_list

    for idx, ii in enumerate(struct_init_name_list):
        if ii == struct_output_name:
            label = idx

    init_data_path_todo = init_data_path_list[label]

    init_data_task_todo = glob.glob(os.path.join(init_data_path_todo, 'task.[0-9]*[0-9]'))
    assert len(init_data_task_todo) > 0, "There is no task in previous calculations path"
    init_data_task_todo.sort()

    task_list = []
    task_num = 0

    if property_type == 'interstitial':
        if os.path.exists(os.path.join(path_to_work, 'element.out')):
            os.remove(os.path.join(path_to_work, 'element.out'))
        fout_element = open(os.path.join(path_to_work, 'element.out'), 'a+')
        fin_element = open(os.path.join(init_data_path_todo, 'element.out'), 'r')

    for ii in init_data_task_todo:
        # get frame number
        task_result = loadfn(os.path.join(ii, 'result_task.json'))
        if reprod_last_frame:
            nframe = 1
        else:
            nframe = len(task_result['energies'])
        if property_type == 'interstitial':
            insert_element = fin_element.readline().split()[0]
        for jj in range(nframe):
            if property_type == 'interstitial':
                print(insert_element, file=fout_element)
            output_task = os.path.join(path_to_work, 'task.%06d' % task_num)
            task_num += 1
            task_list.append(output_task)
            os.makedirs(output_task, exist_ok=True)
            os.chdir(output_task)
            # clear dir
            for kk in ['INCAR', 'POTCAR', 'POSCAR.orig', 'POSCAR', 'conf.lmp', 'in.lammps','STRU']:
                if os.path.exists(kk):
                    os.remove(kk)
            # make conf
            if reprod_last_frame:
                task_result.to('vasp/poscar', 'POSCAR', frame_idx=-1)
            else:
                task_result.to('vasp/poscar', 'POSCAR', frame_idx=jj)
            if inter_param['type'] == 'abacus':
                abacus.poscar2stru("POSCAR",inter_param,"STRU")
                os.remove('POSCAR')

    os.chdir(cwd)

    if property_type == 'interstitial':
        fout_element.close()
        fin_element.close()

    return task_list


def post_repro(init_data_path, init_from_suffix, all_tasks, ptr_data, reprod_last_frame=True):
    ptr_data += "Reproduce: Initial_path Init_E(eV/atom)  Reprod_E(eV/atom)  Difference(eV/atom)\n"
    struct_output_name = all_tasks[0].split('/')[-3]
    property_type = all_tasks[0].split('/')[-2].split('_')[0]
    init_data_path = os.path.join(init_data_path, '*', property_type + '_' + init_from_suffix)
    init_data_path_list = glob.glob(init_data_path)
    init_data_path_list.sort()
    # cwd = os.getcwd()
    struct_init_name_list = []
    for ii in init_data_path_list:
        struct_init_name_list.append(ii.split('/')[-2])

    assert struct_output_name in struct_init_name_list

    for idx, ii in enumerate(struct_init_name_list):
        if ii == struct_output_name:
            label = idx

    init_data_path_todo = init_data_path_list[label]

    init_data_task_todo = glob.glob(os.path.join(init_data_path_todo, 'task.[0-9]*[0-9]'))
    assert len(init_data_task_todo) > 0, "There is no task in previous calculations path"
    init_data_task_todo.sort()

    idid = 0
    init_ener_tot = []
    output_ener_tot = []
    res_data = {}

    for ii in init_data_task_todo:
        init_task_result = loadfn(os.path.join(ii, 'result_task.json'))
        if reprod_last_frame:
            nframe = 1
        else:
            nframe = len(init_task_result['energies'])
        # idid += nframe
        natoms = init_task_result['atom_numbs'][0]
        if reprod_last_frame:
            init_ener = init_task_result['energies'][-1:]
        else:
            init_ener = init_task_result['energies']
        init_ener_tot.extend(list(init_ener))
        output_ener = []
        for jj in range(idid, idid + nframe):
            output_task_result = loadfn(os.path.join(all_tasks[jj], 'result_task.json'))
            output_epa = output_task_result['energies'] / natoms
            output_ener.append(output_epa)
            output_ener_tot.extend(output_task_result['energies'])

            init_epa = init_ener[jj - idid] / natoms
            ptr_data += '%s %7.3f  %7.3f  %7.3f\n' % (ii, init_epa, output_epa, output_epa - init_epa)
        idid += nframe
        output_ener = np.array(output_ener)
        output_ener = np.reshape(output_ener, [-1, 1])
        init_ener = np.reshape(init_ener, [-1, 1]) / natoms
        if reprod_last_frame:
            error_start = 0
        else:
            error_start = 1
            output_ener -= output_ener[-1] - init_ener[-1]
        diff = output_ener - init_ener
        diff = diff[error_start:]
        error = np.linalg.norm(diff) / np.sqrt(np.size(output_ener) - error_start)
        res_data[ii] = {'nframes': len(init_ener), 'error': error}

    if not len(init_ener_tot) == len(output_ener_tot):
        raise RuntimeError("reproduce tasks not equal to init")
    #    for ii in range(len(lmp_ener_tot)):
    #        ptr_data += '%7.3f  %7.3f  %7.3f\n' % (vasp_ener_tot[ii], lmp_ener_tot[ii],
    #                                               lmp_ener_tot[ii] - vasp_ener_tot[ii])
    return res_data, ptr_data
