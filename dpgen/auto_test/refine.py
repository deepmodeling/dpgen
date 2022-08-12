import glob
import os
import re
import dpgen.auto_test.lib.abacus as abacus

def make_refine(init_from_suffix, output_suffix, path_to_work):
    cwd = os.getcwd()
    init_from = re.sub(output_suffix[::-1], init_from_suffix[::-1], path_to_work[::-1], count=1)[::-1]
    if not os.path.exists(init_from):
        raise FileNotFoundError("the initial directory does not exist for refine")

    output = path_to_work
    init_from_task_tot = glob.glob(os.path.join(init_from, 'task.[0-9]*[0-9]'))

    task_num = len(init_from_task_tot)

    task_list = []
    for ii in range(task_num):
        output_task = os.path.join(output, 'task.%06d' % ii)
        os.makedirs(output_task, exist_ok=True)
        os.chdir(output_task)
        for jj in ['INCAR', 'POTCAR', 'POSCAR.orig', 'POSCAR', 'conf.lmp', 'in.lammps','STRU']:
            if os.path.exists(jj):
                os.remove(jj)
        task_list.append(output_task)
        init_from_task = os.path.join(init_from, 'task.%06d' % ii)
        if not os.path.exists(init_from_task):
            raise FileNotFoundError("the initial task directory does not exist for refine")

        if os.path.isfile(os.path.join(init_from_task, 'INPUT')) and os.path.isfile(os.path.join(init_from_task, 'STRU')):
            #if there has INPUT and STRU files in this path, we believe this is a ABACUS job
            CONTCAR = abacus.final_stru(init_from_task)
            POSCAR = 'STRU'
        else:
            CONTCAR = 'CONTCAR'
            POSCAR = 'POSCAR'
        contcar = os.path.join(init_from_task, CONTCAR)
        init_poscar = os.path.join(init_from_task, POSCAR)
        if os.path.exists(contcar):
            os.symlink(os.path.relpath(contcar), POSCAR)
        elif os.path.exists(init_poscar):
            os.symlink(os.path.relpath(init_poscar), POSCAR)
        else:
            raise FileNotFoundError("no %s or %s in the init_from directory" % (CONTCAR,POSCAR))
    os.chdir(cwd)

    return task_list
