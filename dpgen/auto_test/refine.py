import os,re,glob

def make_refine(init_from_suffix, output_suffix, path_to_work, task_num):
    init_from = re.sub(output_suffix, init_from_suffix, path_to_work)
    if not os.path.exists(init_from):
        raise FileNotFoundError("the init_from directory not existed for refine")

    output = path_to_work
    init_from_task_tot = glob.glob(os.path.join(init_from, 'task.[0-9]*[0-9]'))

    if not len(init_from_task_tot) == task_num:
        raise RuntimeError("refine tasks not equal to init")

    task_list = []
    for ii in range(task_num):
        output_task = os.path.join(output, 'task.%06d' % ii)
        os.makedirs(output_task, exist_ok=True)
        os.chdir(output_task)
        for jj in ['INCAR', 'POTCAR', 'POSCAR.orig', 'POSCAR', 'conf.lmp', 'in.lammps']:
            if os.path.exists(jj):
                os.remove(jj)
        task_list.append(output_task)
        init_from_task = os.path.join(init_from, 'task.%06d' % ii)
        if not os.path.exists(init_from_task):
            raise FileNotFoundError("the init_from_task directory not existed for refine")
        contcar = os.path.join(init_from_task, 'CONTCAR')
        init_poscar = os.path.join(init_from_task, 'POSCAR')
        if os.path.exists(contcar):
            os.symlink(os.path.relpath(contcar), 'POSCAR')
        elif os.path.exists(init_poscar):
            os.symlink(os.path.relpath(init_poscar), 'POSCAR')
        else:
            raise FileNotFoundError("no CONTCAR or POSCAR in the init_from directory")

    return task_list
