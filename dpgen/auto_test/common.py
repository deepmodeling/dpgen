import os
from VASP import VASP
from DEEPMD_LMP import DEEPMD_LMP
from MEAM_LMP import MEAM_LMP

def make_interaction(paramters, 
                     path_to_poscar) :
    """
    Make an instance of Interaction
    """
    inter_type = paramters['type']
    if inter_type == 'vasp':
        return VASP(paramters, path_to_poscar)
    elif inter_type == 'deepmd':
        return DEEPMD_LMP(paramters, path_to_poscar)
    elif inter_type == 'meam':
        return MEAM_LMP(paramters, path_to_poscar)
    else:
        raise RuntimeError(f'unknow interaction type {inter_type}')

def make_interaction_files(paramters) :
    """
    Make the forward and backward file of an Interaction
    """
    inter_type = paramters['type']
    if inter_type == 'vasp':
        return VASP.forward_files, VASP.forward_common_files, VASP.backward_files
    elif inter_type == 'deepmd':
        return DEEPMD_LMP.forward_files, DEEPMD_LMP.forward_common_files, DEEPMD_LMP.backward_files
    elif inter_type == 'meam':
        return MEAM_LMP.forward_files, MEAM_LMP.forward_common_files, MEAM_LMP.backward_files
    else:
        raise RuntimeError(f'unknow interaction type {inter_type}')

def make_property(paramters) :
    """
    Make an instance of Property
    """
    prop_type = paramters['type']
    if prop_type == 'eos':
        return EOS(paramters)
    elif prop_type == 'elastic':
        return Elastic(paramters)
    elif prop_type == 'vacancy':
        return Vacancy(paramters)
    else:
        raise RuntimeError(f'unknow property type {prop_type}')

def make_equi(confs, 
              inter_param,
              relax_param):
    # find all POSCARs and their name like mp-xxx
    # ...

    # generate a list of task names like mp-xxx/relaxation
    # ...
    
    # make task directories like mp-xxx/relaxation
    # if mp-xxx/exists then print a warning and exit.
    # ...

    # copy POSCARs to mp-xxx/relaxation
    # ...

    # generate task files
    for ii in task_dirs:
        poscar = os.path.join(ii, 'POSCAR')
        inter = make_interaction(inter_param, poscar)
        inter.make_potential_files(ii)
        inter.make_input_file(ii, 'relaxation', relax_param)


def run_equi(confs, 
             inter_param):
    # find all POSCARs and their name like mp-xxx
    # ...

    # generate a list of task names like mp-xxx/relaxation
    # ...

    # dispatch the tasks
    forward_files, forward_common_files, backward_files = make_interaction_files(inter_param)
    backward_files += logs
    # ...
    pass


def post_equi(confs):
    # find all POSCARs and their name like mp-xxx
    # ...

    # generate a list of task names like mp-xxx/relaxation
    # ...

    # dump the relaxation result.
    for ii in task_dirs:
        poscar = os.path.join(ii, 'POSCAR')
        inter = make_interaction(inter_param, poscar)
        res = inter.post()
        with open('result.json', 'w') as fp:
            json.dump(fp, res, indent=4)
        
    
def make_property(confs, 
                  inter_param,
                  property_list):
    # find all POSCARs and their name like mp-xxx
    # ...
    for ii in conf_dirs:
        for jj in property_list:
            # generate working directory like mp-xxx/eos_00 if jj['type'] == 'eos'
            # handel the exception that the working directory exists
            # ...
            
            # determine the suffix: from scratch or refine
            # ...
            
            property_type = jj['type']
            path_to_equi = os.path.join(ii, 'relaxation')
            path_to_work = os.path.join(ii, property_type + '_' + suffix)
            
            prop = make_property(jj)
            task_list = prop.make_confs(path_to_work, path_to_equi, do_refine)
            
            for kk in task_list:
                poscar = os.path.join(kk, 'POSCAR')                
                inter = make_interaction(inter_param, poscar)
                inter.make_potential_files(kk)
                inter.make_input_file(kk, prop.task_type, prop.task.pararm)


def run_property(confs, 
                 inter_param, 
                 property_list):
    # find all POSCARs and their name like mp-xxx
    # ...
    task_list = []
    for ii in conf_dirs:
        for jj in property_list:
            # determine the suffix: from scratch or refine
            # ...
            property_type = jj['type']
            path_to_work = os.path.join(ii, property_type + '_' + suffix)
           
            tmp_task_list = glob.glob(os.path.join(path_to_work, 'task.[0-9]*[0-9]'))
            tmp_task_list.sort()
            task_list.append(tmp_task_list)
    
    # dispatch the tasks
    forward_files, forward_common_files, backward_files = make_interaction_files(inter_param)
    backward_files += logs
    # ...
    pass

def post_property(confs, 
                  inter_param, 
                  property_list):
    # find all POSCARs and their name like mp-xxx
    # ...
    task_list = []
    for ii in conf_dirs:
        for jj in property_list:
            # determine the suffix: from scratch or refine
            # ...
            property_type = jj['type']
            path_to_work = os.path.join(ii, property_type + '_' + suffix)
            prop = make_property(jj)
            prop.post('result.json', 'result.out', path_to_work)
