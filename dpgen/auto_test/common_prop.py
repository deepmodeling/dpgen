from distutils.version import LooseVersion
import glob
import os
import warnings
from multiprocessing import Pool
import dpgen.auto_test.lib.util as util
from dpgen import dlog
from dpgen.util import sepline
from dpgen.auto_test.EOS import EOS
from dpgen.auto_test.Elastic import Elastic
from dpgen.auto_test.Interstitial import Interstitial
from dpgen.auto_test.Surface import Surface
from dpgen.auto_test.Vacancy import Vacancy
from dpgen.auto_test.calculator import make_calculator
from dpgen.dispatcher.Dispatcher import make_dispatcher
from dpgen.dispatcher.Dispatcher import make_submission
from dpgen.remote.decide_machine import convert_mdata
lammps_task_type = ['deepmd', 'meam', 'eam_fs', 'eam_alloy']


def make_property_instance(paramters):
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
    elif prop_type == 'interstitial':
        return Interstitial(paramters)
    elif prop_type == 'surface':
        return Surface(paramters)
    else:
        raise RuntimeError(f'unknown property type {prop_type}')


def make_property(confs,
                  inter_param,
                  property_list):
    # find all POSCARs and their name like mp-xxx
    # ...
    # conf_dirs = glob.glob(confs)
    # conf_dirs.sort()
    conf_dirs = []
    for conf in confs:
        conf_dirs.extend(glob.glob(conf))
    conf_dirs.sort()
    for ii in conf_dirs:
        sepline(ch=ii, screen=True)
        for jj in property_list:
            if jj.get("skip", False):
                continue
            if 'init_from_suffix' and 'output_suffix' in jj:
                do_refine = True
                suffix = jj['output_suffix']
            elif 'reproduce' in jj and jj['reproduce']:
                do_refine = False
                suffix = 'reprod'
            else:
                do_refine = False
                suffix = '00'
            # generate working directory like mp-xxx/eos_00 if jj['type'] == 'eos'
            # handel the exception that the working directory exists
            # ...

            # determine the suffix: from scratch or refine
            # ...

            property_type = jj['type']
            path_to_equi = os.path.join(ii, 'relaxation', 'relax_task')
            path_to_work = os.path.join(ii, property_type + '_' + suffix)

            if os.path.exists(path_to_work):
                dlog.warning('%s already exists' % path_to_work)
            else:
                os.makedirs(path_to_work)

            prop = make_property_instance(jj)
            task_list = prop.make_confs(path_to_work, path_to_equi, do_refine)

            inter_param_prop = inter_param
            if 'cal_setting' in jj and 'overwrite_interaction' in jj['cal_setting']:
                inter_param_prop = jj['cal_setting']['overwrite_interaction']

            for kk in task_list:
                poscar = os.path.join(kk, 'POSCAR')
                inter = make_calculator(inter_param_prop, poscar)
                inter.make_potential_files(kk)
                dlog.debug(prop.task_type())  ### debug
                inter.make_input_file(kk, prop.task_type(), prop.task_param())

            prop.post_process(task_list)  # generate same KPOINTS file for elastic when doing VASP


def run_property(confs,
                 inter_param,
                 property_list,
                 mdata):
    # find all POSCARs and their name like mp-xxx
    # ...
    # conf_dirs = glob.glob(confs)
    # conf_dirs.sort()
    processes = len(property_list)
    pool = Pool(processes=processes)
    print("Submit job via %d processes" % processes)
    conf_dirs = []
    for conf in confs:
        conf_dirs.extend(glob.glob(conf))
    conf_dirs.sort()
    task_list = []
    work_path_list = []
    for ii in conf_dirs:
        sepline(ch=ii, screen=True)
        for jj in property_list:
            # determine the suffix: from scratch or refine
            # ...
            if jj.get("skip", False):
                continue
            if 'init_from_suffix' and 'output_suffix' in jj:
                suffix = jj['output_suffix']
            elif 'reproduce' in jj and jj['reproduce']:
                suffix = 'reprod'
            else:
                suffix = '00'

            property_type = jj['type']
            path_to_work = os.path.abspath(os.path.join(ii, property_type + '_' + suffix))

            work_path_list.append(path_to_work)
            tmp_task_list = glob.glob(os.path.join(path_to_work, 'task.[0-9]*[0-9]'))
            tmp_task_list.sort()
            task_list.append(tmp_task_list)

            inter_param_prop = inter_param
            if 'cal_setting' in jj and 'overwrite_interaction' in jj['cal_setting']:
                inter_param_prop = jj['cal_setting']['overwrite_interaction']

            # dispatch the tasks
            # POSCAR here is useless
            virtual_calculator = make_calculator(inter_param_prop, "POSCAR")
            forward_files = virtual_calculator.forward_files(property_type)
            forward_common_files = virtual_calculator.forward_common_files(property_type)
            backward_files = virtual_calculator.backward_files(property_type)
            #    backward_files += logs
            # ...
            inter_type = inter_param_prop['type']
            # vasp
            if inter_type == "vasp":
                mdata = convert_mdata(mdata, ["fp"])
            elif inter_type in lammps_task_type:
                mdata = convert_mdata(mdata, ["model_devi"])
            else:
                raise RuntimeError("unknown task %s, something wrong" % inter_type)

            work_path = path_to_work
            all_task = tmp_task_list
            run_tasks = util.collect_task(all_task, inter_type)
            if len(run_tasks) == 0:
                return
            else:
                ret = pool.apply_async(worker, (work_path,
                                                all_task,
                                                forward_common_files,
                                                forward_files,
                                                backward_files,
                                                mdata,
                                                inter_type,
                                                ))
            # run_tasks = [os.path.basename(ii) for ii in all_task]
            # machine, resources, command, group_size = util.get_machine_info(mdata, inter_type)
            # disp = make_dispatcher(machine, resources, work_path, run_tasks, group_size)
            # disp.run_jobs(resources,
            #               command,
            #               work_path,
            #               run_tasks,
            #               group_size,
            #               forward_common_files,
            #               forward_files,
            #               backward_files,
            #               outlog='outlog',
            #               errlog='errlog')
    pool.close()
    pool.join()
    if ret.successful():
        print('finished')


def worker(work_path,
           all_task,
           forward_common_files,
           forward_files,
           backward_files,
           mdata,
           inter_type):
    run_tasks = [os.path.basename(ii) for ii in all_task]
    machine, resources, command, group_size = util.get_machine_info(mdata, inter_type)
    disp = make_dispatcher(machine, resources, work_path, run_tasks, group_size)
    api_version = mdata.get('api_version', '0.9')
    if LooseVersion(api_version) < LooseVersion('1.0'):
        warnings.warn(f"the dpdispatcher will be updated to new version."
            f"And the interface may be changed. Please check the documents for more details")
        disp.run_jobs(resources,
                  command,
                  work_path,
                  run_tasks,
                  group_size,
                  forward_common_files,
                  forward_files,
                  backward_files,
                  outlog='outlog',
                  errlog='errlog')
    elif LooseVersion(api_version) >= LooseVersion('1.0'):
        submission = make_submission(
                mdata_machine=machine,
                mdata_resource=resources,
                commands=[command],
                work_path=work_path,
                run_tasks=run_tasks,
                group_size=group_size,
                forward_common_files=forward_common_files,
                forward_files=forward_files,
                backward_files=backward_files,
                outlog = 'outlog',
                errlog = 'errlog'
            )
        submission.run_submission()

def post_property(confs,
                  #                  inter_param,
                  property_list):
    # find all POSCARs and their name like mp-xxx
    # ...
    #    task_list = []
    # conf_dirs = glob.glob(confs)
    # conf_dirs.sort()
    conf_dirs = []
    for conf in confs:
        conf_dirs.extend(glob.glob(conf))
    conf_dirs.sort()
    for ii in conf_dirs:
        for jj in property_list:
            # determine the suffix: from scratch or refine
            # ...
            if jj.get("skip", False):
                continue
            if 'init_from_suffix' and 'output_suffix' in jj:
                suffix = jj['output_suffix']
            elif 'reproduce' in jj and jj['reproduce']:
                suffix = 'reprod'
            else:
                suffix = '00'
            property_type = jj['type']
            path_to_work = os.path.join(ii, property_type + '_' + suffix)
            prop = make_property_instance(jj)
            prop.compute(os.path.join(path_to_work, 'result.json'), os.path.join(path_to_work, 'result.out'),
                         path_to_work)
