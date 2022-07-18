"""
calypso as model devi engine:
       1. gen_structures
       2. analysis
       3. model devi
"""

import copy
import dpdata
import math
import numpy as np
import os
import random
import re
import glob
import shutil
import sys
from ase.io.vasp import write_vasp
from ase.io.trajectory import Trajectory
from pathlib import Path
from itertools import combinations
from distutils.version import LooseVersion
from dpgen import dlog
from dpgen.generator.lib.utils import make_iter_name
from dpgen.generator.lib.parse_calypso import _parse_calypso_input
from dpgen.dispatcher.Dispatcher import make_dispatcher, make_submission

train_name = '00.train'
model_devi_name = '01.model_devi'
fp_name = '02.fp'
calypso_run_opt_name = 'gen_stru_analy'
calypso_model_devi_name = 'model_devi_results'

def gen_structures(iter_index,jdata,mdata):

    # run calypso
    # vsc means generate elemental, binary and ternary at the same time
    vsc = jdata.get('vsc',False)  # take CALYPSO as confs generator

    model_devi_group_size = mdata['model_devi_group_size']
    model_devi_resources = mdata['model_devi_resources']
    api_version = mdata.get('api_version', '0.9')


    iter_name = make_iter_name(iter_index)
    work_path = os.path.join(iter_name, model_devi_name)
    assert(os.path.isdir(work_path))
    # 
    calypso_run_opt_path = os.path.join(work_path,calypso_run_opt_name)
    calypso_model_devi_path = os.path.join(work_path,calypso_model_devi_name)
    #

    calypso_path = mdata.get('model_devi_calypso_path')
    #calypso_input_path = jdata.get('calypso_input_path')
    popsize = int(_parse_calypso_input('PopSize',calypso_run_opt_path))
    maxstep = int(_parse_calypso_input('MaxStep',calypso_run_opt_path))
    
    all_models = glob.glob(os.path.join(calypso_run_opt_path, 'graph*pb'))
    model_names = [os.path.basename(ii) for ii in all_models]

    deepmdkit_python = mdata.get('model_devi_deepmdkit_python')
    command = "%s calypso_run_opt.py %s 1>> model_devi.log 2>> model_devi.log" % (deepmdkit_python,os.path.abspath(calypso_run_opt_path))
    command += "  ||  %s check_outcar.py %s " % (deepmdkit_python,os.path.abspath(calypso_run_opt_path))
    commands = [command]

    cwd = os.getcwd()
    os.chdir(calypso_run_opt_path)

    forward_files = ['POSCAR', 'calypso_run_opt.py','check_outcar.py','input.dat']
    backward_files = ['OUTCAR','CONTCAR','traj.traj','model_devi.log']

    run_calypso = calypso_path+'/calypso.x | tee log'
    if not vsc:
        Lpickup = _parse_calypso_input('PickUp','.')
        PickUpStep =  _parse_calypso_input('PickStep','.')
        if os.path.exists('tag_pickup_%s'%(str(PickUpStep))):
            dlog.info('caution! tag_pickup_%s exists!'%str(PickUpStep))
            Lpickup = 'F'
        if Lpickup == 'T':
            ftag = open('tag_pickup_%s'%(str(PickUpStep)),'w')
            ftag.close()
            os.remove('step')
            fstep = open('step','w')
            fstep.write('%12s'%str(PickUpStep))
            fstep.close()
        else:
            PickUpStep = 1
            try:
                os.mkdir('opt')
            except:
                pass
                #shutil.rmtree('opt')
                #os.mkdir('opt')


        for ii in range(int(PickUpStep)-1,maxstep+1):
            dlog.info('CALYPSO step %s'%ii)
            if ii == maxstep :  
                #while True:
                #    if len(glob.glob('OUTCAR_*')) == popsize:
                #        break
                os.system('%s'%run_calypso)
                break
            # run calypso

            os.system('%s'%(run_calypso))
            
            for pop in range(ii*int(popsize),(ii+1)*int(popsize)):
                try:
                    os.mkdir('task.%03d'%pop)
                except:
                    shutil.rmtree('task.%03d'%pop)
                    os.mkdir('task.%03d'%pop)
                shutil.copyfile('calypso_run_opt.py',os.path.join('task.%03d'%pop,'calypso_run_opt.py'))
                shutil.copyfile('check_outcar.py',os.path.join('task.%03d'%pop,'check_outcar.py'))
                shutil.copyfile('POSCAR_%s'%str(pop-ii*int(popsize)+1),os.path.join('task.%03d'%(pop),'POSCAR'))
                shutil.copyfile('input.dat',os.path.join('task.%03d'%pop,'input.dat'))
            #for iii in range(1,popsize+1): 
            #    shutil.copyfile('POSCAR_%s'%str(iii),os.path.join('task.%03d'%(iii-1),'POSCAR'))

            all_task = glob.glob( "task.*")
            all_task.sort()

            run_tasks_ = all_task

            run_tasks = [os.path.basename(ii) for ii in run_tasks_]

            if LooseVersion(api_version) < LooseVersion('1.0'):
                warnings.warn(f"the dpdispatcher will be updated to new version."
                    f"And the interface may be changed. Please check the documents for more details")
                dispatcher=make_dispatcher(mdata['model_devi_machine'],mdata['model_devi_resources'],'./', run_tasks, model_devi_group_size)
                dispatcher.run_jobs(mdata['model_devi_resources'],
                                commands,
                                './',
                                run_tasks,
                                model_devi_group_size,
                                model_names,
                                forward_files,
                                backward_files,
                                outlog = 'model_devi.log',
                                errlog = 'model_devi.log')
            elif LooseVersion(api_version) >= LooseVersion('1.0'):
                os.chdir(cwd)
                submission = make_submission(
                    mdata['model_devi_machine'],
                    mdata['model_devi_resources'],
                    commands=commands,
                    work_path=calypso_run_opt_path,
                    run_tasks=run_tasks,
                    group_size=model_devi_group_size,
                    forward_common_files=model_names,
                    forward_files=forward_files,
                    backward_files=backward_files,
                    outlog = 'model_devi.log',
                    errlog = 'model_devi.log')
                submission.run_submission()
                os.chdir(calypso_run_opt_path)
           
            sstep = os.path.join('opt',str(ii))
            os.mkdir(sstep)
            if not os.path.exists('traj'):
                os.mkdir('traj')

            for jjj in range(ii*int(popsize),(ii+1)*int(popsize)):
                # to opt directory
                shutil.copyfile('POSCAR_%s'%str(jjj+1-ii*int(popsize)),os.path.join(sstep,'POSCAR_%s'%str(jjj+1-ii*int(popsize))),)
                shutil.copyfile(os.path.join('task.%03d'%(jjj),'OUTCAR'),os.path.join(sstep,'OUTCAR_%s'%str(jjj+1-ii*int(popsize))),)
                shutil.copyfile(os.path.join('task.%03d'%(jjj),'CONTCAR'),os.path.join(sstep,'CONTCAR_%s'%str(jjj+1-ii*int(popsize))),)
                # to run calypso directory
                shutil.copyfile(os.path.join('task.%03d'%(jjj),'OUTCAR'),'OUTCAR_%s'%str(jjj+1-ii*int(popsize)),)
                shutil.copyfile(os.path.join('task.%03d'%(jjj),'CONTCAR'),'CONTCAR_%s'%str(jjj+1-ii*int(popsize)),)
                # to traj
                shutil.copyfile(os.path.join('task.%03d'%(jjj),'traj.traj'),os.path.join('traj','%s.traj'%str(jjj+1)),)

            if LooseVersion(api_version) < LooseVersion('1.0'):
                os.rename('jr.json','jr_%s.json'%(str(ii)))

            tlist = glob.glob('task.*')
            for t in tlist:
                shutil.rmtree(t)

    else:
        # --------------------------------------------------------------
        # TODO(zhenyu) make this code work for other situation 
        type_map = jdata['type_map']
        how_many_spec = len(type_map)
        if how_many_spec == 1:
            dlog.info('vsc mode can not work in one-element situation' )
            sys.exit()

        comp_temp = list(map(list,list(combinations(type_map,1))))
        for hms in range(2,how_many_spec+1):
            comp_temp.extend(list(map(list,list(combinations(type_map,hms)))))  # comp_temp = [['Mg'],['Al'],['Cu'],['Mg','Al'],['Mg','Cu'],['Al','Cu'],['Mg','Al','Cu']]
        
        component = []
        for comp_temp_ in comp_temp:
            component.append(''.join(comp_temp_))     # component = ['Mg','Al','Cu','MgAl','MgCu','AlCu','MgAlCu']

        print(component)
        #sys.exit()
        calypso_input_path = jdata.get('calypso_input_path')
        
        for idx,com in enumerate(component):
            pwd = os.getcwd()
            os.mkdir(str(idx))
            #shutil.copyfile(os.path.join(cwd,'calypso_input','input.dat.%s'%com),os.path.join(str(idx),'input.dat'))
            shutil.copyfile(os.path.join(calypso_input_path,'input.dat.%s'%com),os.path.join(str(idx),'input.dat'))
            os.chdir(str(idx))
            os.system(run_calypso)
            os.chdir(pwd)
        
        name_list = Path('.').glob('*/POSCAR_*')
        for idx,name in enumerate(name_list):
            shutil.copyfile(name,'POSCAR_%s'%(idx+1))
            try:
                os.mkdir('task.%04d'%(idx+1))
            except:
                shutil.rmtree('task.%04d'%(idx+1))
                os.mkdir('task.%04d'%(idx+1))
            shutil.copyfile('calypso_run_opt.py',os.path.join('task.%04d'%(idx+1),'calypso_run_opt.py'))
            shutil.copyfile('check_outcar.py',os.path.join('task.%04d'%(idx+1),'check_outcar.py'))
            shutil.copyfile('POSCAR_%s'%str(idx+1),os.path.join('task.%04d'%(idx+1),'POSCAR'))
            shutil.copyfile('input.dat',os.path.join('task.%04d'%(idx+1),'input.dat'))

        all_task = glob.glob( "task.*")
        all_task.sort()

        run_tasks_ = all_task

        run_tasks = [os.path.basename(ii) for ii in run_tasks_]

        if LooseVersion(api_version) < LooseVersion('1.0'):
            warnings.warn(f"the dpdispatcher will be updated to new version."
                f"And the interface may be changed. Please check the documents for more details")
            dispatcher=make_dispatcher(mdata['model_devi_machine'],mdata['model_devi_resources'],'./', run_tasks, model_devi_group_size)
            dispatcher.run_jobs(mdata['model_devi_resources'],
                            commands,
                            './',
                            run_tasks,
                            model_devi_group_size,
                            model_names,
                            forward_files,
                            backward_files,
                            outlog = 'model_devi.log',
                            errlog = 'model_devi.log')
        elif LooseVersion(api_version) >= LooseVersion('1.0'):
            os.chdir(cwd)
            submission = make_submission(
                mdata['model_devi_machine'],
                mdata['model_devi_resources'],
                commands=commands,
                work_path=calypso_run_opt_path,
                run_tasks=run_tasks,
                group_size=model_devi_group_size,
                forward_common_files=model_names,
                forward_files=forward_files,
                backward_files=backward_files,
                outlog = 'model_devi.log',
                errlog = 'model_devi.log')
            submission.run_submission()
            os.chdir(calypso_run_opt_path)
        
        os.mkdir('opt')
        if not os.path.exists('traj'):
            os.mkdir('traj')
        for jjj in range(len(all_task)):
            # to opt directory
            shutil.copyfile('POSCAR_%s'%str(jjj+1),os.path.join('opt','POSCAR_%s'%str(jjj+1)),)
            shutil.copyfile(os.path.join('task.%04d'%(jjj+1),'OUTCAR'),os.path.join('opt','OUTCAR_%s'%str(jjj+1)),)
            shutil.copyfile(os.path.join('task.%04d'%(jjj+1),'CONTCAR'),os.path.join('opt','CONTCAR_%s'%str(jjj+1)),)
            # to run calypso directory
            shutil.copyfile(os.path.join('task.%04d'%(jjj+1),'OUTCAR'),'OUTCAR_%s'%str(jjj+1),)
            shutil.copyfile(os.path.join('task.%04d'%(jjj+1),'CONTCAR'),'CONTCAR_%s'%str(jjj+1),)
            # to traj
            shutil.copyfile(os.path.join('task.%04d'%(jjj+1),'traj.traj'),os.path.join('traj','%s.traj'%str(jjj+1)),)

        tlist = glob.glob('task.*')
        for t in tlist:
            shutil.rmtree(t)
        # --------------------------------------------------------------

    os.chdir(cwd)
    os.chdir(work_path)
    f = open('record.calypso','a+')
    f.write('2\n')
    f.close()
    os.chdir(cwd)

def analysis(iter_index,jdata,calypso_run_opt_path,calypso_model_devi_path):

    # Analysis


    #dlog.info('$$$$$$$$$$$$$$$ Analysis Started  $$$$$$$$$$$$$$$$$$')

    ms = dpdata.MultiSystems()

    cwd = os.getcwd()
    iter_name = make_iter_name(iter_index)
    work_path = os.path.join(iter_name, model_devi_name)
    result_path = os.path.join(calypso_run_opt_path,'results')

    # trajs to be model devi
    traj_path = os.path.join(calypso_run_opt_path,'traj')
    traj_list = glob.glob(traj_path+'/*.traj')

    #dlog.info('len(traj_list) %s'%str(len(traj_list)))

    # read confs from traj
    record_traj_num = 0
    for traj_name in traj_list:
        traj_num = os.path.basename(traj_name).split('.')[0]
        trajs_origin = Trajectory(traj_name)

        record_traj_num += len(trajs_origin)
        if len(trajs_origin) >= 20 :
           trajs = [trajs_origin[iii] for iii in [4,9,-10,-5,-1]]
        elif 5<=len(trajs_origin)<20:
           trajs = [trajs_origin[random.randint(1,len(trajs_origin)-1)] for iii in range(4)]
           trajs.append(trajs[-1])
        elif 3<= len(trajs_origin) <5:
           trajs = [trajs_origin[round((len(trajs_origin)-1)/2)] ]
           trajs.append(trajs[-1])
        elif len(trajs_origin) == 2:
           trajs = [trajs_origin[0],trajs_origin[-1] ]
        elif len(trajs_origin) == 1:
           trajs = [trajs_origin[0] ]
        else:
           pass
           
        for idx, traj in enumerate(trajs):
            write_vasp(os.path.join(
                traj_path,'%03d.%03d.poscar' % (
                    int(traj_num), int(idx)
                    )
                ),
                traj)
           
    traj_pos_list = glob.glob(traj_path+'/*.poscar')

    #dlog.info('traj_num %s'%str(len(traj_pos_list)))
    #dlog.info('total_traj_num %s'%str(record_traj_num))

    for npos in traj_pos_list:
        try:
            ms.append(dpdata.System(npos, type_map = jdata['type_map']))
        except Exception as e:
            dlog.info(npos,'failed : ',e)

    if len(ms) == 0:
        print('too little confs, ')
        sys.exit()

    if os.path.exists(os.path.join(result_path,'deepmd')):
        shutil.rmtree(os.path.join(result_path,'deepmd'))
    ms.to_deepmd_raw(os.path.join(result_path,'deepmd'))
    ms.to_deepmd_npy(os.path.join(result_path,'deepmd'))
    

    split_lists = glob.glob(os.path.join(result_path,'deepmd','*'))
    for i,split_list in enumerate(split_lists):
        #ss = dpdata.System(split_list,fmt='deepmd')
        #for j in range(ss.get_nframes()):
        #    ss.to('vasp/poscar',os.path.join(split_list,'%03d.%03d.poscar'%(i,j)),frame_idx=j)
        strus_path = os.path.join(calypso_model_devi_path,'%03d.structures'%i)
        if not os.path.exists(strus_path):
            shutil.copytree(split_list,strus_path)
        else:
            shutil.rmtree(strus_path)
            shutil.copytree(split_list,strus_path)

    os.chdir(cwd)
    os.chdir(work_path)
    f = open('record.calypso','a+')
    f.write('3\n')
    f.close()
    os.chdir(cwd)


def run_calypso_model_devi (iter_index,
                    jdata,
                    mdata) :

    dlog.info('start running CALYPSO')


    iter_name = make_iter_name(iter_index)
    work_path = os.path.join(iter_name, model_devi_name)
    assert(os.path.isdir(work_path))

    calypso_run_opt_path = os.path.join(work_path,calypso_run_opt_name)
    calypso_model_devi_path = os.path.join(work_path,calypso_model_devi_name)

    cwd = os.getcwd()

    record_calypso_path = os.path.join(work_path,'record.calypso')
    while True:
        if not os.path.exists(record_calypso_path):
            f = open(record_calypso_path,'w')
            f.write('1\n')
            lines = '1'
            f.close()
        else:
            f = open(record_calypso_path,'r')
            lines = f.readlines()
            f.close()

        if lines[-1].strip().strip('\n') == '1':
            # Gen Structures
            gen_structures(iter_index,jdata,mdata)

        elif lines[-1].strip().strip('\n') == '2':
            # Analysis & to deepmd/raw
            analysis(iter_index,jdata,calypso_run_opt_path,calypso_model_devi_path)

        elif lines[-1].strip().strip('\n') == '3':
            # Model Devi
            _calypso_run_opt_path = os.path.abspath(calypso_run_opt_path)
            all_models = glob.glob(os.path.join(_calypso_run_opt_path, 'graph*pb'))
            cwd = os.getcwd()
            os.chdir(calypso_model_devi_path)
            args = ' '.join(['calypso_run_model_devi.py', '--all_models',' '.join(all_models),'--type_map',' '.join(jdata.get('type_map'))])
            deepmdkit_python = mdata.get('model_devi_deepmdkit_python')
            os.system(f'{deepmdkit_python} {args} ')
            #Modd(iter_index,calypso_model_devi_path,all_models,jdata)
            os.chdir(cwd)

        elif lines[-1].strip().strip('\n') == '4':
            dlog.info('Model Devi is done.')
            return
