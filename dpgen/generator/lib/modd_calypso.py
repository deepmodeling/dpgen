"""
calypso model devi:
       1. GenStructures
       2. Analysis
       3. Modd
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
from ase.io.vasp import write_vasp
from ase.io.trajectory import Trajectory
from deepmd.infer import calc_model_devi
from deepmd.infer import DeepPot as DP
from dpgen import dlog
from dpgen.generator.lib.utils import make_iter_name
from dpgen.generator.lib.calypso import write_model_devi_out
from dpgen.generator.lib.parse_calypso import _parse_calypso_input,_parse_calypso_dis_mtx
from dpgen.dispatcher.Dispatcher import Dispatcher, _split_tasks, make_dispatcher, make_submission

train_name = '00.train'
model_devi_name = '01.model_devi'
fp_name = '02.fp'
calypso_run_opt_name = 'gen_stru_analy'
calypso_model_devi_name = 'model_devi_results'

def GenStructures(iter_index,jdata,mdata):

    # run calypso
    # vsc means generate elemental, binary and ternary at the same time
    vsc = jdata.get('vsc',False)  # take CALYPSO as confs generator

    model_devi_group_size = mdata['model_devi_group_size']
    model_devi_resources = mdata['model_devi_resources']


    iter_name = make_iter_name(iter_index)
    work_path = os.path.join(iter_name, model_devi_name)
    assert(os.path.isdir(work_path))
    # 
    calypso_run_opt_path = os.path.join(work_path,calypso_run_opt_name)
    calypso_model_devi_path = os.path.join(work_path,calypso_model_devi_name)
    #

    calypso_path = jdata.get('calypso_path')
    #calypso_input_path = jdata.get('calypso_input_path')
    popsize = int(_parse_calypso_input('PopSize',calypso_run_opt_path))
    maxstep = int(_parse_calypso_input('MaxStep',calypso_run_opt_path))
    
    all_models = glob.glob(os.path.join(calypso_run_opt_path, 'graph*pb'))
    model_names = [os.path.basename(ii) for ii in all_models]

    deepmdkit_python = mdata.get('deepmdkit_python')
    command = "%s run_opt.py %s 1>> model_devi.log 2>> model_devi.log" % (deepmdkit_python,os.path.abspath(calypso_run_opt_path))
    command += "  ||  %s check_outcar.py %s " % (deepmdkit_python,os.path.abspath(calypso_run_opt_path))
    commands = [command]

    cwd = os.getcwd()
    os.chdir(calypso_run_opt_path)

    forward_files = ['POSCAR', 'run_opt.py','check_outcar.py']
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
            dlog.info('$$$$$$$$$$$$$$$$$$$$ step %s $$$$$$$$$$$$$$$$$'%ii)
            if ii == maxstep :  
                dlog.info('##################### counting OUTCAR #####################')
                while True:
                    if len(glob.glob('OUTCAR_*')) == popsize:
                        break
                dlog.info('##################### done counting OUTCAR #####################')
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
                shutil.copyfile('run_opt.py',os.path.join('task.%03d'%pop,'run_opt.py'))
                shutil.copyfile('check_outcar.py',os.path.join('task.%03d'%pop,'check_outcar.py'))
                shutil.copyfile('POSCAR_%s'%str(pop-ii*int(popsize)+1),os.path.join('task.%03d'%(pop),'POSCAR'))
            #for iii in range(1,popsize+1): 
            #    shutil.copyfile('POSCAR_%s'%str(iii),os.path.join('task.%03d'%(iii-1),'POSCAR'))
            dlog.info('##################### copy file done #####################')

            all_task = glob.glob( "task.*")
            all_task.sort()

            run_tasks_ = all_task

            run_tasks = [os.path.basename(ii) for ii in run_tasks_]

            dlog.info('$$$$$$$$$$$$$$$ dispatcher running $$$$$$$$$$$$$$$$$$')
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
            dlog.info('$$$$$$$$$$$$$$$ dispatcher $$$$$$$$$$$$$$$$$$')
           
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


            dlog.info('##################### copy file back done #####################')

            os.rename('jr.json','jr_%s.json'%(str(ii)))

            tlist = glob.glob('task.*')
            for t in tlist:
                shutil.rmtree(t)

    else:
        # --------------------------------------------------------------
        component = ['Mg','Al','Cu','MgAl','MgCu','AlCu','MgAlCu']
        
        for idx,com in enumerate(component):
            pwd = os.getcwd()
            os.mkdir(str(idx))
            shutil.copyfile(os.path.join(cwd,'calypso_input','input.dat.%s'%com),os.path.join(str(idx),'input.dat'))
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
            shutil.copyfile('run_opt.py',os.path.join('task.%04d'%(idx+1),'run_opt.py'))
            shutil.copyfile('check_outcar.py',os.path.join('task.%04d'%(idx+1),'check_outcar.py'))
            shutil.copyfile('POSCAR_%s'%str(idx+1),os.path.join('task.%04d'%(idx+1),'POSCAR'))

        dlog.info('##################### copy file done #####################')

        all_task = glob.glob( "task.*")
        all_task.sort()

        run_tasks_ = all_task

        run_tasks = [os.path.basename(ii) for ii in run_tasks_]

        dlog.info('$$$$$$$$$$$$$$$ dispatcher running $$$$$$$$$$$$$$$$$$')
        dispatcher=make_dispatcher(mdata['model_devi_machine'],mdata['model_devi_resources'],'./', run_tasks, model_devi_group_size)
        #print(dispatcher)
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
        dlog.info('$$$$$$$$$$$$$$$ dispatcher $$$$$$$$$$$$$$$$$$')
        
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

        dlog.info('##################### copy file back done #####################')

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

def Analysis(iter_index,jdata,calypso_run_opt_path,calypso_model_devi_path):

    # Analysis


    dlog.info('$$$$$$$$$$$$$$$ Analysis Started  $$$$$$$$$$$$$$$$$$')

    ms = dpdata.MultiSystems()

    cwd = os.getcwd()
    iter_name = make_iter_name(iter_index)
    work_path = os.path.join(iter_name, model_devi_name)
    result_path = os.path.join(calypso_run_opt_path,'results')

    # trajs to be model devi
    traj_path = os.path.join(calypso_run_opt_path,'traj')
    traj_list = glob.glob(traj_path+'/*.traj')

    dlog.info('len(traj_list) %s'%str(len(traj_list)))

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

    dlog.info('traj_num %s'%str(len(traj_pos_list)))
    dlog.info('total_traj_num %s'%str(record_traj_num))

    for npos in traj_pos_list:
        try:
            ms.append(dpdata.System(npos, type_map = jdata['type_map']))
        except Exception as e:
            dlog.info(npos,'failed : ',e)

    if len(ms) == 0:
        print('too little confs, ')
        return

    if os.path.exists(os.path.join(result_path,'deepmd')):
        shutil.rmtree(os.path.join(result_path,'deepmd'))
    ms.to_deepmd_raw(os.path.join(result_path,'deepmd'))
    ms.to_deepmd_npy(os.path.join(result_path,'deepmd'))
    

    split_lists = glob.glob(os.path.join(result_path,'deepmd','*'))
    for i,split_list in enumerate(split_lists):
        ss = dpdata.System(split_list,fmt='deepmd')
        for j in range(ss.get_nframes()):
            ss.to('vasp/poscar',os.path.join(split_list,'%03d.%03d.poscar'%(i,j)),frame_idx=j)
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


def Modd(iter_index,calypso_model_devi_path,all_models,jdata):

    # Model Devi 


    dlog.info('$$$$$$$$$$$$$$$ Model Devi Start    $$$$$$$$$$$$$$$$$$')

    iter_name = make_iter_name(iter_index)
    work_path = os.path.join(iter_name, model_devi_name)

    cwd = os.getcwd()
    Devis = []
    pcount = 0
    strus_lists = glob.glob(os.path.join(calypso_model_devi_path,'*.structures'))
    for num, strus_path in enumerate(strus_lists):
        os.chdir(strus_path)
        os.system('rename .vasp .poscar *')
        os.chdir(cwd)

        structure_list = glob.glob(os.path.join(strus_path,'*.poscar'))

        # every 500 confs in one task dir
        if len(structure_list) == 0:
            continue
        else:
            num_per_task = math.ceil(len(structure_list)/500)
        graphs = [DP(model) for model in all_models]
        for temp in range(num_per_task):
            task_name = os.path.join(calypso_model_devi_path,'task.%03d.%03d'%(num,temp)) 
            put_poscar = os.path.join(task_name,'traj')
            if not os.path.exists(task_name):
                os.mkdir(task_name)
                os.mkdir(put_poscar)
            else:
                shutil.rmtree(task_name)
                os.mkdir(task_name)
                os.mkdir(put_poscar)
            devis = []
            try:
                temp_sl = structure_list[temp*500:(temp+1)*500]
            except Exception as err:
                dlog.info('err %s'%str(err))
                temp_sl = structure_list[temp*500:]
                
            new_index = 0
            for index,sname in enumerate(temp_sl):
                shutil.copyfile(sname,os.path.join(put_poscar,'%s.poscar'%str(index)))
                os.chdir(put_poscar)
                pdata = dpdata.System('%s.poscar'%str(index),type_map = jdata['type_map'])
                nopbc = pdata.nopbc
                coord = pdata.data['coords']
                cell  = pdata.data['cells']
                atom_types = pdata.data['atom_types']
                devi = calc_model_devi(coord,cell,atom_types,graphs,nopbc=nopbc)
                # ------------------------------------------------------------------------------------
                # append min-distance in devi list
                dis_temp = dpdata.System('%s.poscar'%str(index))
                dis = dis_temp.to_ase_structure()[0].get_all_distances(mic=True)
                row,col = np.diag_indices_from(dis)
                dis[row,col] = 10000
                min_dis = np.nanmin(dis)
                devi = np.append(devi[0],min_dis) 
                t = [devi]
                devi = np.array(t)
                # ------------------------------------------------------------------------------------
                temp_d = copy.deepcopy(devi)
                temp_D = copy.deepcopy(devi)
                devis.append(temp_d)
                Devis.append(temp_D)
                devis[index][0][0]  = np.array(index)
                Devis[pcount][0][0] = np.array(pcount)
                pcount += 1
                new_index += 1
                os.chdir(cwd)
            os.chdir(task_name)
            devis = np.vstack(devis)
            write_model_devi_out(devis,'model_devi.out')
            os.chdir(cwd)

    os.chdir(calypso_model_devi_path)
    Devis = np.vstack(Devis)
    write_model_devi_out(Devis,'Model_Devi.out')
    os.chdir(cwd)

    os.chdir(work_path)
    f = open('record.calypso','a+')
    f.write('4\n')
    f.close()
    os.chdir(cwd)

