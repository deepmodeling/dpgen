import argparse
import copy
import dpdata
import glob
import math
import numpy as np
import os
import sys
import shutil
from deepmd.infer import calc_model_devi
from deepmd.infer import DeepPot as DP


def write_model_devi_out(devi, fname):
    assert devi.shape[1] == 8
    #assert devi.shape[1] == 7
    header = '%5s' % 'step'
    for item in 'vf':
        header += '%16s%16s%16s' % (f'max_devi_{item}', f'min_devi_{item}',f'avg_devi_{item}')
    header += '%16s'%str('min_dis')
    np.savetxt(fname,
               devi,
               fmt=['%5d'] + ['%17.6e' for _ in range(7)],
               delimiter='',
               header=header)
    return devi


def Modd(all_models,type_map):

    # Model Devi 

    cwd = os.getcwd()
    graphs = [DP(model) for model in all_models]

    Devis = []
    pcount = 0
    strus_lists = glob.glob(os.path.join(cwd,'*.structures'))
    for num, strus_path in enumerate(strus_lists):

        structures_data = dpdata.System(strus_path,'deepmd/npy',type_map=type_map)

        # every 500 confs in one task dir
        nnum =  structures_data.get_nframes()
        if nnum == 0:
            continue
        else:
            num_per_task = math.ceil(nnum/500)
            

        for temp in range(num_per_task):
            task_name = os.path.join(cwd,'task.%03d.%03d'%(num,temp)) 
            put_poscar = os.path.join(task_name,'traj')
            if not os.path.exists(task_name):
                os.mkdir(task_name)
                os.mkdir(put_poscar)
            else:
                shutil.rmtree(task_name)
                os.mkdir(task_name)
                os.mkdir(put_poscar)
            devis = []
            if (nnum - (temp+1)*500) >= 0:
                temp_sl = range(temp*500,(temp+1)*500)
            else:
                temp_sl = range(temp*500,nnum)
                
            new_index = 0
            for index,frameid in enumerate(temp_sl):
                pdata = structures_data[frameid]
                pdata.to_vasp_poscar(os.path.join(put_poscar,'%s.poscar'%str(index)))
                nopbc = pdata.nopbc
                coord = pdata.data['coords']
                cell  = pdata.data['cells'] if not nopbc else None
                atom_types = pdata.data['atom_types']
                try:
                    devi = calc_model_devi(coord,cell,atom_types,graphs,nopbc=nopbc)
                except TypeError:
                    devi = calc_model_devi(coord,cell,atom_types,graphs)
                # ------------------------------------------------------------------------------------
                # append min-distance in devi list
                dis = pdata.to_ase_structure()[0].get_all_distances(mic=True)
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
            devis = np.vstack(devis)
            write_model_devi_out(devis,os.path.join(task_name, 'model_devi.out'))

    Devis = np.vstack(Devis)
    write_model_devi_out(Devis,os.path.join(cwd,'Model_Devi.out'))

    f = open(os.path.join(os.path.abspath(os.path.join(cwd,os.pardir)),'record.calypso'),'a+')
    f.write('4')
    f.close()

if __name__ == '__main__':

    cwd = os.getcwd()
    model_path = os.path.join(os.path.abspath(os.path.join(cwd,os.pardir)),'gen_stru_analy')
    parser = argparse.ArgumentParser(description='calc model-devi by `all_models` and `type_map`')
    parser.add_argument(
        '--all_models',
        type=str,
        nargs='+',
        default=model_path,
        help='the path of models which will be used to do model-deviation',
    )
    parser.add_argument(
        '--type_map',
        nargs='+',
        help='the type map of models which will be used to do model-deviation',
    )
    args = parser.parse_args()
    #print(vars(args))
    Modd(args.all_models,args.type_map)
    #Modd(sys.argv[1],sys.argv[2])
