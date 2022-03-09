def write_modd():

    ret ="import argparse\n"
    ret+="import copy\n"
    ret+="import dpdata\n"
    ret+="import glob\n"
    ret+="import math\n"
    ret+="import numpy as np\n"
    ret+="import os\n"
    ret+="import sys\n"
    ret+="import shutil\n"
    ret+="\n"
    ret+="\n"
    ret+="def write_model_devi_out(devi, fname):\n"
    ret+="    assert devi.shape[1] == 8\n"
    ret+="    #assert devi.shape[1] == 7\n"
    ret+="    header = '%5s' % 'step'\n"
    ret+="    for item in 'vf':\n"
    ret+="        header += '%16s%16s%16s' % (f'max_devi_{item}', f'min_devi_{item}',f'avg_devi_{item}')\n"
    ret+="    header += '%16s'%str('min_dis')\n"
    ret+="    np.savetxt(fname,\n"
    ret+="               devi,\n"
    ret+="               fmt=['%5d'] + ['%17.6e' for _ in range(7)],\n"
    ret+="               delimiter='',\n"
    ret+="               header=header)\n"
    ret+="    return devi\n"
    ret+="\n"
    ret+="\n"
    ret+="def Modd(all_models,type_map):\n"
    ret+="    from deepmd.infer import calc_model_devi\n"
    ret+="    from deepmd.infer import DeepPot as DP\n"
    ret+="\n"
    ret+="    # Model Devi \n"
    ret+="\n"
    ret+="    cwd = os.getcwd()\n"
    ret+="    graphs = [DP(model) for model in all_models]\n"
    ret+="\n"
    ret+="    Devis = []\n"
    ret+="    pcount = 0\n"
    ret+="    strus_lists = glob.glob(os.path.join(cwd,'*.structures'))\n"
    ret+="    for num, strus_path in enumerate(strus_lists):\n"
    ret+="\n"
    ret+="        structures_data = dpdata.System(strus_path,'deepmd/npy',type_map=type_map)\n"
    ret+="\n"
    ret+="        # every 500 confs in one task dir\n"
    ret+="        nnum =  structures_data.get_nframes()\n"
    ret+="        if nnum == 0:\n"
    ret+="            continue\n"
    ret+="        else:\n"
    ret+="            num_per_task = math.ceil(nnum/500)\n"
    ret+="            \n"
    ret+="\n"
    ret+="        for temp in range(num_per_task):\n"
    ret+="            task_name = os.path.join(cwd,'task.%03d.%03d'%(num,temp)) \n"
    ret+="            put_poscar = os.path.join(task_name,'traj')\n"
    ret+="            if not os.path.exists(task_name):\n"
    ret+="                os.mkdir(task_name)\n"
    ret+="                os.mkdir(put_poscar)\n"
    ret+="            else:\n"
    ret+="                shutil.rmtree(task_name)\n"
    ret+="                os.mkdir(task_name)\n"
    ret+="                os.mkdir(put_poscar)\n"
    ret+="            devis = []\n"
    ret+="            if (nnum - (temp+1)*500) >= 0:\n"
    ret+="                temp_sl = range(temp*500,(temp+1)*500)\n"
    ret+="            else:\n"
    ret+="                temp_sl = range(temp*500,nnum)\n"
    ret+="                \n"
    ret+="            new_index = 0\n"
    ret+="            for index,frameid in enumerate(temp_sl):\n"
    ret+="                pdata = structures_data[frameid]\n"
    ret+="                pdata.to_vasp_poscar(os.path.join(put_poscar,'%s.poscar'%str(index)))\n"
    ret+="                nopbc = pdata.nopbc\n"
    ret+="                coord = pdata.data['coords']\n"
    ret+="                cell  = pdata.data['cells']\n"
    ret+="                atom_types = pdata.data['atom_types']\n"
    ret+="                devi = calc_model_devi(coord,cell,atom_types,graphs,nopbc=nopbc)\n"
    ret+="                # ------------------------------------------------------------------------------------\n"
    ret+="                # append min-distance in devi list\n"
    ret+="                dis = pdata.to_ase_structure()[0].get_all_distances(mic=True)\n"
    ret+="                row,col = np.diag_indices_from(dis)\n"
    ret+="                dis[row,col] = 10000\n"
    ret+="                min_dis = np.nanmin(dis)\n"
    ret+="                devi = np.append(devi[0],min_dis)\n "
    ret+="                t = [devi]\n"
    ret+="                devi = np.array(t)\n"
    ret+="                # ------------------------------------------------------------------------------------\n"
    ret+="                temp_d = copy.deepcopy(devi)\n"
    ret+="                temp_D = copy.deepcopy(devi)\n"
    ret+="                devis.append(temp_d)\n"
    ret+="                Devis.append(temp_D)\n"
    ret+="                devis[index][0][0]  = np.array(index)\n"
    ret+="                Devis[pcount][0][0] = np.array(pcount)\n"
    ret+="                pcount += 1\n"
    ret+="                new_index += 1\n"
    ret+="            devis = np.vstack(devis)\n"
    ret+="            write_model_devi_out(devis,os.path.join(task_name, 'model_devi.out'))\n"
    ret+="\n"
    ret+="    Devis = np.vstack(Devis)\n"
    ret+="    write_model_devi_out(Devis,os.path.join(cwd,'Model_Devi.out'))\n"
    ret+="\n"
    ret+="    f = open(os.path.join(os.path.abspath(os.path.join(cwd,os.pardir)),'record.calypso'),'a+')\n"
    ret+="    f.write('4\n')\n"
    ret+="    f.close()\n"
    ret+="\n"
    ret+="if __name__ == '__main__':\n"
    ret+="\n"
    ret+="    cwd = os.getcwd()\n"
    ret+="    model_path = os.path.join(os.path.abspath(os.path.join(cwd,os.pardir)),'gen_stru_analy')\n"
    ret+="    parser = argparse.ArgumentParser(description='calc model-devi by `all_models` and `type_map`')\n"
    ret+="    parser.add_argument(\n"
    ret+="        '--all_models',\n"
    ret+="        type=str,\n"
    ret+="        nargs='+',\n"
    ret+="        default=model_path,\n"
    ret+="        help='the path of models which will be used to do model-deviation',\n"
    ret+="    )\n"
    ret+="    parser.add_argument(\n"
    ret+="        '--type_map',\n"
    ret+="        nargs='+',\n"
    ret+="        help='the type map of models which will be used to do model-deviation',\n"
    ret+="    )\n"
    ret+="    args = parser.parse_args()\n"
    ret+="    print(vars(args))\n"
    ret+="    #Modd(args.all_models,args.type_map)\n"
    ret+="    #Modd(sys.argv[1],sys.argv[2])\n"
    return ret
