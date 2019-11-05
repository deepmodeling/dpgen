#%%
import os
import argparse
import json
from collections import defaultdict
from itertools import tee

class System(object):
    current_num_of_system = 0
    current_num_of_sub_systems = 0

    @property
    def index_system(self):
        return self._index_system

    @index_system.setter
    def index_system(self,value):
        self._index_system = value
    
    @classmethod
    def register_system(cls):
        cls.current_num_of_system+=1
        return cls.current_num_of_system-1
    
    @classmethod
    def register_sub_system(cls):
        cls.current_num_of_sub_systems+=1
        return cls.current_num_of_sub_systems-1

    def __init__(self, system_prefix=""):
        # print(files_list)
        # if sum(map_relations)>len(files_list):
        #     raise RuntimeError(
        #         "files_list not enough;sum(map_relations):%s>len(files_list):%s, %s" 
        #         % (sum(map_relations),len(files_list),files_list,))
        self.index_system = self.register_system()
        self.sub_system_list = []
        self.system_prefix = system_prefix
        self.current_idx2 = 0
    
    def add_sub_system(self,idx2, files_list):
        idx1 = self.register_sub_system()
        idx2 = self.current_idx2 
        self.sub_system_list.append((idx1, self.index_system, idx2, files_list))
        self.current_idx2 += 1
        
    def get_sub_system(self):
        return self.sub_system_list

    
class Iteration(object):
    current_num_of_itearation = 0
    current_num_of_sub_itearation = 0
        
    @property
    def index_iteration(self):
        return self._index_iteration # pylint: disable=no-member
    
    @index_iteration.setter
    def index_iteration(self, value):
        self._index_sub_iteration = value
    
    @classmethod
    def register_iteration(cls):
        cls.current_num_of_itearation+=1
        return cls.current_num_of_itearation-1
    
    @classmethod
    def register_sub_iteartion(cls):
        cls.current_num_of_sub_itearation +=1
        return cls.current_num_of_sub_itearation-1
        
    def __init__(self, 
                temps,
                nsteps_list=[500, 500, 1000, 1000, 3000, 3000, 6000, 6000],
                sub_iteration_num=8, 
                ensemble='npt', 
                press=[1.0, 10.0, 100.0, 1000.0, 5000.0, 10000.0, 20000.0, 50000.0],
                trj_freq=10):
        if len(nsteps_list) != sub_iteration_num:
            raise RuntimeError(f'{nsteps_list}, {sub_iteration_num}; length does not match')
        self.temps = temps
        self.index_iteration = self.register_iteration()
        self.nsteps_list=nsteps_list
        self.sub_iteration_num=sub_iteration_num
        self.ensemble=ensemble
        self.press = press
        self.trj_freq = trj_freq

    def gen_sub_iter(self, system_list):
        sub_iter_list = []
        for idx2 in range(self.sub_iteration_num):
            iter_dict = {}
            iter_dict['_idx'] = self.register_sub_iteartion()
            iter_dict['ensemble'] = self.ensemble
            iter_dict['nsteps'] = self.nsteps_list[idx2]
            iter_dict['press'] = self.press
            iter_dict['sys_idx'] = [ii[0] for ii in system_list if ii[2]==idx2]
            iter_dict['temps'] = self.temps
            iter_dict['trj_freq'] = self.trj_freq
            sub_iter_list.append(iter_dict)
        return sub_iter_list

def default_map_generator(map_list=[1,1,2,2,2,4,4,4], data_list=None):
    num = 0
    # if len(data_list) < sum(map_list):
    #     raise RuntimeError(f'{data_list} < {map_list};not enough structure to expore, data_list_too_short!')
    if (data_list is None) and ( all(el%10==0 for el in map_list) ):
        for ii in map_list:
            yield [f"{jj:0<5}?" for jj in range(num, num+ii//10)]
            num+=(ii//10)
    elif data_list:
        for ii in map_list:
            yield [data_list[jj] for jj in range(num, num+ii)]
            num += ii
    raise RuntimeError(f"{map_list} length is not enough")
    #   while True:
        # yield [data_list[jj] for jj in range(num, num+ii)]
        # num += ii 

def get_system_list(system_dict,
    map_list=[1,1,2,2,2,4,4,4],
    meta_iter_num=4, 
    sub_iteration_num=8, 
    map_iterator=None,
    file_name="POSCAR"):
    """
    :type map_iterator: Iterable use to generate sys_configs
    :Exmaple [['000000', '000001',], ['00000[2-9]',], ['00001?', '000020',],]
    """
    if sub_iteration_num != len(map_list):
        raise RuntimeError(f"{sub_iteration_num},{map_list};sub_iteration_num does not match the length of map_list")
    
    system_list = []
    for system_prefix,data_list in system_dict.items():
        if map_iterator is None:
            print('12', data_list)
            new_map_iterator = default_map_generator(map_list=map_list, data_list=data_list)
        else:
            origin_one, new_map_iterator = tee(map_iterator) # pylint: disable=unused-variable
        # tee means copy;new_map_generator will become a copy of map_iterator 
        system = System(system_prefix)
        for idx2 in range(sub_iteration_num):
            files_list = [os.path.join(system_prefix, jj) for jj in next(new_map_iterator)]
            system.add_sub_system(idx2=idx2, files_list=files_list)
        system_list.extend(system.get_sub_system())
    return system_list

def scan_files(scan_dir="./" ,file_name="POSCAR", min_allow_files_num=20):
    # will return
    # files_list=[]
    system_dict = defaultdict(list)
    for ii in os.walk(scan_dir):
        if file_name in ii[2]:
            system_prefix = os.path.dirname(ii[0])
            system_suffix = os.path.basename(ii[0])
            system_dict[system_prefix].append(os.path.join(system_suffix, file_name))
    for k,v in list(system_dict.items()):
        if len(v) < min_allow_files_num:
            del system_dict[k]
    return system_dict

# def gen_
    
def default_temps_generator(melt_point, temps_intervel=0.1, num_temps=5):
    temps = [50, ]
    last_temp = 0 
    for ii in range(num_temps-1): # pylint: disable=unused-variable
        last_temp = last_temp + temps_intervel*melt_point
        temps.append(last_temp)
    yield temps
    while True:
        temps = []
        for ii in range(num_temps):
            last_temp = last_temp + temps_intervel*melt_point
            temps.append(last_temp)
        yield temps

def get_model_devi_jobs(melt_point,
    system_list,
    nsteps_list=[500, 500, 1000, 1000, 3000, 3000, 6000, 6000],
    press=[1.0, 10.0, 100.0, 1000.0, 5000.0, 10000.0, 20000.0, 50000.0],
    meta_iter_num=4, 
    sub_iteration_num=8,
    temps_iterator=None,
    ensemble="npt",
    trj_freq=10,
    temps_intervel=0.1,
    num_temps=5):

    if temps_iterator is None:
        temps_iterator = default_temps_generator(melt_point=melt_point, 
            temps_intervel=temps_intervel, num_temps=num_temps)

    if len(nsteps_list) != sub_iteration_num:
        raise RuntimeError(f"{nsteps_list}, {sub_iteration_num};length do not match!")
    model_devi_jobs =[]
    for ii in range(meta_iter_num): # pylint: disable=unused-variable
        temps = next(temps_iterator)
        meta_iter = Iteration(temps=temps,
            nsteps_list=nsteps_list,
            sub_iteration_num=sub_iteration_num,
            ensemble=ensemble,
            press=press,
            trj_freq=trj_freq)
        model_devi_jobs.extend(meta_iter.gen_sub_iter(system_list))
    return model_devi_jobs

def get_sys_configs(system_list):
    sys_configs=[[] for ii in system_list]
    for t in system_list:
        sys_configs[t[0]]=t[3]
    return sys_configs

def get_init_data_sys(scan_dir='./', init_file_name='type.raw'):

    init_data_sys = []
    for t in os.walk(scan_dir):
        if init_file_name in t[2]:
            init_data_sys.append(t[0])
        else: 
            pass
    return init_data_sys


def get_basic_param_json(melt_point,
    out_param_filename='param_basic.json',
    scan_dir="./", 
    file_name='POSCAR',
    init_file_name='type.raw',
    min_allow_files_num=16,
    map_list=[1,1,2,2,2,4,4,4],
    meta_iter_num=4,
    sub_iteration_num=8,
    map_iterator=None,
    nsteps_list=[500, 500, 1000, 1000, 3000, 3000, 6000, 6000],
    press=[1.0, 10.0, 100.0, 1000.0, 5000.0, 10000.0, 20000.0, 50000.0],
    temps_iterator=None,
    ensemble="npt",
    trj_freq=10,
    temps_intervel=0.1,
    num_temps=5,):

    init_data_sys = get_init_data_sys(scan_dir=scan_dir, init_file_name=init_file_name)
    print(f"length of init_data_sys: {len(init_data_sys)} {init_data_sys}")
    system_dict = scan_files(scan_dir, file_name, min_allow_files_num)
    print(f"num of different systems: {len(system_dict)}")
    system_list =get_system_list(system_dict,
        map_list=map_list, 
        meta_iter_num=meta_iter_num, 
        sub_iteration_num=sub_iteration_num, 
        map_iterator=map_iterator,
        file_name=file_name)

    sys_configs = get_sys_configs(system_list)
    print(f"length of sys_configs: {len(sys_configs)}")
    model_devi_jobs = get_model_devi_jobs(melt_point=melt_point,
        system_list=system_list,
        nsteps_list=nsteps_list,
        press=press,
        meta_iter_num=meta_iter_num, 
        sub_iteration_num=sub_iteration_num, 
        temps_iterator=temps_iterator,
        ensemble=ensemble,
        trj_freq=trj_freq,
        temps_intervel=temps_intervel,
        num_temps=num_temps)
    param_dict={
        'init_data_sys': init_data_sys,
        'sys_configs':sys_configs,
        'model_devi_jobs':model_devi_jobs
    }
    with open(out_param_filename, 'w') as p:
        json.dump(param_dict, p, indent=4)

    return param_dict 
def _main():
    parser = argparse.ArgumentParser(description='Collect data from inputs and generate basic param.json')
    parser.add_argument("melt_point", type=float, help="melt_point")
    # parser.addparser.add_argument("JOB_DIR", type=str, help="the directory of the DP-GEN job")
    args = parser.parse_args()
    get_basic_param_json(melt_point=args.melt_point)
  
if __name__=='__main__':
    _main()

def auto_gen_param(args):
    if args.PARAM:
        with open(args.PARAM) as p:
            j = json.load(p)
        melt_point = j['melt_point']
        print('param_basic.json', get_basic_param_json(melt_point=melt_point))
    else:
        raise RuntimeError('must provide melt point or PARAM')

#%%
