import os
import json
from typing import List, TYPE_CHECKING, Tuple, Iterator, Dict, Union

import dpdata
import numpy as np
from dpgen.generator.lib.utils import create_path, expand_idx
from dpgen.generator.model_devi import ModelDeviEngien, Trajectory, Frame
try:
    from gromacs.fileformats.mdp import MDP
except ImportError:
    # JZ: Do not excuse users using other eigiens
    #dlog.info("GromacsWrapper>=0.8.0 is needed for DP-GEN + Gromacs.")
    pass


@ModelDeviEngien.register("gromacs")
class GromacsEngien(ModelDeviEngien):
    def make_input(self, iter_index: int, sys_index: int, directories: Iterator[str], conf_name: str, models: List[str]):
        # JZ: doesn't support dpdata here??
        self._make_model_devi_native_gromacs(iter_index, directories, conf_name, models)
    
    def get_running_parameters(self, work_path: str) -> Dict[str, Union[str, List[str]]]:
        with open (os.path.join(work_path, 'cur_job.json'), 'r') as fp:
            cur_job = json.load (fp)
        gromacs_settings = self.jdata.get("gromacs_settings", {})
        mdp_filename = gromacs_settings.get("mdp_filename", "md.mdp")
        topol_filename = gromacs_settings.get("topol_filename", "processed.top")
        conf_filename = gromacs_settings.get("conf_filename", "conf.gro")
        index_filename = gromacs_settings.get("index_filename", "index.raw")
        deffnm = gromacs_settings.get("deffnm", "deepmd")
        maxwarn = gromacs_settings.get("maxwarn", 1)
        nsteps = cur_job["nsteps"]
        lmp_exec = self.mdata['lmp_command']
        # Angus: lmp_exec name should be changed to model_devi_exec.
        # We should also change make_dispatcher
        # For now, I will use this name for gromacs command

        command = "%s grompp -f %s -p %s -c %s -o %s -maxwarn %d" % (lmp_exec, mdp_filename, topol_filename, conf_filename, deffnm, maxwarn)
        command += "&& %s mdrun -deffnm %s -nsteps %d" %(lmp_exec, deffnm, nsteps) 
        
        forward_files = [mdp_filename, topol_filename, conf_filename, index_filename,  "input.json" ]
        backward_files = ["%s.tpr" % deffnm, "%s.log" %deffnm , 'model_devi.out', 'model_devi.log']
        return {
            "command": command,
            "forward_files": forward_files,
            "backward_files": backward_files,
            "common_files": [],
        }

    def extract_trajectory(self, directory) -> 'Trajectory':
        return GromacsTrajectory(self, directory)

    def _make_model_devi_native_gromacs(self, iter_index, task_paths, cc, models):
        model_devi_jobs = self.jdata['model_devi_jobs']
        cur_job = model_devi_jobs[iter_index]
        dt = cur_job.get("dt", None)
        if dt is not None:
            model_devi_dt = dt
        else:
            model_devi_dt = self.jdata['model_devi_dt']
        nsteps = cur_job.get("nsteps", None)
        if nsteps is None:
            raise RuntimeError("nsteps is None, you should set nsteps in model_devi_jobs!")
        # Currently Gromacs engine is not supported for different temperatures!
        # If you want to change temperatures, you should change it in mdp files.
    
        sys_idx = expand_idx(cur_job['sys_idx'])
        if (len(sys_idx) != len(list(set(sys_idx)))) :
            raise RuntimeError("system index should be uniq")

        task_path = next(task_paths)
        # dlog.info(task_path)
        create_path(task_path)
        #create_path(os.path.join(task_path, 'traj'))
        #loc_conf_name = 'conf.lmp'
        gromacs_settings = self.jdata.get("gromacs_settings" , "")
        for key,file in gromacs_settings.items():
            if key != "traj_filename" and key != "mdp_filename":
                # cc is abspath; use relpath instead
                os.symlink(os.path.relpath(os.path.join(cc,file), task_path), os.path.join(task_path, file))
        
        # input.json for DP-Gromacs
        with open(os.path.join(cc, "input.json")) as f:
            input_json = json.load(f)
        input_json["graph_file"] = models[0]
        with open(os.path.join(task_path,'input.json'), 'w') as _outfile:
            json.dump(input_json, _outfile, indent = 4)

        # trj_freq
        trj_freq = cur_job.get("trj_freq", 10)
        mdp = MDP()
        mdp.read(os.path.join(cc, gromacs_settings['mdp_filename']))
        mdp['nstcomm'] = trj_freq
        mdp['nstxout'] = trj_freq
        mdp['nstlog'] = trj_freq
        mdp['nstenergy'] = trj_freq
        # dt
        mdp['dt'] = dt
        mdp.write(os.path.join(task_path, gromacs_settings['mdp_filename']))

        cwd_ = os.getcwd()
        os.chdir(task_path)
        job = {}
        
        job["model_devi_dt"] =  model_devi_dt
        job["nsteps"] = nsteps
        with open('job.json', 'w') as _outfile:
            json.dump(job, _outfile, indent = 4)
        os.chdir(cwd_)


class GromacsTrajectory(Trajectory):
    def get_model_deviations(self) -> np.ndarray:
        all_conf = np.loadtxt(os.path.join(self.directory, 'model_devi.out'))
        self.time = all_conf[:, 0]
        return all_conf[:, 4]

    def get_frame(self, idx: int) -> "GromacsFrame":
        return GromacsFrame(self, idx)


class GromacsFrame(Frame):
    def read_frame(self) -> dpdata.System:
        fidx = self.idx[0]
        time = self.trajectory.time[fidx]
        conf_name = os.path.join(self.trajectory.directory, 'traj', '%d.gromacstrj' % time)
        type_map = self.trajectory.engien.jdata['type_map']
        system = dpdata.System(conf_name, type_map=type_map, fmt='gromacs/gro')
        return system
