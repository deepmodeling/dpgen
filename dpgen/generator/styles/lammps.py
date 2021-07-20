"""This file is aimmed to implement LAMMPS engien."""
from __future__ import annotations
import os
import json
import shutil
import copy
from typing import List, TYPE_CHECKING, Tuple, Iterator, Dict
from distutils.version import LooseVersion

import dpdata
import numpy as np
from dpgen.generator.lib.utils import create_path, _get_param_alias, expand_idx
from dpgen.generator.model_devi import ModelDeviEngien, Trajectory, Frame
from dpgen.generator.lib.lammps import make_lammps_input
from dpgen.generator.lib.gaussian import take_cluster


@ModelDeviEngien.register("lammps")
class LAMMPSEngien(ModelDeviEngien):
    def make_input(self, iter_index: int, sys_index: int, directories: Iterator[str], conf_name: str, models: List[str]):
        # the default system format is vasp/poscar
        fmt = self.jdata.get('sys_format', 'vasp/poscar')
        system = dpdata.System(conf_name, fmt = fmt, type_map = self.jdata['type_map'])
        if self.jdata.get('model_devi_nopbc', False):
            system.remove_pbc()
        new_conf_name = os.path.splitext(conf_name)[0] + ".lmp"
        system.to_lammps_lmp(new_conf_name)
        model_devi_jobs = self.jdata['model_devi_jobs']
        cur_job = model_devi_jobs[iter_index]
        input_mode = "native"
        if "template" in cur_job:
            input_mode = "revise_template"
        if input_mode == "native":
            _make_model_devi_native(iter_index, directories, self.jdata, self.mdata, new_conf_name, models)
        elif input_mode == "revise_template":
            _make_model_devi_revmat(iter_index, directories, self.jdata, self.mdata, new_conf_name, models, sys_index)
        else:
            raise RuntimeError('unknown model_devi input mode', input_mode)

    def get_running_parameters(self, work_path: str) -> Dict[str]:
        use_plm = self.jdata.get('model_devi_plumed', False)
        use_plm_path = self.jdata.get('model_devi_plumed_path', False)
        lmp_exec = self.mdata['lmp_command']
        command = "{ if [ ! -f dpgen.restart.10000 ]; then %s -i input.lammps -v restart 0; else %s -i input.lammps -v restart 1; fi }" % (lmp_exec, lmp_exec)
        command = "/bin/sh -c '%s'" % command
        forward_files = ['conf.lmp', 'input.lammps', 'traj']
        backward_files = ['model_devi.out', 'model_devi.log', 'traj']
        if use_plm:
            forward_files += ['input.plumed']
           # backward_files += ['output.plumed']
            backward_files += ['output.plumed','COLVAR','dump.0.xyz']
            if use_plm_path:
                forward_files += ['plmpath.pdb']
        return {
            "command": command,
            "forward_files": forward_files,
            "backward_files": backward_files,
            "common_files": [],
        }

    def extract_trajectory(self, directory) -> 'Trajectory':
        return LAMMPSTrajectory(self, directory)

class LAMMPSTrajectory(Trajectory):
    def get_model_deviations(self) -> np.ndarray:
        all_conf = np.loadtxt(os.path.join(self.directory, 'model_devi.out'))
        # in current LAMMPS DeePMD plugin, it should has 7 colums for frame info
        # the fifth column (index 4) is max f
        # the atomic deviation starts from index 7
        self.time = all_conf[:, 0]
        if all_conf.shape[1] == 7:
            # 1D
            return all_conf[:, 4]
        elif all_conf.shape[1] > 7:
            # 2D
            return all_conf[:, 7:]
        else:
            raise RuntimeError("The format of model_devi.out is not supported")
    
    def get_frame(self, idx: int) -> "LAMMPSFrame":
        return LAMMPSFrame(self, idx)


class LAMMPSFrame(Frame):
    def read_frame(self) -> dpdata.System:
        fidx = self.idx[0]
        time = self.trajectory.time[fidx]
        conf_name = os.path.join(self.trajectory.directory, 'traj', '%d.lammpstrj' % time)
        type_map = self.trajectory.engien.jdata['type_map']
        system = dpdata.System(conf_name, type_map=type_map, fmt='lammps/dump')
        if len(self.idx) > 1:
            # cluster
            cidx = self.idx[1]
            return take_cluster(conf_name, cidx, self.trajectory.engien.jdata)
        # else: frame
        return system


def _make_model_devi_native(iter_index, task_paths, jdata, mdata, conf_name, models):
    model_devi_jobs = jdata['model_devi_jobs']
    cur_job = model_devi_jobs[iter_index]
    ensemble, nsteps, trj_freq, temps, press, pka_e, dt = parse_cur_job(cur_job)
    if dt is not None :
        model_devi_dt = dt
    sys_idx = expand_idx(cur_job['sys_idx'])
    if (len(sys_idx) != len(list(set(sys_idx)))) :
        raise RuntimeError("system index should be uniq")

    use_ele_temp = jdata.get('use_ele_temp', 0)
    model_devi_dt = jdata['model_devi_dt']
    model_devi_neidelay = None
    if 'model_devi_neidelay' in jdata :
        model_devi_neidelay = jdata['model_devi_neidelay']
    model_devi_taut = 0.1
    if 'model_devi_taut' in jdata :
        model_devi_taut = jdata['model_devi_taut']
    model_devi_taup = 0.5
    if 'model_devi_taup' in jdata :
        model_devi_taup = jdata['model_devi_taup']
    mass_map = jdata['mass_map']
    nopbc = jdata.get('model_devi_nopbc', False)

    for tt_ in temps:
        if use_ele_temp:
            if type(tt_) == list:
                tt = tt_[0]
                if use_ele_temp == 1:
                    te_f = tt_[1]
                    te_a = None
                else:
                    te_f = None
                    te_a = tt_[1]
            else:
                assert(type(tt_) == float or type(tt_) == int)
                tt = float(tt_)
                if use_ele_temp == 1:
                    te_f = tt
                    te_a = None
                else:
                    te_f = None
                    te_a = tt
        else :
            tt = tt_
            te_f = None
            te_a = None
        for pp in press:
            task_path = next(task_paths)
            # dlog.info(task_path)
            create_path(task_path)
            create_path(os.path.join(task_path, 'traj'))
            loc_conf_name = 'conf.lmp'
            # conf_name is absolute path; generate relative path to conf_name
            os.symlink(os.path.relpath(conf_name, task_path), os.path.join(task_path, loc_conf_name) )
            cwd_ = os.getcwd()
            os.chdir(task_path)
            deepmd_version = mdata.get('deepmd_version', '1')
            file_c = make_lammps_input(ensemble,
                                        loc_conf_name,
                                        models,
                                        nsteps,
                                        model_devi_dt,
                                        model_devi_neidelay,
                                        trj_freq,
                                        mass_map,
                                        tt,
                                        jdata = jdata,
                                        tau_t = model_devi_taut,
                                        pres = pp,
                                        tau_p = model_devi_taup,
                                        pka_e = pka_e,
                                        ele_temp_f = te_f,
                                        ele_temp_a = te_a,
                                        nopbc = nopbc,
                                        deepmd_version = deepmd_version)
            job = {}
            job["ensemble"] = ensemble
            job["press"] = pp
            job["temps"] = tt
            if te_f is not None:
                job["ele_temp"] = te_f
            if te_a is not None:
                job["ele_temp"] = te_a
            job["model_devi_dt"] =  model_devi_dt
            with open('job.json', 'w') as _outfile:
                json.dump(job, _outfile, indent = 4)
            os.chdir(cwd_)
            with open(os.path.join(task_path, 'input.lammps'), 'w') as fp :
                fp.write(file_c)


def _make_model_devi_revmat(iter_index, task_paths, jdata, mdata, conf_name, models, sys_index):
    model_devi_jobs = jdata['model_devi_jobs']
    cur_job = model_devi_jobs[iter_index]
    sys_idx = expand_idx(cur_job['sys_idx'])
    if (len(sys_idx) != len(list(set(sys_idx)))) :
        raise RuntimeError("system index should be uniq")
    use_plm = jdata.get('model_devi_plumed', False)
    use_plm_path = jdata.get('model_devi_plumed_path', False)
    trj_freq = _get_param_alias(cur_job, ['t_freq', 'trj_freq','traj_freq'])

    rev_keys, rev_mat, num_lmp = parse_cur_job_revmat(cur_job, use_plm = use_plm)
    lmp_templ = cur_job['template']['lmp']
    lmp_templ = os.path.abspath(lmp_templ)
    if use_plm:
        plm_templ = cur_job['template']['plm']
        plm_templ = os.path.abspath(plm_templ)
        if use_plm_path:
            plm_path_templ = cur_job['template']['plm_path']
            plm_path_templ = os.path.abspath(plm_path_templ)

    deepmd_version = mdata.get('deepmd_version', '1')

    sys_rev = cur_job.get('sys_rev_mat', None)
    total_rev_keys = rev_keys
    total_rev_mat = rev_mat
    total_num_lmp = num_lmp
    if sys_rev is not None:
        total_rev_mat = []
        sys_rev_keys, sys_rev_mat, sys_num_lmp = parse_cur_job_sys_revmat(cur_job,
                                                                            sys_idx=sys_idx[sys_index],
                                                                            use_plm=use_plm)
        _lmp_keys = rev_keys[:num_lmp] + sys_rev_keys[:sys_num_lmp]
        if use_plm:
            _plm_keys = rev_keys[num_lmp:] + sys_rev_keys[sys_num_lmp:]
            _lmp_keys += _plm_keys
        total_rev_keys = _lmp_keys
        total_num_lmp = num_lmp + sys_num_lmp
        for pub in rev_mat:
            for pri in sys_rev_mat:
                _lmp_mat = pub[:num_lmp] + pri[:sys_num_lmp]
                if use_plm:
                    _plm_mat = pub[num_lmp:] + pri[sys_num_lmp:]
                    _lmp_mat += _plm_mat
                total_rev_mat.append(_lmp_mat)
    for ii in range(len(total_rev_mat)):
        total_rev_item = total_rev_mat[ii]
        task_path = next(task_paths)
        # create task path
        create_path(task_path)
        create_path(os.path.join(task_path, 'traj'))
        # link conf
        loc_conf_name = 'conf.lmp'
        os.symlink(os.path.relpath(conf_name, task_path), os.path.join(task_path, loc_conf_name) )
        cwd_ = os.getcwd()
        # chdir to task path
        os.chdir(task_path)
        shutil.copyfile(lmp_templ, 'input.lammps')
        # revise input of lammps
        with open('input.lammps') as fp:
            lmp_lines = fp.readlines()
        lmp_lines = revise_lmp_input_model(lmp_lines, models, trj_freq, deepmd_version = deepmd_version)
        lmp_lines = revise_lmp_input_dump(lmp_lines, trj_freq)
        lmp_lines = revise_by_keys(
            lmp_lines, total_rev_keys[:total_num_lmp], total_rev_item[:total_num_lmp]
        )
        # revise input of plumed
        if use_plm:
            lmp_lines = revise_lmp_input_plm(lmp_lines, 'input.plumed')
            shutil.copyfile(plm_templ, 'input.plumed')
            with open('input.plumed') as fp:
                plm_lines = fp.readlines()
            # allow using the same list as lmp
            # user should not use the same key name for plm
            plm_lines = revise_by_keys(
                plm_lines, total_rev_keys, total_rev_item
            )
            with open('input.plumed', 'w') as fp:
                fp.write(''.join(plm_lines))
            if use_plm_path:
                shutil.copyfile(plm_path_templ, 'plmpath.pdb')
        # dump input of lammps
        with open('input.lammps', 'w') as fp:
            fp.write(''.join(lmp_lines))
        with open('job.json', 'w') as fp:
            job = {}
            for ii,jj in zip(total_rev_keys, total_rev_item) : job[ii] = jj
            json.dump(job, fp, indent = 4)
        os.chdir(cwd_)


def parse_cur_job(cur_job) :
    ensemble = _get_param_alias(cur_job, ['ens', 'ensemble'])
    temps = [-1]
    press = [-1]
    if 'npt' in ensemble :
        temps = _get_param_alias(cur_job, ['Ts','temps'])
        press = _get_param_alias(cur_job, ['Ps','press'])
    elif 'nvt' == ensemble or 'nve' == ensemble:
        temps = _get_param_alias(cur_job, ['Ts','temps'])
    nsteps = _get_param_alias(cur_job, ['nsteps'])
    trj_freq = _get_param_alias(cur_job, ['t_freq', 'trj_freq','traj_freq'])
    if 'pka_e' in cur_job :
        pka_e = _get_param_alias(cur_job, ['pka_e'])
    else :
        pka_e = None
    if 'dt' in cur_job :
        dt = _get_param_alias(cur_job, ['dt'])
    else :
        dt = None
    return ensemble, nsteps, trj_freq, temps, press, pka_e, dt


def parse_cur_job_revmat(cur_job, use_plm = False):
    templates = [cur_job['template']['lmp']]
    if use_plm :
        templates.append(cur_job['template']['plm'])
    revise_keys = []
    revise_values = []
    if 'rev_mat' not in cur_job.keys():
        cur_job['rev_mat'] = {}
    if 'lmp' not in cur_job['rev_mat'].keys():
        cur_job['rev_mat']['lmp'] = {}
    for ii in cur_job['rev_mat']['lmp'].keys():
        revise_keys.append(ii)
        revise_values.append(cur_job['rev_mat']['lmp'][ii])
    n_lmp_keys = len(revise_keys)
    if use_plm:
        if 'plm' not in cur_job['rev_mat'].keys():
            cur_job['rev_mat']['plm'] = {}
        for ii in cur_job['rev_mat']['plm'].keys():
            revise_keys.append(ii)
            revise_values.append(cur_job['rev_mat']['plm'][ii])
    revise_matrix = expand_matrix_values(revise_values)
    return revise_keys, revise_matrix, n_lmp_keys


def parse_cur_job_sys_revmat(cur_job, sys_idx, use_plm=False):
    templates = [cur_job['template']['lmp']]
    if use_plm:
        templates.append(cur_job['template']['plm'])
    sys_revise_keys = []
    sys_revise_values = []
    if 'sys_rev_mat' not in cur_job.keys():
        cur_job['sys_rev_mat'] = {}
    local_rev = cur_job['sys_rev_mat'].get(str(sys_idx), {})
    if 'lmp' not in local_rev.keys():
        local_rev['lmp'] = {}
    for ii in local_rev['lmp'].keys():
        sys_revise_keys.append(ii)
        sys_revise_values.append(local_rev['lmp'][ii])
    n_sys_lmp_keys = len(sys_revise_keys)
    if use_plm:
        if 'plm' not in local_rev.keys():
            local_rev['plm'] = {}
        for ii in local_rev['plm'].keys():
            sys_revise_keys.append(ii)
            sys_revise_values.append(local_rev['plm'][ii])
    sys_revise_matrix = expand_matrix_values(sys_revise_values)
    return sys_revise_keys, sys_revise_matrix, n_sys_lmp_keys

def revise_lmp_input_model(lmp_lines, task_model_list, trj_freq, deepmd_version = '1'):
    idx = find_only_one_key(lmp_lines, ['pair_style', 'deepmd'])
    graph_list = ' '.join(task_model_list)
    if LooseVersion(deepmd_version) < LooseVersion('1'):
        lmp_lines[idx] = "pair_style      deepmd %s %d model_devi.out\n" % (graph_list, trj_freq)
    else:
        lmp_lines[idx] = "pair_style      deepmd %s out_freq %d out_file model_devi.out\n" % (graph_list, trj_freq)
    return lmp_lines


def revise_lmp_input_dump(lmp_lines, trj_freq):
    idx = find_only_one_key(lmp_lines, ['dump', 'dpgen_dump'])
    lmp_lines[idx] = "dump            dpgen_dump all custom %d traj/*.lammpstrj id type x y z\n" % trj_freq
    return lmp_lines


def revise_lmp_input_plm(lmp_lines, in_plm, out_plm = 'output.plumed'):
    idx = find_only_one_key(lmp_lines, ['fix', 'dpgen_plm'])
    lmp_lines[idx] = "fix            dpgen_plm all plumed plumedfile %s outfile %s\n" % (in_plm, out_plm)
    return lmp_lines


def revise_by_keys(lmp_lines, keys, values):
    for kk,vv in zip(keys, values):
        for ii in range(len(lmp_lines)):
            lmp_lines[ii] = lmp_lines[ii].replace(kk, str(vv))
    return lmp_lines

def expand_matrix_values(target_list, cur_idx = 0):
    nvar = len(target_list)
    if cur_idx == nvar :
        return [[]]
    else :
        res = []
        prev = expand_matrix_values(target_list, cur_idx+1)
        for ii in target_list[cur_idx]:
            tmp = copy.deepcopy(prev)
            for jj in tmp:
               jj.insert(0, ii)
               res.append(jj)
        return res


def find_only_one_key(lmp_lines, key):
    found = []
    for idx in range(len(lmp_lines)):
        words = lmp_lines[idx].split()
        nkey = len(key)
        if len(words) >= nkey and words[:nkey] == key :
            found.append(idx)
    if len(found) > 1:
        raise RuntimeError('found %d keywords %s' % (len(found), key))
    if len(found) == 0:
        raise RuntimeError('failed to find keyword %s' % (key))
    return found[0]