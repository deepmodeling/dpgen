import os
import warnings
import dpgen.auto_test.lib.lammps as lammps
from dpgen import dlog
from monty.serialization import loadfn, dumpfn
from dpgen.auto_test.Task import Task
from dpgen.auto_test.lib.lammps import inter_deepmd, inter_meam, inter_eam_fs, inter_eam_alloy

supported_inter = ["deepmd", 'meam', 'eam_fs', 'eam_alloy']


class Lammps(Task):
    def __init__(self,
                 inter_parameter,
                 path_to_poscar):
        self.inter = inter_parameter
        self.inter_type = inter_parameter['type']
        self.type_map = inter_parameter['type_map']
        self.in_lammps = inter_parameter.get('in_lammps', 'auto')
        if self.inter_type == 'meam':
            self.model = list(map(os.path.abspath, inter_parameter['model']))
        else:
            self.model = os.path.abspath(inter_parameter['model'])
        self.path_to_poscar = path_to_poscar
        assert self.inter_type in supported_inter
        self.set_inter_type_func()

    def set_inter_type_func(self):

        if self.inter_type == "deepmd":
            self.inter_func = inter_deepmd

        elif self.inter_type == 'meam':
            self.inter_func = inter_meam

        elif self.inter_type == 'eam_fs':
            self.inter_func = inter_eam_fs

        else:
            self.inter_func = inter_eam_alloy

    def set_model_param(self):

        if self.inter_type == "deepmd":
            model_name = os.path.basename(self.model)
            deepmd_version = self.inter.get("deepmd_version", "1.2.0")
            self.model_param = {'model_name': [model_name],
                                'param_type': self.type_map,
                                'deepmd_version': deepmd_version}
        elif self.inter_type == 'meam':
            model_name = list(map(os.path.basename, self.model))
            self.model_param = {'model_name': [model_name],
                                'param_type': self.type_map}
        else:
            model_name = os.path.basename(self.model)
            self.model_param = {'model_name': [model_name],
                                'param_type': self.type_map}

    def make_potential_files(self,
                             output_dir):
        cwd = os.getcwd()
        if self.inter_type == 'meam':
            model_lib = os.path.basename(self.model[0])
            model_file = os.path.basename(self.model[1])
            os.chdir(os.path.join(output_dir, '../'))
            if os.path.islink(model_lib):
                link_lib = os.readlink(model_lib)
                if not os.path.abspath(link_lib) == self.model[0]:
                    os.remove(model_lib)
                    os.symlink(os.path.relpath(self.model[0]), model_lib)
            else:
                os.symlink(os.path.relpath(self.model[0]), model_lib)

            if os.path.islink(model_file):
                link_file = os.readlink(model_file)
                if not os.path.abspath(link_file) == self.model[1]:
                    os.remove(model_file)
                    os.symlink(os.path.relpath(self.model[1]), model_file)
            else:
                os.symlink(os.path.relpath(self.model[1]), model_file)
            os.chdir(output_dir)
            if not os.path.islink(model_lib):
                os.symlink(os.path.join('..', model_lib), model_lib)
            elif not os.path.join('..', model_lib) == os.readlink(model_lib):
                os.remove(model_lib)
                os.symlink(os.path.join('..', model_lib), model_lib)

            if not os.path.islink(model_file):
                os.symlink(os.path.join('..', model_file), model_file)
            elif not os.path.join('..', model_file) == os.readlink(model_file):
                os.remove(model_file)
                os.symlink(os.path.join('..', model_file), model_file)
            os.chdir(cwd)

        else:
            model_file = os.path.basename(self.model)
            os.chdir(os.path.join(output_dir, '../'))
            if os.path.islink(model_file):
                link_file = os.readlink(model_file)
                if not os.path.abspath(link_file) == self.model:
                    os.remove(model_file)
                    os.symlink(os.path.relpath(self.model), model_file)
            else:
                os.symlink(os.path.relpath(self.model), model_file)
            os.chdir(output_dir)
            if not os.path.islink(model_file):
                os.symlink(os.path.join('..', model_file), model_file)
            elif not os.path.join('..', model_file) == os.readlink(model_file):
                os.remove(model_file)
                os.symlink(os.path.join('..', model_file), model_file)
            os.chdir(cwd)

        dumpfn(self.inter, os.path.join(output_dir, 'inter.json'), indent=4)

    def make_input_file(self,
                        output_dir,
                        task_type,
                        task_param):
        lammps.cvt_lammps_conf(os.path.join(output_dir, 'POSCAR'), os.path.join(output_dir, 'conf.lmp'), lammps.element_list(self.type_map))

        # dumpfn(task_param, os.path.join(output_dir, 'task.json'), indent=4)

        etol = 1e-12
        ftol = 1e-6
        maxiter = 5000
        maxeval = 500000
        B0 = 70
        bp = 0
        ntypes = len(self.type_map)

        cal_type = task_param['cal_type']
        cal_setting = task_param['cal_setting']

        self.set_model_param()

        # deal with user input in.lammps for relaxation
        if os.path.isfile(self.in_lammps) and task_type == 'relaxation':
            with open(self.in_lammps, 'r') as fin:
                fc = fin.read()
        # user input in.lammps for property calculation
        if 'input_prop' in cal_setting and os.path.isfile(cal_setting['input_prop']):
            with open(os.path.abspath(cal_setting['input_prop']), 'r') as fin:
                fc = fin.read()

        else:
            if 'etol' in cal_setting:
                dlog.info("%s setting etol to %s" % (self.make_input_file.__name__, cal_setting['etol']))
                etol = cal_setting['etol']
            if 'ftol' in cal_setting:
                dlog.info("%s setting ftol to %s" % (self.make_input_file.__name__, cal_setting['ftol']))
                ftol = cal_setting['ftol']
            if 'maxiter' in cal_setting:
                dlog.info("%s setting maxiter to %s" % (self.make_input_file.__name__, cal_setting['maxiter']))
                maxiter = cal_setting['maxiter']
            if 'maxeval' in cal_setting:
                dlog.info("%s setting maxeval to %s" % (self.make_input_file.__name__, cal_setting['maxeval']))
                maxeval = cal_setting['maxeval']

            if cal_type == 'relaxation':
                relax_pos = cal_setting['relax_pos']
                relax_shape = cal_setting['relax_shape']
                relax_vol = cal_setting['relax_vol']

                if [relax_pos, relax_shape, relax_vol] == [True, False, False]:
                    fc = lammps.make_lammps_equi('conf.lmp', self.type_map, self.inter_func, self.model_param,
                                                 etol, ftol, maxiter, maxeval, False)
                elif [relax_pos, relax_shape, relax_vol] == [True, True, True]:
                    fc = lammps.make_lammps_equi('conf.lmp', self.type_map, self.inter_func, self.model_param,
                                                 etol, ftol, maxiter, maxeval, True)
                elif [relax_pos, relax_shape, relax_vol] == [True, True, False] and not task_type == 'eos':
                    if 'scale2equi' in task_param:
                        scale2equi = task_param['scale2equi']
                        fc = lammps.make_lammps_press_relax('conf.lmp', self.type_map, scale2equi[int(output_dir[-6:])],
                                                            self.inter_func,
                                                            self.model_param, B0, bp, etol, ftol, maxiter, maxeval)
                    else:
                        fc = lammps.make_lammps_equi('conf.lmp', self.type_map, self.inter_func, self.model_param,
                                                     etol, ftol, maxiter, maxeval, True)
                elif [relax_pos, relax_shape, relax_vol] == [True, True, False] and task_type == 'eos':
                    task_param['cal_setting']['relax_shape'] = False
                    fc = lammps.make_lammps_equi('conf.lmp', self.type_map, self.inter_func, self.model_param,
                                                 etol, ftol, maxiter, maxeval, False)
                elif [relax_pos, relax_shape, relax_vol] == [False, False, False]:
                    fc = lammps.make_lammps_eval('conf.lmp', self.type_map, self.inter_func, self.model_param)

                else:
                    raise RuntimeError("not supported calculation setting for LAMMPS")

            elif cal_type == 'static':
                fc = lammps.make_lammps_eval('conf.lmp', self.type_map, self.inter_func, self.model_param)

            else:
                raise RuntimeError("not supported calculation type for LAMMPS")

        dumpfn(task_param, os.path.join(output_dir, 'task.json'), indent=4)

        in_lammps_not_link_list = ['eos']
        if task_type not in in_lammps_not_link_list:
            with open(os.path.join(output_dir, '../in.lammps'), 'w') as fp:
                fp.write(fc)
            cwd = os.getcwd()
            os.chdir(output_dir)
            if not (os.path.islink('in.lammps') or os.path.isfile('in.lammps')):
                os.symlink('../in.lammps', 'in.lammps')
            else:
                os.remove('in.lammps')
                os.symlink('../in.lammps', 'in.lammps')
            os.chdir(cwd)
        else:
            with open(os.path.join(output_dir, 'in.lammps'), 'w') as fp:
                fp.write(fc)

    def compute(self,
                output_dir):
        log_lammps = os.path.join(output_dir, 'log.lammps')
        dump_lammps = os.path.join(output_dir, 'dump.relax')
        if not os.path.isfile(log_lammps):
            warnings.warn("cannot find log.lammps in " + output_dir + " skip")
            return None
        if not os.path.isfile(dump_lammps):
            warnings.warn("cannot find dump.relax in " + output_dir + " skip")
            return None
        else:
            box = []
            coord = []
            vol = []
            energy = []
            force = []
            virial = []
            stress = []
            with open(dump_lammps, 'r') as fin:
                dump = fin.read().split('\n')
            dumptime = []
            for idx, ii in enumerate(dump):
                if ii == 'ITEM: TIMESTEP':
                    box.append([])
                    coord.append([])
                    force.append([])
                    dumptime.append(int(dump[idx + 1]))
                    natom = int(dump[idx + 3])
                    xlo_bound = float(dump[idx + 5].split()[0])
                    xhi_bound = float(dump[idx + 5].split()[1])
                    xy = float(dump[idx + 5].split()[2])
                    ylo_bound = float(dump[idx + 6].split()[0])
                    yhi_bound = float(dump[idx + 6].split()[1])
                    xz = float(dump[idx + 6].split()[2])
                    zlo = float(dump[idx + 7].split()[0])
                    zhi = float(dump[idx + 7].split()[1])
                    yz = float(dump[idx + 7].split()[2])
                    xx = xhi_bound - max([0, xy, xz, xy + xz]) - (xlo_bound - min([0, xy, xz, xy + xz]))
                    yy = yhi_bound - max([0, yz]) - (ylo_bound - min([0, yz]))
                    zz = zhi - zlo
                    box[-1].append([xx, 0.0, 0.0])
                    box[-1].append([xy, yy, 0.0])
                    box[-1].append([xz, yz, zz])
                    vol.append(xx * yy * zz)
                    type_list = []
                    for jj in range(natom):
                        type_list.append(int(dump[idx + 9 + jj].split()[1]) - 1)
                        if 'xs ys zs' in dump[idx + 8]:
                            a_x = float(dump[idx + 9 + jj].split()[2]) * xx + float(
                                dump[idx + 9 + jj].split()[3]) * xy \
                                  + float(dump[idx + 9 + jj].split()[4]) * xz
                            a_y = float(dump[idx + 9 + jj].split()[3]) * yy + float(
                                dump[idx + 9 + jj].split()[4]) * yz
                            a_z = float(dump[idx + 9 + jj].split()[4]) * zz
                        else:
                            a_x = float(dump[idx + 9 + jj].split()[2])
                            a_y = float(dump[idx + 9 + jj].split()[3])
                            a_z = float(dump[idx + 9 + jj].split()[4])
                        coord[-1].append([a_x, a_y, a_z])
                        fx = float(dump[idx + 9 + jj].split()[5])
                        fy = float(dump[idx + 9 + jj].split()[6])
                        fz = float(dump[idx + 9 + jj].split()[7])
                        force[-1].append([fx, fy, fz])

            with open(log_lammps, 'r') as fp:
                if 'Total wall time:' not in fp.read():
                    warnings.warn("lammps not finished " + log_lammps + " skip")
                    return None
                else:
                    fp.seek(0)
                    lines = fp.read().split('\n')
                    idid = -1
                    for ii in dumptime:
                        idid += 1
                        for jj in lines:
                            line = jj.split()
                            if len(line) and str(ii) == line[0]:
                                try:
                                    [float(kk) for kk in line]
                                except:
                                    continue
                                stress.append([])
                                virial.append([])
                                energy.append(float(line[1]))
                                # virials = stress * vol * 1e5 *1e-30 * 1e19/1.6021766208
                                stress[-1].append([float(line[2]) / 1000.0, float(line[5]) / 1000.0, float(line[6]) / 1000.0])
                                stress[-1].append([float(line[5]) / 1000.0, float(line[3]) / 1000.0, float(line[7]) / 1000.0])
                                stress[-1].append([float(line[6]) / 1000.0, float(line[7]) / 1000.0, float(line[4]) / 1000.0])
                                stress_to_virial = vol[idid] * 1e5 * 1e-30 * 1e19 / 1.6021766208
                                virial[-1].append([float(line[2]) * stress_to_virial, float(line[5]) * stress_to_virial,
                                                   float(line[6]) * stress_to_virial])
                                virial[-1].append([float(line[5]) * stress_to_virial, float(line[3]) * stress_to_virial,
                                                   float(line[7]) * stress_to_virial])
                                virial[-1].append([float(line[6]) * stress_to_virial, float(line[7]) * stress_to_virial,
                                                   float(line[4]) * stress_to_virial])
                                break

            _tmp = self.type_map
            dlog.debug(_tmp)
            type_map_list = lammps.element_list(self.type_map)

            type_map_idx = list(range(len(type_map_list)))
            atom_numbs = []
            for ii in type_map_idx:
                count = 0
                for jj in type_list:
                    if jj == ii:
                        count += 1
                atom_numbs.append(count)

            # d_dump = dpdata.System(dump_lammps, fmt='lammps/dump', type_map=type_map_list)
            # d_dump.to('vasp/poscar', contcar, frame_idx=-1)

            result_dict = {"@module": "dpdata.system", "@class": "LabeledSystem", "data": {"atom_numbs": atom_numbs,
                                                                                           "atom_names": type_map_list,
                                                                                           "atom_types": {
                                                                                               "@module": "numpy",
                                                                                               "@class": "array",
                                                                                               "dtype": "int64",
                                                                                               "data": type_list},
                                                                                           "orig": {"@module": "numpy",
                                                                                                    "@class": "array",
                                                                                                    "dtype": "int64",
                                                                                                    "data": [0, 0, 0]},
                                                                                           "cells": {"@module": "numpy",
                                                                                                     "@class": "array",
                                                                                                     "dtype": "float64",
                                                                                                     "data": box},
                                                                                           "coords": {
                                                                                               "@module": "numpy",
                                                                                               "@class": "array",
                                                                                               "dtype": "float64",
                                                                                               "data": coord},
                                                                                           "energies": {
                                                                                               "@module": "numpy",
                                                                                               "@class": "array",
                                                                                               "dtype": "float64",
                                                                                               "data": energy},
                                                                                           "forces": {
                                                                                               "@module": "numpy",
                                                                                               "@class": "array",
                                                                                               "dtype": "float64",
                                                                                               "data": force},
                                                                                           "virials": {
                                                                                               "@module": "numpy",
                                                                                               "@class": "array",
                                                                                               "dtype": "float64",
                                                                                               "data": virial},
                                                                                           "stress": {
                                                                                               "@module": "numpy",
                                                                                               "@class": "array",
                                                                                               "dtype": "float64",
                                                                                               "data": stress}}}

            contcar = os.path.join(output_dir, 'CONTCAR')
            dumpfn(result_dict, contcar, indent=4)
            d_dump = loadfn(contcar)
            d_dump.to('vasp/poscar', contcar, frame_idx=-1)

            return result_dict

    def forward_files(self, property_type='relaxation'):
        if self.inter_type == 'meam':
            return ['conf.lmp', 'in.lammps'] + list(map(os.path.basename, self.model))
        else:
            return ['conf.lmp', 'in.lammps', os.path.basename(self.model)]

    def forward_common_files(self, property_type='relaxation'):
        if property_type not in ['eos']:
            if self.inter_type == 'meam':
                return ['in.lammps'] + list(map(os.path.basename, self.model))
            else:
                return ['in.lammps', os.path.basename(self.model)]
        else:
            if self.inter_type == 'meam':
                return list(map(os.path.basename, self.model))
            else:
                return [os.path.basename(self.model)]

    def backward_files(self, property_type='relaxation'):
        return ['log.lammps', 'outlog', 'dump.relax']
