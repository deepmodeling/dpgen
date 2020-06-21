import os
import warnings
import json
from Task import Task
import dpdata
import dpgen.auto_test.lib.lammps as lammps


class EAM_FS_LMP(Task):
    def __init__(self,
                 inter_parameter,
                 path_to_poscar):
        self.inter = inter_parameter
        self.model = os.path.abspath(inter_parameter['model'])
        self.type_map = inter_parameter['type_map']
        self.path_to_poscar = path_to_poscar

    def make_potential_files(self,
                             output_dir):
        model_file = os.path.basename(self.model)
        cwd = os.getcwd()
        os.chdir(output_dir)
        if os.path.islink(model_file):
            link_file = os.readlink(model_file)
            if not os.path.abspath(link_file) == self.model:
                os.remove(model_file)
                os.symlink(os.path.relpath(self.model), model_file)
        else:
            os.symlink(os.path.relpath(self.model), model_file)
        os.chdir(cwd)

        with open(os.path.join(output_dir, 'inter.json'), 'w') as fp:
            json.dump(self.inter, fp, indent=4)

    def make_input_file(self,
                        output_dir,
                        task_type,
                        task_param):
        lammps.cvt_lammps_conf(os.path.join(output_dir, 'POSCAR'), 'conf.lmp')
        with open(os.path.join(output_dir, 'task.json'), 'w') as fp:
            json.dump(task_param, fp, indent=4)
        # lines in lammps.in related to model
        #         line_model = "pair_style      eam/fs \n"
        #         line_model += "pair_coeff      * * %s " % (os.path.basename(self.model))
        #         for ii in self.type_map:
        #             line_model += ii + ' '
        #         line_model += '\n'

        etol = 1e-12
        ftol = 1e-6
        maxiter = 5000
        maxeval = 500000
        change_box = True
        B0 = 70
        bp = 0
        scale2equi = 1
        ntypes = len(self.type_map)
        reprod_opt = False
        static = False

        if 'etol' in task_param:
            etol = task_param['etol']
        if 'ftol' in task_param:
            ftol = task_param['ftol']
        if 'maxiter' in task_param:
            maxiter = task_param['maxiter']
        if 'maxeval' in task_param:
            maxeval = task_param['maxeval']
        if 'change_box' in task_param:
            change_box = task_param['change_box']
        if 'scale2equi' in task_param:
            scale2equi = task_param['scale2equi']
        if 'reprod_opt' in task_param:
            reprod_opt = task_param['reprod_opt']
        if 'static-opt' in task_param:
            static = task_param['static-opt']

        model_name = os.path.basename(self.model)

        model_param = {'model_name': model_name,
                       'param_type': self.type_map}

        fc = ''
        if task_type == 'relaxation' \
                or (task_type == 'eos' and not change_box) \
                or (task_type == 'surface' and not static):
            fc = lammps.make_lammps_equi('conf.lmp', ntypes, lammps.inter_eam_fs, model_param,
                                         etol, ftol, maxiter, maxeval, change_box)

        if task_type == 'static' \
                or (task_type == 'surface' and static):
            fc = lammps.make_lammps_eval('conf.lmp', ntypes, lammps.inter_eam_fs, model_param)

        if task_type == 'elastic':
            fc = lammps.make_lammps_elastic('conf.lmp', ntypes, lammps.inter_eam_fs, model_param,
                                            etol, ftol, maxiter, maxeval)

        if task_type == 'vacancy' \
                or (task_type == 'eos' and change_box) \
                or (task_type == 'interstitial'):
            fc = lammps.make_lammps_press_relax('conf.lmp', ntypes, scale2equi, lammps.inter_eam_fs,
                                                model_param, B0, bp, etol, ftol, maxiter, maxeval)

        if reprod_opt:
            fc = lammps.make_lammps_eval('conf.lmp', ntypes, lammps.inter_eam_fs, model_param)

        with open(os.path.join(output_dir, 'in.lammps'), 'w') as fp:
            fp.write(fc)

    def compute(self,
                output_dir):
        log_lammps = os.path.join(output_dir, 'log.lammps')
        if not os.path.isfile(log_lammps):
            warnings.warn("cannot find log.lammps in " + output_dir + " skip")
            return None
        else:
            with open(log_lammps, 'r') as fp:
                if 'Total wall time:' not in fp.read():
                    warnings.warn("lammps not finished " + log_lammps + " skip")
                    return None
                else:
                    fp.seek(0)
                    lines = fp.read().split('\n')
                    for ii in lines:
                        if ("Total number of atoms" in ii) and (not 'print' in ii):
                            natoms = float(ii.split('=')[1].split()[0])
                        if ("Final energy per atoms" in ii) and (not 'print' in ii):
                            epa = float(ii.split('=')[1].split()[0])

                    dump = os.path.join(output_dir, 'dump.relax')
                    contcar = os.path.join(output_dir, 'CONTCAR')
                    d_dump = dpdata.System(dump, fmt='lammps/dump')
                    d_dump.to('vasp/poscar', contcar, frame_idx=-1)
                    force = d_dump['forces']

                    result_dict = {"energy": natoms * epa, "force": list(force[-1].reshape(natoms * 3))}
                    return result_dict

    def forward_files(self):
        return ['conf.lmp', 'in.lammps', os.path.basename(self.model)]

    def forward_common_files(self):
        return ['in.lammps', os.path.basename(self.model)]

    def backward_files(self):
        return ['log.lammps', 'eam_fs.out', 'dump.relax']
