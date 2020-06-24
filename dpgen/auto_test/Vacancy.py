from dpgen.auto_test.Property import Property
from dpgen.auto_test.refine import make_refine
from dpgen.auto_test.reproduce import make_repro
from pymatgen.core.structure import Structure
from pymatgen.analysis.defects.generators import VacancyGenerator
import dpgen.auto_test.lib.vasp as vasp
import dpgen.auto_test.lib.lammps as lammps
from monty.serialization import loadfn,dumpfn
from dpgen import dlog
import numpy as np
import os,json


class Vacancy (Property) :
    def __init__ (self,
                  parameter) :
        self.parameter = parameter
        default_supercell = [1, 1, 1]
        self.supercell = parameter.get('supercell', default_supercell)
        self.reprod = parameter.get('reprod-opt', False)
        self.cal_type = parameter.get('cal_type', 'relaxation')
        default_cal_setting = {"relax_pos": True,
                               "relax_shape": True,
                               "relax_vol": True}
        self.cal_setting = parameter.get('cal_setting', default_cal_setting)

    def make_confs(self,
                   path_to_work,
                   path_to_equi,
                   refine=False):
        path_to_work = os.path.abspath(path_to_work)
        if os.path.exists(path_to_work):
            dlog.warning('%s already exists' % path_to_work)
        else:
            os.makedirs(path_to_work)
        path_to_equi = os.path.abspath(path_to_equi)
        task_list = []
        cwd = os.getcwd()

        print('gen vacancy with supercell ' + str(self.supercell))

        equi_contcar = os.path.join(path_to_equi, 'CONTCAR')
        if not os.path.exists(equi_contcar):
            raise RuntimeError("please do relaxation first")

        ss = Structure.from_file(equi_contcar)
        vds = VacancyGenerator(ss)
        dss = []
        for jj in vds:
            dss.append(jj.generate_defect_structure(self.supercell))

        if refine:
            task_list = make_refine(self.parameter['init_from_suffix'],
                                    self.parameter['output_suffix'],
                                    path_to_work,
                                    len(dss))
            for ii in task_list:
                os.chdir(ii)
                np.savetxt('supercell.out', self.supercell, fmt='%d')
            os.chdir(cwd)

        if self.reprod:
            self.cal_type = 'static'
            if 'vasp_lmp_path' not in self.parameter:
                raise RuntimeError("please provide the vasp_lmp_path for reproduction")
            vasp_lmp_path = os.path.abspath(self.parameter['vasp_lmp_path'])
            task_list = make_repro(vasp_lmp_path, path_to_work)
            os.chdir(cwd)

        else:
            os.chdir(path_to_work)
            if os.path.isfile('POSCAR'):
                os.remove('POSCAR')
            os.symlink(os.path.relpath(equi_contcar), 'POSCAR')
#           task_poscar = os.path.join(output, 'POSCAR')

            for ii in range(len(dss)):
                output_task = os.path.join(path_to_work, 'task.%06d' % ii)
                os.makedirs(output_task, exist_ok=True)
                os.chdir(output_task)
                for jj in ['INCAR', 'POTCAR', 'POSCAR', 'conf.lmp', 'in.lammps']:
                    if os.path.exists(jj):
                        os.remove(jj)
                task_list.append(output_task)
                dss[ii].to('POSCAR', 'POSCAR')
                np.savetxt('supercell.out', self.supercell, fmt='%d')
            os.chdir(cwd)
        return task_list

    def task_type(self):
        return self.parameter['type']

    def task_param(self):
        return self.parameter

    def _compute_lower(self,
                       output_file,
                       all_tasks,
                       all_res):
        output_file = os.path.abspath(output_file)
        res_data = {}
        ptr_data = os.path.dirname(output_file)+'\n'


        if not self.reprod:
            ptr_data += "Structure: \tVac_E(eV)  E(eV) equi_E(eV)\n"
            for ii in all_tasks:
                with open(os.path.join(ii, 'inter.json')) as fp:
                    idata = json.load(fp)
                inter_type = idata['type']
                equi_path = os.path.abspath(os.path.join(os.path.dirname(output_file), '../relaxation'))
                structure_dir = os.path.basename(ii)

                if inter_type == 'vasp':
                    equi_outcar = os.path.join(equi_path, 'OUTCAR')
                    equi_natoms, equi_epa, equi_vpa = vasp.get_nev(equi_outcar)
                elif inter_type in ['deepmd', 'meam', 'eam_fs', 'eam_alloy']:
                    equi_log = os.path.join(equi_path, 'log.lammps')
                    equi_natoms, equi_epa, equi_vpa = lammps.get_nev(equi_log)

                else:
                    raise RuntimeError('interaction type not supported')

                natoms = len(all_res[ii]['force']) / 3
                epa = all_res[ii]['energy'] / natoms
                evac = epa * natoms - equi_epa * natoms
                ptr_data += "%s: %7.3f  %7.3f %7.3f \n" % (structure_dir, evac, epa * natoms, equi_epa * natoms)
                res_data[structure_dir] = [evac, epa * natoms, equi_epa * natoms]

        else:
            if 'vasp_lmp_path' not in self.parameter:
                raise RuntimeError("please provide the vasp_lmp_path for reproduction")
            vasp_lmp_path = os.path.abspath(self.parameter['vasp_lmp_path'])
            res_data, ptr_data = reproduce.post_repro(vasp_lmp_path,all_tasks,ptr_data)

        with open(output_file, 'w') as fp:
            json.dump(res_data, fp, indent=4)

        return res_data, ptr_data
