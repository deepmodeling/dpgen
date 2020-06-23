from dpgen.auto_test.Property import Property
from dpgen.auto_test.refine import make_refine
from dpgen.auto_test  import reproduce
import dpgen.auto_test.lib.vasp as vasp
import dpgen.auto_test.lib.lammps as lammps
from pymatgen.core.structure import Structure
from pymatgen.core.surface import generate_all_slabs
import numpy as np
import os,json


class Surface(Property):
    def __init__(self,
                 parameter):
        self.parameter = parameter
        self.min_slab_size = parameter['min_slab_size']
        self.min_vacuum_size = parameter['min_vacuum_size']
        self.pert_xz = parameter['pert_xz']
        default_max_miller = 2
        self.miller = parameter.get('max_miller', default_max_miller)
        self.static = parameter.get('static-opt', False)
        self.relax = parameter.get('change_box', False)
        self.reprod = parameter.get('reprod-opt', False)

    def make_confs(self,
                   path_to_work,
                   path_to_equi,
                   refine=False):
        path_to_work = os.path.abspath(path_to_work)
        path_to_equi = os.path.abspath(path_to_equi)
        task_list = []
        cwd = os.getcwd()

        equi_contcar = os.path.join(path_to_equi, 'CONTCAR')
        ptypes = vasp.get_poscar_types(equi_contcar)
        if not os.path.exists(equi_contcar):
            raise RuntimeError("please do relaxation first")
        # gen structure
        ss = Structure.from_file(equi_contcar)
        # gen slabs
        all_slabs = generate_all_slabs(ss, self.miller, self.min_slab_size, self.min_vacuum_size)


        if refine:
            task_list = make_refine(self.parameter['init_from_suffix'],
                                    self.parameter['output_suffix'],
                                    path_to_work,
                                    len(all_slabs))
            # record miller
            for ii in range(len(task_list)):
                os.chdir(task_list[ii])
                np.savetxt('miller.out', all_slabs[ii].miller_index, fmt='%d')
            os.chdir(cwd)

        if self.reprod:
            if 'vasp_lmp_path' not in self.parameter:
                raise RuntimeError("please provide the vasp_lmp_path for reproduction")
            vasp_lmp_path = os.path.abspath(self.parameter['vasp_lmp_path'])
            task_list = reproduce.make_repro(vasp_lmp_path, path_to_work)
            os.chdir(cwd)

        else:
            os.chdir(path_to_work)
            if os.path.isfile('POSCAR'):
                os.remove('POSCAR')
            os.symlink(os.path.relpath(equi_contcar), 'POSCAR')
            #           task_poscar = os.path.join(output, 'POSCAR')
            for ii in range(len(all_slabs)):
                output_task = os.path.join(path_to_work, 'task.%06d' % ii)
                os.makedirs(output_task, exist_ok=True)
                os.chdir(output_task)
                for jj in ['INCAR', 'POTCAR', 'POSCAR', 'conf.lmp', 'in.lammps']:
                    if os.path.exists(jj):
                        os.remove(jj)
                task_list.append(output_task)
                print("# %03d generate " % ii, output_task, " \t %d atoms" % len(all_slabs[ii].sites))
                # make confs
                all_slabs[ii].to('POSCAR', 'POSCAR.tmp')
                vasp.regulate_poscar('POSCAR.tmp', 'POSCAR')
                vasp.sort_poscar('POSCAR', 'POSCAR', ptypes)
                vasp.perturb_xz('POSCAR', 'POSCAR', self.pert_xz)
                # record miller
                np.savetxt('miller.out', all_slabs[ii].miller_index, fmt='%d')
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
        ptr_data = os.path.dirname(output_file) + '\n'


        if not self.reprod:
            ptr_data += "Miller_Indices: \tSurf_E(J/m^2) EpA(eV) equi_EpA(eV)\n"
            for ii in all_tasks:
                with open(os.path.join(ii, 'inter.json')) as fp:
                    idata = json.load(fp)
                inter_type = idata['type']
                equi_path = os.path.abspath(os.path.join(os.path.dirname(output_file), '../relaxation'))
                structure_dir = os.path.basename(ii)
                if inter_type == 'vasp':
                    equi_outcar = os.path.join(equi_path, 'OUTCAR')
                    equi_natoms, equi_epa, equi_vpa = vasp.get_nev(equi_outcar)
                    outcar = os.path.join(ii, 'OUTCAR')
                    natoms, epa, vpa = vasp.get_nev(outcar)
                    if self.static:
                        e0 = np.array(vasp.get_energies(outcar)) / natoms
                        epa = e0[0]
                    boxes = vasp.get_boxes(outcar)
                    AA = np.linalg.norm(np.cross(boxes[0][0], boxes[0][1]))

                elif inter_type in ['deepmd', 'meam', 'eam_fs', 'eam_alloy']:
                    equi_log = os.path.join(equi_path, 'log.lammps')
                    equi_natoms, equi_epa, equi_vpa = lammps.get_nev(equi_log)
                    lmp_log = os.path.join(ii, 'log.lammps')
                    natoms, epa, vpa = lammps.get_nev(lmp_log)
                    AA = lammps.get_base_area(lmp_log)

                else:
                    raise RuntimeError('interaction type not supported')

                Cf = 1.60217657e-16 / (1e-20 * 2) * 0.001
                evac = (epa * natoms - equi_epa * natoms) / AA * Cf
                miller_index = []
                with open(os.path.join(ii,'miller.out'),'r') as fin:
                    ss = int(fin.readline().split()[0])
                    miller_index.append(ss)

                ptr_data += "%s: \t%7.3f    %8.3f %8.3f\n" % (miller_index, evac, epa, equi_epa)
                res_data[miller_index] = [evac, epa, equi_epa]

        else:
            if 'vasp_lmp_path' not in self.parameter:
                raise RuntimeError("please provide the vasp_lmp_path for reproduction")
            vasp_lmp_path = os.path.abspath(self.parameter['vasp_lmp_path'])
            res_data, ptr_data = reproduce.post_repro(vasp_lmp_path, all_tasks, ptr_data)

        with open(output_file, 'w') as fp:
            json.dump(res_data, fp, indent=4)

        return res_data, ptr_data
