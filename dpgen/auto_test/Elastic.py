from dpgen import dlog
from dpgen.auto_test.Property import Property
from dpgen.auto_test.refine import make_refine
import dpgen.auto_test.lib.vasp as vasp
import dpgen.auto_test.lib.lammps as lammps
from pymatgen.core.structure import Structure
from pymatgen.analysis.elasticity.strain import DeformedStructureSet, Strain
from pymatgen.analysis.elasticity.stress import Stress
from pymatgen.analysis.elasticity.elastic import ElasticTensor
from monty.serialization import loadfn, dumpfn
import numpy as np
import os


class Elastic(Property):
    def __init__(self,
                 parameter):
        self.parameter = parameter
        default_norm_def = 2e-3
        default_shear_def = 5e-3
        self.norm_deform = parameter.get('norm_deform', default_norm_def)
        self.shear_deform = parameter.get('shear_deform', default_shear_def)
        self.cal_type = parameter.get('cal_type', 'relaxation')
        default_cal_setting = {"relax_pos": True,
                               "relax_shape": False,
                               "relax_vol": False}
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

        norm_def = self.norm_deform
        shear_def = self.shear_deform
        norm_strains = [-norm_def, -0.5 * norm_def, 0.5 * norm_def, norm_def]
        shear_strains = [-shear_def, -0.5 * shear_def, 0.5 * shear_def, shear_def]
        print('gen with norm ' + str(norm_strains))
        print('gen with shear ' + str(shear_strains))

        equi_contcar = os.path.join(path_to_equi, 'CONTCAR')
        if not os.path.exists(equi_contcar):
            raise RuntimeError("please do relaxation first")

        ss = Structure.from_file(equi_contcar)
        dfm_ss = DeformedStructureSet(ss,
                                      symmetry=False,
                                      norm_strains=norm_strains,
                                      shear_strains=shear_strains)
        n_dfm = len(dfm_ss)

        os.chdir(path_to_work)
        if os.path.isfile('POSCAR'):
            os.remove('POSCAR')
        os.symlink(os.path.relpath(equi_contcar), 'POSCAR')
        #           task_poscar = os.path.join(output, 'POSCAR')
        # stress
        equi_outcar = os.path.join(path_to_equi, 'OUTCAR')
        equi_log = os.path.join(path_to_equi, 'log.lammps')
        if os.path.exists(equi_outcar):
            stress = vasp.get_stress(equi_outcar)
            np.savetxt('equi.stress.out', stress)
        elif os.path.exists(equi_log):
            stress = lammps.get_stress(equi_log)
            np.savetxt('equi.stress.out', stress)
        os.chdir(cwd)

        if refine:
            task_list = make_refine(self.parameter['init_from_suffix'],
                                    self.parameter['output_suffix'],
                                    path_to_work,
                                    n_dfm)
            os.chdir(cwd)
        else:
            for ii in range(n_dfm):
                output_task = os.path.join(path_to_work, 'task.%06d' % ii)
                os.makedirs(output_task, exist_ok=True)
                os.chdir(output_task)
                for jj in ['INCAR', 'POTCAR', 'POSCAR', 'conf.lmp', 'in.lammps']:
                    if os.path.exists(jj):
                        os.remove(jj)
                task_list.append(output_task)
                dfm_ss.deformed_structures[ii].to('POSCAR', 'POSCAR')
                # record strain
                df = Strain.from_deformation(dfm_ss.deformations[ii])
                dumpfn(df.as_dict(), 'strain.json', indent=4)
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
        equi_stress = Stress(np.loadtxt(os.path.join(os.path.dirname(output_file), 'equi.stress.out')))
        lst_strain = []
        lst_stress = []
        for ii in all_tasks:
            idata = loadfn(os.path.join(ii, 'inter.json'))
            inter_type = idata['type']
            strain = loadfn(os.path.join(ii, 'strain.json'))
            if inter_type == 'vasp':
                stress = vasp.get_stress(os.path.join(ii, 'OUTCAR'))
                # convert from pressure in kB to stress
                stress *= -1000
                lst_strain.append(strain)
                lst_stress.append(Stress(stress))
            elif inter_type in ['deepmd', 'meam', 'eam_fs', 'eam_alloy']:
                stress = lammps.get_stress(os.path.join(ii, 'log.lammps'))
                # convert from pressure to stress
                stress = -stress
                lst_strain.append(strain)
                lst_stress.append(Stress(stress))
        et = ElasticTensor.from_independent_strains(lst_strain, lst_stress, eq_stress=equi_stress, vasp=False)
        res_data['elastic_tensor'] = []
        for ii in range(6):
            for jj in range(6):
                res_data['elastic_tensor'].append(et.voigt[ii][jj] / 1e4)
                ptr_data += "%7.2f " % (et.voigt[ii][jj] / 1e4)
            ptr_data += '\n'

        BV = et.k_voigt / 1e4
        GV = et.g_voigt / 1e4
        EV = 9 * BV * GV / (3 * BV + GV)
        uV = 0.5 * (3 * BV - 2 * GV) / (3 * BV + GV)

        res_data['BV'] = BV
        res_data['GV'] = GV
        res_data['EV'] = EV
        res_data['uV'] = uV
        ptr_data += "# Bulk   Modulus BV = %.2f GPa\n" % BV
        ptr_data += "# Shear  Modulus GV = %.2f GPa\n" % GV
        ptr_data += "# Youngs Modulus EV = %.2f GPa\n" % EV
        ptr_data += "# Poission Ratio uV = %.2f " % uV

        dumpfn(res_data, output_file, indent=4)

        return res_data, ptr_data
