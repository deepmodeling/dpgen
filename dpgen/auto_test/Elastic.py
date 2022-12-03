import glob
import os
from shutil import copyfile
import re

from monty.serialization import loadfn, dumpfn
from pymatgen.analysis.elasticity.elastic import ElasticTensor
from pymatgen.analysis.elasticity.strain import DeformedStructureSet, Strain
from pymatgen.analysis.elasticity.stress import Stress
from pymatgen.core.structure import Structure
from pymatgen.io.vasp import Incar, Kpoints

import dpgen.auto_test.lib.vasp as vasp
from dpgen import dlog
from dpgen.auto_test.Property import Property
from dpgen.auto_test.refine import make_refine
from dpgen.generator.lib.vasp import incar_upper

import dpgen.auto_test.lib.abacus as abacus
import dpgen.generator.lib.abacus_scf as abacus_scf

class Elastic(Property):
    def __init__(self,
                 parameter,inter_param=None):
        if not ('init_from_suffix' in parameter and 'output_suffix' in parameter):
            default_norm_def = 1e-2
            default_shear_def = 1e-2
            parameter['norm_deform'] = parameter.get('norm_deform', default_norm_def)
            self.norm_deform = parameter['norm_deform']
            parameter['shear_deform'] = parameter.get('shear_deform', default_shear_def)
            self.shear_deform = parameter['shear_deform']
        parameter['cal_type'] = parameter.get('cal_type', 'relaxation')
        self.cal_type = parameter['cal_type']
        default_cal_setting = {"relax_pos": True,
                               "relax_shape": False,
                               "relax_vol": False}
        if 'cal_setting' not in parameter:
            parameter['cal_setting'] = default_cal_setting
        else:
            if "relax_pos" not in parameter['cal_setting']:
                parameter['cal_setting']['relax_pos'] = default_cal_setting['relax_pos']
            if "relax_shape" not in parameter['cal_setting']:
                parameter['cal_setting']['relax_shape'] = default_cal_setting['relax_shape']
            if "relax_vol" not in parameter['cal_setting']:
                parameter['cal_setting']['relax_vol'] = default_cal_setting['relax_vol']
        self.cal_setting = parameter['cal_setting']
        # parameter['reproduce'] = False
        # self.reprod = parameter['reproduce']
        self.parameter = parameter
        self.inter_param = inter_param if inter_param != None else {'type': 'vasp'}

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

        if 'start_confs_path' in self.parameter and os.path.exists(self.parameter['start_confs_path']):
            init_path_list = glob.glob(os.path.join(self.parameter['start_confs_path'], '*'))
            struct_init_name_list = []
            for ii in init_path_list:
                struct_init_name_list.append(ii.split('/')[-1])
            struct_output_name = path_to_work.split('/')[-2]
            assert struct_output_name in struct_init_name_list
            path_to_equi = os.path.abspath(os.path.join(self.parameter['start_confs_path'],
                                                        struct_output_name, 'relaxation', 'relax_task'))

        task_list = []
        cwd = os.getcwd()

        if self.inter_param['type'] == 'abacus':
            CONTCAR = abacus.final_stru(path_to_equi)
            POSCAR = 'STRU'
        else:
            CONTCAR = 'CONTCAR'
            POSCAR = 'POSCAR'

        equi_contcar = os.path.join(path_to_equi, CONTCAR)

        os.chdir(path_to_work)
        if os.path.isfile(POSCAR):
            os.remove(POSCAR)
        if os.path.islink(POSCAR):
            os.remove(POSCAR)
        os.symlink(os.path.relpath(equi_contcar), POSCAR)
        #           task_poscar = os.path.join(output, 'POSCAR')

        # stress, deal with unsupported stress in dpdata
        #with open(os.path.join(path_to_equi, 'result.json')) as fin:
        #    equi_result = json.load(fin)
        #equi_stress = np.array(equi_result['stress']['data'])[-1]
        equi_result = loadfn(os.path.join(path_to_equi, 'result.json'))
        equi_stress = equi_result['stress'][-1]
        dumpfn(equi_stress, 'equi.stress.json', indent=4)
        os.chdir(cwd)

        if refine:
            print('elastic refine starts')
            task_list = make_refine(self.parameter['init_from_suffix'],
                                    self.parameter['output_suffix'],
                                    path_to_work)

            # record strain
            # df = Strain.from_deformation(dfm_ss.deformations[idid])
            # dumpfn(df.as_dict(), 'strain.json', indent=4)
            init_from_path = re.sub(self.parameter['output_suffix'][::-1],
                                    self.parameter['init_from_suffix'][::-1],
                                    path_to_work[::-1], count=1)[::-1]
            task_list_basename = list(map(os.path.basename, task_list))

            for ii in task_list_basename:
                init_from_task = os.path.join(init_from_path, ii)
                output_task = os.path.join(path_to_work, ii)
                os.chdir(output_task)
                if os.path.isfile('strain.json'):
                    os.remove('strain.json')
                copyfile(os.path.join(init_from_task, 'strain.json'), 'strain.json')
                #os.symlink(os.path.relpath(
                #    os.path.join((re.sub(self.parameter['output_suffix'], self.parameter['init_from_suffix'], ii)),
                #                 'strain.json')),
                #           'strain.json')
            os.chdir(cwd)
        else:
            norm_def = self.norm_deform
            shear_def = self.shear_deform
            norm_strains = [-norm_def, -0.5 * norm_def, 0.5 * norm_def, norm_def]
            shear_strains = [-shear_def, -0.5 * shear_def, 0.5 * shear_def, shear_def]

            if not os.path.exists(equi_contcar):
                raise RuntimeError("please do relaxation first")

            if self.inter_param['type'] == 'abacus':
                ss = abacus.stru2Structure(equi_contcar)
            else:
                ss = Structure.from_file(equi_contcar)
            dfm_ss = DeformedStructureSet(ss,
                                          symmetry=False,
                                          norm_strains=norm_strains,
                                          shear_strains=shear_strains)
            n_dfm = len(dfm_ss)

            print('gen with norm ' + str(norm_strains))
            print('gen with shear ' + str(shear_strains))
            for ii in range(n_dfm):
                output_task = os.path.join(path_to_work, 'task.%06d' % ii)
                os.makedirs(output_task, exist_ok=True)
                os.chdir(output_task)
                for jj in ['INCAR', 'POTCAR', 'POSCAR', 'conf.lmp', 'in.lammps','STRU']:
                    if os.path.exists(jj):
                        os.remove(jj)
                task_list.append(output_task)
                dfm_ss.deformed_structures[ii].to('POSCAR', 'POSCAR')
                if self.inter_param['type'] == 'abacus':
                    abacus.poscar2stru("POSCAR",self.inter_param,"STRU")
                    os.remove('POSCAR')
                # record strain
                df = Strain.from_deformation(dfm_ss.deformations[ii])
                dumpfn(df.as_dict(), 'strain.json', indent=4)
            os.chdir(cwd)
        return task_list

    def post_process(self, task_list):
        if self.inter_param['type'] == 'abacus':
            POSCAR = 'STRU'
            INCAR = 'INPUT'
            KPOINTS = 'KPT'
        else:
            POSCAR = 'POSCAR'
            INCAR = 'INCAR'
            KPOINTS = 'KPOINTS'

        cwd = os.getcwd()
        poscar_start = os.path.abspath(os.path.join(task_list[0], '..', POSCAR))
        os.chdir(os.path.join(task_list[0], '..'))
        if os.path.isfile(os.path.join(task_list[0], INCAR)):
            if self.inter_param['type'] == 'abacus':
                input_aba = abacus_scf.get_abacus_input_parameters('INPUT')
                if 'kspacing' in input_aba:
                    kspacing = float(input_aba['kspacing'])
                    kpt = abacus.make_kspacing_kpt(poscar_start,kspacing)
                    kpt += [0,0,0]
                    abacus.write_kpt('KPT',kpt)
                    del input_aba['kspacing']
                    os.remove('INPUT')
                    abacus.write_input('INPUT',input_aba)
                else:
                    os.rename(os.path.join(task_list[0], 'KPT'),'./KPT')
            else:
                incar = incar_upper(Incar.from_file(os.path.join(task_list[0], 'INCAR')))
                kspacing = incar.get('KSPACING')
                kgamma = incar.get('KGAMMA', False)
                ret = vasp.make_kspacing_kpoints(poscar_start, kspacing, kgamma)
                kp = Kpoints.from_string(ret)
                if os.path.isfile('KPOINTS'):
                    os.remove('KPOINTS')
                kp.write_file("KPOINTS")

            os.chdir(cwd)
            kpoints_universal = os.path.abspath(os.path.join(task_list[0], '..', KPOINTS))
            for ii in task_list:
                if os.path.isfile(os.path.join(ii, KPOINTS)):
                    os.remove(os.path.join(ii, KPOINTS))
                if os.path.islink(os.path.join(ii, KPOINTS)):
                    os.remove(os.path.join(ii, KPOINTS))
                os.chdir(ii)
                os.symlink(os.path.relpath(kpoints_universal), KPOINTS)

        os.chdir(cwd)

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
        equi_stress = Stress(loadfn(os.path.join(os.path.dirname(output_file), 'equi.stress.json')))
        equi_stress *= -1000
        lst_strain = []
        lst_stress = []
        for ii in all_tasks:
            strain = loadfn(os.path.join(ii, 'strain.json'))
            # stress, deal with unsupported stress in dpdata
            #with open(os.path.join(ii, 'result_task.json')) as fin:
            #    task_result = json.load(fin)
            #stress = np.array(task_result['stress']['data'])[-1]
            stress = loadfn(os.path.join(ii, 'result_task.json'))['stress'][-1]
            lst_strain.append(strain)
            lst_stress.append(Stress(stress * -1000))

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
        ptr_data += "# Poission Ratio uV = %.2f\n " % uV

        dumpfn(res_data, output_file, indent=4)

        return res_data, ptr_data
