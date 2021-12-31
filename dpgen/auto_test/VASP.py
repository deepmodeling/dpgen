import os
from dpgen import dlog
from dpgen.util import sepline
import dpgen.auto_test.lib.vasp as vasp
from dpgen.auto_test.Task import Task
from dpgen.generator.lib.vasp import incar_upper
from dpdata import LabeledSystem
from monty.serialization import dumpfn
from pymatgen.io.vasp import Incar, Kpoints
from pymatgen.core.structure import Structure


class VASP(Task):
    def __init__(self,
                 inter_parameter,
                 path_to_poscar):
        self.inter = inter_parameter
        self.inter_type = inter_parameter['type']
        self.incar = inter_parameter['incar']
        self.potcar_prefix = inter_parameter.get('potcar_prefix', '')
        self.potcars = inter_parameter['potcars']
        self.path_to_poscar = path_to_poscar

    def make_potential_files(self,
                             output_dir):
        potcar_not_link_list = ['vacancy', 'interstitial']
        task_type = output_dir.split('/')[-2].split('_')[0]

        ele_pot_list = [key for key in self.potcars.keys()]
        poscar = os.path.abspath(os.path.join(output_dir, 'POSCAR'))
        pos_str = Structure.from_file(poscar)
        ele_pos_list_tmp = list(ii.as_dict()['element'] for ii in pos_str.species)

        ele_pos_list = [ele_pos_list_tmp[0]]
        for ii in range(1, len(ele_pos_list_tmp)):
            if not ele_pos_list_tmp[ii] == ele_pos_list_tmp[ii - 1]:
                ele_pos_list.append(ele_pos_list_tmp[ii])

        if task_type in potcar_not_link_list:
            with open(os.path.join(output_dir, 'POTCAR'), 'w') as fp:
                for ii in ele_pos_list:
                    for jj in ele_pot_list:
                        if ii == jj:
                            with open(os.path.join(self.potcar_prefix, self.potcars[jj]), 'r') as fin:
                                for line in fin:
                                    print(line.strip('\n'), file=fp)

        else:
            if not os.path.isfile(os.path.join(output_dir, '../POTCAR')):
                with open(os.path.join(output_dir, '../POTCAR'), 'w') as fp:
                    for ii in ele_pos_list:
                        for jj in ele_pot_list:
                            if ii == jj:
                                with open(os.path.join(self.potcar_prefix, self.potcars[jj]), 'r') as fin:
                                    for line in fin:
                                        print(line.strip('\n'), file=fp)
            cwd = os.getcwd()
            os.chdir(output_dir)
            if not os.path.islink('POTCAR'):
                os.symlink('../POTCAR', 'POTCAR')
            elif not '../POTCAR' == os.readlink('POTCAR'):
                os.remove('POTCAR')
                os.symlink('../POTCAR', 'POTCAR')
            os.chdir(cwd)

        dumpfn(self.inter, os.path.join(output_dir, 'inter.json'), indent=4)

    def make_input_file(self,
                        output_dir,
                        task_type,
                        task_param):
        sepline(ch=output_dir)
        dumpfn(task_param, os.path.join(output_dir, 'task.json'), indent=4)

        assert (os.path.exists(self.incar)), 'no INCAR file for relaxation'
        relax_incar_path = os.path.abspath(self.incar)
        incar_relax = incar_upper(Incar.from_file(relax_incar_path))

        # deal with relaxation

        cal_type = task_param['cal_type']
        cal_setting = task_param['cal_setting']

        # user input INCAR for property calculation
        if 'input_prop' in cal_setting and os.path.isfile(cal_setting['input_prop']):
            incar_prop = os.path.abspath(cal_setting['input_prop'])
            incar = incar_upper(Incar.from_file(incar_prop))

        # revise INCAR based on the INCAR provided in the "interaction"
        else:
            incar = incar_relax
            if cal_type == 'relaxation':
                relax_pos = cal_setting['relax_pos']
                relax_shape = cal_setting['relax_shape']
                relax_vol = cal_setting['relax_vol']
                if [relax_pos, relax_shape, relax_vol] == [True, False, False]:
                    isif = 2
                elif [relax_pos, relax_shape, relax_vol] == [True, True, True]:
                    isif = 3
                elif [relax_pos, relax_shape, relax_vol] == [True, True, False]:
                    isif = 4
                elif [relax_pos, relax_shape, relax_vol] == [False, True, False]:
                    isif = 5
                elif [relax_pos, relax_shape, relax_vol] == [False, True, True]:
                    isif = 6
                elif [relax_pos, relax_shape, relax_vol] == [False, False, True]:
                    isif = 7
                elif [relax_pos, relax_shape, relax_vol] == [False, False, False]:
                    nsw = 0
                    isif = 2
                    if not ('NSW' in incar and incar.get('NSW') == nsw):
                        dlog.info("%s setting NSW to %d" % (self.make_input_file.__name__, nsw))
                        incar['NSW'] = nsw
                else:
                    raise RuntimeError("not supported calculation setting for VASP")

                if not ('ISIF' in incar and incar.get('ISIF') == isif):
                    dlog.info("%s setting ISIF to %d" % (self.make_input_file.__name__, isif))
                    incar['ISIF'] = isif

            elif cal_type == 'static':
                nsw = 0
                if not ('NSW' in incar and incar.get('NSW') == nsw):
                    dlog.info("%s setting NSW to %d" % (self.make_input_file.__name__, nsw))
                    incar['NSW'] = nsw

            else:
                raise RuntimeError("not supported calculation type for VASP")

            if 'ediff' in cal_setting:
                dlog.info("%s setting EDIFF to %s" % (self.make_input_file.__name__, cal_setting['ediff']))
                incar['EDIFF'] = cal_setting['ediff']

            if 'ediffg' in cal_setting:
                dlog.info("%s setting EDIFFG to %s" % (self.make_input_file.__name__, cal_setting['ediffg']))
                incar['EDIFFG'] = cal_setting['ediffg']

            if 'encut' in cal_setting:
                dlog.info("%s setting ENCUT to %s" % (self.make_input_file.__name__, cal_setting['encut']))
                incar['ENCUT'] = cal_setting['encut']

            if 'kspacing' in cal_setting:
                dlog.info("%s setting KSPACING to %s" % (self.make_input_file.__name__, cal_setting['kspacing']))
                incar['KSPACING'] = cal_setting['kspacing']

            if 'kgamma' in cal_setting:
                dlog.info("%s setting KGAMMA to %s" % (self.make_input_file.__name__, cal_setting['kgamma']))
                incar['KGAMMA'] = cal_setting['kgamma']

        try:
            kspacing = incar.get('KSPACING')
        except KeyError:
            raise RuntimeError("KSPACING must be given in INCAR")

        if 'KGAMMA' in incar:
            kgamma = incar.get('KGAMMA')
        else:
            kgamma = False

        incar.write_file(os.path.join(output_dir, '../INCAR'))
        cwd = os.getcwd()
        os.chdir(output_dir)
        if not os.path.islink('INCAR'):
            os.symlink('../INCAR', 'INCAR')
        elif not '../INCAR' == os.readlink('INCAR'):
            os.remove('INCAR')
            os.symlink('../INCAR', 'INCAR')
        os.chdir(cwd)
        ret = vasp.make_kspacing_kpoints(self.path_to_poscar, kspacing, kgamma)
        kp = Kpoints.from_string(ret)
        kp.write_file(os.path.join(output_dir, "KPOINTS"))

    def compute(self,
                output_dir):
        outcar = os.path.join(output_dir, 'OUTCAR')
        if not os.path.isfile(outcar):
            dlog.warning("cannot find OUTCAR in " + output_dir + " skip")
            return None
        else:
            ls = LabeledSystem(outcar)
            stress = []
            with open(outcar, 'r') as fin:
                lines = fin.read().split('\n')
            for line in lines:
                if 'in kB' in line:
                    stress_xx = float(line.split()[2])
                    stress_yy = float(line.split()[3])
                    stress_zz = float(line.split()[4])
                    stress_xy = float(line.split()[5])
                    stress_yz = float(line.split()[6])
                    stress_zx = float(line.split()[7])
                    stress.append([])
                    stress[-1].append([stress_xx, stress_xy, stress_zx])
                    stress[-1].append([stress_xy, stress_yy, stress_yz])
                    stress[-1].append([stress_zx, stress_yz, stress_zz])

            outcar_dict = ls.as_dict()
            outcar_dict['data']['stress'] = {"@module": "numpy", "@class": "array", "dtype": "float64", "data": stress}

            return outcar_dict

    def forward_files(self, property_type='relaxation'):
        return ['INCAR', 'POSCAR', 'KPOINTS', 'POTCAR']

    def forward_common_files(self, property_type='relaxation'):
        potcar_not_link_list = ['vacancy', 'interstitial']
        if property_type == 'elastic':
            return ['INCAR', 'KPOINTS', 'POTCAR']
        elif property_type in potcar_not_link_list:
            return ['INCAR']
        else:
            return ['INCAR', 'POTCAR']

    def backward_files(self, property_type='relaxation'):
        return ['OUTCAR', 'outlog', 'CONTCAR', 'OSZICAR', 'XDATCAR']
