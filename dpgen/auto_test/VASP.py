import os
import json
import warnings
from Task import Task
from dpgen.generator.lib.vasp import incar_upper
from pymatgen.io.vasp import Incar, Kpoints
from dpgen import dlog
import dpgen.auto_test.lib.vasp as vasp


class VASP(Task):
    def __init__(self,
                 inter_parameter,
                 path_to_poscar):
        self.inter = inter_parameter
        self.incar = inter_parameter['incar']
        self.potcars = inter_parameter['potcars']
        default_potcar_prefix = ''
        self.potcar_prefix = inter_parameter.get('potcar_prefix',default_potcar_prefix)
        self.path_to_poscar = path_to_poscar

    def make_potential_files(self,
                             output_dir):
        with open(os.path.join(output_dir, 'POTCAR'), 'w') as fp:
            for ii in self.potcars:
                with open(os.path.join(self.potcar_prefix, self.potcars[ii]), 'r') as fin:
                    for line in fin:
                        print(line.strip('\n'), file=fp)

        with open(os.path.join(output_dir, 'inter.json'), 'w') as fp:
            json.dump(self.inter, fp, indent=4)

    def make_input_file(self,
                        output_dir,
                        task_type,
                        task_param):
        with open(os.path.join(output_dir, 'task.json'), 'w') as fp:
            json.dump(task_param, fp, indent=4)

        assert (os.path.exists(self.incar)), 'no INCAR file for relaxation'
        relax_incar_path = os.path.abspath(self.incar)
        incar = incar_upper(Incar.from_file(relax_incar_path))

        if 'ISIF' in incar:
            isif = incar.get('ISIF')
        else:
            isif = 3

        if 'NSW' in incar:
            nsw = incar.get('NSW')
        else:
            nsw = 200

        try:
            kspacing = incar.get('KSPACING')
        except KeyError:
            raise RuntimeError("KSPACING must be given in INCAR")

        if 'KGAMMA' in incar:
            kgamma = incar.get('KGAMMA')
        else:
            kgamma = False


        if task_type in ['relaxation', 'vacancy', 'interstitial']:
            isif = 3


        if task_type == 'eos':
            if 'change_box' in task_param and not task_param['change_box']:
                isif = 2
            else:
                isif = 4

        if task_type == 'elastic':
            isif = 2

        if task_type == 'surface':
            if 'static-opt' in task_param and task_param['static-opt']:
                nsw = 0
            elif 'change_box' in task_param and task_param['change_box']:
                isif = 4
            else:
                isif = 2

        if task_type == 'static' \
                or ('reprod_opt' in task_param and task_param['reprod_opt']):
            nsw = 0

        if not ('ISIF' in incar and incar.get('ISIF') == isif):
            dlog.info("%s:%s setting ISIF to %d" % (__file__, self.make_input_file.__name__, isif))
            incar['ISIF'] = isif

        if not ('NSW' in incar and incar.get('NSW') == nsw):
            dlog.info("%s:%s setting NSW to %d" % (__file__, self.make_input_file.__name__, nsw))
            incar['NSW'] = nsw


        if 'ediff' in task_param:
            dlog.info("%s:%s setting ediff to %s" % (__file__, self.make_input_file.__name__, task_param['ediff']))
            incar['EDIFF'] = task_param['ediff']

        if 'ediffg' in task_param:
            dlog.info("%s:%s setting ediff to %s" % (__file__, self.make_input_file.__name__, task_param['ediffg']))
            incar['EDIFFG'] = task_param['ediffg']

        if 'encut' in task_param:
            dlog.info("%s:%s setting ediff to %s" % (__file__, self.make_input_file.__name__, task_param['encut']))
            incar['ENCUT'] = task_param['encut']

        if 'kspacing' in task_param:
            dlog.info("%s:%s setting ediff to %s" % (__file__, self.make_input_file.__name__, task_param['kspacing']))
            incar['KSPACING'] = task_param['kspacing']

        if 'kgamma' in task_param:
            dlog.info("%s:%s setting ediff to %s" % (__file__, self.make_input_file.__name__, task_param['kgamma']))
            incar['KGAMMA'] = task_param['kgamma']

        fc = incar.get_string()
        # write incar
        with open(os.path.join(output_dir, 'INCAR'), 'w') as fp:
            fp.write(fc)

        ret = vasp.make_kspacing_kpoints(self.path_to_poscar, kspacing, kgamma)
        kp = Kpoints.from_string(ret)
        kp.write_file(os.path.join(output_dir, "KPOINTS"))

    def compute(self,
                output_dir):
        outcar = os.path.join(output_dir, 'OUTCAR')
        if not os.path.isfile(outcar):
            warnings.warn("cannot find OUTCAR in " + output_dir + " skip")
            return None
        else:
            force = []
            position = []
            energy = []
            with open(outcar, 'r') as fp:
                if 'Elapsed time (sec):' not in fp.read():
                    warnings.warn("incomplete job " + outcar+ " skip")
                    return None
                else:
                    fp.seek(0)
                    for line in fp:
                        if 'TOTAL-FORCE' in line:
                            position.append([])
                            force.append([])
                            fp.readline()
                            while True:
                                ss = fp.readline().split()
                                if len(ss) != 6:
                                    break
                                position[-1].append(float(ss[0]))
                                force[-1].append(float(ss[3]))
                                position[-1].append(float(ss[1]))
                                force[-1].append(float(ss[4]))
                                position[-1].append(float(ss[2]))
                                force[-1].append(float(ss[5]))
                        elif 'free  energy   TOTEN' in line:
                            energy.append(float(line.split()[4]))
        if len(force) > 0 and len(energy) > 0:
            result_dict = {"energy": energy[-1], "force": force[-1]}
            return result_dict

    def forward_files(self):
        return ['INCAR', 'POSCAR', 'POTCAR']

    def forward_common_files(self):
        return ['INCAR', 'POTCAR']

    def backward_files(self):
        return ['OUTCAR', 'VASP.out', 'CONTCAR', 'OSZICAR']
