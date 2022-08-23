import os
from dpgen import dlog
from dpgen.util import sepline
import dpgen.auto_test.lib.abacus as abacus
import dpgen.generator.lib.abacus_scf as abacus_scf
from dpgen.auto_test.Task import Task

from dpdata import LabeledSystem
from monty.serialization import dumpfn
import numpy as np


class ABACUS(Task):
    def __init__(self,
                 inter_parameter,
                 path_to_poscar):
        self.inter = inter_parameter
        self.inter_type = inter_parameter['type']
        self.incar = inter_parameter.get('incar',{})
        self.potcar_prefix = inter_parameter.get('potcar_prefix', '')
        self.potcars = inter_parameter.get('potcars',None)
        self.orbfile = inter_parameter.get('orb_files',None)
        self.deepks  = inter_parameter.get('deepks_desc',None)
        self.path_to_poscar = path_to_poscar
        self.if_define_orb_file = False if self.orbfile == None else True

    def make_potential_files(self,
                             output_dir):
        stru = os.path.abspath(os.path.join(output_dir, 'STRU'))
        if not os.path.isfile(stru):
            raise FileNotFoundError("No file %s" % stru)
        stru_data = abacus_scf.get_abacus_STRU(stru)
        atom_names = stru_data['atom_names']
        orb_files  = stru_data['orb_files']
        pp_files   = stru_data["pp_files"]
        dpks_descriptor = stru_data['dpks_descriptor'] 

        if os.path.islink(os.path.join(output_dir, 'STRU')):
            stru_path,tmpf = os.path.split(os.readlink(os.path.join(output_dir, 'STRU')))
        else:
            stru_path = output_dir

        if pp_files == None:
            raise RuntimeError("No pseudopotential information in STRU file")   

        pp_dir = os.path.abspath(self.potcar_prefix)
        cwd = os.getcwd()
        os.chdir(output_dir)
        if not os.path.isdir("./pp_orb"): os.mkdir("./pp_orb")
        for i in range(len(atom_names)):
            pp_orb_file = [[pp_files[i],self.potcars]]
            if orb_files != None:
                pp_orb_file.append([orb_files[i],self.orbfile])
            elif self.orbfile != None:
                assert(atom_names[i] in self.orbfile),"orb_file of %s is not defined" % atom_names[i]
                pp_orb_file.append([self.orbfile[atom_names[i]],self.orbfile])

            if dpks_descriptor != None:
                pp_orb_file.append([dpks_descriptor[i],self.deepks])
            elif self.deepks != None:
                pp_orb_file.append([self.deepks,self.deepks])

            for tmpf,tmpdict in pp_orb_file:
                atom = atom_names[i]
                if os.path.isfile(os.path.join(stru_path,tmpf)):
                    linked_file = os.path.join(stru_path,tmpf)
                elif tmpdict != None and os.path.isfile(os.path.join(pp_dir,tmpdict[atom])):
                    linked_file = os.path.join(pp_dir,tmpdict[atom])
                else:
                    raise RuntimeError("Can not find file %s" % tmpf.split('/')[-1]) 
                target_file = os.path.join("./pp_orb/",tmpf.split('/')[-1])
                if os.path.isfile(target_file):
                    os.remove(target_file)
                os.symlink(linked_file, target_file)

        os.chdir(cwd)

        dumpfn(self.inter, os.path.join(output_dir, 'inter.json'), indent=4)

    def modify_input(self,incar,x,y):
        if x in incar and incar[x] != y:
            dlog.info("setting %s to %s" % (x,y))
        incar[x] = y

    def make_input_file(self,
                        output_dir,
                        task_type,
                        task_param):
        sepline(ch=output_dir)
        dumpfn(task_param, os.path.join(output_dir, 'task.json'), indent=4)

        assert (os.path.exists(self.incar)), 'no INPUT file for relaxation'
        relax_incar_path = os.path.abspath(self.incar)
        incar_relax = abacus_scf.get_abacus_input_parameters(relax_incar_path)

        # deal with relaxation

        cal_type = task_param['cal_type']
        cal_setting = task_param['cal_setting']

        # user input INCAR for property calculation
        if 'input_prop' in cal_setting and os.path.isfile(cal_setting['input_prop']):
            incar_prop = os.path.abspath(cal_setting['input_prop'])
            incar = abacus_scf.get_abacus_input_parameters(incar_prop)
            dlog.info("Detected 'input_prop' in 'relaxation', use %s as INPUT, and ignore 'cal_setting'" % incar_prop)

        # revise INCAR based on the INCAR provided in the "interaction"
        else:
            incar = incar_relax
            for key in cal_setting:
                if key in ['relax_pos','relax_shape','relax_vol','K_POINTS','']:continue
                if key[0] == '_' : continue
                if 'interaction' in key.lower():continue
                incar[key.lower()] = cal_setting[key]

            fix_atom = [False,False,False]
            if cal_type == 'relaxation':
                relax_pos = cal_setting['relax_pos']
                relax_shape = cal_setting['relax_shape']
                relax_vol = cal_setting['relax_vol']
                if [relax_pos, relax_shape, relax_vol] == [True, False, False]:
                    self.modify_input(incar,'calculation','relax')
                elif [relax_pos, relax_shape, relax_vol] == [True, True, True]:
                    self.modify_input(incar,'calculation','cell-relax')
                elif [relax_pos, relax_shape, relax_vol] == [True, True, False]:
                    self.modify_input(incar,'calculation','cell-relax')
                    self.modify_input(incar,'fixed_axes','volume')
                elif [relax_pos, relax_shape, relax_vol] == [False, True, False]:
                    self.modify_input(incar,'calculation','cell-relax')
                    self.modify_input(incar,'fixed_axes','volume')
                    fix_atom = [True,True,True]
                elif [relax_pos, relax_shape, relax_vol] == [False, True, True]:
                    self.modify_input(incar,'calculation','cell-relax')
                    fix_atom = [True,True,True]
                elif [relax_pos, relax_shape, relax_vol] == [False, False, True]:
                    raise RuntimeError("relax volume but fix shape is not supported for ABACUS")
                elif [relax_pos, relax_shape, relax_vol] == [False, False, False]:
                    self.modify_input(incar,'calculation','scf')
                else:
                    raise RuntimeError("not supported calculation setting for ABACUS")

            elif cal_type == 'static':
                self.modify_input(incar,'calculation','scf')

            else:
                raise RuntimeError("not supported calculation type for ABACUS")

            #modify STRU file base on the value of fix_atom
            abacus.stru_fix_atom(os.path.join(output_dir, 'STRU'),fix_atom)

        if 'basis_type' not in incar:
            dlog.info("'basis_type' is not defined, set to be 'pw'!")
            self.modify_input(incar,'basis_type','pw')
        if 'ntype' not in incar:
            raise RuntimeError("ntype is not defined in INPUT")
        if 'lcao' in incar['basis_type'].lower() and not self.if_define_orb_file:
            mess = "The basis_type is %s, but not define orbital file!!!" % incar['basis_type']
            raise RuntimeError(mess)
        abacus.write_input(os.path.join(output_dir, '../INPUT'),incar)
        cwd = os.getcwd()
        os.chdir(output_dir)
        if not os.path.islink('INPUT'):
            os.symlink('../INPUT', 'INPUT')
        elif not '../INPUT' == os.readlink('INPUT'):
            os.remove('INPUT')
            os.symlink('../INPUT', 'INPUT')
        os.chdir(cwd)

        if 'kspacing' in incar:
            kspacing = float(incar['kspacing'])
            if os.path.isfile(os.path.join(output_dir, 'STRU')):
                kpt = abacus.make_kspacing_kpt(os.path.join(output_dir, 'STRU'),kspacing)
                kpt += [0,0,0]
            else:
                kpt = [1,1,1,0,0,0]
        elif 'K_POINTS' in cal_setting:
            kpt = cal_setting['K_POINTS']
        else:
            mess  = "K point information is not defined\n"
            mess += "You can set key word 'kspacing' (unit in 1/bohr) as a float value in INPUT\n"
            mess += "or set key word 'K_POINTS' as a list in 'cal_setting', e.g. [1,2,3,0,0,0]\n"
            raise RuntimeError(mess)
        abacus.write_kpt(os.path.join(output_dir, 'KPT'),kpt)

    def compute(self,
                output_dir):
        if not os.path.isfile(os.path.join(output_dir,'INPUT')):
            dlog.warning("cannot find INPUT in " + output_dir + " skip")   
            return None     
        ls = LabeledSystem(output_dir,fmt='abacus/relax')
        outcar_dict = ls.as_dict()
        return outcar_dict

    def forward_files(self, property_type='relaxation'):
        return ['INPUT', 'STRU', 'KPT', 'pp_orb']

    def forward_common_files(self, property_type='relaxation'):
        return []

    def backward_files(self, property_type='relaxation'):
        return []
