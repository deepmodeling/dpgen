import os,sys,shutil
import unittest
import json
import numpy as np
import tarfile
from glob import glob

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
__package__ = 'database'
from .context import dpgen
from .context import Entry
from .context import VaspInput,DPPotcar
from .context import parsing_vasp
from dpdata import System,LabeledSystem
from monty.shutil import remove 
from monty.serialization import loadfn,dumpfn
from pymatgen.io.vasp import Potcar,Poscar,Incar,Kpoints
from .context import setUpModule

iter_pat="02.fp/task.007.00000*"
init_pat="al.bcc.02x02x02/02.md/sys-0016/scale-1.000/00000*"

def tar_file(path,outdir='.'):
    tar = tarfile.open(path)
    names = tar.getnames()
    for name in names:
      tar.extract(name,path=outdir)
    tar.close()


class Test(unittest.TestCase):
    def setUp(self):
        self.cwd=os.getcwd()
        self.r_init_path=os.path.join(self.cwd,'init')
        self.data=os.path.join(self.cwd,'data')
        self.r_iter_path=os.path.join(self.cwd,'iter.000000')
        for path in [self.r_init_path, self.r_iter_path]+[self.data]:
           if os.path.isdir(path) :
              shutil.rmtree(path)
           tar_file(path+'.tar.gz')
           assert os.path.isdir(path)
        self.ref_init_input=loadfn(os.path.join(self.cwd,'data/init_input.json'))
        self.ref_entries=loadfn(os.path.join(self.cwd,'data/entries.json'))
        self.init_path=sorted(glob(os.path.join(self.r_init_path,init_pat)))
        self.iter_path=sorted(glob(os.path.join(self.r_iter_path,iter_pat)))
        with open("param_Al.json", "r") as fr:
          jdata = json.load(fr)
        self.config_info_dict = jdata["config_info_dict"]
        self.skip_init = jdata["skip_init"]
        self.output = jdata["output"]

    def testDPPotcar(self):
            
        refd={'@module': 'dpgen.database.vasp',
              '@class': 'DPPotcar',
              'symbols': ['Al'],
              'elements': ['Al'],
              'hashs': '',
              'functional': 'PBE'}
        try:
           Potcar(['Al'])
           #ps  TITEL  = PAW_PBE Al 04Jan2001
           refd.update({'hashs':['9aafba2c552fad8414179cae2e888e67']})
        except Exception:
           pass

        for f in self.init_path+self.iter_path:
            fpp=os.path.join(f,'POTCAR')
            pp=DPPotcar.from_file(fpp)
            self.assertEqual( pp.elements, refd['elements'])
            self.assertEqual( pp.symbols, refd['symbols'])
            self.assertEqual( pp.hashs, refd['hashs'])
            self.assertEqual( pp.functional, refd['functional'])
            self.assertEqual( pp.as_dict(), refd)
       
    def testVaspInput(self):
        for f in self.init_path:
            vi=VaspInput.from_directory(f)
            self.assertEqual(vi['INCAR'],self.ref_init_input['INCAR'])
            self.assertEqual(str(vi['POTCAR']),str(self.ref_init_input['POTCAR']))
            self.assertEqual(vi['POSCAR'].structure,self.ref_init_input['POSCAR'].structure)

    def testEntry(self):
        entries=[]
        for i,f in enumerate(self.iter_path):
            vi=VaspInput.from_directory(f)
            ls=LabeledSystem(os.path.join(f,'OUTCAR'))
            attrib=loadfn(os.path.join(f,'job.json'))
            comp=vi['POSCAR'].structure.composition
            entry=Entry(comp,
                       'vasp',
                       vi.as_dict(),
                       ls.as_dict(),
                       entry_id='pku-'+str(i),
                       attribute=attrib)
            entries.append(entry)
        self.assertEqual( len(entries), len(self.ref_entries))
        ret0=entries[0]
        r0=self.ref_entries[0]
        self.assertEqual(
                     Incar.from_dict(ret0.inputs['INCAR']),
                     Incar.from_dict(r0.inputs['INCAR'])
                    )
        self.assertEqual(
                     str(r0.inputs['KPOINTS']),
                     str(Kpoints.from_dict(ret0.inputs['KPOINTS']))
                    )

        self.assertEqual(
                        ret0.inputs['POTCAR'],
                        r0.inputs['POTCAR'].as_dict()
                        )
        self.assertEqual(
                        Poscar.from_dict(ret0.inputs['POSCAR']).structure,
                        r0.inputs['POSCAR'].structure
                        )
        self.assertEqual(ret0.entry_id,'pku-0')

    def testParsingVasp(self):
        parsing_vasp(self.cwd, self.config_info_dict, self.skip_init,self.output, id_prefix=dpgen.SHORT_CMD )
        #try:
        #   Potcar(['Al'])
        #   ref=os.path.join(self.cwd,'data/all_data_pp.json')
        #except Exception:
        #   ref=os.path.join(self.cwd,'data/all_data.json')
        #Potcar(['Al'])
        ref=os.path.join(self.cwd,'data/all_data_pp.json')

        ret=os.path.join(self.cwd,'dpgen_db.json')

        retd=loadfn(ret)
        retd=sorted(retd,key= lambda x: int(x.entry_id.split('_')[-1]))

        refd=loadfn(ref)
        refd=sorted(refd,key= lambda x: int(x.entry_id.split('_')[-1]))
        self.assertEqual(len(retd),len(refd)) 
        for i,j in zip(retd,refd):
            self.assertEqual(i.entry_id,j.entry_id)
            self.assertEqual(i.calculator,j.calculator)
            self.assertEqual(len(i.data),len(j.data))
            self.assertEqual(i.number_element , j.number_element )
            self.assertEqual(i.entry_id , j.entry_id )
            self.assertEqual(len(i.composition),len(j.composition))
            self.assertEqual(len(i.attribute),len(j.attribute))
        os.remove(os.path.join(self.cwd,'dpgen_db.json'))
        

    def tearDown(self):
        for path in [self.r_init_path, self.r_iter_path, self.data]:
           if os.path.isdir(path) :
              shutil.rmtree(path)
        if os.path.isfile("dpgen.log"):
          os.remove("dpgen.log")
        if os.path.isfile("record.database"):
          os.remove("record.database")
        
