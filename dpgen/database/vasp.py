#/usr/bin/env python
# coding: utf-8
# Copyright (c) PThe Dpmodeling Team.


import os
import warnings
from monty.io import zopen
from monty.os.path import zpath
from monty.json import MSONable, MontyDecoder
from pymatgen.io.vasp import Potcar,Incar,Kpoints,Poscar,PotcarSingle 
"""
Classes for reading/manipulating/writing VASP input files. All major VASP input
files.
"""

class DPPotcar(MSONable):
    def __init__(self,symbols=None,functional="PBE",pp_file=None,pp_lists=None):
        if pp_lists and pp_file is None:
           for pp in pp_lists: 
               assert isinstance(pp,PotcarSingle)  
           self.potcars=pp_lists  
        elif  pp_file and pp_list is None:
           self.potcars=Potcar.from_file(pp_file)
        elif pp_file and pp_list:
           self.potcars=Potcar.from_file(pp_file)
        else:
           try:
              self.potcars=Potcar(symbols=symbols, functional=functional)
           except Exception:
              warnings.warn ("""Inproperly configure of POTCAR !""")
              self.potcars=None
         
        if self.potcars is not None:
           self.symbols  = [pp.symbol for pp in self.potcars]
           self.functional = list(set([pp.functional for pp in self.potcars]))[0]
           self.hashs  = [pp.get_potcar_hash() for pp in self.potcars]
        else:
           self.symbols=symbols
           self.functional=functional
           self.hashs = ''
        self.elements = self._get_elements()
   
    def __repr__(self):
        return str(self)

    def __str__(self):
        if self.potcars is not None:
           return str(self.potcars)
        else:
           ret  ="Functional: %s\n"%self.functional
           ret +=" ".join(self.symbols)+"\n"
           return ret

    def _get_elements(self):
        elements=[]
        for el in self.symbols:
            if '_' in el:
               elements.append(el.split('_')[0])
            else:
               elements.append(el)
        return elements

    @classmethod
    def from_dict(cls,d):
        return cls(symbols=d['symbols'],functional=d['functional'])

    def as_dict(self):
        d={}
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        d['symbols']=self.symbols
        d['elements']=self.elements
        d['hashs']=self.hashs
        d['functional']=self.functional
        return d

    @classmethod
    def from_file(cls,filename):
        try:
             potcars=Potcar.from_file(filename)
             return cls(pp_lists=potcars)
        except Exception:
             with open(filename,'r') as f:
                  content=f.readlines()
             functional=content[0].strip().split(':')[-1].strip()
             symbols=content[1].strip().split()
             return cls(symbols=symbols,functional=functional)
          

    def write_file(self,filename):
        with open(filename,'w') as f:
             f.write(str(self))


class VaspInput(dict, MSONable):
    """
    Class to contain a set of vasp input objects corresponding to a run.

    Args:
        incar: Incar object.
        kpoints: Kpoints object.
        poscar: Poscar object.
        potcar: Potcar object.
        optional_files: Other input files supplied as a dict of {
            filename: object}. The object should follow standard pymatgen
            conventions in implementing a as_dict() and from_dict method.
    """

    def __init__(self, incar, poscar, potcar, kpoints=None, optional_files=None,
                 **kwargs):
        super().__init__(**kwargs)
        
        self.update({'INCAR': incar,
                     'POSCAR': poscar,
                     'POTCAR': potcar})
        if kpoints:
           self.update({'KPOINTS': kpoints})
        if optional_files is not None:
            self.update(optional_files)

    def __str__(self):
        output = []
        for k, v in self.items():
            output.append(k)
            output.append(str(v))
            output.append("")
        return "\n".join(output)
 
    def as_dict(self):
        d = {k: v.as_dict() for k, v in self.items()}
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        return d

    @classmethod
    def from_dict(cls, d):
        dec = MontyDecoder()
        sub_d = {"optional_files": {}}
        for k, v in d.items():
            if k in ["INCAR", "POSCAR", "POTCAR", "KPOINTS"]:
                sub_d[k.lower()] = dec.process_decoded(v)
            elif k not in ["@module", "@class"]:
                sub_d["optional_files"][k] = dec.process_decoded(v)
        return cls(**sub_d)

    def write_input(self, output_dir=".", make_dir_if_not_present=True):
        """
        Write VASP input to a directory.

        Args:
            output_dir (str): Directory to write to. Defaults to current
                directory (".").
            make_dir_if_not_present (bool): Create the directory if not
                present. Defaults to True.
        """
        if make_dir_if_not_present and not os.path.exists(output_dir):
            os.makedirs(output_dir)
        for k, v in self.items():
            with zopen(os.path.join(output_dir, k), "wt") as f:
                f.write(v.__str__())

    @staticmethod
    def from_directory(input_dir, optional_files=None):
        """
        Read in a set of VASP input from a directory. Note that only the
        standard INCAR, POSCAR, POTCAR and KPOINTS files are read unless
        optional_filenames is specified.

        Args:
            input_dir (str): Directory to read VASP input from.
            optional_files (dict): Optional files to read in as well as a
                dict of {filename: Object type}. Object type must have a
                static method from_file.
        """
        sub_d = {}
        try:
            for fname, ftype in [("INCAR", Incar), ("KPOINTS", Kpoints),
                                 ("POSCAR", Poscar), ("POTCAR", DPPotcar)]:
                fullzpath = zpath(os.path.join(input_dir, fname))
                sub_d[fname.lower()] = ftype.from_file(fullzpath)
        except Exception:
            for fname, ftype in [("INCAR", Incar), 
                                 ("POSCAR", Poscar), ("POTCAR", DPPotcar)]:
                fullzpath = zpath(os.path.join(input_dir, fname))
                sub_d[fname.lower()] = ftype.from_file(fullzpath)

        sub_d["optional_files"] = {}
        if optional_files is not None:
            for fname, ftype in optional_files.items():
                sub_d["optional_files"][fname] = \
                    ftype.from_file(os.path.join(input_dir, fname))
        return VaspInput(**sub_d)

