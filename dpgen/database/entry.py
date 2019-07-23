#/usr/bin/env python
# coding: utf-8
# Copyright (c) The Dpmodeling Team.

import json
import warnings
from uuid import uuid4
from dpdata import System,LabeledSystem
from dpgen.database.vasp import VaspInput
from pymatgen.core.composition import Composition
from monty.json import MontyEncoder, MontyDecoder, MSONable

"""
This module implements equivalents of the basic Entry objects, which
is the basic entity that can be used to perform many analyses. Entries
contain calculated information, typically from VASP or other electronic
structure codes. For example, Entries can be used as inputs for DeepMD-Kit.
"""


class Entry(MSONable):
    """
    An lightweight Entry object containing key computed data
    for storing purpose. 

    """

    def __init__(self, composition, calculator, inputs,
                 data, entry_id=None, attribute=None, tag=None):
        """
        Initializes a Entry.

        Args:
            composition (Composition): Composition of the entry. For
                flexibility, this can take the form of all the typical input
                taken by a Composition, including a {symbol: amt} dict,
                a string formula, and others.
            inputs (dict): An dict of parameters associated with
                the entry. Defaults to None.
            data (dict): An dict of any additional data associated
                with the entry. Defaults to None.
            entry_id (obj): An optional id to uniquely identify the entry.
            attribute: Optional attribute of the entry. This can be used to
                specify that the entry is a newly found compound, or to specify
                a particular label for the entry, or else ... Used for further
                analysis and plotting purposes. An attribute can be anything
                but must be MSONable.
        """
        self.composition = Composition(composition)
        self.calculator  = calculator
        self.inputs = inputs
        self.data = data 
        self.entry_id = entry_id
        self.name = self.composition.reduced_formula
        self.attribute = attribute
        self.tag = tag

    #def __eq__(self,other):
    #    if not self.composition == other.composition:
    #       return False
    #    if not self.calculator == other.calculator:
    #       return False
    #    if not self.inputs == other.inputs:
    #       return False
    #    if not self.data == other.data:
    #       return False
    #    if not self.name == other.name:
    #       return False
    #    if not self.attribute  == other.attribute:
    #       return False
    #    if not self.tag == other.tag:
    #       return False
    #    return True

    @property
    def number_element(self):
        return len(self.composition)

    def __repr__(self):
        output = ["Entry {} - {}".format(self.entry_id, self.composition.formula),
                  "calculator: {}".format(self.calculator) 
                  ]
        return "\n".join(output)

    def __str__(self):
        return self.__repr__()

    @classmethod
    def from_dict(cls, d):
        dec = MontyDecoder()
        return cls(d["composition"], d["calculator"], 
                   inputs={k: dec.process_decoded(v)
                               for k, v in d.get("inputs", {}).items()},
                   data={k: dec.process_decoded(v)
                         for k, v in d.get("data", {}).items()},
                   entry_id=d.get("entry_id", None),
                   attribute=d["attribute"] if "attribute" in d else None,
                   tag=d["tag"] if "tag" in d else None
                   )

    def as_dict(self):
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "composition": self.composition.as_dict(),
                "calculator": self.calculator,
                "inputs": json.loads(json.dumps(self.inputs,
                                                    cls=MontyEncoder)),
                "data": json.loads(json.dumps(self.data, cls=MontyEncoder)),
                "entry_id": self.entry_id,
                "attribute": self.attribute,
                "tag": self.tag}
