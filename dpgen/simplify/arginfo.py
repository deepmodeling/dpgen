from typing import List
from dargs import Argument, Variant

from dpgen.arginfo import general_mdata_arginfo
from dpgen.generator.arginfo import (
    basic_args,
    data_args,
    training_args,
    fp_style_vasp_args,
    fp_style_gaussian_args,
)


def general_simplify_arginfo() -> Argument:
    """General simplify arginfo.

    Returns
    -------
    Argument
        arginfo
    """
    doc_labeled = "If true, the initial data is labeled."
    doc_pick_data = "(List of) Path to the directory with the pick data with the deepmd/npy or the HDF5 file with deepmd/hdf5 format. Systems are detected recursively."
    doc_init_pick_number = "The number of initial pick data."
    doc_iter_pick_number = "The number of pick data in each iteration."
    doc_model_devi_f_trust_lo = "The lower bound of forces for the selection for the model deviation."
    doc_model_devi_f_trust_hi = "The higher bound of forces for the selection for the model deviation."

    return [
        Argument("labeled", bool, optional=True, default=False, doc=doc_labeled),
        Argument("pick_data", [str, list], doc=doc_pick_data),
        Argument("init_pick_number", int, doc=doc_init_pick_number),
        Argument("iter_pick_number", int, doc=doc_iter_pick_number),
        Argument("model_devi_f_trust_lo", float, optional=False, doc=doc_model_devi_f_trust_lo),
        Argument("model_devi_f_trust_hi", float, optional=False, doc=doc_model_devi_f_trust_hi),
    ]


def fp_style_variant_type_args() -> Variant:
    """Generate variant for fp style variant type.
    
    Returns
    -------
    Variant
        variant for fp style
    """
    doc_fp_style = 'Software for First Principles, if `labeled` is false. Options include “vasp”, “gaussian” up to now.'
    doc_fp_style_none = 'No fp.'
    doc_fp_style_vasp = 'VASP.'
    doc_fp_style_gaussian = 'Gaussian. The command should be set as `g16 < input`.'

    return Variant("fp_style", [
        Argument("none", dict, doc=doc_fp_style_none),
        # simplify use the same fp method as run
        Argument("vasp", dict, fp_style_vasp_args(), doc=doc_fp_style_vasp),
        Argument("gaussian", dict, fp_style_gaussian_args(),
                 doc=doc_fp_style_gaussian),
    ],
        optional=True,
        default_tag="none",
        doc=doc_fp_style)


def fp_args() -> List[Argument]:
    """Generate arginfo for fp.

    Returns
    -------
    List[Argument]
        arginfo
    """
    doc_fp_task_max = 'Maximum of structures to be calculated in 02.fp of each iteration.'
    doc_fp_task_min = 'Minimum of structures to be calculated in 02.fp of each iteration.'
    doc_fp_accurate_threshold = 'If the accurate ratio is larger than this number, no fp calculation will be performed, i.e. fp_task_max = 0.'
    doc_fp_accurate_soft_threshold = 'If the accurate ratio is between this number and fp_accurate_threshold, the fp_task_max linearly decays to zero.'

    return [
        Argument("fp_task_max", int, optional=True, doc=doc_fp_task_max),
        Argument("fp_task_min", int, optional=True, doc=doc_fp_task_min),
        Argument("fp_accurate_threshold", float,
                 optional=True, doc=doc_fp_accurate_threshold),
        Argument("fp_accurate_soft_threshold", float,
                 optional=True, doc=doc_fp_accurate_soft_threshold),
    ]


def simplify_jdata_arginfo() -> Argument:
    """Generate arginfo for dpgen simplify jdata.

    Returns
    -------
    Argument
        arginfo
    """
    doc_run_jdata = "Parameters for simplify.json, the first argument of `dpgen simplify`."
    return Argument("simplify_jdata",
                    dict,
                    sub_fields=[
                        *basic_args(),
                        # TODO: we may remove sys_configs; it is required in train method
                        *data_args(),
                        *general_simplify_arginfo(),
                        # simplify use the same training method as run
                        *training_args(),
                        *fp_args(),
                    ],
                    sub_variants=[
                        fp_style_variant_type_args(),
                    ],
                    doc=doc_run_jdata,
                    )


def simplify_mdata_arginfo() -> Argument:
    """Generate arginfo for dpgen simplify mdata.

    Returns
    -------
    Argument
        arginfo
    """
    return general_mdata_arginfo("simplify_mdata", ("train", "model_devi", "fp"))
