from dargs import Argument

from dpgen.arginfo import general_mdata_arginfo

def simplify_mdata_arginfo() -> Argument:
    """Generate arginfo for dpgen simplify mdata.
    
    Returns
    -------
    Argument
        arginfo
    """
    return general_mdata_arginfo("simplify_mdata", ("train", "model_devi", "fp"))
