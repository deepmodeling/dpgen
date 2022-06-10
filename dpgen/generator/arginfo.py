from dargs import Argument

from dpgen.arginfo import general_mdata_arginfo

def run_mdata_arginfo() -> Argument:
    """Generate arginfo for dpgen run mdata.
    
    Returns
    -------
    Argument
        arginfo
    """
    return general_mdata_arginfo("run_mdata", ("train", "model_devi", "fp"))
