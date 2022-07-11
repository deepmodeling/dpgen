from typing import Tuple

from dargs import Argument

from dpgen.dispatcher.Dispatcher import mdata_arginfo


def general_mdata_arginfo(name: str, tasks: Tuple[str]) -> Argument:
    """Generate arginfo for general mdata.

    Parameters
    ----------
    name : str
        mdata name
    tasks : tuple[str]
        tuple of task keys, e.g. ("train", "model_devi", "fp")
    
    Returns
    -------
    Argument
        arginfo
    """
    
    doc_api_version = "Please set to 1.0"
    doc_run_mdata = "machine.json file"
    arg_api_version = Argument("api_version", str, optional=False, doc=doc_api_version)

    sub_fields = [arg_api_version]
    doc_mdata = "Parameters of command, machine, and resources for %s"
    for task in tasks:
        sub_fields.append(Argument(
            task, dict, optional=False, sub_fields=mdata_arginfo(),
            doc=doc_mdata % task,
        ))
    return Argument(name, dict, sub_fields=sub_fields, doc=doc_run_mdata)
