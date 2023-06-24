from typing import Callable, Tuple

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
    doc_deepmd_version = "DeePMD-kit version, e.g. 2.1.3"
    doc_run_mdata = "machine.json file"
    arg_api_version = Argument(
        "api_version", str, default="1.0", optional=True, doc=doc_api_version
    )
    arg_deepmd_version = Argument(
        "deepmd_version", str, optional=True, default="2", doc=doc_deepmd_version
    )

    sub_fields = [arg_api_version, arg_deepmd_version]
    doc_mdata = "Parameters of command, machine, and resources for %s"
    for task in tasks:
        sub_fields.append(
            Argument(
                task,
                dict,
                optional=False,
                sub_fields=mdata_arginfo(),
                doc=doc_mdata % task,
            )
        )
    return Argument(name, dict, sub_fields=sub_fields, doc=doc_run_mdata)


def check_nd_list(dimesion: int = 2) -> Callable:
    """Return a method to check if the input is a nd list.

    Parameters
    ----------
    dimesion : int, default=2
        dimension of the array

    Returns
    -------
    callable
        check function
    """
    
    def check(value, dimension=dimesion):
        if value is None:
            # do not check null
            return True
        if dimension:
            if not isinstance(value, list):
                return False
        if dimension > 1:
            if not all(check(v, dimension=dimesion-1) for v in value):
                return False
        return True
    return check

errmsg_nd_list = "Must be a %d-dimension list."
