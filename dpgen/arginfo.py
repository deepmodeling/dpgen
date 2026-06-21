from dargs import Argument


def _is_mdata_task_config(value):
    """Return whether a machine task section matches runtime accepted shapes.

    Parameters
    ----------
    value : object
        Candidate value for one machine task section, such as ``train``.

    Returns
    -------
    bool
        Whether the value is a dict or a non-empty list of dicts.
    """
    if isinstance(value, dict):
        return True
    if isinstance(value, list):
        return len(value) > 0 and all(isinstance(item, dict) for item in value)
    return False


from dpgen.dispatcher.Dispatcher import mdata_arginfo


def general_mdata_arginfo(name: str, tasks: tuple[str]) -> Argument:
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
                (dict, list),
                optional=False,
                sub_fields=mdata_arginfo(),
                extra_check=_is_mdata_task_config,
                extra_check_errmsg=(
                    "expected a dict or a list of dicts, matching "
                    "convert_mdata() runtime semantics"
                ),
                doc=doc_mdata % task,
            )
        )
    return Argument(name, dict, sub_fields=sub_fields, doc=doc_run_mdata)
