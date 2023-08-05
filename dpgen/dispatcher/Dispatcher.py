import os
from distutils.version import LooseVersion

# import dargs
from dargs.dargs import Argument
from dpdispatcher import Machine, Resources, Submission, Task


def make_submission(
    mdata_machine,
    mdata_resources,
    commands,
    work_path,
    run_tasks,
    group_size,
    forward_common_files,
    forward_files,
    backward_files,
    outlog,
    errlog,
):
    if mdata_machine.get("local_root", "./") != "./":
        raise RuntimeError("local_root must be './' in dpgen's machine.json.")

    abs_local_root = os.path.abspath("./")

    abs_mdata_machine = mdata_machine.copy()
    abs_mdata_machine["local_root"] = abs_local_root

    machine = Machine.load_from_dict(abs_mdata_machine)
    resources = Resources.load_from_dict(mdata_resources)

    command = "&&".join(commands)

    task_list = []
    for ii in run_tasks:
        task = Task(
            command=command,
            task_work_path=ii,
            forward_files=forward_files,
            backward_files=backward_files,
            outlog=outlog,
            errlog=errlog,
        )
        task_list.append(task)

    submission = Submission(
        work_base=work_path,
        machine=machine,
        resources=resources,
        task_list=task_list,
        forward_common_files=forward_common_files,
        backward_common_files=[],
    )
    return submission


def mdata_arginfo() -> list[Argument]:
    """This method generates arginfo for a single mdata.

    A submission requires the following keys: command, machine,
    and resources.

    Returns
    -------
    list[Argument]
        arginfo
    """
    doc_command = "Command of a program."
    doc_user_forward_files = "Files to be forwarded to the remote machine."
    doc_user_backward_files = "Files to be backwarded from the remote machine."
    command_arginfo = Argument("command", str, optional=False, doc=doc_command)
    machine_arginfo = Machine.arginfo()
    machine_arginfo.name = "machine"
    resources_arginfo = Resources.arginfo()
    resources_arginfo.name = "resources"
    user_forward_files_arginfo = Argument(
        "user_forward_files", list, optional=True, doc=doc_user_forward_files
    )
    user_backward_files_arginfo = Argument(
        "user_backward_files", list, optional=True, doc=doc_user_backward_files
    )

    return [
        command_arginfo,
        machine_arginfo,
        resources_arginfo,
        user_forward_files_arginfo,
        user_backward_files_arginfo,
    ]


def make_submission_compat(
    machine: dict,
    resources: dict,
    commands: list[str],
    work_path: str,
    run_tasks: list[str],
    group_size: int,
    forward_common_files: list[str],
    forward_files: list[str],
    backward_files: list[str],
    outlog: str = "log",
    errlog: str = "err",
    api_version: str = "1.0",
) -> None:
    """Make submission with compatibility of both dispatcher API v0 and v1.

    If `api_version` is less than 1.0, raise RuntimeError. If
    `api_version` is large than 1.0, use `make_submission`.

    Parameters
    ----------
    machine : dict
        machine dict
    resources : dict
        resource dict
    commands : list[str]
        list of commands
    work_path : str
        working directory
    run_tasks : list[str]
        list of paths to running tasks
    group_size : int
        group size
    forward_common_files : list[str]
        forwarded common files shared for all tasks
    forward_files : list[str]
        forwarded files for each task
    backward_files : list[str]
        backwarded files for each task
    outlog : str, default=log
        path to log from stdout
    errlog : str, default=err
        path to log from stderr
    api_version : str, default=1.0
        API version. 1.0 is required
    """
    if LooseVersion(api_version) < LooseVersion("1.0"):
        raise RuntimeError(
            "API version %s has been removed. Please upgrade to 1.0." % api_version
        )

    elif LooseVersion(api_version) >= LooseVersion("1.0"):
        submission = make_submission(
            machine,
            resources,
            commands=commands,
            work_path=work_path,
            run_tasks=run_tasks,
            group_size=group_size,
            forward_common_files=forward_common_files,
            forward_files=forward_files,
            backward_files=backward_files,
            outlog=outlog,
            errlog=errlog,
        )
        submission.run_submission()
