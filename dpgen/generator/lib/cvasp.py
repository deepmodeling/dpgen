#!/usr/bin/env python3

import argparse

from custodian import Custodian
from custodian.vasp.handlers import (
    FrozenJobErrorHandler,
    NonConvergingErrorHandler,
    PositiveEnergyErrorHandler,
    StdErrHandler,
    UnconvergedErrorHandler,
    VaspErrorHandler,
    WalltimeHandler,
)
from custodian.vasp.jobs import VaspJob as cvj
from custodian.vasp.validators import VaspFilesValidator, VasprunXMLValidator

handlers = [
    VaspErrorHandler(),
    FrozenJobErrorHandler(),
    StdErrHandler(),
    NonConvergingErrorHandler(),
    WalltimeHandler(),
    PositiveEnergyErrorHandler(),
    UnconvergedErrorHandler(),
]
validators = [VaspFilesValidator(), VasprunXMLValidator()]


def runvasp(
    cmd,
    opt=False,
    max_errors=3,
    backup=False,
    auto_gamma=False,
    auto_npar=False,
    ediffg=-0.05,
):
    """Cmd example:
    cmd=['mpirun', '-np', '32' , '-machinefile', 'hosts','vasp_std'].
    """
    if opt:
        jobs = cvj.full_opt_run(
            cmd,
            auto_npar=auto_npar,
            ediffg=ediffg,
            backup=backup,
            auto_gamma=auto_gamma,
        )
    else:
        jobs = [cvj(cmd, auto_npar=auto_npar, backup=backup, auto_gamma=auto_gamma)]
    c = Custodian(handlers, jobs, validators=validators, max_errors=max_errors)
    c.run()


def __main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "CMD",
        type=str,
        help="""The command for runing vasp, e.g.,
                              'mpirun -np 32 /path/vasp_std' or
                              'srun /path/vasp_std'
                             """,
    )
    parser.add_argument(
        "MAXERR", type=int, help="The maximum error time for runing vasp"
    )
    args = parser.parse_args()
    cmd = args.CMD.split()
    runvasp(cmd=cmd, max_errors=args.MAXERR)


if __name__ == "__main__":
    __main()
    # vasp="/sharedext4/vasp/vasp.5.4.4/bin/vasp_std"
    # runvasp(cmd=['srun', vasp],max_errors=3,backup=True)
    # runvasp(cmd=['mpirun', '-np', ncpu, fp_cmd],max_errors=max_errors)
