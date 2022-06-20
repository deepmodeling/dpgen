from dargs import Argument, ArgumentEncoder

from dpgen.arginfo import general_mdata_arginfo


def init_bulk_mdata_arginfo() -> Argument:
    """Generate arginfo for dpgen init_bulk mdata.
    
    Returns
    -------
    Argument
        arginfo
    """
    return general_mdata_arginfo("init_bulk_mdata", ("fp",))


def init_surf_mdata_arginfo() -> Argument:
    """Generate arginfo for dpgen init_surf mdata.
    
    Returns
    -------
    Argument
        arginfo
    """
    return general_mdata_arginfo("init_surf_mdata", ("fp",))


def init_reaction_mdata_arginfo() -> Argument:
    """Generate arginfo for dpgen init_reaction mdata.
    
    Returns
    -------
    Argument
        arginfo
    """
    return general_mdata_arginfo("init_reaction_mdata", ("reaxff", "build", "fp"))


def init_reaction_jdata_arginfo() -> Argument:
    """Generate arginfo for dpgen init_reaction jdata.
    
    Returns
    -------
    Argument
        dpgen init_reaction jdata arginfo
    """
    doc_init_reaction = "Generate initial data for reactive systems for small gas-phase molecules, from a ReaxFF NVT MD trajectory."
    doc_type_map = "Type map, which should match types in the initial data. e.g. [\"C\", \"H\", \"O\"]"
    doc_reaxff = "Parameters for ReaxFF NVT MD."
    doc_data = "Path to initial LAMMPS data file. The atom_style should be charge."
    doc_ff = "Path to ReaxFF force field file. Available in the lammps/potentials directory."
    doc_control = "Path to ReaxFF control file."
    doc_temp = "Target Temperature for the NVT MD simulation. Unit: K."
    doc_dt = "Real time for every time step. Unit: fs."
    doc_tau_t = "Time to determine how rapidly the temperature. Unit: fs."
    doc_dump_frep = "Frequency of time steps to collect trajectory."
    doc_nstep = "Total steps to run the ReaxFF MD simulation."
    doc_cutoff = "Cutoff radius to take clusters from the trajectory. Note that only a complete molecule or free radical will be taken."
    doc_dataset_size = "Collected dataset size for each bond type."
    doc_qmkeywords = "Gaussian keywords for first-principle calculations. e.g. force mn15/6-31g**. Note that \"force\" job is necessary to collect data."

    return Argument("init_reaction_jdata", dict, [
        Argument("type_map", list, doc=doc_type_map),
        Argument("reaxff", dict, [
            Argument("data", str, doc=doc_data),
            Argument("ff", str, doc=doc_ff),
            Argument("control", str, doc=doc_control),
            Argument("temp", [float, int], doc=doc_temp),
            Argument("dt", [float, int], doc=doc_dt),
            Argument("tau_t", [float, int], doc=doc_tau_t),
            Argument("dump_freq", int, doc=doc_dump_frep),
            Argument("nstep", int, doc=doc_nstep),
        ], doc=doc_reaxff),
        Argument("cutoff", float, doc=doc_cutoff),
        Argument("dataset_size", int, doc=doc_dataset_size),
        Argument("qmkeywords", str, doc=doc_qmkeywords),
    ], doc=doc_init_reaction)
