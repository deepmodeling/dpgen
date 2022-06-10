from dargs import Argument

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
