from dargs import Argument, ArgumentEncoder, Variant
from typing import Dict, List

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

def init_bulk_vasp_args() -> List[Argument]:
    return []

def init_bulk_abacus_args() -> List[Argument]:
    doc_relax_kpt = 'Path of `KPT` file for relaxation in stage 1. Only useful if `init_fp_style` is "ABACUS".'
    doc_md_kpt = 'Path of `KPT` file for MD simulations in stage 3. Only useful if `init_fp_style` is "ABACUS".'
    doc_atom_masses = 'List of atomic masses of elements. The order should be the same as `Elements`. Only useful if `init_fp_style` is "ABACUS".'
    return [
        Argument("relax_kpt", str, optional=True, doc=doc_relax_kpt),
        Argument("md_kpt", str, optional=True, doc=doc_md_kpt),
        Argument("atom_masses", list, optional=True, doc=doc_atom_masses),
    ]
    

def init_bulk_variant_type_args() -> List[Variant]:
    doc_init_fp_style = "First-principle software. If this key is absent."
    return [Variant("init_fp_style", [
            Argument("VASP", dict, init_bulk_vasp_args(), doc="No more parameters is needed to be added."),
            Argument("ABACUS", dict, init_bulk_abacus_args(), doc="ABACUS"),
        ], default_tag="VASP", optional=True, doc=doc_init_fp_style)]


def init_bulk_jdata_arginfo() -> Argument:
    """Generate arginfo for dpgen init_bulk jdata.
    
    Returns
    -------
    Argument
        dpgen init_bulk jdata arginfo
    """
    doc_init_bulk = "Generate initial data for bulk systems."
    doc_stages = "Stages for `init_bulk`."
    doc_elements = "Atom types."
    doc_potcars = "Path of POTCAR."
    doc_cell_type = "Specifying which typical structure to be generated. **Options** include fcc, hcp, bcc, sc, diamond."
    doc_super_cell = "Size of supercell."
    doc_from_poscar = "Deciding whether to use a given poscar as the beginning of relaxation. If it's true, keys (`cell_type`, `latt`) will be aborted. Otherwise, these two keys are **necessary**."
    doc_from_poscar_path = "Path of POSCAR for VASP or STRU for ABACUS. **Necessary** if `from_poscar` is true."
    doc_relax_incar = "Path of INCAR for VASP or INPUT for ABACUS for relaxation in VASP. **Necessary** if `stages` include 1."
    doc_md_incar = "Path of INCAR for VASP or INPUT for ABACUS for MD in VASP. **Necessary** if `stages` include 3."
    doc_scale = "Scales for isotropic transforming cells."
    doc_skip_relax = "If it's true, you may directly run stage 2 (perturb and scale) using an unrelaxed POSCAR."
    doc_pert_numb = "Number of perturbations for each scaled (key `scale`) POSCAR."
    doc_pert_box = "Anisotropic Perturbation for cells (independent changes of lengths of three box vectors as well as angel among) in decimal formats. 9 elements of the 3x3 perturbation matrix will be randomly sampled from a uniform distribution (default) in the range [-pert_box, pert_box]. Such a perturbation matrix adds the identity matrix gives the actual transformation matrix for this perturbation operation."
    doc_pert_atom = "Perturbation of atom coordinates (Angstrom). Random perturbations are performed on three coordinates of each atom by adding values randomly sampled from a uniform distribution in the range [-pert_atom, pert_atom]."
    doc_md_nstep = "Steps of AIMD in stage 3. If it's not equal to settings via `NSW` in `md_incar`, DP-GEN will follow `NSW`."
    doc_coll_ndata = "Maximal number of collected data."
    doc_type_map = "The indices of elements in deepmd formats will be set in this order."

    return Argument("init_bulk_jdata", dict, [
        Argument("stages", list, optional=False, doc=doc_stages),
        Argument("elements", list, optional=False, doc=doc_elements),
        Argument("potcars", list, optional=True, doc=doc_potcars),
        Argument("cell_type", str, optional=True, doc=doc_cell_type),
        Argument("super_cell", list, optional=False, doc=doc_super_cell),
        Argument("from_poscar", bool, optional=True, default=False, doc=doc_from_poscar),
        Argument("from_poscar_path", str, optional=True, doc=doc_from_poscar_path),
        Argument("relax_incar", str, optional=True, doc=doc_relax_incar),
        Argument("md_incar", str, optional=True, doc=doc_md_incar),
        Argument("scale", list, optional=False, doc=doc_scale),
        Argument("skip_relax", bool, optional=False, doc=doc_skip_relax),
        Argument("pert_numb", int, optional=False, doc=doc_pert_numb),
        Argument("pert_box", float, optional=False, doc=doc_pert_box),
        Argument("pert_atom", float, optional=False, doc=doc_pert_atom),
        Argument("md_nstep", int, optional=False, doc=doc_md_nstep),
        Argument("coll_ndata", int, optional=False, doc=doc_coll_ndata),
        Argument("type_map", list, optional=True, doc=doc_type_map),
    ], sub_variants=init_bulk_variant_type_args(),
    doc=doc_init_bulk) 

def init_surf_jdata_arginfo() -> Argument:
    """Generate arginfo for dpgen init_surf jdata.
    
    Returns
    -------
    Argument
        dpgen init_surf jdata arginfo
    """
    doc_init_surf = "Generate initial data for surface systems."
    doc_stages = "Stages for `init_surf`."
    doc_elements = "Atom types."
    doc_potcars = "Path of POTCAR."
    doc_cell_type = "Specifying which typical structure to be generated. **Options** include fcc, hcp, bcc, sc, diamond."
    doc_super_cell = "Size of supercell."
    doc_from_poscar = "Deciding whether to use a given poscar as the beginning of relaxation. If it's true, keys (`cell_type`, `latt`) will be aborted. Otherwise, these two keys are **necessary**."
    doc_from_poscar_path = "Path of POSCAR for VASP or STRU for ABACUS. **Necessary** if `from_poscar` is true."
    doc_latt = "Lattice constant for single cell."
    doc_layer_numb = "Number of atom layers constructing the slab."
    doc_z_min = "Thickness of slab without vacuum (Angstrom). If `layer_numb` is set, `z_min` will be ignored."
    doc_vacuum_max = "Maximal thickness of vacuum (Angstrom)."
    doc_vacuum_min = "Minimal thickness of vacuum (Angstrom). Default value is 2 times atomic radius."
    doc_vacuum_resol = "Interval of thickness of vacuum. If size of `vacuum_resol` is 1, the interval is fixed to its value. If size of `vacuum_resol` is 2, the interval is `vacuum_resol[0]` before `mid_point`, otherwise `vacuum_resol[1]` after `mid_point`."
    doc_vacuum_numb = "The total number of vacuum layers **Necessary** if vacuum_resol is empty."
    doc_mid_point = "The mid point separating head region and tail region. **Necessary** if the size of vacuum_resol is 2 or 0."
    doc_head_ratio  = "Ratio of vacuum layers in the nearby region with denser intervals(head region). **Necessary** if vacuum_resol is empty."
    doc_millers = "Miller indices."
    doc_relax_incar = "Path of INCAR for relaxation in VASP. **Necessary** if `stages` include 1."
    doc_scale = "Scales for isotropic transforming cells."
    doc_skip_relax = "If it's true, you may directly run stage 2 (perturb and scale) using an unrelaxed POSCAR."
    doc_pert_numb = "Number of perturbations for each scaled (key `scale`) POSCAR."
    doc_pert_box = "Anisotropic Perturbation for cells (independent changes of lengths of three box vectors as well as angel among) in decimal formats. 9 elements of the 3x3 perturbation matrix will be randomly sampled from a uniform distribution (default) in the range [-pert_box, pert_box]. Such a perturbation matrix adds the identity matrix gives the actual transformation matrix for this perturbation operation."
    doc_pert_atom = "Perturbation of atom coordinates (Angstrom). Random perturbations are performed on three coordinates of each atom by adding values randomly sampled from a uniform distribution in the range [-pert_atom, pert_atom]."
    doc_coll_ndata = "Maximal number of collected data."
    return Argument("init_surf_jdata", dict, [
        Argument("stages", list, optional=False, doc=doc_stages),
        Argument("elements", list, optional=False, doc=doc_elements),
        Argument("potcars", list, optional=True, doc=doc_potcars),
        Argument("cell_type", str, optional=True, doc=doc_cell_type),
        Argument("super_cell", list, optional=False, doc=doc_super_cell),
        Argument("from_poscar", bool, optional=True, default=False, doc=doc_from_poscar),
        Argument("from_poscar_path", str, optional=True, doc=doc_from_poscar_path),
        Argument("latt", float, optional=False, doc=doc_latt),
        Argument("layer_numb", int, optional=True, doc=doc_layer_numb),
        Argument("z_min", int, optional=True, doc=doc_z_min),
        Argument("vacuum_max", float, optional=False, doc=doc_vacuum_max),
        Argument("vacuum_min", float, optional=True, doc=doc_vacuum_min),
        Argument("vacuum_resol", list, optional=False, doc=doc_vacuum_resol),
        Argument("vacuum_numb", int, optional=True, doc=doc_vacuum_numb),
        Argument("mid_point", float, optional=True, doc=doc_mid_point),
        Argument("head_ratio", float, optional=True, doc=doc_head_ratio),
        Argument("millers", list, optional=False, doc=doc_millers),
        Argument("relax_incar", str, optional=True, doc=doc_relax_incar),
        Argument("scale", list, optional=False, doc=doc_scale),
        Argument("skip_relax", bool, optional=False, doc=doc_skip_relax),
        Argument("pert_numb", int, optional=False, doc=doc_pert_numb),
        Argument("pert_box", float, optional=False, doc=doc_pert_box),
        Argument("pert_atom", float, optional=False, doc=doc_pert_atom),
        Argument("coll_ndata", int, optional=False, doc=doc_coll_ndata),
    ], doc=doc_init_surf)

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
    doc_qmkeywords = "Gaussian keywords for first-principle calculations. e.g. force mn15/6-31g** Geom=PrintInputOrient. Note that \"force\" job is necessary to collect data. Geom=PrintInputOrient should be used when there are more than 50 atoms in a cluster."

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
