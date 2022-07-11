from typing import Dict, List
from dargs import Argument, Variant

from dpgen.arginfo import general_mdata_arginfo


def run_mdata_arginfo() -> Argument:
    """Generate arginfo for dpgen run mdata.
    
    Returns
    -------
    Argument
        arginfo
    """
    return general_mdata_arginfo("run_mdata", ("train", "model_devi", "fp"))

# basics
def basic_args() -> List[Argument]:
    doc_type_map = 'Atom types.'
    doc_mass_map = 'Standard atom weights.'
    doc_use_ele_temp = 'Currently only support fp_style vasp. \n\n\
- 0: no electron temperature. \n\n\
- 1: eletron temperature as frame parameter. \n\n\
- 2: electron temperature as atom parameter.'

    return [
        Argument("type_map", list, optional=False, doc=doc_type_map),
        Argument("mass_map", list, optional=False, doc=doc_mass_map),
        Argument("use_ele_temp", int, optional=True,
                 default=0, doc=doc_use_ele_temp),
    ]


def data_args() -> List[Argument]:
    doc_init_data_prefix = 'Prefix of initial data directories.'
    doc_init_data_sys = 'Directories of initial data. You may use either absolute or relative path here. Systems will be detected recursively in the directories.'
    doc_sys_format = 'Format of initial data.'
    doc_init_batch_size = 'Each number is the batch_size of corresponding system for training in init_data_sys. One recommended rule for setting the sys_batch_size and init_batch_size is that batch_size mutiply number of atoms ot the stucture should be larger than 32. If set to auto, batch size will be 32 divided by number of atoms.'
    doc_sys_configs_prefix = 'Prefix of sys_configs.'
    doc_sys_configs = 'Containing directories of structures to be explored in iterations.Wildcard characters are supported here.'
    doc_sys_batch_size = 'Each number is the batch_size for training of corresponding system in sys_configs. If set to auto, batch size will be 32 divided by number of atoms.'

    return [
        Argument("init_data_prefix", str, optional=True,
                 doc=doc_init_data_prefix),
        Argument("init_data_sys", list,
                 optional=False, doc=doc_init_data_sys),
        Argument("sys_format", str, optional=True, default='vasp/poscar', doc=doc_sys_format),
        Argument("init_batch_size", str, optional=True,
                 doc=doc_init_batch_size),
        Argument("sys_configs_prefix", str, optional=True,
                 doc=doc_sys_configs_prefix),
        Argument("sys_configs", list, optional=False, doc=doc_sys_configs),
        Argument("sys_batch_size", list, optional=True,
                 doc=doc_sys_batch_size),
    ]

# Training


def training_args() -> List[Argument]:
    doc_numb_models = 'Number of models to be trained in 00.train. 4 is recommend.'
    doc_training_iter0_model_path = 'The model used to init the first iter training. Number of element should be equal to numb_models.'
    doc_training_init_model = 'Iteration > 0, the model parameters will be initilized from the model trained at the previous iteration. Iteration == 0, the model parameters will be initialized from training_iter0_model_path.'
    doc_default_training_param = 'Training parameters for deepmd-kit in 00.train. You can find instructions from here: (https://github.com/deepmodeling/deepmd-kit).'
    doc_dp_compress = 'Use dp compress to compress the model.'

    return [
        Argument("numb_models", int, optional=False, doc=doc_numb_models),
        Argument("training_iter0_model_path", list, optional=True,
                 doc=doc_training_iter0_model_path),
        Argument("training_init_model", bool, optional=True,
                 doc=doc_training_init_model),
        Argument("default_training_param", dict, optional=False,
                 doc=doc_default_training_param),
        Argument("dp_compress", bool, optional=True,
                 default=False, doc=doc_dp_compress),
    ]


# Exploration
def model_devi_jobs_args() -> List[Argument]:
    # this may be not correct
    doc_sys_idx = 'Systems to be selected as the initial structure of MD and be explored. The index corresponds exactly to the sys_configs.'
    doc_temps = 'Temperature (K) in MD.'
    doc_press = 'Pressure (Bar) in MD.'
    doc_trj_freq = 'Frequecy of trajectory saved in MD.'
    doc_nsteps = 'Running steps of MD.'
    doc_ensembles = 'Determining which ensemble used in MD, options include “npt” and “nvt”.'
    doc_neidelay = 'delay building until this many steps since last build.'
    doc_taut = 'Coupling time of thermostat (ps).'
    doc_taup = 'Coupling time of barostat (ps).'
    doc_model_devi_f_trust_lo = 'Lower bound of forces for the selection. If dict, should be set for each index in sys_idx, respectively.'
    doc_model_devi_f_trust_hi = 'Upper bound of forces for the selection. If dict, should be set for each index in sys_idx, respectively.'
    doc_model_devi_v_trust_lo = 'Lower bound of virial for the selection. If dict, should be set for each index in sys_idx, respectively. Should be used with DeePMD-kit v2.x.'
    doc_model_devi_v_trust_hi = 'Upper bound of virial for the selection. If dict, should be set for each index in sys_idx, respectively. Should be used with DeePMD-kit v2.x.'

    args = [
        Argument("sys_idx", list, optional=False, doc=doc_sys_idx),
        Argument("temps", list, optional=False, doc=doc_temps),
        Argument("press", list, optional=False, doc=doc_press),
        Argument("trj_freq", int, optional=False, doc=doc_trj_freq),
        Argument("nsteps", int, optional=False, doc=doc_nsteps),
        Argument("ensembles", str, optional=False, doc=doc_ensembles),
        Argument("neidelay", int, optional=True, doc=doc_neidelay),
        Argument("taut", float, optional=True, doc=doc_taut),
        Argument("taup", float, optional=True, doc=doc_taup),
        Argument("model_devi_f_trust_lo", [
                 float, dict], optional=False, doc=doc_model_devi_f_trust_lo),
        Argument("model_devi_f_trust_hi", [
                 float, dict], optional=False, doc=doc_model_devi_f_trust_hi),
        Argument("model_devi_v_trust_lo", [
                 float, dict], optional=False, doc=doc_model_devi_v_trust_lo),
        Argument("model_devi_v_trust_hi", [
                 float, dict], optional=False, doc=doc_model_devi_v_trust_hi),
    ]

    doc_model_devi_jobs = 'Settings for exploration in 01.model_devi. Each dict in the list corresponds to one iteration. The index of model_devi_jobs exactly accord with index of iterations'
    return Argument("model_devi_jobs", list, args, [], repeat=True, doc=doc_model_devi_jobs)


def model_devi_lmp_args() -> List[Argument]:
    doc_model_devi_dt = 'Timestep for MD. 0.002 is recommend.'
    doc_model_devi_skip = 'Number of structures skipped for fp in each MD.'
    doc_model_devi_f_trust_lo = 'Lower bound of forces for the selection. If list or dict, should be set for each index in sys_configs, respectively.'
    doc_model_devi_f_trust_hi = 'Upper bound of forces for the selection. If list or dict, should be set for each index in sys_configs, respectively.'
    doc_model_devi_v_trust_lo = 'Lower bound of virial for the selection. If list or dict, should be set for each index in sys_configs, respectively. Should be used with DeePMD-kit v2.x.'
    doc_model_devi_v_trust_hi = 'Upper bound of virial for the selection. If list or dict, should be set for each index in sys_configs, respectively. Should be used with DeePMD-kit v2.x.'
    doc_model_devi_adapt_trust_lo = 'Adaptively determines the lower trust levels of force and virial. This option should be used together with model_devi_numb_candi_f, model_devi_numb_candi_v and optionally with model_devi_perc_candi_f and model_devi_perc_candi_v. dpgen will make two sets:\n\n\
- 1. From the frames with force model deviation lower than model_devi_f_trust_hi, select max(model_devi_numb_candi_f, model_devi_perc_candi_f*n_frames) frames with largest force model deviation. \n\n\
- 2. From the frames with virial model deviation lower than model_devi_v_trust_hi, select max(model_devi_numb_candi_v, model_devi_perc_candi_v*n_frames) frames with largest virial model deviation. \n\n\
The union of the two sets is made as candidate dataset.'
    doc_model_devi_numb_candi_f = 'See model_devi_adapt_trust_lo.'
    doc_model_devi_numb_candi_v = 'See model_devi_adapt_trust_lo.'
    doc_model_devi_perc_candi_f = 'See model_devi_adapt_trust_lo.'
    doc_model_devi_perc_candi_v = 'See model_devi_adapt_trust_lo.'
    doc_model_devi_f_avg_relative = 'Normalized the force model deviations by the RMS force magnitude along the trajectory. This key should not be used with use_relative.'
    doc_model_devi_clean_traj = 'If type of model_devi_clean_traj is bool type then it denote whether to clean traj folders in MD since they are too large. If it is Int type, then the most recent n iterations of traj folders will be retained, others will be removed.'
    doc_model_devi_nopbc = 'Assume open boundary condition in MD simulations.'
    doc_model_devi_activation_func = 'Set activation functions for models, length of the list should be the same as numb_models, and two elements in the list of string respectively assign activation functions to the embedding and fitting nets within each model. Backward compatibility: the orginal "list of String" format is still supported, where embedding and fitting nets of one model use the same activation function, and the length of the list should be the same as numb_models.'

    return [
        model_devi_jobs_args(),
        Argument("model_devi_dt", float,
                 optional=False, doc=doc_model_devi_dt),
        Argument("model_devi_skip", int, optional=False,
                 doc=doc_model_devi_skip),
        Argument("model_devi_f_trust_lo", [
                 float, list, dict], optional=False, doc=doc_model_devi_f_trust_lo),
        Argument("model_devi_f_trust_hi", [
                 float, list, dict], optional=False, doc=doc_model_devi_f_trust_hi),
        Argument("model_devi_v_trust_lo", [
                 float, list, dict], optional=False, doc=doc_model_devi_v_trust_lo),
        Argument("model_devi_v_trust_hi", [
                 float, list, dict], optional=False, doc=doc_model_devi_v_trust_hi),
        Argument("model_devi_adapt_trust_lo", bool, optional=True,
                 doc=doc_model_devi_adapt_trust_lo),
        Argument("model_devi_numb_candi_f", int, optional=True,
                 doc=doc_model_devi_numb_candi_f),
        Argument("model_devi_numb_candi_v", int, optional=True,
                 doc=doc_model_devi_numb_candi_v),
        Argument("model_devi_perc_candi_f", float,
                 optional=True, doc=doc_model_devi_perc_candi_f),
        Argument("model_devi_perc_candi_v", float,
                 optional=True, doc=doc_model_devi_perc_candi_v),
        Argument("model_devi_f_avg_relative", bool, optional=True,
                 doc=doc_model_devi_f_avg_relative),
        Argument("model_devi_clean_traj", [
                 bool, int], optional=False, doc=doc_model_devi_clean_traj),
        Argument("model_devi_nopbc", bool, optional=False,
                 doc=doc_model_devi_nopbc),
        Argument("model_devi_activation_func", list, optional=True,
                 doc=doc_model_devi_activation_func),
    ]


def model_devi_args() -> List[Variant]:
    doc_model_devi_engine = "Engine for the model deviation task."
    return [Variant("model_devi_engine", [
            Argument("lammps", dict, model_devi_lmp_args(), doc="LAMMPS"),
        ], default_tag="lammps", optional=True, doc=doc_model_devi_engine)]


# Labeling
# vasp
def fp_style_vasp_args() -> List[Argument]:
    doc_fp_pp_path = 'Directory of psuedo-potential file to be used for 02.fp exists.'
    doc_fp_pp_files = 'Psuedo-potential file to be used for 02.fp. Note that the order of elements should correspond to the order in type_map.'
    doc_fp_incar = 'Input file for VASP. INCAR must specify KSPACING and KGAMMA.'
    doc_fp_aniso_kspacing = 'Set anisotropic kspacing. Usually useful for 1-D or 2-D materials. Only support VASP. If it is setting the KSPACING key in INCAR will be ignored.'
    doc_cvasp = 'If cvasp is true, DP-GEN will use Custodian to help control VASP calculation.'

    return [
        Argument("fp_pp_path", str, optional=False, doc=doc_fp_pp_path),
        Argument("fp_pp_files", list, optional=False, doc=doc_fp_pp_files),
        Argument("fp_incar", str, optional=False, doc=doc_fp_incar),
        Argument("fp_aniso_kspacing", list, optional=False,
                 doc=doc_fp_aniso_kspacing),
        Argument("cvasp", bool, optional=False, doc=doc_cvasp),
    ]


# gaussian
def fp_style_gaussian_args() -> List[Argument]:
    doc_keywords = 'Keywords for Gaussian input.'
    doc_multiplicity = 'Spin multiplicity for Gaussian input. If set to auto, the spin multiplicity will be detected automatically. If set to frag, the "fragment=N" method will be used.'
    doc_nproc = 'The number of processors for Gaussian input.'

    args = [
        Argument("keywords", [str or list],
                 optional=False, doc=doc_keywords),
        Argument("multiplicity", [int or str],
                 optional=False, doc=doc_multiplicity),
        Argument("nproc", int, optional=False, doc=doc_nproc),
    ]

    doc_use_clusters = 'If set to true, clusters will be taken instead of the whole system. This option does not work with DeePMD-kit 0.x.'
    doc_cluster_cutoff = 'The cutoff radius of clusters if use_clusters is set to true.'
    doc_fp_params_gaussian = 'Parameters for Gaussian calculation.'

    return [
        Argument("use_clusters", bool, optional=True, default=False, doc=doc_use_clusters),
        Argument("cluster_cutoff", float,
                 optional=True, doc=doc_cluster_cutoff),
        Argument("fp_params", dict, args, [],
                 optional=False, doc=doc_fp_params_gaussian),
    ]

# siesta
def fp_style_siesta_args() -> List[Argument]:
    doc_ecut = 'Define the plane wave cutoff for grid.'
    doc_ediff = 'Tolerance of Density Matrix.'
    doc_kspacing = 'Sample factor in Brillouin zones.'
    doc_mixingweight = 'Proportion a of output Density Matrix to be used for the input Density Matrix of next SCF cycle (linear mixing).'
    doc_NumberPulay = 'Controls the Pulay convergence accelerator.'

    args = [
        Argument("ecut", int, optional=False, doc=doc_ecut),
        Argument("ediff", float, optional=False, doc=doc_ediff),
        Argument("kspacing", float, optional=False, doc=doc_kspacing),
        Argument("mixingweight", float, optional=False, doc=doc_mixingweight),
        Argument("NumberPulay", int, optional=False, doc=doc_NumberPulay),
    ]

    doc_use_clusters = 'If set to true, clusters will be taken instead of the whole system. This option does not work with DeePMD-kit 0.x.'
    doc_cluster_cutoff = 'The cutoff radius of clusters if use_clusters is set to true.'
    doc_fp_params_siesta = 'Parameters for siesta calculation.'

    return [
        Argument("use_clusters", bool, optional=False, doc=doc_use_clusters),
        Argument("cluster_cutoff", float,
                 optional=False, doc=doc_cluster_cutoff),
        Argument("fp_params", dict, args, [],
                 optional=False, doc=doc_fp_params_siesta),
    ]

# cp2k


def fp_style_cp2k_args() -> List[Argument]:
    doc_user_fp_params = 'Parameters for cp2k calculation. find detail in manual.cp2k.org. only the kind section must be set before use. we assume that you have basic knowledge for cp2k input.'
    doc_external_input_path = 'Conflict with key:user_fp_params, use the template input provided by user, some rules should be followed, read the following text in detail.'

    return [
        Argument("user_fp_params", dict, optional=False,
                 doc=doc_user_fp_params),
        Argument("external_input_path", str, optional=False,
                 doc=doc_external_input_path),
    ]


def fp_style_variant_type_args() -> Variant:
    doc_fp_style = 'Software for First Principles.'

    return Variant("fp_style", [Argument("vasp", dict, fp_style_vasp_args()),
                                Argument("gaussian", dict,
                                         fp_style_gaussian_args()),
                                Argument("siesta", dict,
                                         fp_style_siesta_args()),
                                Argument("cp2k", dict, fp_style_cp2k_args())],
                   optional=False,
                   doc=doc_fp_style)


def fp_args() -> List[Argument]:
    doc_fp_task_max = 'Maximum of structures to be calculated in 02.fp of each iteration.'
    doc_fp_task_min = 'Minimum of structures to be calculated in 02.fp of each iteration.'
    doc_fp_accurate_threshold = 'If the accurate ratio is larger than this number, no fp calculation will be performed, i.e. fp_task_max = 0.'
    doc_fp_accurate_soft_threshold = 'If the accurate ratio is between this number and fp_accurate_threshold, the fp_task_max linearly decays to zero.'
    doc_fp_cluster_vacuum = 'If the vacuum size is smaller than this value, this cluster will not be choosen for labeling.'

    return [
        Argument("fp_task_max", int, optional=False, doc=doc_fp_task_max),
        Argument("fp_task_min", int, optional=False, doc=doc_fp_task_min),
        Argument("fp_accurate_threshold", float,
                 optional=True, doc=doc_fp_accurate_threshold),
        Argument("fp_accurate_soft_threshold", float,
                 optional=True, doc=doc_fp_accurate_soft_threshold),
        Argument("fp_cluster_vacuum", float,
                 optional=True, doc=doc_fp_cluster_vacuum),
    ]


def run_jdata_arginfo() -> Argument:
    """Argument information for dpgen run mdata.
    
    Returns
    -------
    Argument
        argument information
    """
    doc_run_jdata = "param.json file"
    return Argument("run_jdata",
                    dict,
                    sub_fields=basic_args() + data_args() + training_args() + fp_args(),
                    sub_variants=model_devi_args() + [fp_style_variant_type_args()],
                    doc=doc_run_jdata)
