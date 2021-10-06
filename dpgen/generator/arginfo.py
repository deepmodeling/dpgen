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
    doc_mass_map = 'Standard atomic weights (default: "auto"). You can also set a list of atomic masses corresponding to type_map. Tips: at present the default value will not be applied automatically, so you need to set "mass_map" to "auto" manually in param.json.'
    doc_use_ele_temp = 'Currently only support fp_style vasp. \n\n\
- 0: no electron temperature. \n\n\
- 1: eletron temperature as frame parameter. \n\n\
- 2: electron temperature as atom parameter.'

    return [
        Argument("type_map", list, optional=False, doc=doc_type_map),
        Argument("mass_map", [str, list], optional=True, default="auto", doc=doc_mass_map),
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
        Argument("init_batch_size", [list, str], optional=True,
                 doc=doc_init_batch_size),
        Argument("sys_configs_prefix", str, optional=True,
                 doc=doc_sys_configs_prefix),
        Argument("sys_configs", list, optional=False, doc=doc_sys_configs),
        Argument("sys_batch_size", list, optional=True,
                 doc=doc_sys_batch_size),
    ]

# Training


def training_args() -> List[Argument]:
    """Traning arguments.
    
    Returns
    -------
    list[dargs.Argument]
        List of training arguments.
    """
    doc_numb_models = 'Number of models to be trained in 00.train. 4 is recommend.'
    doc_training_iter0_model_path = 'The model used to init the first iter training. Number of element should be equal to numb_models.'
    doc_training_init_model = 'Iteration > 0, the model parameters will be initilized from the model trained at the previous iteration. Iteration == 0, the model parameters will be initialized from training_iter0_model_path.'
    doc_default_training_param = 'Training parameters for deepmd-kit in 00.train. You can find instructions from here: (https://github.com/deepmodeling/deepmd-kit).'
    doc_dp_compress = 'Use dp compress to compress the model.'
    doc_training_reuse_iter = "The minimal index of iteration that continues training models from old models of last iteration."
    doc_reusing = " This option is only adopted when continuing training models from old models. This option will override default parameters."
    doc_training_reuse_old_ratio = "The probability proportion of old data during training." + doc_reusing
    doc_training_reuse_numb_steps = "Number of training batch." + doc_reusing
    doc_training_reuse_start_lr = "The learning rate the start of the training." + doc_reusing
    doc_training_reuse_start_pref_e = "The prefactor of energy loss at the start of the training." + doc_reusing
    doc_training_reuse_start_pref_f = "The prefactor of force loss at the start of the training." + doc_reusing
    doc_model_devi_activation_func = "The activation function in the model. The shape of list should be (N_models, 2), where 2 represents the embedding and fitting network. This option will override default parameters."

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
        Argument("training_reuse_iter", [None, int], optional=True, doc=doc_training_reuse_iter),
        Argument("training_reuse_old_ratio", [None, float], optional=True, doc=doc_training_reuse_old_ratio),
        Argument("training_reuse_numb_steps", [None, int], alias=["training_reuse_stop_batch"], optional=True, default=400000, doc=doc_training_reuse_numb_steps),
        Argument("training_reuse_start_lr", [None, float], optional=True, default=1e-4, doc=doc_training_reuse_start_lr),
        Argument("training_reuse_start_pref_e", [None, float, int], optional=True, default=0.1, doc=doc_training_reuse_start_pref_e),
        Argument("training_reuse_start_pref_f", [None, float, int], optional=True, default=100, doc=doc_training_reuse_start_pref_f),
        Argument("model_devi_activation_func", [None, list], optional=True, doc=doc_model_devi_activation_func),
    ]


# Exploration
def model_devi_jobs_args() -> List[Argument]:
    # this may be not correct
    doc_sys_idx = 'Systems to be selected as the initial structure of MD and be explored. The index corresponds exactly to the sys_configs.'
    doc_temps = 'Temperature (K) in MD.'
    doc_press = 'Pressure (Bar) in MD. Required when ensemble is npt.'
    doc_trj_freq = 'Frequecy of trajectory saved in MD.'
    doc_nsteps = 'Running steps of MD.'
    doc_ensemble = 'Determining which ensemble used in MD, options include “npt” and “nvt”.'
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
        Argument("press", list, optional=True, doc=doc_press),
        Argument("trj_freq", int, optional=False, doc=doc_trj_freq),
        Argument("nsteps", int, optional=False, doc=doc_nsteps),
        Argument("ensemble", str, optional=False, doc=doc_ensemble),
        Argument("neidelay", int, optional=True, doc=doc_neidelay),
        Argument("taut", float, optional=True, doc=doc_taut),
        Argument("taup", float, optional=True, doc=doc_taup),
        Argument("model_devi_f_trust_lo", [
                 float, dict], optional=True, doc=doc_model_devi_f_trust_lo),
        Argument("model_devi_f_trust_hi", [
                 float, dict], optional=True, doc=doc_model_devi_f_trust_hi),
        Argument("model_devi_v_trust_lo", [
                 float, dict], optional=True, doc=doc_model_devi_v_trust_lo),
        Argument("model_devi_v_trust_hi", [
                 float, dict], optional=True, doc=doc_model_devi_v_trust_hi),
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
    doc_shuffle_poscar = 'Shuffle atoms of each frame before running simulations. The purpose is to sample the element occupation of alloys.'
    doc_use_relative = 'Calculate relative force model deviation.'
    doc_epsilon = 'The level parameter for computing the relative force model deviation.'
    doc_use_relative_v = 'Calculate relative virial model deviation.'
    doc_epsilon_v = 'The level parameter for computing the relative virial model deviation.'

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
                 float, list, dict], optional=True, default=1e10, doc=doc_model_devi_v_trust_lo),
        Argument("model_devi_v_trust_hi", [
                 float, list, dict], optional=True, default=1e10, doc=doc_model_devi_v_trust_hi),
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
        Argument("model_devi_nopbc", bool, optional=True, default=False,
                 doc=doc_model_devi_nopbc),
        Argument("shuffle_poscar", bool, optional=True, default=False, doc=doc_shuffle_poscar),
        Argument("use_relative", bool, optional=True, default=False, doc=doc_use_relative),
        Argument("epsilon", float, optional=True, doc=doc_epsilon),
        Argument("use_relative_v", bool, optional=True, default=False, doc=doc_use_relative_v),
        Argument("epsilon_v", float, optional=True, doc=doc_epsilon_v),
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
    doc_ratio_failed = 'Check the ratio of unsuccessfully terminated jobs. If too many FP tasks are not converged, RuntimeError will be raised.'
    doc_fp_skip_bad_box = 'Skip the configurations that are obviously unreasonable before 02.fp'

    return [
        Argument("fp_pp_path", str, optional=False, doc=doc_fp_pp_path),
        Argument("fp_pp_files", list, optional=False, doc=doc_fp_pp_files),
        Argument("fp_incar", str, optional=False, doc=doc_fp_incar),
        Argument("fp_aniso_kspacing", list, optional=True,
                 doc=doc_fp_aniso_kspacing),
        Argument("cvasp", bool, optional=True, doc=doc_cvasp),
        Argument("ratio_failed", float, optional=True,
                 doc=doc_ratio_failed),
        Argument("fp_skip_bad_box", str, optional=True,
                 doc=doc_fp_skip_bad_box),
    ]

# abacus
def fp_style_abacus_args() -> List[Argument]:
    doc_fp_pp_path = 'Directory of psuedo-potential or numerical orbital files to be used for 02.fp exists.'
    doc_fp_pp_files = 'Psuedo-potential file to be used for 02.fp. Note that the order of elements should correspond to the order in type_map.'
    doc_fp_orb_files = 'numerical orbital file to be used for 02.fp when using LCAO basis. Note that the order of elements should correspond to the order in type_map.'
    doc_fp_incar = 'Input file for ABACUS. This is optinal but priority over user_fp_params, one can also setting the key and value of INPUT in user_fp_params.'
    doc_fp_kpt_file = 'KPT file for ABACUS.'
    doc_fp_dpks_descriptor = 'DeePKS descriptor file name. The file should be in pseudopotential directory.'
    doc_user_fp_params = 'Set the key and value of INPUT.'
    doc_k_points = 'Monkhorst-Pack k-grids setting for generating KPT file of ABACUS'

    return [
        Argument("fp_pp_path", str, optional=False, doc=doc_fp_pp_path),
        Argument("fp_pp_files", list, optional=False, doc=doc_fp_pp_files),
        Argument("fp_orb_files", list, optional=True, doc=doc_fp_orb_files),
        Argument("fp_incar", str, optional=True, doc=doc_fp_incar),
        Argument("fp_kpt_file", str, optional=True, doc=doc_fp_kpt_file),
        Argument("fp_dpks_descriptor", str, optional=True, doc=doc_fp_dpks_descriptor),
        Argument("user_fp_params", dict, optional=True, doc=doc_user_fp_params),
        Argument("k_points", list, optional=True, doc=doc_k_points),
    ]



# gaussian
def fp_style_gaussian_args() -> List[Argument]:
    """Gaussian fp style arguments.
    
    Returns
    -------
    list[dargs.Argument]
        list of Gaussian fp style arguments
    """
    doc_keywords = 'Keywords for Gaussian input, e.g. force b3lyp/6-31g**. If a list, run multiple steps.'
    doc_multiplicity = ('Spin multiplicity for Gaussian input. If `auto`, multiplicity will be detected automatically, '
                        'with the following rules: when fragment_guesses=True, multiplicity will +1 for each radical, '
                        'and +2 for each oxygen molecule; when fragment_guesses=False, multiplicity will be 1 or 2, '
                        'but +2 for each oxygen molecule.')
    doc_nproc = 'The number of processors for Gaussian input.'
    doc_charge = 'Molecule charge. Only used when charge is not provided by the system.'
    doc_fragment_guesses = 'Initial guess generated from fragment guesses. If True, `multiplicity` should be `auto`.'
    doc_basis_set = 'Custom basis set.'
    doc_keywords_high_multiplicity = ('Keywords for points with multiple raicals. `multiplicity` should be `auto`. '
                                      'If not set, fallback to normal keywords.')

    args = [
        Argument("keywords", [str, list],
                 optional=False, doc=doc_keywords),
        Argument("multiplicity", [int, str],
                 optional=True, default="auto", doc=doc_multiplicity),
        Argument("nproc", int, optional=False, doc=doc_nproc),
        Argument("charge", int, optional=True, default=0, doc=doc_charge),
        Argument("fragment_guesses", bool, optional=True, default=False, doc=doc_fragment_guesses),
        Argument("basis_set", str, optional=True, doc=doc_basis_set),
        Argument("keywords_high_multiplicity", str, optional=True, doc=doc_keywords_high_multiplicity),
    ]

    doc_use_clusters = 'If set to true, clusters will be taken instead of the whole system.'
    doc_cluster_cutoff = ('The soft cutoff radius of clusters if `use_clusters` is set to true. Molecules will be taken '
                          'as whole even if part of atoms is out of the cluster. Use `cluster_cutoff_hard` to only '
                          'take atoms within the hard cutoff radius.')
    doc_cluster_cutoff_hard = ('The hard cutoff radius of clusters if `use_clusters` is set to true. Outside the hard cutoff radius, '
                               'atoms will not be taken even if they are in a molecule where some atoms are within the cutoff radius.')
    doc_cluster_minify = ('If enabled, when an atom within the soft cutoff radius connects a single bond with '
                          'a non-hydrogen atom out of the soft cutoff radius, the outer atom will be replaced by a '
                          'hydrogen atom. When the outer atom is a hydrogen atom, the outer atom will be '
                          'kept. In this case, other atoms out of the soft cutoff radius will be removed.')
    doc_fp_params_gaussian = 'Parameters for Gaussian calculation.'

    return [
        Argument("use_clusters", bool, optional=True, default=False, doc=doc_use_clusters),
        Argument("cluster_cutoff", float,
                 optional=True, doc=doc_cluster_cutoff),
        Argument("cluster_cutoff_hard", float, optional=True, doc=doc_cluster_cutoff_hard),
        Argument("cluster_minify", bool, optional=True, default=False, doc=doc_cluster_minify),
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
    doc_fp_pp_path = 'Directory of psuedo-potential or numerical orbital files to be used for 02.fp exists.'
    doc_fp_pp_files = 'Psuedo-potential file to be used for 02.fp. Note that the order of elements should correspond to the order in type_map.'

    args = [
        Argument("ecut", int, optional=False, doc=doc_ecut),
        Argument("ediff", float, optional=False, doc=doc_ediff),
        Argument("kspacing", float, optional=False, doc=doc_kspacing),
        Argument("mixingWeight", float, optional=False, doc=doc_mixingweight),
        Argument("NumberPulay", int, optional=False, doc=doc_NumberPulay),
    ]

    doc_use_clusters = 'If set to true, clusters will be taken instead of the whole system. This option does not work with DeePMD-kit 0.x.'
    doc_cluster_cutoff = 'The cutoff radius of clusters if use_clusters is set to true.'
    doc_fp_params_siesta = 'Parameters for siesta calculation.'

    return [
        Argument("use_clusters", bool, optional=True, doc=doc_use_clusters),
        Argument("cluster_cutoff", float,
                 optional=True, doc=doc_cluster_cutoff),
        Argument("fp_params", dict, args, [],
                 optional=False, doc=doc_fp_params_siesta),
        Argument("fp_pp_path", str, optional=False, doc=doc_fp_pp_path),
        Argument("fp_pp_files", list, optional=False, doc=doc_fp_pp_files),
    ]

# cp2k


def fp_style_cp2k_args() -> List[Argument]:
    doc_user_fp_params = 'Parameters for cp2k calculation. find detail in manual.cp2k.org. only the kind section must be set before use. we assume that you have basic knowledge for cp2k input.'
    doc_external_input_path = 'Conflict with key:user_fp_params, use the template input provided by user, some rules should be followed, read the following text in detail.'
    doc_ratio_failed = 'Check the ratio of unsuccessfully terminated jobs. If too many FP tasks are not converged, RuntimeError will be raised.'

    return [
        Argument("user_fp_params", dict, optional=True,
                 doc=doc_user_fp_params, alias=["fp_params"]),
        Argument("external_input_path", str, optional=True,
                 doc=doc_external_input_path),
        Argument("ratio_failed", float, optional=True,
                 doc=doc_ratio_failed),
    ]


def fp_style_variant_type_args() -> Variant:
    doc_fp_style = 'Software for First Principles.'

    return Variant("fp_style", [Argument("vasp", dict, fp_style_vasp_args()),
                                Argument("gaussian", dict,
                                         fp_style_gaussian_args()),
                                Argument("siesta", dict,
                                         fp_style_siesta_args()),
                                Argument("cp2k", dict, fp_style_cp2k_args()),
                                Argument("abacus", dict, fp_style_abacus_args())],
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
