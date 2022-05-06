type_map: 
    | type: ``list``
    | argument path: ``type_map``

    Atom types.

mass_map: 
    | type: ``list``
    | argument path: ``mass_map``

    Standard atom weights.

use_ele_temp: 
    | type: ``int``
    | argument path: ``use_ele_temp``

    Currently only support fp_style vasp. 

    - 0: no electron temperature. 

    - 1: eletron temperature as frame parameter. 

    - 2: electron temperature as atom parameter.

init_data_prefix: 
    | type: ``str``, optional
    | argument path: ``init_data_prefix``

    Prefix of initial data directories.

init_data_prefix: 
    | type: ``list``
    | argument path: ``init_data_prefix``

    Directories of initial data. You may use either absolute or relative path here.

sys_format: 
    | type: ``str``
    | argument path: ``sys_format``

    Format of initial data. It will be vasp/poscar if not set.

init_multi_systems: 
    | type: ``bool``, optional
    | argument path: ``init_multi_systems``

    If set to true, init_data_sys directories should contain sub-directories of various systems. DP-GEN will regard all of these sub-directories as inital data systems.

init_batch_size: 
    | type: ``str``, optional
    | argument path: ``init_batch_size``

    Each number is the batch_size of corresponding system for training in init_data_sys. One recommended rule for setting the sys_batch_size and init_batch_size is that batch_size mutiply number of atoms ot the stucture should be larger than 32. If set to auto, batch size will be 32 divided by number of atoms.

sys_configs_prefix: 
    | type: ``str``, optional
    | argument path: ``sys_configs_prefix``

    Prefix of sys_configs.

sys_configs: 
    | type: ``str``
    | argument path: ``sys_configs``

    Containing directories of structures to be explored in iterations.Wildcard characters are supported here.

sys_batch_size: 
    | type: ``list``, optional
    | argument path: ``sys_batch_size``

    Each number is the batch_size for training of corresponding system in sys_configs. If set to auto, batch size will be 32 divided by number of atoms.

numb_models: 
    | type: ``int``
    | argument path: ``numb_models``

    Number of models to be trained in 00.train. 4 is recommend.

training_iter0_model_path: 
    | type: ``list``, optional
    | argument path: ``training_iter0_model_path``

    The model used to init the first iter training. Number of element should be equal to numb_models.

training_init_model: 
    | type: ``bool``, optional
    | argument path: ``training_init_model``

    Iteration > 0, the model parameters will be initilized from the model trained at the previous iteration. Iteration == 0, the model parameters will be initialized from training_iter0_model_path.

default_training_param: 
    | type: ``dict``
    | argument path: ``default_training_param``

    Training parameters for deepmd-kit in 00.train. You can find instructions from here: (https://github.com/deepmodeling/deepmd-kit).

dp_compress: 
    | type: ``bool``, optional, default: ``False``
    | argument path: ``dp_compress``

    Use dp compress to compress the model.

model_devi_dt: 
    | type: ``float``
    | argument path: ``model_devi_dt``

    Timestep for MD. 0.002 is recommend.

model_devi_skip: 
    | type: ``int``
    | argument path: ``model_devi_skip``

    Number of structures skipped for fp in each MD.

model_devi_f_trust_lo: 
    | type: ``list`` | ``float``
    | argument path: ``model_devi_f_trust_lo``

    Lower bound of forces for the selection. If list, should be set for each index in sys_configs, respectively.

model_devi_f_trust_hi: 
    | type: ``list`` | ``float``
    | argument path: ``model_devi_f_trust_hi``

    Upper bound of forces for the selection. If list, should be set for each index in sys_configs, respectively.

model_devi_v_trust_lo: 
    | type: ``list`` | ``float``
    | argument path: ``model_devi_v_trust_lo``

    Lower bound of virial for the selection. If list, should be set for each index in sys_configs, respectively. Should be used with DeePMD-kit v2.x.

model_devi_v_trust_hi: 
    | type: ``list`` | ``float``
    | argument path: ``model_devi_v_trust_hi``

    Upper bound of virial for the selection. If list, should be set for each index in sys_configs, respectively. Should be used with DeePMD-kit v2.x.

model_devi_adapt_trust_lo: 
    | type: ``bool``, optional
    | argument path: ``model_devi_adapt_trust_lo``

    Adaptively determines the lower trust levels of force and virial. This option should be used together with model_devi_numb_candi_f, model_devi_numb_candi_v and optionally with model_devi_perc_candi_f and model_devi_perc_candi_v. dpgen will make two sets:

    - 1. From the frames with force model deviation lower than model_devi_f_trust_hi, select max(model_devi_numb_candi_f, model_devi_perc_candi_f*n_frames) frames with largest force model deviation. 

    - 2. From the frames with virial model deviation lower than model_devi_v_trust_hi, select max(model_devi_numb_candi_v, model_devi_perc_candi_v*n_frames) frames with largest virial model deviation. 

    The union of the two sets is made as candidate dataset.

model_devi_numb_candi_f: 
    | type: ``int``, optional
    | argument path: ``model_devi_numb_candi_f``

    See model_devi_adapt_trust_lo.

model_devi_numb_candi_v: 
    | type: ``int``, optional
    | argument path: ``model_devi_numb_candi_v``

    See model_devi_adapt_trust_lo.

model_devi_perc_candi_f: 
    | type: ``float``, optional
    | argument path: ``model_devi_perc_candi_f``

    See model_devi_adapt_trust_lo.

model_devi_perc_candi_v: 
    | type: ``float``, optional
    | argument path: ``model_devi_perc_candi_v``

    See model_devi_adapt_trust_lo.

model_devi_f_avg_relative: 
    | type: ``bool``, optional
    | argument path: ``model_devi_f_avg_relative``

    Normalized the force model deviations by the RMS force magnitude along the trajectory. This key should not be used with use_relative.

model_devi_clean_traj: 
    | type: ``bool`` | ``int``
    | argument path: ``model_devi_clean_traj``

    If type of model_devi_clean_traj is bool type then it denote whether to clean traj folders in MD since they are too large. If it is Int type, then the most recent n iterations of traj folders will be retained, others will be removed.

model_devi_nopbc: 
    | type: ``bool``
    | argument path: ``model_devi_nopbc``

    Assume open boundary condition in MD simulations.

model_devi_activation_func: 
    | type: ``list``, optional
    | argument path: ``model_devi_activation_func``

    Set activation functions for models, length of the list should be the same as numb_models, and two elements in the list of string respectively assign activation functions to the embedding and fitting nets within each model. Backward compatibility: the orginal "list of String" format is still supported, where embedding and fitting nets of one model use the same activation function, and the length of the list should be the same as numb_models.

model_devi_jobs: 
    | type: ``dict`` | ``list``
    | argument path: ``model_devi_jobs``

    Settings for exploration in 01.model_devi. Each dict in the list corresponds to one iteration. The index of model_devi_jobs exactly accord with index of iterations

    sys_idx: 
        | type: ``list``
        | argument path: ``model_devi_jobs/sys_idx``

        Systems to be selected as the initial structure of MD and be explored. The index corresponds exactly to the sys_configs.

    temps: 
        | type: ``list``
        | argument path: ``model_devi_jobs/temps``

        Temperature (K) in MD.

    press: 
        | type: ``list``
        | argument path: ``model_devi_jobs/press``

        Pressure (Bar) in MD.

    trj_freq: 
        | type: ``int``
        | argument path: ``model_devi_jobs/trj_freq``

        Frequecy of trajectory saved in MD.

    nsteps: 
        | type: ``int``
        | argument path: ``model_devi_jobs/nsteps``

        Running steps of MD.

    ensembles: 
        | type: ``str``
        | argument path: ``model_devi_jobs/ensembles``

        Determining which ensemble used in MD, options include “npt” and “nvt”.

    neidelay: 
        | type: ``int``, optional
        | argument path: ``model_devi_jobs/neidelay``

        delay building until this many steps since last build.

    taut: 
        | type: ``float`` | ``str``, optional, default: ``log``
        | argument path: ``model_devi_jobs/taut``

        Coupling time of thermostat (ps).

    taup: 
        | type: ``float`` | ``str``, optional, default: ``log``
        | argument path: ``model_devi_jobs/taup``

        Coupling time of barostat (ps).

fp_style: 
    | type: ``dict``
    | argument path: ``fp_style``

    Assume open boundary condition in MD simulations.


    Depending on the value of *fp_style*, different sub args are accepted. 

    fp_style:
        | type: ``str`` (flag key)
        | argument path: ``fp_style/fp_style`` 
        | possible choices: vasp, gaussian, siesta, cp2k

        The code used for fp tasks.


    When *fp_style* is set to ``vasp``: 

    fp_pp_path: 
        | type: ``str``
        | argument path: ``fp_style[vasp]/fp_pp_path``

        Directory of psuedo-potential file to be used for 02.fp exists.

    fp_pp_files: 
        | type: ``list``
        | argument path: ``fp_style[vasp]/fp_pp_files``

        Psuedo-potential file to be used for 02.fp. Note that the order of elements should correspond to the order in type_map.

    fp_incar: 
        | type: ``str``
        | argument path: ``fp_style[vasp]/fp_incar``

        Input file for VASP. INCAR must specify KSPACING and KGAMMA.

    fp_aniso_kspacing: 
        | type: ``list``
        | argument path: ``fp_style[vasp]/fp_aniso_kspacing``

        Set anisotropic kspacing. Usually useful for 1-D or 2-D materials. Only support VASP. If it is setting the KSPACING key in INCAR will be ignored.

    cvasp: 
        | type: ``bool``
        | argument path: ``fp_style[vasp]/cvasp``

        If cvasp is true, DP-GEN will use Custodian to help control VASP calculation.


    When *fp_style* is set to ``gaussian``: 

    use_clusters: 
        | type: ``bool``
        | argument path: ``fp_style[gaussian]/use_clusters``

        If set to true, clusters will be taken instead of the whole system. This option does not work with DeePMD-kit 0.x.

    cluster_cutoff: 
        | type: ``float``
        | argument path: ``fp_style[gaussian]/cluster_cutoff``

        The cutoff radius of clusters if use_clusters is set to true.

    fp_params: 
        | type: ``dict``
        | argument path: ``fp_style[gaussian]/fp_params``

        Parameters for Gaussian calculation.

        doc_keywords: 
            | type: ``str``
            | argument path: ``fp_style[gaussian]/fp_params/doc_keywords``

            Keywords for Gaussian input.

        multiplicity: 
            | type: ``int``
            | argument path: ``fp_style[gaussian]/fp_params/multiplicity``

            Spin multiplicity for Gaussian input. If set to auto, the spin multiplicity will be detected automatically. If set to frag, the "fragment=N" method will be used.

        nproc: 
            | type: ``int``
            | argument path: ``fp_style[gaussian]/fp_params/nproc``

            The number of processors for Gaussian input.


    When *fp_style* is set to ``siesta``: 

    use_clusters: 
        | type: ``bool``
        | argument path: ``fp_style[siesta]/use_clusters``

        If set to true, clusters will be taken instead of the whole system. This option does not work with DeePMD-kit 0.x.

    cluster_cutoff: 
        | type: ``float``
        | argument path: ``fp_style[siesta]/cluster_cutoff``

        The cutoff radius of clusters if use_clusters is set to true.

    fp_params: 
        | type: ``dict``
        | argument path: ``fp_style[siesta]/fp_params``

        Parameters for siesta calculation.

        ecut: 
            | type: ``int``
            | argument path: ``fp_style[siesta]/fp_params/ecut``

            Define the plane wave cutoff for grid.

        ediff: 
            | type: ``float``
            | argument path: ``fp_style[siesta]/fp_params/ediff``

            Tolerance of Density Matrix.

        kspacing: 
            | type: ``float``
            | argument path: ``fp_style[siesta]/fp_params/kspacing``

            Sample factor in Brillouin zones.

        mixingweight: 
            | type: ``float``
            | argument path: ``fp_style[siesta]/fp_params/mixingweight``

            Proportion a of output Density Matrix to be used for the input Density Matrix of next SCF cycle (linear mixing).

        NumberPulay: 
            | type: ``int``
            | argument path: ``fp_style[siesta]/fp_params/NumberPulay``

            Controls the Pulay convergence accelerator.


    When *fp_style* is set to ``cp2k``: 

    user_fp_params: 
        | type: ``dict``
        | argument path: ``fp_style[cp2k]/user_fp_params``

        Parameters for cp2k calculation. find detail in manual.cp2k.org. only the kind section must be set before use. we assume that you have basic knowledge for cp2k input.

    external_input_path: 
        | type: ``str``
        | argument path: ``fp_style[cp2k]/external_input_path``

        Conflict with key:user_fp_params, use the template input provided by user, some rules should be followed, read the following text in detail.

fp_task_max: 
    | type: ``int``
    | argument path: ``fp_task_max``

    Maximum of structures to be calculated in 02.fp of each iteration.

fp_task_min: 
    | type: ``int``
    | argument path: ``fp_task_min``

    Minimum of structures to be calculated in 02.fp of each iteration.

fp_accurate_threshold: 
    | type: ``float``, optional
    | argument path: ``fp_accurate_threshold``

    If the accurate ratio is larger than this number, no fp calculation will be performed, i.e. fp_task_max = 0.

fp_accurate_soft_threshold: 
    | type: ``float``, optional
    | argument path: ``fp_accurate_soft_threshold``

    If the accurate ratio is between this number and fp_accurate_threshold, the fp_task_max linearly decays to zero.

fp_cluster_vacuum: 
    | type: ``float``, optional
    | argument path: ``fp_cluster_vacuum``

    If the vacuum size is smaller than this value, this cluster will not be choosen for labeling.

