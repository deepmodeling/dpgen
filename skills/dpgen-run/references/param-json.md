# DP-GEN Run param.json Guide

Load this reference when creating, patching, explaining, or reviewing `param.json`.

The current authoritative schema is [dpgen/generator/arginfo.py](../../../dpgen/generator/arginfo.py).

## System and data fields

Collect or preserve:

- `type_map`: ordered element symbols
- `mass_map`: optional; defaults to `"auto"`
- `init_data_prefix`: optional prefix for `init_data_sys`
- `init_data_sys`: DeepMD NumPy-system directories
- `sys_configs_prefix`: optional prefix for `sys_configs`
- `sys_configs`: list of systems, each represented by a list of structure paths

`init_data_sys` entries must resolve to directories containing `type_map.raw`, `type.raw`, and `set.000/`. The ordering in every `type_map.raw` must match `type_map`.

`sys_configs` is a list of lists. The outer list separates systems; each inner list contains structure files for that system.

## Training fields

Required core inputs include:

- `numb_models`: ensemble size; four is a common recommendation, not a schema default
- `default_training_param.model`
- `default_training_param.learning_rate`
- `default_training_param.loss`
- `default_training_param.training`

DeePMD-kit 1.x uses `training.systems`. DeePMD-kit 2.x and 3.x use `training.training_data.systems`. Leave the version-appropriate systems value unset or empty because DP-GEN fills it from `init_data_sys`. Verify `deepmd_version` and installed DeePMD-kit before using version-specific features.

### Descriptor selection

Do not change a user-selected descriptor silently.

| Descriptor    | Typical use                                   | Notes                                 |
| ------------- | --------------------------------------------- | ------------------------------------- |
| `se_e2_a`     | simple systems and conservative compatibility | explicit `sel` list                   |
| `se_atten`    | multicomponent systems                        | supports `sel: "auto"`                |
| `se_atten_v2` | modern attention workflow                     | spell exactly; supports `sel: "auto"` |

DP-GEN exposes `train_backend` values `"tensorflow"` and `"pytorch"`, but this alone does not guarantee every descriptor or DeePMD-kit feature works in the installed stack. Verify backend, `deepmd_version`, and DeePMD-kit compatibility. Prefer dpgen2 when a requested DPA-2/DPA-3 or PyTorch-only workflow is not verified with this DP-GEN installation.

## Exploration fields

Collect or preserve:

- `model_devi_dt`: MD timestep in ps
- `model_devi_skip`: initial frames to skip
- `model_devi_f_trust_lo`
- `model_devi_f_trust_hi`
- `model_devi_clean_traj`
- `model_devi_jobs`

Each `model_devi_jobs` item normally defines `sys_idx`, `temps`, `press`, `trj_freq`, `nsteps`, and `ensemble`.

Do not set `model_devi_engine` for ordinary LAMMPS exploration; LAMMPS is the default. Set it only for a verified alternative engine.

### Trust-window interpretation

- below `model_devi_f_trust_lo`: accurate, normally not selected
- between low and high: candidate for labeling
- above `model_devi_f_trust_hi`: failed or unreliable, normally discarded

Values such as 0.05-0.10 eV/Angstrom for the low threshold and 0.15-0.30 eV/Angstrom for the high threshold are starting ranges, not universal defaults. Tighter thresholds increase labeling cost and require system-specific justification.

## First-principles fields

`fp_style` must be a backend accepted by the current schema, for example `"vasp"`, `"cp2k"`, `"abacus"`, `"gaussian"`, or `"pwscf"`. The value `"none"` is not valid for current `dpgen run`.

Always collect:

- `fp_style`
- `fp_task_max`
- `fp_task_min`

Backend-specific inputs include:

- VASP: `fp_pp_path`, `fp_pp_files`, `fp_incar`
- CP2K: `user_fp_params`
- ABACUS: `user_fp_params`, `fp_pp_path`, `fp_pp_files`, `fp_orb_files`
- Gaussian: `fp_params`
- PWSCF: `user_fp_params`

CP2K is supported natively. For multiple KIND sections in `user_fp_params`, use `"_": ["H", "O"]` with parallel `POTENTIAL` and `BASIS_SET` arrays.

## Logical construction order

1. elements and masses
1. initial data and exploration structures
1. model ensemble and DeePMD training input
1. model-deviation controls and job schedule
1. first-principles backend and task limits

Do not assume a field is optional from an old example. Validate against the current schema.

## Checked repository examples

- [CP2K methane parameters](../../../examples/run/dp2.x-lammps-cp2k/param_CH4_deepmd-kit-2.0.1.json)
- [VASP CH4 parameters](../../../examples/run/dp2.x-lammps-vasp/CH4/param_CH4_deepmd-kit-2.x.json)

External references:

- https://docs.deepmodeling.com/projects/dpgen/en/latest/run/param.html
- https://docs.deepmodeling.com/projects/deepmd/en/latest/model/index.html
