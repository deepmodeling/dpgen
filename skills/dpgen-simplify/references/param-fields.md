# Simplify Parameter Notes

This file gives concise notes for the most important fields in `param.json` for `dpgen simplify`.

For a complete official-style example template, see:

- `assets/param.example.qm7.from-official-docs.json`

## Core dataset fields

### `pick_data`

Path to the candidate dataset to simplify.

Use this field to point to the existing DeepMD-compatible data source.

### `sys_configs`

Configuration discovery pattern for systems.

Often provided as nested lists of path patterns.

### `init_data_prefix`

Prefix used by DP-GEN when resolving initial data systems.

### `init_data_sys`

List of initial system indices for training.

Can be empty when starting fully from `pick_data`.

### `sys_batch_size`

Batch size policy per system, often set to `["auto"]`.

### `type_map`

Ordered list of element symbols.

This ordering should stay consistent with the dataset and the DeePMD training configuration.

### `mass_map`

Atomic masses corresponding to `type_map`.

Keep the order consistent with `type_map`.

## Simplify selection fields

### `init_pick_number`

Number of structures picked at the initial stage.

### `iter_pick_number`

Number of structures picked in each later simplify iteration.

### `model_devi_f_trust_lo`

Lower bound of the force-deviation trust window.

### `model_devi_f_trust_hi`

Upper bound of the force-deviation trust window.

In general, these two values define the force-deviation region used during structure selection.

### `model_devi_e_trust_lo`

Lower bound of the energy-deviation trust window (optional).

### `model_devi_e_trust_hi`

Upper bound of the energy-deviation trust window (optional).

## Training-related fields

### `numb_models`

Number of DeePMD models trained as an ensemble.

### `train_backend`

Training backend (for example `pytorch`) when needed by your environment.

### `labeled`

Whether the input candidate data is already labeled.

### `default_training_param`

Training template used during simplify.

This usually contains:

- model descriptor
- fitting net
- learning rate
- loss
- training steps

## FP-related fields

### `fp_style`

Backend used for first-principles calculations.

Set to `none` if no FP stage is intended.
If data is already labeled and no re-labeling is requested, use `none`.

### `fp_task_max`

Maximum number of FP tasks.

### `fp_task_min`

Minimum number of FP tasks.

### `fp_params`

Backend-specific parameters for the chosen `fp_style`.

### `fp_pp_path`

Path to pseudopotential or backend support files if needed.

### `fp_pp_files`

List of files needed by the FP backend.

## Practical advice

When modifying `param.json`:

1. keep `type_map` consistent
1. do not silently switch descriptor families
1. do not silently change the FP backend
1. make threshold changes explicit and traceable
1. if doing experiments, derive all variants from one base template
