---
name: dpgen-run
description: Prepare, explain, validate, and run DP-GEN concurrent learning workflows for training deep potential models via iterative exploration. Use when the user wants to generate or modify `param.json` and `machine.json` for `dpgen run`, configure training/exploration/labeling iterations, select descriptor types, set trust levels, define model_devi_jobs, or inspect run outputs.
license: LGPL-3.0-or-later
metadata:
  author: MatMaster
  version: 0.1.0
  repository: https://github.com/deepmodeling/dpgen
---

# DP-GEN Run (Concurrent Learning)

Use this skill when the user wants to prepare, explain, validate, or execute the `dpgen run` concurrent learning workflow.

This skill is for the main DP-GEN iterative loop: train an ensemble of deep potential models, explore configuration space via LAMMPS MD, select uncertain structures, label them with first-principles calculations, and feed new data back into training.

## Requirements

Preparation and validation require Python and DP-GEN. Actual workflow execution
also requires DeePMD-kit, a compatible exploration engine such as LAMMPS, and
the software selected by `fp_style`. Scheduler stages must activate their
runtime explicitly in `resources.source_list`.

## Core Rule (Critical)

DP-GEN run always uses **two parameter classes** and therefore **two JSON files**:

- **Workflow parameters** -> `param.json`
- **Execution / machine parameters** -> `machine.json`

Run exactly:

```bash
dpgen run param.json machine.json
```

Environment boundary rule:

- Outer layer: run `dpgen run param.json machine.json` in an activated environment where `dpgen -h` works.
- Inner layer: for scheduler stages, explicitly activate runtime in `resources.source_list` on the server side.

## Critical Pitfalls (Must Embed)

These are verified failure modes discovered through testing. Treat as hard rules:

1. **`se_atten_v2` spelling** — The descriptor type is `se_atten_v2` (double-t in "atten"). Writing `se_attn_v2` will silently fail or error. Always verify the exact string.

1. **`remote_root` is context-dependent** — require it for remote or scheduler
   contexts whose installed dpdispatcher schema needs a remote working root.
   Lazy/local contexts can omit it. Preserve an existing working omission and
   validate against the installed dpdispatcher version.

1. **`batch_type` must be registered** — canonical names and lowercase aliases
   can both be valid. Preserve an existing working value and validate new values
   against the installed dpdispatcher version.

1. **CP2K `user_fp_params` KIND format** — When multiple element kinds are needed, use `"_": ["H", "O"]` with parallel arrays for `POTENTIAL` and `BASIS_SET`. This maps to repeated KIND sections.

1. **`type_map.raw` ordering** — The `type_map.raw` in your `init_data_sys` directories must exactly match the `type_map` array in `param.json`. Mismatches cause silent data corruption.

1. **DeePMD-kit input format** — `default_training_param` uses `model`, `learning_rate`, `loss`, and `training` blocks. DeePMD-kit 2.x–3.x place systems under `training.training_data`, while 1.x uses `training.systems`; dpgen fills the version-appropriate systems field. Check the configured `deepmd_version` and installed DeePMD-kit before generating the rest of the training input or using 3.x-specific features.

1. **`model_devi_engine` defaults to LAMMPS** — No need to set `model_devi_engine` explicitly when using LAMMPS. Only set it for alternative engines.

1. **`fp_style: "cp2k"` is native** — dpgen v0.12+ supports CP2K natively via `user_fp_params` dict. No external plugins required.

1. **`sys_configs` format** — List of lists: outer list = systems, inner list = file paths (POSCAR/xyz) for that system. Example: `[["path/to/POSCAR_1"], ["path/to/POSCAR_2"]]`.

1. **`init_data_sys` format** — List of paths to directories in deepmd/npy format. Each directory must contain `type_map.raw`, `type.raw`, and `set.000/` with numpy arrays.

## Agent responsibilities

When using this skill, the agent should:

1. confirm that the task is a `dpgen run` concurrent learning workflow
1. check whether existing configs or templates are already available
1. collect only the missing training, exploration, FP, and machine inputs
1. select appropriate descriptor type based on system complexity
1. generate or patch `param.json`
1. generate or patch `machine.json`
1. explain important parameters in plain language when asked
1. validate the workflow before execution
1. provide the exact command for running
1. after execution, summarize outputs and next inspection targets

## Working policy

### 0. Treat execution as separately authorized

Prepare, explain, and validate by default. Execute `dpgen run` only when the
user has explicitly requested execution and, after validation, confirmed the
exact command. Supplying configuration files, passing validation, or showing the
command does not by itself authorize an HPC or first-principles workload.

### 1. Ask only for missing inputs

Do not ask the user for everything if part of the configuration is already available.

If the user already provides:

- a partial `param.json`
- a partial `machine.json`
- a known training template
- a known cluster template

then patch those files instead of rebuilding everything from scratch.

### 2. Preserve the user's scientific choices

Do not silently change:

- descriptor family or type
- fitting net structure
- fp backend
- trust thresholds
- `type_map` ordering
- `model_devi_jobs` schedule
- ensemble or temperature settings

If a value looks scientifically questionable, explain the concern instead of silently replacing it.

### 3. Descriptor selection guidance

dpgen exposes `train_backend` options for `"tensorflow"` and `"pytorch"`, but that option alone does not guarantee every PyTorch descriptor or every DeePMD-kit 3.x feature works in the user's installed stack. If backend or descriptor compatibility is unclear, verify it against the repo's `train_backend` options, the configured `deepmd_version`, and the installed DeePMD-kit version before generating configs.

When the user has not specified a descriptor type, recommend based on:

| Descriptor    | Use case                                     | Notes                                                                               |
| ------------- | -------------------------------------------- | ----------------------------------------------------------------------------------- |
| `se_e2_a`     | Simple systems, fast training, well-tested   | Classic two-body embedding, requires explicit `sel` list                            |
| `se_atten`    | Multi-component systems, moderate complexity | Attention-based, supports `sel: "auto"`                                             |
| `se_atten_v2` | Modern default, best accuracy/cost balance   | **Spell exactly `se_atten_v2`** (double-t). Supports `sel: "auto"`, `attn_layer: 0` |

For new projects, prefer `se_atten_v2` unless the user has specific needs. If the user requires PyTorch-only models or DPA-2/DPA-3, prefer **dpgen2** unless the current dpgen and DeePMD-kit installation are verified to support the requested model.

### 4. Keep local and scheduler execution explicit

If the user wants local execution, produce local-friendly commands.

If the user wants scheduler execution, produce scheduler-friendly commands and keep queue, partition, and resource requests explicit.

Do not invent scheduler module names or executable paths.

### 5. Do not invent environment activation commands

If the user already has a working activation command such as:

- `conda activate ...`
- `module load ...`
- `source ...`

reuse it exactly.

If execution is requested and the activation method is unknown, ask the user for the precise activation command.

Do not guess conda environment names, module names, or site-specific paths.

### 5.1 Outer launcher policy

Use an activated DP-GEN environment and verify with:

```bash
dpgen -h
```

Do not start run from a shell where `dpgen` is unavailable.

### 5.2 Outer vs inner runtime boundaries (critical)

Treat run execution as two separate environment layers:

1. Outer layer: the shell that launches `dpgen run param.json machine.json` (must have `dpgen` in PATH)
1. Inner layer: stage tasks dispatched by DP-GEN (`train` / `model_devi` / `fp`) on server/runtime side

Even if the outer layer is correct, inner stage tasks still need explicit runtime setup in `machine.json`.
Do not assume the outer shell environment will be inherited by dispatched stage jobs.
For scheduler-style execution, `resources.source_list` must explicitly activate the required runtime environment.

### 6. Prefer reproducible output layout

When generating a run workflow, keep files organized and predictable.

Recommended structure:

```text
project/
├── param.json
├── machine.json
├── init_data/
│   └── system_000/
│       ├── type_map.raw
│       ├── type.raw
│       └── set.000/
├── assets/
│   └── structures/
│       └── POSCAR_*
├── cp2k_basis_pp_file/    (if using CP2K)
│   ├── BASIS_MOLOPT
│   └── GTH_POTENTIALS
└── iter.*/                (created by dpgen)
```

## Minimum required inputs

Collect the following information before generating files.

### System information

- `type_map` — ordered element symbols
- `mass_map` — optional atomic masses matching `type_map` order; defaults to `"auto"`
- `init_data_prefix` — optional prefix for `init_data_sys` paths
- `init_data_sys` — list of paths to initial training data (deepmd/npy format)
- `sys_configs_prefix` — optional prefix for `sys_configs` paths
- `sys_configs` — list of lists of structure file paths

### Training setup

- `numb_models` — number of ensemble models (4 is recommended)
- `default_training_param` — full DeePMD training input:
  - `model.descriptor` — descriptor type and settings
  - `model.fitting_net` — fitting network settings
  - `learning_rate` — learning rate schedule
  - `loss` — loss function weights
  - `training` — training steps and output settings

### Exploration setup (model_devi)

- `model_devi_dt` — MD timestep in ps (e.g. 0.0005 = 0.5 fs)
- `model_devi_skip` — number of initial frames to skip
- `model_devi_f_trust_lo` — lower force deviation threshold (eV/Å)
- `model_devi_f_trust_hi` — upper force deviation threshold (eV/Å)
- `model_devi_clean_traj` — whether to clean trajectory files after selection
- `model_devi_jobs` — list of exploration job specifications:
  - `sys_idx` — which systems to explore
  - `temps` — temperatures (K)
  - `press` — pressures (bar)
  - `trj_freq` — trajectory save frequency
  - `nsteps` — number of MD steps
  - `ensemble` — `"nvt"` or `"npt"`

### FP setup

- `fp_style` — a labeling backend registered in the [run argument schema](../../dpgen/generator/arginfo.py), such as `"vasp"`, `"cp2k"`, `"abacus"`, `"gaussian"`, or `"pwscf"`; `"none"` is not a valid current choice
- `fp_task_max` — maximum number of FP tasks per iteration
- `fp_task_min` — minimum number to trigger FP
- Backend-specific settings:
  - VASP: `fp_pp_path`, `fp_pp_files`, and `fp_incar`
  - CP2K: `user_fp_params` (nested dict representing cp2k input)
  - ABACUS: `user_fp_params`, `fp_pp_path`, `fp_pp_files`, `fp_orb_files`
  - Gaussian: `fp_params` (keywords, nproc, multiplicity)
  - PWSCF: `user_fp_params`

### Execution setup

For each stage `train`, `model_devi`, and `fp`, collect or preserve:

- `command`
- `machine.batch_type`
- `machine.context_type`
- `machine.local_root`
- `machine.remote_root` when required by the selected context; preserve a
  working omission for lazy/local execution
- `resources.number_node`
- `resources.cpu_per_node`
- `resources.gpu_per_node`
- `resources.group_size`
- `resources.source_list` (required for scheduler jobs; use it to activate environment explicitly)
- any explicit queue / partition / custom scheduler flags if the user already uses them

Choose a runtime profile first, then fill the matching template:

- server-local scheduler:
  [existing scheduler example](../../examples/machine/DeePMD-kit-1.x/machine-lsf-slurm-cp2k.json)
- pure local shell testing:
  [existing local example](../../examples/machine/DeePMD-kit-1.x/machine-local.json)

## How to build `param.json`

Construct `param.json` around these logical blocks:

1. element and mass definitions (`type_map`, `mass_map`)
1. data source configuration (`init_data_prefix`, `init_data_sys`, `sys_configs_prefix`, `sys_configs`)
1. model ensemble count (`numb_models`)
1. default DeePMD training parameters (`default_training_param`)
1. exploration settings (`model_devi_dt`, `model_devi_skip`, `model_devi_f_trust_lo/hi`, `model_devi_clean_traj`)
1. exploration job schedule (`model_devi_jobs`)
1. FP backend settings (`fp_style`, `fp_task_max`, `fp_task_min`, backend-specific params)

Key fields always required:

- `type_map`
- `init_data_sys`
- `sys_configs`
- `numb_models`
- `default_training_param`
- `model_devi_dt`
- `model_devi_skip`
- `model_devi_f_trust_lo`
- `model_devi_f_trust_hi`
- `model_devi_jobs`
- `fp_style`
- `fp_task_max`
- `fp_task_min`

Trust level guidance:

- `model_devi_f_trust_lo`: structures below this are "accurate" — not selected
- `model_devi_f_trust_hi`: structures above this are "failed" — discarded
- Between lo and hi: "candidate" — selected for FP labeling
- Typical starting values: lo=0.05–0.10, hi=0.15–0.30 (eV/Å)
- Tighter thresholds = more FP cost, better accuracy
- System-specific tuning is recommended after initial iterations

Official reference examples:

- CP2K:
  [methane example](../../examples/run/dp2.x-lammps-cp2k/param_CH4_deepmd-kit-2.0.1.json)
- VASP:
  [CH4 example](../../examples/run/dp2.x-lammps-vasp/CH4/param_CH4_deepmd-kit-2.x.json)

## How to build `machine.json`

Construct `machine.json` with separate stage blocks for:

- `train`
- `model_devi`
- `fp`

For each stage, keep the following explicit:

- `command` — the executable (`dp` for train, `lmp` for model_devi, FP command for fp)
- machine or context configuration
- resources
- queue or partition if needed
- cpu and gpu counts
- custom scheduler flags
- environment activation commands

Do not merge all stages into one vague machine block.

Stage-specific commands:

- `train`: `"dp"` (DeePMD-kit training)
- `model_devi`: `"lmp"` (LAMMPS with DeePMD plugin)
- `fp`: backend-specific (`"vasp_std"`, `"cp2k.popt"`, `"abacus"`, `"pw.x"`, `"g16"`)

## Validation before run

Before execution, validate the workflow in this order:

- Step 1: confirm outer-layer `dpgen` is available:

  ```bash
  dpgen -h
  ```

- Step 2: validate JSON syntax:

  ```bash
  python -m json.tool param.json
  python -m json.tool machine.json
  ```

- Step 3: verify `init_data_sys` directories exist and contain proper deepmd/npy format:

  Resolve each `init_data_sys` entry against `init_data_prefix` when a prefix
  is set; do not assume the literal directory `init_data/`. Each NumPy system
  must contain `type_map.raw`, `type.raw`, and `set.000/`.

- Step 4: verify `type_map.raw` content matches `param.json` `type_map` ordering.

- Step 5: verify `sys_configs` structure files exist.

- Step 6: verify stage commands match the selected software stack.

- Step 7: for CP2K, verify basis set and potential files are accessible.

- Step 8: only then run:

  ```bash
  dpgen run param.json machine.json
  ```

## Output contract

Always provide:

1. final absolute paths to `param.json` and `machine.json`
1. the exact run command (`dpgen run param.json machine.json`)
1. a short pre-run checklist
1. any unresolved required fields
1. if execution was performed, the main output locations and next files to inspect

## Guardrails

- Never merge workflow and machine parameters into one file.
- Never run `dpgen run` before both JSON files are present.
- Never hardcode personal cluster, account, queue, or path settings as universal defaults.
- Never silently change the user's scientific choices.
- Keep `type_map` ordering consistent with `init_data_sys` type_map.raw files.
- If required inputs are missing, stop and ask instead of guessing.
- Always spell `se_atten_v2` correctly (not `se_attn_v2`).
- Do not use `fp_style: "none"`; stop before `dpgen run` when only preparation or validation was authorized.
- Include `remote_root` when the selected context requires it; do not add it
  solely to rewrite a working lazy/local configuration.
- Use a `batch_type` accepted by the installed dpdispatcher version and
  preserve an existing working canonical name or alias.
- For CP2K KIND sections with multiple elements, use `"_": ["elem1", "elem2"]` array format.
- Do not assume outer-shell activation is inherited by stage jobs; for scheduler execution, require explicit `source_list` per stage.
- If the user already has working templates, patch them rather than overwriting them blindly.
- Do not set `model_devi_engine` unless using a non-LAMMPS engine (it defaults to LAMMPS).
- In `default_training_param`, leave `training.training_data.systems` unset or empty for DeePMD-kit 2.x–3.x, or `training.systems` unset or empty for 1.x — `dpgen/generator/run.py` fills it from `init_data_sys` automatically.

## Repository examples and references

Use these checked-in repository examples as starting points:

- [local machine example](../../examples/machine/DeePMD-kit-1.x/machine-local.json)
- [CP2K parameter example](../../examples/run/dp2.x-lammps-cp2k/param_CH4_deepmd-kit-2.0.1.json)
- [VASP parameter example](../../examples/run/dp2.x-lammps-vasp/CH4/param_CH4_deepmd-kit-2.x.json)
- [scheduler machine example](../../examples/machine/DeePMD-kit-1.x/machine-lsf-slurm-cp2k.json)

External references:

- DP-GEN run overview: https://docs.deepmodeling.com/projects/dpgen/en/latest/run/index.html
- run parameter definitions: https://docs.deepmodeling.com/projects/dpgen/en/latest/run/param.html
- run machine definitions: https://docs.deepmodeling.com/projects/dpgen/en/latest/run/mdata.html
- DeePMD-kit model documentation: https://docs.deepmodeling.com/projects/deepmd/en/latest/model/index.html
