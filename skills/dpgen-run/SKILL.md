---
name: dpgen-run
description: Prepare, explain, validate, and run DP-GEN concurrent learning workflows for training deep potential models via iterative exploration. Use when the user wants to generate or modify `param.json` and `machine.json` for `dpgen run`, configure training/exploration/labeling iterations, select descriptor types, set trust levels, define model_devi_jobs, or inspect run outputs.
compatibility: Requires a runnable environment with Python and an activated DP-GEN runtime where `dpgen` is available in PATH for the outer run command. Real execution also requires DeePMD-kit, LAMMPS (with DeePMD plugin for model_devi), and any backend-specific software required by the selected `fp_style`. For scheduler execution, each stage environment must be explicitly activated in `resources.source_list`.
license: LGPL-3.0-or-later
metadata:
  author: MatMaster
  version: 0.1.0
  repository: https://github.com/deepmodeling/dpgen
---

# DP-GEN Run (Concurrent Learning)

Use this skill when the user wants to prepare, explain, validate, or execute the `dpgen run` concurrent learning workflow.

This skill is for the main DP-GEN iterative loop: train an ensemble of deep potential models, explore configuration space via LAMMPS MD, select uncertain structures, label them with first-principles calculations, and feed new data back into training.

## Core Rule (Critical)

DP-GEN run always uses **two parameter classes** and therefore **two JSON files**:

- **Workflow parameters** -> `param.json`
- **Execution / machine parameters** -> `machine.json`

Run exactly:

```bash
dpgen run param.json machine.json
```

Environment boundary rule:

- Outer layer: run `dpgen run param.json machine.json` in an activated environment where `dpgen --version` works.
- Inner layer: for scheduler stages, explicitly activate runtime in `resources.source_list` on the server side.

## Critical Pitfalls (Must Embed)

These are verified failure modes discovered through testing. Treat as hard rules:

1. **`se_atten_v2` spelling** — The descriptor type is `se_atten_v2` (double-t in "atten"). Writing `se_attn_v2` will silently fail or error. Always verify the exact string.

1. **`remote_root` is mandatory** — dpdispatcher requires `remote_root` even for local Shell execution. Set it to a writable path like `/tmp/dpgen_train`. Omitting it causes runtime errors.

1. **`batch_type: "Shell"` is case-sensitive** — Must be capitalized `"Shell"`, not `"shell"` or `"SHELL"`. Same for `"Slurm"`, `"PBS"`, etc.

1. **CP2K `user_fp_params` KIND format** — When multiple element kinds are needed, use `"_": ["H", "O"]` with parallel arrays for `POTENTIAL` and `BASIS_SET`. This maps to repeated KIND sections.

1. **`type_map.raw` ordering** — The `type_map.raw` in your `init_data_sys` directories must exactly match the `type_map` array in `param.json`. Mismatches cause silent data corruption.

1. **DeePMD-kit 3.x format** — `default_training_param` uses nested structure: `model.descriptor`, `model.fitting_net`, `learning_rate`, `loss`, `training`. The `training` block includes `training_data.systems` (left empty — dpgen fills it).

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

dpgen includes `train_backend` options for both `"tensorflow"` and `"pytorch"`, so backend support may include PyTorch depending on the installed DeePMD-kit stack. If backend or descriptor compatibility is unclear, verify it against the repo's supported `train_backend` values (`"tensorflow"`, `"pytorch"`) and the user's installed software before generating configs.

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
dpgen --version
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
- `mass_map` — atomic masses matching `type_map` order
- `init_data_prefix` — prefix for init_data_sys paths
- `init_data_sys` — list of paths to initial training data (deepmd/npy format)
- `sys_configs_prefix` — prefix for sys_configs paths
- `sys_configs` — list of lists of structure file paths

### Training setup

- `numb_models` — number of ensemble models (default: 4)
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

- `fp_style` — backend: `"vasp"`, `"cp2k"`, `"abacus"`, `"gaussian"`, `"pwscf"`, or `"none"`
- `fp_task_max` — maximum number of FP tasks per iteration
- `fp_task_min` — minimum number to trigger FP
- Backend-specific settings:
  - VASP: `fp_pp_path`, `fp_pp_files`, `fp_incar` or `fp_params`
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
- `machine.remote_root` **(mandatory — even for local Shell)**
- `resources.number_node`
- `resources.cpu_per_node`
- `resources.gpu_per_node`
- `resources.group_size`
- `resources.source_list` (required for scheduler jobs; use it to activate environment explicitly)
- any explicit queue / partition / custom scheduler flags if the user already uses them

Choose a runtime profile first, then fill the matching template:

- server-local Slurm: `assets/machine.template.server-local-slurm.json`
- pure local shell testing: `assets/machine.template.local-shell.json`

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
- `mass_map`
- `init_data_sys`
- `sys_configs`
- `numb_models`
- `default_training_param`
- `model_devi_dt`
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

- `assets/param.example.water-cp2k.json`
- `assets/param.example.water-vasp.json`

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

1. confirm outer-layer `dpgen` is available:

```bash
dpgen --version
```

2. validate JSON syntax:

```bash
python -m json.tool param.json
python -m json.tool machine.json
```

3. verify `init_data_sys` directories exist and contain proper deepmd/npy format:

```bash
# Each directory must have type_map.raw, type.raw, set.000/
ls init_data/*/type_map.raw
```

4. verify `type_map.raw` content matches `param.json` `type_map` ordering
1. verify `sys_configs` structure files exist
1. verify stage commands match the selected software stack
1. for CP2K: verify basis set and potential files are accessible
1. only then run:

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
- Always include `remote_root` in machine config — it is mandatory even for local execution.
- Always capitalize `batch_type` values: `"Shell"`, `"Slurm"`, `"PBS"`.
- For CP2K KIND sections with multiple elements, use `"_": ["elem1", "elem2"]` array format.
- Do not assume outer-shell activation is inherited by stage jobs; for scheduler execution, require explicit `source_list` per stage.
- If the user already has working templates, patch them rather than overwriting them blindly.
- Do not set `model_devi_engine` unless using a non-LAMMPS engine (it defaults to LAMMPS).
- In `default_training_param`, leave `training.training_data.systems` unset or empty — `dpgen/generator/run.py` fills it from `init_data_sys` automatically.

## References and bundled files

Use these bundled files:

- `assets/param.template.json`
- `assets/param.example.water-cp2k.json`
- `assets/param.example.water-vasp.json`
- `assets/machine.template.json`
- `assets/machine.template.server-local-slurm.json`
- `assets/machine.template.local-shell.json`
- `references/param-fields.md`
- `references/machine-fields.md`
- `references/workflow-notes.md`
- `references/descriptor-types.md`

External references:

- DP-GEN run overview: https://docs.deepmodeling.com/projects/dpgen/en/latest/run/index.html
- run parameter definitions: https://docs.deepmodeling.com/projects/dpgen/en/latest/run/param.html
- run machine definitions: https://docs.deepmodeling.com/projects/dpgen/en/latest/run/mdata.html
- DeePMD-kit model documentation: https://docs.deepmodeling.com/projects/deepmd/en/latest/model/index.html
