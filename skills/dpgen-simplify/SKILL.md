---
name: dpgen-simplify
description: Prepare, explain, validate, and run DP-GEN simplify workflows for reducing repeated or redundant DeepMD datasets. Use when the user wants to generate or modify `param.json` and `machine.json`, run `dpgen simplify param.json machine.json`, organize repeated simplify experiments, or inspect simplify outputs.
compatibility: Requires a runnable environment with Python and an activated DP-GEN runtime where `dpgen` is available in PATH for the outer simplify command. Real execution also requires DeePMD-kit and any backend-specific software required by the selected `fp_style`. For scheduler execution, each stage environment must be explicitly activated in `resources.source_list`.
license: LGPL-3.0-or-later
metadata:
  author: hyb1109
  version: 0.2.0
  repository: https://github.com/deepmodeling/dpgen
---

# DP-GEN Simplify

Use this skill when the user wants to prepare, explain, validate, or execute the `dpgen simplify` workflow.

This skill is for dataset simplification workflows where the user already has candidate data in DeepMD-compatible format and wants to reduce repeated or redundant structures through iterative selection.

## Core Rule (Critical)

DP-GEN simplify always uses **two parameter classes** and therefore **two JSON files**:

- **Workflow parameters** -> `param.json`
- **Execution / machine parameters** -> `machine.json`

Run exactly:

```bash
dpgen simplify param.json machine.json
```

Environment boundary rule:

- Outer layer: run `dpgen simplify param.json machine.json` in an activated environment where `dpgen --version` works.
- Inner layer: for scheduler stages, explicitly activate runtime in `resources.source_list` on the server side.

## Agent responsibilities

When using this skill, the agent should:

1. confirm that the task is a simplify workflow
1. check whether existing configs or templates are already available
1. collect only the missing dataset, training, FP, and machine inputs
1. generate or patch `param.json`
1. generate or patch `machine.json`
1. explain important simplify parameters in plain language when asked
1. validate the workflow before execution
1. provide the exact command for running simplify
1. if requested, help structure repeated experiments
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

- descriptor family
- fitting net structure
- fp backend
- trust thresholds
- `type_map` ordering

If a value looks scientifically questionable, explain the concern instead of silently replacing it.

### 3. Keep local and scheduler execution explicit

If the user wants local execution, produce local-friendly commands.

If the user wants scheduler execution, produce scheduler-friendly commands and keep queue, partition, and resource requests explicit.

Do not invent scheduler module names or executable paths.

### 4. Do not invent environment activation commands

If the user already has a working activation command such as:

- `conda activate ...`
- `module load ...`
- `source ...`

reuse it exactly.

If execution is requested and the activation method is unknown, ask the user for the precise activation command.

Do not guess conda environment names, module names, or site-specific paths.

### 4.1 Outer launcher policy

Use an activated DP-GEN environment and verify with:

```bash
dpgen --version
```

Do not start simplify from a shell where `dpgen` is unavailable.

### 4.2 Outer vs inner runtime boundaries (critical)

Treat simplify execution as two separate environment layers:

1. Outer layer: the shell that launches `dpgen simplify param.json machine.json` (must have `dpgen` in PATH)
1. Inner layer: stage tasks dispatched by DP-GEN (`train` / `model_devi` / `fp`) on server/runtime side

Even if the outer layer is correct, inner stage tasks still need explicit runtime setup in `machine.json`.
Do not assume the outer shell environment will be inherited by dispatched stage jobs.
For scheduler-style execution, `resources.source_list` must explicitly activate the required runtime environment.

### 5. Prefer reproducible output layout

When generating a simplify workflow, keep files organized and predictable.

Recommended structure:

```text
project/
├── param.json
├── machine.json
├── run.sh
├── logs/
└── summary/
```

For repeated experiments:

```text
project/
├── base/
├── exp_01/
├── exp_02/
├── exp_03/
└── summary/
```

## Minimum required inputs

Collect the following information before generating files.

### Dataset information

- `pick_data`
- `sys_configs`
- `init_data_prefix`
- `init_data_sys`
- `sys_batch_size`
- dataset format
- `type_map`
- `mass_map` if needed
- `labeled`

### Simplify controls

- `init_pick_number`
- `iter_pick_number`
- `model_devi_f_trust_lo`
- `model_devi_f_trust_hi`
- `model_devi_e_trust_lo` / `model_devi_e_trust_hi` if energy trust is used
- `numb_models` if not already specified

### Training setup

- `train_backend` if required by environment (for example `pytorch`)
- `default_training_param`
  - descriptor settings
  - fitting network settings
  - learning rate settings
  - loss settings
  - training step settings

### FP setup

- `fp_style`
- If data is already labeled (energy/force/virial available) and no re-labeling is requested, set `fp_style` to `none`.
- if `fp_style != "none"`, collect matching FP runtime settings such as:
  - `fp_task_max`
  - `fp_task_min`
  - `fp_params`
  - pseudopotential or backend file paths if required

### Execution setup

For each stage `train`, `model_devi`, and `fp`, collect or preserve:

- `command`
- `machine.batch_type`
- `machine.context_type`
- `machine.local_root`
- `machine.remote_root`
- `resources.number_node`
- `resources.cpu_per_node`
- `resources.gpu_per_node`
- `resources.group_size`
- `resources.source_list` (required for scheduler jobs; use it to activate environment explicitly)
- any explicit queue / partition / custom scheduler flags if the user already uses them

Choose a runtime profile first, then fill the matching template:

- server-local Slurm: `assets/machine.template.server-local-slurm.json`
- local machine -> remote Slurm via SSH: `assets/machine.template.ssh-remote-slurm.json`
- pure local shell testing: `assets/machine.template.local-shell.json`

## How to build `param.json`

Construct `param.json` around these logical blocks:

1. element and mass definitions
1. data source and batch settings
1. model ensemble count
1. default DeePMD training parameters
1. FP backend settings
1. simplify pick settings
1. trust thresholds

Key fields usually include:

- `type_map`
- `mass_map`
- `pick_data`
- `init_data_prefix`
- `init_data_sys`
- `sys_batch_size`
- `numb_models`
- `default_training_param`
- `fp_style`
- `shuffle_poscar`
- `fp_task_max`
- `fp_task_min`
- `fp_pp_path`
- `fp_pp_files`
- `fp_params`
- `init_pick_number`
- `iter_pick_number`
- `model_devi_f_trust_lo`
- `model_devi_f_trust_hi`

If the user is doing grid experiments, keep a base template and derive variants from it.

Official reference example (QM7-style, adapted with path placeholders):

- `assets/param.example.qm7.from-official-docs.json`

## How to build `machine.json`

Construct `machine.json` with separate stage blocks for:

- `train`
- `model_devi`
- `fp`

For each stage, keep the following explicit:

- `command`
- machine or context configuration
- resources
- queue or partition if needed
- cpu and gpu counts
- custom scheduler flags
- environment activation commands

Do not merge all stages into one vague machine block.

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

3. verify required dataset paths exist
1. verify stage commands match the selected software stack
1. if `fp_style` is `none`, do not require FP-specific backend settings
1. only then run:

```bash
dpgen simplify param.json machine.json
```

## Output contract

Always provide:

1. final absolute paths to `param.json` and `machine.json`
1. the exact simplify command to run (`dpgen simplify param.json machine.json`)
1. a short pre-run checklist
1. any unresolved required fields
1. if execution was performed, the main output locations and next files to inspect

## Guardrails

- Never merge workflow and machine parameters into one file.
- Never run `dpgen simplify` before both JSON files are present.
- Never hardcode personal cluster, account, queue, or path settings as universal defaults.
- Never silently change the user's scientific choices.
- Keep `type_map` ordering consistent with dataset typing.
- If required inputs are missing, stop and ask instead of guessing.
- If `fp_style` is `none`, skip FP-specific prompts and keep FP-specific settings disabled or unset.
- If data is already labeled and the user does not request new labels, enforce `fp_style = "none"` and do not require active FP runtime fields.
- Do not assume outer-shell activation is inherited by stage jobs; for scheduler execution, require explicit `source_list` per stage.
- If the user already has working templates, patch them rather than overwriting them blindly.

## References and bundled files

Use these bundled files:

- `assets/param.template.json`
- `assets/param.example.qm7.from-official-docs.json`
- `assets/machine.template.json`
- `assets/machine.template.server-local-slurm.json`
- `assets/machine.template.ssh-remote-slurm.json`
- `assets/machine.template.local-shell.json`
- `references/param-fields.md`
- `references/machine-fields.md`
- `references/workflow-notes.md`

External references:

- DP-GEN simplify overview: https://docs.deepmodeling.com/projects/dpgen/en/latest/simplify/simplify.html
- simplify parameter definitions: https://docs.deepmodeling.com/projects/dpgen/en/latest/simplify/simplify-jdata.html
- simplify machine definitions: https://docs.deepmodeling.com/projects/dpgen/en/latest/simplify/simplify-mdata.html
