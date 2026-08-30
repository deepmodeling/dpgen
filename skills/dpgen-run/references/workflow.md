# DP-GEN Run Workflow Policy

Load this reference for end-to-end preparation, task scoping, project layout, or result reporting.

## Scope and requirements

`dpgen run` implements the concurrent-learning loop:

1. train an ensemble of Deep Potential models
1. explore configuration space, normally with LAMMPS MD
1. select configurations in the model-deviation candidate window
1. label selected structures with the configured first-principles backend
1. add labeled data to the next training iteration

Preparation and validation require Python and DP-GEN. Execution additionally requires DeePMD-kit, a compatible exploration engine, the software selected by `fp_style`, and any scheduler runtime.

The workflow always uses:

- `param.json` for scientific and workflow parameters
- `machine.json` for commands, contexts, and resources
- `dpgen run param.json machine.json` as the launcher

## Working policy

### Inspect before asking

Look for existing configuration files, training inputs, dataset metadata, and machine templates. Ask only for values that cannot be discovered. Patch existing working files rather than rebuilding them without need.

### Preserve scientific choices

Do not silently change:

- descriptor family or settings
- fitting network
- training backend
- first-principles backend
- trust thresholds
- `type_map` ordering
- `model_devi_jobs` schedule
- ensemble size, temperatures, pressures, or MD ensembles

Explain a concern and request direction before changing a scientific choice.

### Preserve site-specific choices

Reuse known activation commands, executable names, queues, partitions, accounts, paths, and scheduler flags exactly. Never guess conda environments, modules, personal paths, or cluster policy.

### Keep execution separately authorized

Preparing files, validating them, or displaying a command does not authorize an HPC, MD, training, or first-principles workload. Run only when the user explicitly requests execution and confirms the exact validated command.

## Recommended layout

```text
project/
|-- param.json
|-- machine.json
|-- init_data/
|   `-- system_000/
|       |-- type_map.raw
|       |-- type.raw
|       `-- set.000/
|-- assets/
|   `-- structures/
|-- cp2k_basis_pp_file/   # only when CP2K needs it
`-- iter.*/              # created by DP-GEN
```

Keep repeated experiments in separate, clearly named directories derived from one reviewed base configuration.

## Result contract

Before execution, report:

1. absolute paths to `param.json` and `machine.json`
1. exact launcher command
1. validation results
1. unresolved required inputs
1. expected cost-bearing stages

After execution, report:

1. run and current iteration status
1. failed or pending stages
1. main output and log locations
1. candidate and labeled counts when available
1. the next files that need inspection

## General guardrails

- Never run before both JSON files exist and pass validation.
- Keep `type_map` consistent from data through training and labeling.
- Do not overwrite working user templates blindly.
- Do not assume outer-shell activation reaches dispatched jobs.
- Stop rather than guess a missing scientific or site-specific value.

Official overview: https://docs.deepmodeling.com/projects/dpgen/en/latest/run/index.html
