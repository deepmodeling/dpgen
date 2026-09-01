---
name: dpgen-run
description: Prepare, explain, validate, and run DP-GEN concurrent-learning workflows. Use for dpgen run param.json or machine.json, training and exploration settings, labeling backends, trust levels, model_devi_jobs, execution, or output inspection.
license: LGPL-3.0-or-later
metadata:
  author: MatMaster
  version: 0.2.1
  repository: https://github.com/deepmodeling/dpgen
---

# DP-GEN Run

Use this skill for the iterative DP-GEN loop: train an ensemble, explore with MD, select uncertain structures, label them, and add the labels to training data.

## Core contract

- Keep workflow settings in `param.json` and execution settings in `machine.json`.
- The exact launcher is `dpgen run param.json machine.json`.
- Treat the launcher shell and dispatched stage environments as separate layers.
- Prepare and validate by default. Execute only when the user explicitly requests it and confirms the exact command after validation.

## Workflow

1. Confirm this is a `dpgen run` task and inspect existing files first.
1. Collect only missing inputs; patch working configurations instead of rebuilding them.
1. Load only the references needed for the current task.
1. Preserve scientific and site-specific choices unless the user approves changes.
1. Validate files, data, commands, and environments before proposing execution.
1. Report exact paths, the command, unresolved fields, and post-run inspection targets.

## Progressive references

Read only the relevant files:

- [Workflow policy](references/workflow.md): scope, requirements, missing-input policy, layout, and reporting.
- [param.json guide](references/param-json.md): systems, training, exploration, FP fields, descriptors, thresholds, and examples.
- [machine.json guide](references/machine-json.md): stage commands, contexts, resources, scheduler environments, and examples.
- [Validation and execution](references/validation-and-run.md): schema checks, data checks, authorization, launch, and output inspection.

## Non-negotiable guardrails

- Never merge the two JSON files or invent cluster paths, queues, modules, or activation commands.
- Preserve descriptor, fitting network, FP backend, thresholds, schedules, ensemble settings, and `type_map` ordering.
- Spell `se_atten_v2` exactly; do not use `fp_style: "none"` for current `dpgen run`.
- Stop and ask when required scientific or execution inputs cannot be discovered safely.
