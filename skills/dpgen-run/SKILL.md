---
name: dpgen-run
description: Prepare, explain, validate, and run DP-GEN concurrent-learning workflows using param.json and machine.json.
license: LGPL-3.0-or-later
metadata:
  author: MatMaster
  version: 0.3.0
  repository: https://github.com/deepmodeling/dpgen
---

# DP-GEN Run

Use this skill for the iterative train, explore, label, and retrain workflow.

## Contract

- Keep scientific settings in `param.json` and execution settings in `machine.json`.
- Run exactly `dpgen run param.json machine.json`.
- Prepare and validate by default; execute only after explicit confirmation.

## Workflow

1. Inspect existing data and configurations before asking for inputs.
1. Patch working files and preserve scientific/site-specific choices.
1. Load only the references needed for the task.
1. Validate schemas, generated task inputs, stage state, and monitoring signals.
1. Report paths, command, per-system evidence, risks, and next checks.

## References

- [Workflow guidance](references/workflow.md): scope, layout, policy, and reporting.
- [Configuration, validation, and execution](references/validation-and-run.md): both JSON files, generated-input checks, and safe launch/restart checks.
- [Monitoring and troubleshooting](references/monitoring.md): accuracy trends, stage progress, and evidence-based diagnosis.

## Guardrails

- Never merge the JSON files or invent paths, queues, modules, or commands.
- Preserve descriptor, backend, thresholds, schedules, ensemble, and `type_map` choices.
- Check `se_atten_v2` against the installed backend; `fp_style: "none"` is invalid.
- Stop and ask when required scientific or execution inputs are undiscoverable.
