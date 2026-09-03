# DP-GEN Run Workflow Guidance

Load this reference for scope, preparation policy, project layout, or reporting.

## Loop

`dpgen run` repeats five activities:

1. train an ensemble of Deep Potential models;
1. explore configurations, normally with LAMMPS MD;
1. select structures in the model-deviation candidate window;
1. label selected structures with the configured first-principles backend;
1. add labels to the next training iteration.

Use `param.json` for scientific/workflow settings, `machine.json` for commands,
contexts, and resources, and `dpgen run param.json machine.json` as the launcher.
Preparation needs Python and DP-GEN. Execution also needs DeePMD-kit, the
exploration engine, the selected FP software, and any scheduler runtime.

## Preparation policy

Inspect existing configuration files, training inputs, dataset metadata, and
machine templates first. Ask only for values that cannot be discovered. Patch
working files instead of rebuilding them.

Preserve descriptor and fitting settings, training and FP backends, thresholds,
`type_map` ordering, `model_devi_jobs`, ensemble size, temperatures, pressures,
and MD ensembles. Explain concerns and request direction before changing them.

Reuse known activation commands, executables, queues, partitions, accounts,
paths, and scheduler flags exactly. Never guess site policy or conda/module setup.

## Generic layout

```text
project/
|-- param.json
|-- machine.json
|-- init_data/
|   `-- raw_xxx/          # DeepMD NumPy system(s)
|-- assets/               # structures and runtime support files
`-- iter.*/               # created by DP-GEN
```

Keep repeated experiments in separate, clearly named directories derived from a
reviewed base configuration. Backend-specific files belong under `assets/` or
the path explicitly required by the selected backend.

## Reporting contract

Before execution, report absolute JSON paths, the exact command, validation
results, unresolved inputs, and cost-bearing stages. Execute only after the user
confirms that exact validated command.

After execution, report the current iteration and stage status, failed/pending
tasks, main logs and outputs, candidate and labeled counts, accuracy trends, and
the next files to inspect. Use [monitoring guidance](monitoring.md) when an
iteration stalls or accuracy fails to improve.

## Guardrails

- Never run before both JSON files exist and pass validation.
- Keep `type_map` consistent through data, training, and labeling.
- Do not overwrite working templates blindly or assume outer activation reaches jobs.
- Stop rather than guess a missing scientific or site-specific value.

Official overview: https://docs.deepmodeling.com/projects/dpgen/en/latest/run/index.html
