# DP-GEN Configuration, Validation, and Execution

Load this reference when creating, reviewing, validating, launching, or
recovering a `dpgen run` configuration.

## 1. `param.json` essentials

Preserve or collect:

- `type_map` (ordered symbols), optional `mass_map`, `init_data_prefix`, and `init_data_sys`;
- `sys_configs_prefix` and `sys_configs` (a list of systems, each a list of structures);
- `numb_models` and complete `default_training_param.model`, `learning_rate`, `loss`, and `training`;
- `model_devi_dt`, `model_devi_skip`, trust thresholds, cleanup flag, and `model_devi_jobs`;
- `fp_style`, `fp_task_min`, `fp_task_max`, and backend-specific FP inputs.

DeepMD-kit 1.x uses `training.systems`; 2.x/3.x use
`training.training_data.systems`. Leave the version-appropriate value empty so
DP-GEN can fill it from `init_data_sys`. Verify `deepmd_version` and the
installed stack before using version-specific descriptors or features.

`train_backend` accepts `tensorflow` and `pytorch`, but compatibility still
depends on the installed DeePMD-kit and descriptor. Verify descriptor spelling
against that backend; current releases register `se_atten_v2`.
Ordinary LAMMPS exploration uses the default
`model_devi_engine`; set an alternative only when verified. Current `dpgen run`
does not accept `fp_style: "none"`.

Each `init_data_sys` entry may be a DeepMD NumPy system directory, a parent
directory containing multiple systems, or an HDF5 file. Resolve entries with
DP-GEN's `init_data_prefix` and discovery rules, then validate each NumPy leaf:
it needs `type.raw` and one loadable `set.*` directory. `type_map.raw` is
optional; when present, its symbols must be included in `param.json.type_map`.
Subsets and different orders are valid because DeePMD remaps by element name;
without it, `type.raw` must already use the configured model order.

## 2. `machine.json` essentials

Keep separate `train`, `model_devi`, and `fp` blocks. Each block supplies a
`command`, `machine` (`batch_type`, `context_type`, `local_root`, and any
required `remote_root`/profile), and `resources` (`number_node`, CPU/GPU counts,
`group_size`, queue/partition flags, and scheduler `source_list`). See the
[repository machine schema walkthrough](../../../doc/run/example-of-machine.md)
for generic local, scheduler, and SSH shapes.

The launcher shell and dispatched tasks are different environments. Scheduler
jobs must activate their own runtime through `source_list`; outer activation is
not inherited reliably. Keep the training command as `dp`; DP-GEN appends
backend flags such as `--pt` or `--jax` from `train_backend`. Validate
commands, context types, batch aliases, roots, and resources against installed
DPDispatcher rather than normalizing by assumption.

## 3. Validate generated inputs before launch

```bash
dpgen -h
python -m json.tool param.json
python -m json.tool machine.json
```

Normalize `param.json` with `run_jdata_arginfo()` and convert `machine.json` with
`convert_mdata()`. Resolve all data/structure paths, compare type maps, check
executables and scheduler access, and verify FP inputs and cost limits.

Before submission, inspect each stage's generated backend artifact. For training,
check task `input.json` for model count, seeds, checkpoints, data systems, and
final steps; for LAMMPS exploration, check `input.lammps` and `job.json` for
`sys_idx`, temperatures/pressures, and effective MD settings; for first-principles
labeling, inspect backend inputs such as VASP `INCAR`. When reuse or iteration
overrides exist, assert both default and effective values.

Cross-check the intended iteration against `record.dpgen`, stage directories,
DPDispatcher work base/submission metadata, and logs. A stage is not complete
because a directory or job exists.

## 4. Restart and recovery

Treat a restart from `record.dpgen` as cost-bearing. Keep one controller for a
work directory. Before editing state, stop it, back up the state file, record
stage/job IDs and old/new settings, and query the scheduler.

A changed active-stage command, context, remote root, queue, resources,
`group_size`, flags, activation, task list, or common files can create a new
DPDispatcher submission identity. It will not import completion state from the
old identity. Distinguish recovery of the old submission from an intentional
rerun with a new remote root; never resubmit solely because a connection timed
out.

Positive recovery evidence is the old submission hash, unchanged job IDs, and a
recovery message; absence of a new-submission message is required.
`record.dpgen` is a workflow pointer, not a scheduler completion record.
Apply later-iteration machine changes only after the current cost-bearing stage
advances.

## 5. Confirm and inspect

Show the exact files, validation summary, unresolved risks, and
`dpgen run param.json machine.json`. Execute only after explicit confirmation.
After launch, inspect `iter.*`, stage logs, failure/pending states, selected and
labeled counts, and [monitoring outputs](monitoring.md).
External references: https://docs.deepmodeling.com/projects/dpgen/en/latest/run/index.html
