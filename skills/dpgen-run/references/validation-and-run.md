# DP-GEN Run Validation and Execution

Load this reference before claiming a configuration is valid or before executing `dpgen run`.

## 1. Verify the outer environment

```bash
dpgen -h
```

Stop if the launcher is unavailable. Actual execution also requires DeePMD-kit, the exploration engine, the selected first-principles software, and scheduler access where applicable.

## 2. Check JSON syntax

```bash
python -m json.tool param.json
python -m json.tool machine.json
```

## 3. Validate against installed schemas

Use the installed DP-GEN and dpdispatcher code, not documentation alone.

```python
import json

from dpgen.generator.arginfo import run_jdata_arginfo
from dpgen.remote.decide_machine import convert_mdata
from dpgen.util import normalize

with open("param.json") as stream:
    normalize(run_jdata_arginfo(), json.load(stream), strict_check=False)

with open("machine.json") as stream:
    convert_mdata(json.load(stream))
```

Schema normalization is necessary but not sufficient; it does not prove paths, executables, scheduler permissions, or scientific choices are correct.

## 4. Validate data and structures

Resolve every `init_data_sys` entry against `init_data_prefix` when set. Verify each NumPy system contains:

- `type_map.raw`
- `type.raw`
- `set.000/` with the expected arrays

Compare every `type_map.raw` line-for-line with `param.json.type_map`. Verify every structure referenced through `sys_configs_prefix` and `sys_configs` exists.

## 5. Validate scientific stages

Confirm:

- DeePMD-kit version, backend, descriptor, and training-input layout agree
- `train_backend` and the machine `train.command` select the same backend; for DeePMD-kit 3.x PyTorch, verify `dp --pt train -h` and use `dp --pt` when required
- each exploration job references valid systems and MD settings
- force trust thresholds are ordered and scientifically intentional
- FP inputs and support files match `fp_style`
- CP2K basis and potential files are accessible when CP2K is selected
- `fp_task_min` and `fp_task_max` reflect intended cost limits

## 6. Validate execution stages

For `train`, `model_devi`, and `fp`, verify:

- commands exist in the dispatched environment
- context and batch types are accepted by installed dpdispatcher
- required local and remote roots exist
- scheduler resources and flags are valid
- `resources.source_list` activates the correct inner runtime

A successful outer `dpgen -h` does not validate dispatched environments.

## 7. Restart and recovery safety

Treat a restart from an existing `record.dpgen` as a potentially cost-bearing
operation, especially after changing `machine.json`.

Current DPDispatcher submission identity includes the resolved machine,
resources, common files, work path, and generated jobs. Job identity includes
the resources and tasks. Consequently, changes such as command, context,
remote root, queue, node count, `group_size`, scheduler flags, environment
activation, or `strategy.ratio_unfinished` can produce a new submission or
job hash. SSH contexts place work beneath a submission-hash-specific remote
directory, so a new identity does not recover the old submission JSON, job IDs,
or completion tags and may submit the stage as fresh work.

### What the machine-file checksum means

DPDispatcher does not hash the raw bytes of `machine.json`. The checksum is
useful provenance, but restart behavior follows the resolved current-stage
submission identity:

| Restart situation | Expected behavior |
| --- | --- |
| `machine.json` is byte-identical, the DP-GEN/DPDispatcher versions are unchanged, and the same stage generates the same tasks and common-file lists | The same submission hash is expected, so DPDispatcher can recover the existing submission JSON and scheduler job IDs. |
| The file checksum changed, but the resolved block for the active stage did not (for example, only `fp` changed while `train` is active) | The active-stage submission hash may remain unchanged and recover normally. The raw checksum alone cannot decide this. |
| A resolved machine/resource field for the active stage changed, or the generated tasks/common-file lists changed | A new submission hash is expected. DP-GEN re-enters the stage named by `record.dpgen`, but DPDispatcher treats its generated task list as a fresh submission rather than importing completion state from the old identity. |

An unchanged machine file is therefore not sufficient if `param.json`, task
generation, common files, or serialization behavior changed. Conversely, a
changed file checksum does not necessarily restart the active stage.

After restart, require positive evidence of recovery: the log should report
`Find old submission; recover submission from json file` with the previous
submission hash, the scheduler job IDs should be unchanged, and no new
`was submitted` message should appear. `record.dpgen` is a workflow-stage
pointer; it is not a scheduler completion record and does not cancel an older
active job when a new submission identity is created.

`ratio_unfinished` controls how many tasks may remain unfinished while
evaluating one submission. It does not select only unfinished tasks from a
different submission identity and does not migrate completion state between
hashes.

Before restarting:

1. Stop or rule out another active `dpgen run` controller for the same work
   directory.
1. Record the current `record.dpgen`, scheduler job IDs, submission records,
   and the old and proposed stage configuration. A file checksum is useful for
   provenance, but inspect the resolved stage values to judge identity changes.
1. Query the scheduler before retrying an ambiguous submission. Do not infer
   that a timeout or disconnected SSH session means nothing was submitted.
1. If the old submission is active, either let it finish under the exact old
   configuration or deliberately cancel/migrate it after preserving job IDs
   and logs. Do not start a second controller against the edited stage.
1. If enough outputs already exist, verify task counts, completion markers,
   parser-readable results, and the configured failure policy before moving to
   post-processing. Prefer completing retrieval with the old configuration.
   Advance or repair `record.dpgen` only when the installed DP-GEN stage
   semantics are verified and the recovery is documented.
1. If a new submission is intentional, independently enumerate its exact task
   set and confirm it contains only the work meant to be rerun. Do not assume
   DPDispatcher will infer this from existing local outputs.

Apply machine changes intended for later iterations only after the current
cost-bearing stage has advanced beyond them.

## 8. Confirm before launch

Show the exact files, validation summary, unresolved risks, and command:

```bash
dpgen run param.json machine.json
```

Execute only after the user explicitly requests the run and confirms this exact validated command. If the user requested only preparation or validation, stop here.

## 9. Inspect outputs

After launch, report the current `iter.*` directory and the status of training, model-deviation, and FP tasks. Summarize failures from stage logs and report selected and labeled structure counts when available.
