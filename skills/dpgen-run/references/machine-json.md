# DP-GEN Run machine.json Guide

Load this reference when creating, patching, explaining, or reviewing `machine.json`.

## Runtime boundary

Treat execution as two environments:

1. the outer shell that launches `dpgen run param.json machine.json`
1. the dispatched `train`, `model_devi`, and `fp` tasks

The outer shell must resolve `dpgen`. Scheduler jobs must activate their own software in `resources.source_list`; outer activation is not inherited reliably.

## Stage blocks

Keep separate `train`, `model_devi`, and `fp` blocks. For each stage, collect or preserve:

- `command`
- `machine.batch_type`
- `machine.context_type`
- `machine.local_root`
- `machine.remote_root` when the selected context requires it
- `resources.number_node`
- `resources.cpu_per_node`
- `resources.gpu_per_node`
- `resources.group_size`
- `resources.source_list` for scheduler runtime activation
- queue, partition, account, and custom scheduler flags when applicable

`train_backend` does not rewrite the machine command. `dp` selects the DeePMD installation's default backend; a DeePMD-kit 3.x PyTorch installation may require `dp --pt`. Verify the installed entry point with `dp train -h` or `dp --pt train -h`, then make the `train` command agree with `train_backend`.

Typical later-stage commands are `lmp` for LAMMPS exploration and a backend executable such as `vasp_std`, `cp2k.popt`, `abacus`, `pw.x`, or `g16` for labeling. Preserve known working commands.

## Context and batch compatibility

`remote_root` is context-dependent. Require it for remote or scheduler contexts whose installed dpdispatcher schema needs a remote working root. Preserve a working omission for lazy or local execution and validate against the installed dpdispatcher version.

`batch_type` must be registered by the installed dpdispatcher. Canonical names and lowercase aliases can both be valid. Preserve a working value; validate newly introduced values rather than normalizing them by assumption.

## Starting profiles

Choose a profile that matches where the outer command runs:

- server-local scheduler: [scheduler example](../../../examples/machine/DeePMD-kit-1.x/machine-lsf-slurm-cp2k.json)
- pure local shell: [local example](../../../examples/machine/DeePMD-kit-1.x/machine-local.json)

Patch the closest existing working configuration. Do not transplant site-specific paths or scheduler settings blindly.

## Guardrails

- Never combine all stages into one vague block.
- Never invent executables, modules, activation commands, queues, accounts, or paths.
- Keep CPU, GPU, node, and grouping requests explicit.
- Require explicit `source_list` activation for scheduler stages.
- Preserve working local omissions and installed-version aliases.
- Ensure the training command selects the same backend as `train_backend` and all commands match the scientific stack selected in `param.json`.

External machine reference: https://docs.deepmodeling.com/projects/dpgen/en/latest/run/mdata.html
