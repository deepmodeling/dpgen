# Simplify Machine File Notes

This file gives concise notes for `machine.json` used by `dpgen simplify`.

## General rule

Keep three stage blocks separate:

- `train`
- `model_devi`
- `fp`

Do not collapse them into one ambiguous runtime block.

## Runtime profiles

Use one of these profiles based on where `dpgen` is launched:

1. Server-local Slurm (already logged into cluster login node)

   - `context_type = "LocalContext"`
   - `batch_type = "Slurm"`
   - template: `assets/machine.template.server-local-slurm.json`

1. Local workstation -> remote Slurm cluster

   - `context_type = "SSHContext"`
   - `batch_type = "Slurm"`
   - requires `remote_profile`
   - template: `assets/machine.template.ssh-remote-slurm.json`

1. Local single-machine shell testing

   - `context_type = "LazyLocalContext"`
   - `batch_type = "Shell"`
   - template: `assets/machine.template.local-shell.json`

If your current workflow is "on server, submit Slurm jobs", use profile 1.

## Each stage should make these explicit

- `command`
- machine or context information
- resources
- queue / partition if needed
- cpu count
- gpu count
- environment activation commands
- custom scheduler flags if needed

Important boundary: outer `dpgen simplify` environment and inner stage-job environments are different layers.
Outer layer must be an activated DP-GEN environment (`dpgen --version` passes).
Do not assume outer activation is inherited by stage jobs. For scheduler profiles, set `resources.source_list` explicitly.

## `train`

Used for model training.

Typical concerns:

- DeePMD environment
- GPU availability
- training queue
- environment activation

## `model_devi`

Used for model deviation evaluation.

Typical concerns:

- DeePMD runtime
- consistency with training environment
- output and log handling

## `fp`

Used for first-principles calculations.

Typical concerns:

- backend executable
- pseudopotential / basis support files
- scheduler settings
- backend-specific environment

If `fp_style` is `none`, keep this stage disabled/unset and do not require active FP executable settings.

## Practical advice

When building `machine.json`:

1. do not invent executable names
1. do not invent scheduler module names
1. keep environment activation explicit
1. keep queue and resource requests explicit
1. if the user already has a working template, patch it instead of rewriting everything
