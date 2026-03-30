# Simplify Machine File Notes

This file gives concise notes for `machine.json` used by `dpgen simplify`.

## General rule

Keep three stage blocks separate:

- `train`
- `model_devi`
- `fp`

Do not collapse them into one ambiguous runtime block.

## Each stage should make these explicit

- `command`
- machine or context information
- resources
- queue / partition if needed
- cpu count
- gpu count
- environment activation commands
- custom scheduler flags if needed

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

## Practical advice

When building `machine.json`:

1. do not invent executable names
1. do not invent scheduler module names
1. keep environment activation explicit
1. keep queue and resource requests explicit
1. if the user already has a working template, patch it instead of rewriting everything
