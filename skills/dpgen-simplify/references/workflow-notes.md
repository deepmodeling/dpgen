# Simplify Workflow Notes

## Standard command

```bash
dpgen simplify param.json machine.json
```

Fallback launcher when local `dpgen` is unavailable:

```bash
uvx --from dpgen dpgen simplify param.json machine.json
```

## Recommended workflow

1. confirm this is really a simplify task
1. inspect dataset source
1. define simplify thresholds
1. build or patch `param.json`
1. build or patch `machine.json`
1. validate the workflow
1. execute if requested
1. summarize outputs and next checks

## Validation checklist

Run these before execution:

```bash
dpgen --version
uvx --from dpgen dpgen --version  # fallback if local dpgen is missing
python -m json.tool param.json
python -m json.tool machine.json
```

Also confirm:

- `pick_data` paths exist
- stage commands match the intended software stack
- FP-specific settings are present only when `fp_style != "none"`

## Recommended repeated-experiment structure

For multiple simplify experiments, use one base template and derive variants.

Example:

```text
project/
├── base/
├── exp_lo020_hi040/
├── exp_lo025_hi045/
├── exp_lo030_hi050/
└── summary/
```

## What to summarize after a run

At minimum, report:

- run status
- stage status
- output locations
- simplify thresholds
- picked counts
- next files to inspect
