# CH4 workflow on the Lebesgue platform

This directory is a self-contained methane example for preparing initial data
and running DP-GEN with DeePMD-kit, LAMMPS, VASP, and the Lebesgue platform.

Before running the example:

1. provide `POTCAR_H` and `POTCAR_C` in this directory;
1. replace the Lebesgue account placeholders in `lebesgue_v2_machine.json`;
1. verify the program, image, and machine-type settings available to your
   Lebesgue account.

Generate the initial data with:

```bash
dpgen init_bulk init.json lebesgue_v2_machine.json
```

Then run the iterative workflow from the generated data directory with:

```bash
dpgen run param_CH4_deepmd-kit-2.0.1.json lebesgue_v2_machine.json
```

The machine configuration uses placeholder credentials and resource settings.
Consult the current DPDispatcher Bohrium/Lebesgue documentation before
submitting production jobs.
