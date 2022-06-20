# init-reaction

`dpgen init-reaction` is a workflow to initilize data for reactive systems of small gas-phase molecules. The workflow was introduced in the "Initialization" section of [Energy & Fuels, 2021, 35 (1), 762–769](https://10.1021/acs.energyfuels.0c03211).

To start the workflow, one needs a box containing reactive systems. The following packages are required for each of the step:
- Exploring: [LAMMPS](https://github.com/lammps/lammps)
- Sampling: [MDDatasetBuilder](https://github.com/tongzhugroup/mddatasetbuilder)
- Labeling: [Gaussian](https://gaussian.com/)

The Exploring step uses LAMMPS [pair_style reaxff](https://docs.lammps.org/latest/pair_reaxff.html) to run a short ReaxMD NVT MD simulation. In the Sampling step, molecular clusters are taken and k-means clustering algorithm is applied to remove the redundancy, which is described in [Nature Communications, 11, 5713 (2020)](10.1038/s41467-020-19497-z). The Labeling step calculates energies and forces using the Gaussian package.

For detailed parameters, see [parametes](init-reaction-jdata.rst) and [machine parameters](init-reaction-mdata.rst).

The genereated data can be used to continue DP-GEN concurrent learning workflow. Read [Energy & Fuels, 2021, 35 (1), 762–769](https://10.1021/acs.energyfuels.0c03211) for details.
