# DP-GEN: A concurrent learning platform for the generation of reliable deep learning based potential energy models

[![GitHub release](https://img.shields.io/github/release/deepmodeling/dpgen.svg?maxAge=86400)](https://github.com/deepmodeling/dpgen/releases/)
[![doi:10.1016/j.cpc.2020.107206](https://img.shields.io/badge/DOI-10.1016%2Fj.cpc.2020.107206-blue)](https://doi.org/10.1016/j.cpc.2020.107206)
![Citations](https://citations.njzjz.win/10.1016/j.cpc.2020.107206)
[![conda install](https://img.shields.io/conda/dn/conda-forge/dpgen?label=conda%20install)](https://anaconda.org/conda-forge/dpgen)
[![pip install](https://img.shields.io/pypi/dm/dpgen?label=pip%20install)](https://pypi.org/project/dpgen)

DP-GEN (Deep Generator) is a software written in Python, delicately designed to generate a deep learning based model of interatomic potential energy and force field. DP-GEN is dependent on [DeePMD-kit](https://github.com/deepmodeling/deepmd-kit/). With highly scalable interface with common softwares for molecular simulation, DP-GEN is capable to  automatically prepare scripts and maintain job queues on HPC machines (High Performance Cluster) and analyze results.

If you use this software in any publication, please cite:

Yuzhi Zhang, Haidi Wang, Weijie Chen, Jinzhe Zeng, Linfeng Zhang, Han Wang, and Weinan E, DP-GEN: A concurrent learning platform for the generation of reliable deep learning based potential energy models, Computer Physics Communications, 2020, 253, 107206.

## Highlighted features
+ **Accurate and efficient**: DP-GEN is capable to sample more than tens of million structures and select only a few for first principles calculation. DP-GEN will finally obtain a uniformly accurate model.
+ **User-friendly and automatic**: Users may install and run DP-GEN easily. Once successfully running, DP-GEN can dispatch and handle all jobs on HPCs, and thus there's no need for any personal effort.
+ **Highly scalable**: With modularized code structures, users and developers can easily extend DP-GEN for their most relevant needs. DP-GEN currently supports for HPC systems ([Slurm](https://slurm.schedmd.com/), [PBS](https://www.openpbs.org/), LSF and cloud machines), Deep Potential interface with DeePMD-kit, MD interface with [LAMMPS](https://www.lammps.org/), [Gromacs](http://www.gromacs.org/), [AMBER](https://ambermd.org/), Calypso and *ab-initio* calculation interface with [VASP](https://www.vasp.at/), [PWSCF](https://www.quantum-espresso.org/), [CP2K](https://www.cp2k.org/), [SIESTA](https://departments.icmab.es/leem/siesta/), [Gaussian](https://gaussian.com/), Abacus, [PWmat](http://www.pwmat.com/), etc. We're sincerely welcome and embraced to users' contributions, with more possibilities and cases to use DP-GEN.

## Download and Install

DP-GEN only supports Python 3.8 and above.

One can download the source code of dpgen by
```bash
git clone https://github.com/deepmodeling/dpgen.git
```
then you may install DP-GEN easily by:
```bash
cd dpgen
pip install --user .
```
With this command, the dpgen executable is install to `$HOME/.local/bin/dpgen`. You may want to export the `PATH` by
```bash
export PATH=$HOME/.local/bin:$PATH
```
To test if the installation is successful, you may execute
```bash
dpgen -h
```

## Workflows and usage

DP-GEN contains the following workflows:

* [`dpgen run`](https://docs.deepmodeling.com/projects/dpgen/en/latest/run/): Main process of Deep Generator.
* [Init](https://docs.deepmodeling.com/projects/dpgen/en/latest/init/): Generating initial data.
  * `dpgen init_bulk`: Generating initial data for bulk systems.
  * `dpgen init_surf`: Generating initial data for surface systems.
  * `dpgen init_reaction`: Generating initial data for reactive systems.
* [`dpgen simplify`](https://docs.deepmodeling.com/projects/dpgen/en/latest/simplify/): Reducing the amount of existing dataset.
* [`dpgen autotest`](https://docs.deepmodeling.com/projects/dpgen/en/latest/autotest/): Autotest for Deep Potential.

For detailed usage and parameters, read [DP-GEN documentation](https://docs.deepmodeling.com/projects/dpgen/).

## Tutorials and examples

* [Tutorials](https://tutorials.deepmodeling.com/en/latest/Tutorials/DP-GEN/): basic tutorials for DP-GEN.
* [Examples](examples): input files in [JSON](https://docs.python.org/3/library/json.html) format.
* [Publications](https://deepmodeling.com/blog/papers/dpgen/): Published research articles using DP-GEN.
* [User guide](https://docs.deepmodeling.com/projects/dpgen/en/latest/user-guide/): frequently asked questions listed in troubleshooting.

## License
The project dpgen is licensed under [GNU LGPLv3.0](./LICENSE).

## Contributing

DP-GEN is maintained by [DeepModeling's developers](https://docs.deepmodeling.com/projects/dpgen/en/latest/credits.html). Contributors are always welcome.
