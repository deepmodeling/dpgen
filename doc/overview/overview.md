# Overview

## About DP-GEN

[![GitHub release](https://img.shields.io/github/release/deepmodeling/dpgen.svg?maxAge=86400)](https://github.com/deepmodeling/dpgen/releases/)
[![doi:10.1016/j.cpc.2020.107206](https://img.shields.io/badge/DOI-10.1016%2Fj.cpc.2020.107206-blue)](https://doi.org/10.1016/j.cpc.2020.107206)
![Citations](https://citations.njzjz.win/10.1016/j.cpc.2020.107206)
[![conda install](https://img.shields.io/conda/dn/conda-forge/dpgen?label=conda%20install)](https://anaconda.org/conda-forge/dpgen)
[![pip install](https://img.shields.io/pypi/dm/dpgen?label=pip%20install)](https://pypi.org/project/dpgen)

DP-GEN (Deep Generator)  is a software written in Python, delicately designed to generate a deep learning based model of interatomic potential energy and force field. DP-GEN is dependent on [DeepMD-kit](https://github.com/deepmodeling/deepmd-kit/blob/master/README.md). With highly scalable interface with common softwares for molecular simulation, DP-GEN is capable to  automatically prepare scripts and maintain job queues on HPC machines (High Performance Cluster) and analyze results.

If you use this software in any publication, please cite:

Yuzhi Zhang, Haidi Wang, Weijie Chen, Jinzhe Zeng, Linfeng Zhang, Han Wang, and Weinan E, DP-GEN: A concurrent learning platform for the generation of reliable deep learning based potential energy models, Computer Physics Communications, 2020, 107206.

### Highlighted features
+ **Accurate and efficient**: DP-GEN is capable to sample more than tens of million structures and select only a few for first principles calculation. DP-GEN will finally obtain a uniformly accurate model.
+ **User-friendly and automatic**: Users may install and run DP-GEN easily. Once succusefully running, DP-GEN can dispatch and handle all jobs on HPCs, and thus there's no need for any personal effort.
+ **Highly scalable**: With modularized code structures, users and developers can easily extend DP-GEN for their most relevant needs. DP-GEN currently supports for HPC systems (Slurm, PBS, LSF and cloud machines ), Deep Potential interface with DeePMD-kit, MD interface with [LAMMPS](https://www.lammps.org/), [Gromacs](http://www.gromacs.org/)   and *ab-initio* calculation interface with VASP, PWSCF, CP2K, SIESTA and Gaussian, Abacus, PWMAT, etc . We're sincerely welcome and embraced to users' contributions, with more possibilities and cases to use DP-GEN.

## Download and install

DP-GEN only supports Python 3.9 and above. You can use one of the following methods to install DP-GEN:

- Install via pip: `pip install dpgen`
- Install via conda: `conda install -c conda-forge dpgen``
- Install from source code: `git clone https://github.com/deepmodeling/dpgen && pip install ./dpgen`

To test if the installation is successful, you may execute

```bash
dpgen -h
```

## Use DP-GEN

A quick-start on using DPGEN can be found [here](https://tutorials.deepmodeling.com/en/latest/Tutorials/DP-GEN/index.html). You can follow the [Handson tutorial](https://tutorials.deepmodeling.com/en/latest/Tutorials/DP-GEN/learnDoc/DP-GEN_handson.html), it is friendly to new users.


## Case Studies

- [Practical-Guidelines-for-DP](https://tutorials.deepmodeling.com/en/latest/CaseStudies/Practical-Guidelines-for-DP/index.html)

Before starting a new Deep Potential (DP) project, we suggest people (especially those who are newbies) read the following context first to get some insights into what tools we can use, what kinds of risks and difficulties we may meet, and how we can advance a new DP project smoothly.

- [Convergence-Test](https://tutorials.deepmodeling.com/en/latest/CaseStudies/Convergence-Test/index.html)

to ensure the data quality, the reliability of the final model, as well as the feasibility of the project, a convergence test should be done first.

- [Gas-phase](https://tutorials.deepmodeling.com/en/latest/CaseStudies/Gas-phase/index.html)

In this tutorial, we will take the simulation of methane combustion as an example and introduce the procedure of DP-based MD simulation.

- [Mg-Y_alloy](https://tutorials.deepmodeling.com/en/latest/CaseStudies/Mg-Y_alloy/index.html)

 We will briefly analyze the candidate configurational space of a metallic system by taking Mg-based Mg-Y binary alloy as an example. The task is divided into steps during the DP-GEN process.

- [Transfer-learning](https://tutorials.deepmodeling.com/en/latest/CaseStudies/Transfer-learning/index.html)

 This tutorial will introduce how to implement potential energy surface (PES) transfer-learning by using the DP-GEN software. In DP-GEN (version > 0.8.0), the “simplify” module is designed for this purpose.

## License
The project dpgen is licensed under [GNU LGPLv3.0](https://github.com/deepmodeling/dpgen/blob/master/LICENSE)
