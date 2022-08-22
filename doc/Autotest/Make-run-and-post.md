There are three operations in auto test package, namely `make`, `run`, and `post`. Here we take `eos` property as an example for property type.

### Make
The `INCAR`, `POSCAR`, `POTCAR` input files for VASP or `in.lammps`, `conf.lmp`, and the interatomic potential files for LAMMPS will be generated in the directory `confs/mp-*/relaxation/relax_task` for relaxation or `confs/mp-*/eos_00/task.[0-9]*[0-9]` for EOS. The `machine.json` file is not needed for `make`. Example: 
```bash
dpgen autotest make relaxation.json 
```

### Run
The jobs would be dispatched according to the parameter in `machine.json` file and the calculation results would be sent back. Example:
```bash
dpgen autotest run relaxation.json machine.json
```

### Post
The post process of calculation results would be performed. `result.json` in json format will be generated in `confs/mp-*/relaxation/relax_task` for relaxation and `result.json` in json format and `result.out` in txt format in `confs/mp-*/eos_00` for EOS. The `machine.json` file is also not needed for `post`. Example:
```bash
dpgen autotest post relaxation.json 
```
