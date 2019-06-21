<span style="font-size:larger;">dpgen Manual</span>
========

# Table of contents
- [About dpgen](#About-dpgen)
- [genenrator](#generator)
- [auto_test](#auto_test)
  - [How to write `param.json`](##How-to-write-`param.json`)
  - [How to write `machine.json`](##How-to-write-`machine.json`)
- [Troubleshooting](#Troubleshooting)
# About dpgen
The deep potential generator
# generator
# auto_test
At this step, we assume that you have prepared some graph files like `graph.*.pb` and the particular pseudopotential `POTCAR`.
Then you only need one command to achieve automatic testing of physical properties.
```
python run.py param.json machine.json
```
## How to write `param.json`
We take Cu as an example to show the parameter settings of `param.json`.

The first part is the fundamental setting for particular alloy system. 
```
    "_comment": "models",
    "potcar_map" : {
	"Cu" : "/somewhere/POTCAR",
        "Zr" : "/elsewhere/POTCAR"
    },
    "conf_dir":"the folder of configuration",
    "task_type":"deepmd",
    "task":"eos",
```
You need to add the specified paths of necessary `POTCAR` files in "potcar_map". The different `POTCAR` paths are separated by commas.
Then you also need to add the folder path of particular configuration, which contains `POSCAR` file. For your convenience, we recommend that you use `gen_confs.py` to generate configurations by the following command.
```
python gen_confs.py Cu
```
It will store the various configurations of the given element or alloy in **confs** folder.

So far, there are  3 test types for "task_type" (i.e. "vasp", "deepmd" and "meam") and 7 test items for "task" as follows.

>00.**equi**:(default task) the equilibrium state, return energy and volume per atom

>01.**eos**: the equation of state

>02.**elastic**: the elasticity like Young's module

>03.**vacancy**: the vacancy defect 

>04.**interstitial**: the interstitial defect

>05.**surf**: the surface energy

>06.**phonon**:(beta version) the phonon specturm

## How to write `machine.json`

# Troubleshooting
