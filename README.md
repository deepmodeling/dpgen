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
        "Mg" : "/elsewhere/POTCAR"
    },
    "conf_dir":"the folder of configuration",
    "task_type":"deepmd",
    "task":"eos",
```
You need to add the specified paths of necessary `POTCAR` files in "potcar_map". The different `POTCAR` paths are separated by commas.

## How to write `machine.json`

# Troubleshooting
