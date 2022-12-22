# Init_surf

You may prepare initial data for surface systems with VASP by:

```bash
dpgen init_surf PARAM [MACHINE]
```
The MACHINE configure file is optional. If this parameter exists, then the optimization
tasks or MD tasks will be submitted automatically according to MACHINE.json. That is to say, if one only wants to prepare `surf-xxx/sys-xxx` folders for the second stage but wants to skip relaxation, `dpgen init_surf PARAM` should be used (without `MACHINE`).
"stages" and "skip_relax" in `PARAM` should be set as: 
```json
  "stages": [1,2],
  "skip_relax": true,
```

Basically `init_surf` can be divided into two parts , denoted as {dargs:argument}`stages <init_surf_jdata/stages>` in `PARAM`:
1. Build specific surface in folder `00.place_ele`
2. Pertub and scale in folder `01.scale_pert`

All stages must be **in order**.

Generally, `init_surf` does not run AIMD but only generates a lot of configurations. Compared with `init_bulk`, which runs DFT calculations twice, init_surf does once. Usually, we do `init_bulk`, run many rounds of DP-GEN iterations, collect enough data for the bulk system, and do `init_surf` after that. At this point, the lattice constant has been determined, and the lattice constant required for the initial configuration of `init_surf` can be used directly. These configurations made by `init_surf` are prepared for `01.model_devi`. Candidates will do DFT calculation in `02.fp`. 

- Generate vacuum layers

According to [the source code of pert_scaled](https://github.com/deepmodeling/dpgen/blob/8dea29ef125f66be9641afe5ac4970433a9c9ce1/dpgen/data/surf.py#L484), init_surf will generate a series of surface structures with specified separations between the sample layer and its periodic image. There are two ways to specify the interval in generating the vacuum layers: 1) to set the interval value and 2) to set the number of intervals.

You can use {dargs:argument}`layer_numb <init_surf_jdata/layer_numb>` (the number of layers of the slab) or {dargs:argument}`z_min <init_surf_jdata/z_min>` (the total thickness) to specify the thickness of the atoms below. Then `vacuum_*` parameters specify the vacuum layers above. `dpgen init_surf` will make a series of structures with the thickness of vacuum layers from {dargs:argument}`vacuum_min <init_surf_jdata/vacuum_min>` to {dargs:argument}`vacuum_max <init_surf_jdata/vacuum_max>`. The number of vacuum layers is controlled by the parameter {dargs:argument}`vacuum_resol <init_surf_jdata/vacuum_resol>`. 

The layers will be generated even when the size of {dargs:argument}`vacuum_resol <init_surf_jdata/vacuum_resol>` is 1. When the size of {dargs:argument}`vacuum_resol <init_surf_jdata/vacuum_resol>` is 2 or it is empty, the whole interval range is divided into the nearby region with denser intervals (head region) and the far-away region with sparser intervals (tail region), which are divided by {dargs:argument}`mid_point <init_surf_jdata/mid_point>`. 

When the size of {dargs:argument}`vacuum_resol <init_surf_jdata/vacuum_resol>` is 2, two elements respectively decide the number of intervals in head region and tail region.

When {dargs:argument}`vacuum_resol <init_surf_jdata/vacuum_resol>` is empty, the number of intervals in the head region = vacuum_num * head_ratio. {dargs:argument}`vacuum_num <init_surf_jdata/vacuum_numb>` and {dargs:argument}`head_ratio <init_surf_jdata/head_ratio>` are both keys in `param.json`.

- Attach files in the task path

One can use the machine parameter `forward_files` to upload other files besides POSCAR, INCAR, and POTCAR. For example, "vdw_kernal.bindat" for each task. 

See [the document of task parameters](https://docs.deepmodeling.com/projects/dpdispatcher/en/latest/task.html#argument:task/forward_files).

Following is an example for `PARAM`, which generates data from a typical structure fcc.
```json
{
  "stages": [
    1,
    2
  ],
  "cell_type": "fcc",
  "latt": 4.034,
  "super_cell": [
    2,
    2,
    2
  ],
  "layer_numb": 3,
  "vacuum_max": 9.0,
  "vacuum_resol": [
    0.5,
    1
  ],
  "mid_point": 4.0,
  "millers": [
    [
      1,
      0,
      0
    ],
    [
      1,
      1,
      0
    ],
    [
      1,
      1,
      1
    ]
  ],
  "elements": [
    "Al"
  ],
  "potcars": [
    "....../POTCAR"
  ],
  "relax_incar": "....../INCAR_metal_rlx_low",
  "scale": [
    1.0
  ],
  "skip_relax": true,
  "pert_numb": 2,
  "pert_box": 0.03,
  "pert_atom": 0.01,
  "_comment": "that's all"
}
```

Another example is `from_poscar` method. Here you need to specify the POSCAR file. 

```json
{
  "stages": [
    1,
    2
  ],
  "cell_type": "fcc",
  "from_poscar":	true,
  "from_poscar_path":	"POSCAR",
  "super_cell": [
    1,
    1,
    1
  ],
  "layer_numb": 3,
  "vacuum_max": 5.0,
  "vacuum_resol": [0.5,2],
  "mid_point": 2.0,
  "millers": [
    [
      1,
      0,
      0
    ]
  ],
  "elements": [
    "Al"
  ],
  "potcars": [
    "./POTCAR"
  ],
  "relax_incar" : "INCAR_metal_rlx_low",
  "scale": [
    1.0
  ],
  "skip_relax": true,
  "pert_numb": 5,
  "pert_box": 0.03,
  "pert_atom": 0.01,
  "coll_ndata": 5000,
  "_comment": "that's all"
}
```
