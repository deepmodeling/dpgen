# Init_surf

You may prepare initial data for surface systems with VASP by:

```bash
dpgen init_surf PARAM [MACHINE]
```
The MACHINE configure file is optional. If this parameter exists, then the optimization
tasks or MD tasks will be submitted automatically according to MACHINE.json.


Basically `init_surf` can be divided into two parts , denoted as `stages` in `PARAM`:
1. Build specific surface in folder `00.place_ele`
2. Pertub and scale in folder `01.scale_pert`

All stages must be **in order**.


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
  "vacuum_max": 9,
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
  "vacuum_max": 5,
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
