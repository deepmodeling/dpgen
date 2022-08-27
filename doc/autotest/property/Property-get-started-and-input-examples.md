## Property get started and input examples

Here we take deepmd for example and the input file for other task types is similar.

```json
{
    "structures":       ["confs/std-*"],
    "interaction": {
        "type":          "deepmd",
        "model":         "frozen_model.pb",
        "type_map":     {"Al": 0}
    },
    "properties": [
        {
         "type":         "eos",
         "vol_start":    0.9,
         "vol_end":      1.1,
         "vol_step":     0.01
        },
        {
         "type":         "elastic",
         "norm_deform":  1e-2,
         "shear_deform": 1e-2
        },
        {
         "type":             "vacancy",
         "supercell":        [3, 3, 3],
         "start_confs_path": "../vasp/confs"
        },
        {
         "type":         "interstitial",
         "supercell":   [3, 3, 3],
         "insert_ele":  ["Al"],
         "conf_filters":{"min_dist": 1.5},
         "cal_setting": {"input_prop": "lammps_input/lammps_high"}
        },
        {
         "type":           "surface",
         "min_slab_size":  10,
         "min_vacuum_size":11,
         "max_miller":     2,
         "cal_type":       "static"
        },
        {
         "type": "gamma",
         "lattice_type": "fcc",
         "miller_index": [1, 1, 1],
         "displace_direction": [1, 1, 0],
         "supercell_size": [1, 1, 10],
         "min_vacuum_size": 10,
         "add_fix": ["true", "true", "false"],
         "n_steps": 20
        }
        ]
}
```
Universal key words for properties

Key words | data structure | example | description
---|---|---|---
**type** | String | "eos" | property type
skip | Boolean | true | whether to skip current property or not
start_confs_path | String | "../vasp/confs" | start from the equilibrium configuration in other path only for the current property type
cal_setting["input_prop"] | String | "lammps_input/lammps_high" |input commands file 
cal_setting["overwrite_interaction"] | Dict | | overwrite the interaction in the `interaction` part only for the current property type

other parameters in `cal_setting` and `cal_type` in `relaxation` also apply in `property`.

Key words for **EOS**

Key words | data structure | example | description
---|---|---|---
**vol_start** | Float | 0.9 | the starting volume related to the equilibrium structure
**vol_end** | Float | 1.1 | the biggest volume related to the equilibrium structure
**vol_step** | Float | 0.01 | the volume increment related to the equilibrium structure
vol_abs | Boolean | false | whether to treat vol_start, vol_end and vol_step as absolute volume or not (as relative volume), default = false

Key words for **Elastic**

Key words | data structure | example | description
---|---|---|---
norm_deform | Float | 1e-2 | deformation in xx, yy, zz, default = 1e-2
shear_deform | Float | 1e-2 | deformation in other directions, default = 1e-2

Key words for **Vacancy**

Key words | data structure | example | description
---|---|---|---
supercell | List of Int | [3,3,3] | the supercell to be constructed, default = [1,1,1]

Key words for **Interstitial**

Key words | data structure | example | description
---|---|---|---
**insert_ele** | List of String | ["Al"] | the element to be inserted
supercell | List of Int | [3,3,3] | the supercell to be constructed, default = [1,1,1]
conf_filters | Dict | "min_dist": 1.5 | filter out the undesirable configuration
bcc_self | Boolean | false | whether to do the self-interstitial calculations for bcc structures, default = false

Key words for **Surface**

Key words | data structure | example | description
---|---|---|---
**min_slab_size** | Int | 10 | minimum size of slab thickness
**min_vacuum_size** | Int | 11 | minimum size of vacuum width
pert_xz | Float | 0.01 | perturbation through xz direction used to compute surface energy, default = 0.01
max_miller | Int | 2 | the maximum miller index, default = 2

Key words for **Gamma**

Key words | data structure | example | description
---|---|---|---
**lattice_type** | String | "fcc" | "bcc" or "fcc" at this stage
**miller_index** | List of Int | [1,1,1] | slip plane for gamma-line calculation
**displace_direction** | List of Int | [1,1,0] | slip direction for gamma-line calculation
supercell_size | List of Int | [1,1,10] | the supercell to be constructed, default = [1,1,5]
min_vacuum_size | Int or Float | 10 | minimum size of vacuum width, default = 20
add_fix | List of String | ['true','true','false'] | whether to fix atoms in the direction, default = ['true','true','false'] (standard method)
n_steps | Int | 20 | Number of points for gamma-line calculation, default = 10