## Relaxation-post

Take `deepmd` post for example:
```bash
dpgen autotest post relaxation.json
tree confs/std-fcc/relaxation/
```
the output will be:
```
confs/std-fcc/relaxation/
|-- frozen_model.pb -> ../../../frozen_model.pb
|-- in.lammps
|-- jr.json
`-- relax_task
    |-- conf.lmp
    |-- CONTCAR
    |-- dump.relax
    |-- frozen_model.pb -> ../frozen_model.pb
    |-- in.lammps -> ../in.lammps
    |-- inter.json
    |-- log.lammps
    |-- outlog
    |-- POSCAR -> ../../POSCAR
    |-- result.json
    `-- task.json
```
`result.json` stores the box cell, coordinates, energy, force, virial,... information of each frame in the relaxation trajectory and `CONTCAR` is the final equilibrium configuration.

`result.json`:

```json
{
    "@module": "dpdata.system",
    "@class": "LabeledSystem",
    "data": {
        "atom_numbs": [
            1
        ],
        "atom_names": [
            "Al"
        ],
        "atom_types": {
            "@module": "numpy",
            "@class": "array",
            "dtype": "int64",
            "data": [
                0
            ]
        },
        "orig": {
            "@module": "numpy",
            "@class": "array",
            "dtype": "int64",
            "data": [
                0,
                0,
                0
            ]
        },
        "cells": {
            "@module": "numpy",
            "@class": "array",
            "dtype": "float64",
            "data": [
                [
                    [
                        2.8637824638,
                        0.0,
                        0.0
                    ],
                    [
                        1.4318912319,
                        2.4801083646,
                        0.0
                    ],
                    [
                        1.4318912319,
                        0.8267027882,
                        2.3382685902
                    ]
                ],
                [
                    [
                        2.8549207998018438,
                        0.0,
                        0.0
                    ],
                    [
                        1.4274603999009239,
                        2.472433938457684,
                        0.0
                    ],
                    [
                        1.4274603999009212,
                        0.8241446461525599,
                        2.331033071844216
                    ]
                ],
                [
                    [
                        2.854920788303194,
                        0.0,
                        0.0
                    ],
                    [
                        1.427460394144466,
                        2.472433928487206,
                        0.0
                    ],
                    [
                        1.427460394154763,
                        0.8241446428350139,
                        2.331033062460779
                    ]
                ]
            ]
        },
        "coords": {
            "@module": "numpy",
            "@class": "array",
            "dtype": "float64",
            "data": [
                [
                    [
                        0.0,
                        0.0,
                        0.0
                    ]
                ],
                [
                    [
                        5.709841595683707e-25,
                        -4.3367974740910857e-19,
                        0.0
                    ]
                ],
                [
                    [
                        -8.673606219968035e-19,
                        8.673619637565944e-19,
                        8.673610853102186e-19
                    ]
                ]
            ]
        },
        "energies": {
            "@module": "numpy",
            "@class": "array",
            "dtype": "float64",
            "data": [
                -3.745029,
                -3.7453815,
                -3.7453815
            ]
        },
        "forces": {
            "@module": "numpy",
            "@class": "array",
            "dtype": "float64",
            "data": [
                [
                    [
                        0.0,
                        -6.93889e-18,
                        -3.46945e-18
                    ]
                ],
                [
                    [
                        1.38778e-17,
                        6.93889e-18,
                        -1.73472e-17
                    ]
                ],
                [
                    [
                        1.38778e-17,
                        1.73472e-17,
                        -4.51028e-17
                    ]
                ]
            ]
        },
        "virials": {
            "@module": "numpy",
            "@class": "array",
            "dtype": "float64",
            "data": [
                [
                    [
                        -0.07534992071654338,
                        1.2156615579052586e-17,
                        1.3904892126132796e-17
                    ],
                    [
                        1.2156615579052586e-17,
                        -0.07534992071654338,
                        4.61571024026576e-12
                    ],
                    [
                        1.3904892126132796e-17,
                        4.61571024026576e-12,
                        -0.07534992071654338
                    ]
                ],
                [
                    [
                        -9.978994290457664e-08,
                        -3.396452753975288e-15,
                        8.785831629151552e-16
                    ],
                    [
                        -3.396452753975288e-15,
                        -9.991375413666671e-08,
                        5.4790751628409565e-12
                    ],
                    [
                        8.785831629151552e-16,
                        5.4790751628409565e-12,
                        -9.973497959053003e-08
                    ]
                ],
                [
                    [
                        1.506940521266962e-11,
                        1.1152016233536118e-11,
                        -8.231900529157644e-12
                    ],
                    [
                        1.1152016233536118e-11,
                        -6.517665029355618e-11,
                        -6.33706710415926e-12
                    ],
                    [
                        -8.231900529157644e-12,
                        -6.33706710415926e-12,
                        5.0011471096530724e-11
                    ]
                ]
            ]
        },
        "stress": {
            "@module": "numpy",
            "@class": "array",
            "dtype": "float64",
            "data": [
                [
                    [
                        -7.2692250000000005,
                        1.1727839e-15,
                        1.3414452e-15
                    ],
                    [
                        1.1727839e-15,
                        -7.2692250000000005,
                        4.4529093000000003e-10
                    ],
                    [
                        1.3414452e-15,
                        4.4529093000000003e-10,
                        -7.2692250000000005
                    ]
                ],
                [
                    [
                        -9.71695e-06,
                        -3.3072633e-13,
                        8.5551193e-14
                    ],
                    [
                        -3.3072633e-13,
                        -9.729006000000001e-06,
                        5.3351969e-10
                    ],
                    [
                        8.5551193e-14,
                        5.3351969e-10,
                        -9.711598e-06
                    ]
                ],
                [
                    [
                        1.4673689e-09,
                        1.0859169e-09,
                        -8.0157343e-10
                    ],
                    [
                        1.0859169e-09,
                        -6.3465139e-09,
                        -6.1706584e-10
                    ],
                    [
                        -8.0157343e-10,
                        -6.1706584e-10,
                        4.8698191e-09
                    ]
                ]
            ]
        }
    }
}
```
