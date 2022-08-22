```bash
dpgen autotest post reproduce.json
```

the output will be:

`result.out:`

```
/root/auto_test_example/deepmd/confs/std-fcc/interstitial_reprod
Reproduce: Initial_path Init_E(eV/atom)  Reprod_E(eV/atom)  Difference(eV/atom)
.../vasp/confs/std-fcc/interstitial_00/task.000000  -3.020   -3.240   -0.220
.../vasp/confs/std-fcc/interstitial_00/task.000000  -3.539   -3.541   -0.002
.../vasp/confs/std-fcc/interstitial_00/task.000000  -3.582   -3.582   -0.001
.../vasp/confs/std-fcc/interstitial_00/task.000000  -3.582   -3.581    0.001
.../vasp/confs/std-fcc/interstitial_00/task.000000  -3.594   -3.593    0.001
.../vasp/confs/std-fcc/interstitial_00/task.000000  -3.594   -3.594    0.001
.../vasp/confs/std-fcc/interstitial_00/task.000000  -3.598   -3.597    0.001
.../vasp/confs/std-fcc/interstitial_00/task.000000  -3.600   -3.600    0.001
.../vasp/confs/std-fcc/interstitial_00/task.000000  -3.600   -3.600    0.001
.../vasp/confs/std-fcc/interstitial_00/task.000000  -3.601   -3.600    0.001
.../vasp/confs/std-fcc/interstitial_00/task.000000  -3.602   -3.601    0.001
.../vasp/confs/std-fcc/interstitial_00/task.000000  -3.603   -3.602    0.001
.../vasp/confs/std-fcc/interstitial_00/task.000000  -3.603   -3.602    0.001
.../vasp/confs/std-fcc/interstitial_00/task.000000  -3.603   -3.602    0.001
.../vasp/confs/std-fcc/interstitial_00/task.000000  -3.603   -3.602    0.001
.../vasp/confs/std-fcc/interstitial_00/task.000000  -3.603   -3.602    0.001
.../vasp/confs/std-fcc/interstitial_00/task.000000  -3.603   -3.602    0.001
.../vasp/confs/std-fcc/interstitial_00/task.000000  -3.603   -3.602    0.001
.../vasp/confs/std-fcc/interstitial_00/task.000001  -3.345   -3.372   -0.027
.../vasp/confs/std-fcc/interstitial_00/task.000001  -3.546   -3.556   -0.009
.../vasp/confs/std-fcc/interstitial_00/task.000001  -3.587   -3.593   -0.007
.../vasp/confs/std-fcc/interstitial_00/task.000001  -3.593   -3.599   -0.006
.../vasp/confs/std-fcc/interstitial_00/task.000001  -3.600   -3.606   -0.006
.../vasp/confs/std-fcc/interstitial_00/task.000001  -3.600   -3.606   -0.006
.../vasp/confs/std-fcc/interstitial_00/task.000001  -3.624   -3.631   -0.006
.../vasp/confs/std-fcc/interstitial_00/task.000001  -3.634   -3.640   -0.007
.../vasp/confs/std-fcc/interstitial_00/task.000001  -3.637   -3.644   -0.007
.../vasp/confs/std-fcc/interstitial_00/task.000001  -3.637   -3.644   -0.007
.../vasp/confs/std-fcc/interstitial_00/task.000001  -3.638   -3.645   -0.007
.../vasp/confs/std-fcc/interstitial_00/task.000001  -3.638   -3.645   -0.007
.../vasp/confs/std-fcc/interstitial_00/task.000001  -3.639   -3.646   -0.007
.../vasp/confs/std-fcc/interstitial_00/task.000001  -3.639   -3.646   -0.007
.../vasp/confs/std-fcc/interstitial_00/task.000001  -3.639   -3.646   -0.007
.../vasp/confs/std-fcc/interstitial_00/task.000001  -3.639   -3.646   -0.007
.../vasp/confs/std-fcc/interstitial_00/task.000001  -3.639   -3.646   -0.007
.../vasp/confs/std-fcc/interstitial_00/task.000001  -3.639   -3.646   -0.007
.../vasp/confs/std-fcc/interstitial_00/task.000001  -3.639   -3.646   -0.007
.../vasp/confs/std-fcc/interstitial_00/task.000001  -3.639   -3.646   -0.007
.../vasp/confs/std-fcc/interstitial_00/task.000001  -3.639   -3.646   -0.007
```

the comparison of the initial and reproduced results as well as the absolute path of the initial data is recorded.

`result.json:`

```json
{
    "/root/auto_test_example/vasp/confs/std-fcc/interstitial_00/task.000000": {
        "nframes": 18,
        "error": 0.0009738182472213228
    },
    "/root/auto_test_example/vasp/confs/std-fcc/interstitial_00/task.000001": {
        "nframes": 21,
        "error": 0.0006417039154057605
    }
}
```

the error analysis corresponding to the initial data is recorded and the error of the first frame is disregarded when all the frames are considered in reproduce.