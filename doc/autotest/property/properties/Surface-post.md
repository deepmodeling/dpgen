## Surface-post

For `Surface`, we need to calculate the energy difference between a crystal structure with and without a surface with a certain miller index divided by the surface area.

The examples of the output files `result.json` in json format and `result.out` in txt format are given below.

#### result.json
```txt
{
    [1,1,1]: [1.200  -297.099 -298.299],
    [1,2,1]: [1.400  -296.899 -298.299],
    ...
}
```

#### result.out
```txt
conf_dir: confs/mp-100/surface_00
Miller_Indices:     Surf_E(J/m^2) EpA(eV) equi_EpA(eV)
[1,1,1]:       1.200  -297.099 -298.299
[1,2,1]:       1.400  -296.899 -298.299
```
