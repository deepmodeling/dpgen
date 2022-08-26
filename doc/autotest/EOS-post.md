## EOS-post

The post processing of EOS would go to every directory in `confs/mp-*/eos_00` and do the post processing. Let's suppose we are now in `confs/mp-100/eos_00` and there are `task.000000, task.000001,..., task.000039` in this directory. By reading `inter.json` file in every task directory, the task type can be determined and the energy and force information of every task can further be obtained. By appending the `dict` of energy and force into a list, an example of the list with 1 atom is given as:
```txt
[
    {"energy": E1, "force": [fx1, fy1, fz1]},
    {"energy": E2, "force": [fx2, fy2, fz2]},
    ...
    {"energy": E40, "force": [fx40, fy40, fz40]}
]
```
Then the volume can be calculated from the task id and the corresponding energy can be obtained from the list above. Finally, there would be `result.json` in json format and `result.out` in txt format in `confs/mp-100/eos_00` containing the EOS results.

An example of `result.json` is give as:
```txt
{
    10.00: -3.0245,
    10.50: -3.0216,
         ...
    29.50: -7.9740
}
```

An example of `result.out` is given below:

```txt
conf_dir: confs/mp-100/eos_00
 VpA(A^3)  EpA(eV)
 10.000   -3.0245
 10.500   -3.0216
  ...        ...
 29.500   -7.9740
```

 
