## EOS post

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
    "14.808453313267595": -3.7194474,
    "14.972991683415014": -3.7242038,
        ...
    "17.934682346068534": -3.7087655
}
```

An example of `result.out` is given below:

```txt
onf_dir: /root/auto_test_example/deepmd/confs/std-fcc/eos_00
 VpA(A^3)  EpA(eV)
 14.808   -3.7194
 14.973   -3.7242
   ...      ...
 17.935   -3.7088
```
