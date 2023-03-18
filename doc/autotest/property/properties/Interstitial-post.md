## Interstitial post

For `Interstitial`, we need to calculate the energy difference between a crystal structure with and without atom added in.
The examples of the output files `result.json` in json format and `result.out` in txt format are given below.

#### result.json
```json
{
    "Al-[3, 3, 3]-task.000000": [
        4.022952000000004,
        -100.84773,
        -104.870682
    ],
    "Al-[3, 3, 3]-task.000001": [
        2.7829520000000088,
        -102.08773,
        -104.870682
    ]
}
```

#### result.out
```txt
/root/auto_test_example/deepmd/confs/std-fcc/interstitial_00
Insert_ele-Struct: Inter_E(eV)  E(eV) equi_E(eV)
Al-[3, 3, 3]-task.000000:   4.023  -100.848 -104.871
Al-[3, 3, 3]-task.000001:   2.783  -102.088 -104.871
```
