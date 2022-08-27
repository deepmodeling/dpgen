## Elastic-post

The `ElasticTensor` module in [pymatgen.analysis.elasticity.elastic](https://pymatgen.org/pymatgen.analysis.elasticity.elastic.html) is used to get the elastic tensor, Bv, and Gv. The mechanical properties of a crystal structure would be written in `result.json` in json format and `result.out` in txt format. The example of the output file is give below.

#### result.json

```json
{
    "elastic_tensor": [
    130.24, 114.43, 96.89, 0.00, 0.00, 0.00, 112.91, 132.64, 97.92, 0.00, 0.00, -0.00, 96.10, 96.10, 263.48, 
    0.00, 0.00, 0.00, 0.00, -0.00, -0.00, 33.50, 0.00, 0.00, 0.00, -0.00, -0.00, 0.00, 33.29, 0.00, 0.00, 0.00, 
   -0.00, -0.00, -0.00, 18.70],
    "BV": 126.75,
    "GV": 31.57,
    "EV": 87.46,
    "uv": 0.38
}
```

The order of `elastic_tensor` is C11, C12, ..., C16, C21, C22, ..., C26, ..., C66 and the unit of Bv, Gv, Ev, and uv is GPa.

#### result.out

```txt
conf_dir: confs/mp-100/elastic_00
 130.24  114.43   96.89    0.00    0.00    0.00
 112.91  132.64   97.92    0.00    0.00   -0.00
  96.10   96.10  263.48    0.00    0.00    0.00
   0.00   -0.00   -0.00   33.50    0.00    0.00
   0.00   -0.00   -0.00    0.00   33.29    0.00
   0.00    0.00   -0.00   -0.00   -0.00   18.70
# Bulk   Modulus BV = 126.75 GPa
# Shear  Modulus GV = 31.57 GPa
# Youngs Modulus EV = 87.46 GPa
# Poission Ratio uV = 0.38
```
