## Init_bulk

You may prepare initial data for bulk systems with VASP by:

```bash
dpgen init_bulk PARAM [MACHINE]
```
The MACHINE configure file is optional. If this parameter exists, then the optimization
tasks or MD tasks will be submitted automatically according to MACHINE.json.

Basically `init_bulk` can be divided into four parts , denoted as `stages` in `PARAM`:
1. Relax in folder `00.place_ele`
2. Perturb and scale in folder `01.scale_pert`
3. Run a short AIMD in folder `02.md`
4. Collect data in folder `02.md`.

All stages must be **in order**. One doesn't need to run all stages. For example, you may run stage 1 and 2, generating supercells as starting point of exploration in `dpgen run`.

If MACHINE is None, there should be only one stage in stages. Corresponding tasks will be generated, but user's intervention should be involved in, to manually run the scripts.

Following is an example for `PARAM`, which generates data from a typical structure hcp.
```json
{
    "stages" : [1,2,3,4],
    "cell_type":    "hcp",
    "latt":     4.479,
    "super_cell":   [2, 2, 2],
    "elements":     ["Mg"],
    "potcars":      ["....../POTCAR"],
    "relax_incar": "....../INCAR_metal_rlx",
    "md_incar" : "....../INCAR_metal_md",
    "scale":        [1.00],
    "skip_relax":   false,
    "pert_numb":    2,
    "md_nstep" : 5,
    "pert_box":     0.03,
    "pert_atom":    0.01,
    "coll_ndata":   5000,
    "type_map" : [ "Mg", "Al"],
    "_comment":     "that's all"
}
```

If you want to specify a structure as starting point for `init_bulk`, you may set in `PARAM` as follows.

```json
"from_poscar":	true,
"from_poscar_path":	"....../C_mp-47_conventional.POSCAR",
```
`init_bulk` supports both VASP and ABACUS for first-principle calculation. You can choose the software by specifying the key {dargs:argument}`init_fp_style <init_bulk_jdata/init_fp_style>`. If {dargs:argument}`init_fp_style <init_bulk_jdata/init_fp_style>` is not specified, the default software will be VASP. 

When using ABACUS for {dargs:argument}`init_fp_style <init_bulk_jdata/init_fp_style>`, the keys of the paths of `INPUT` files for relaxation and MD simulations are the same as `INCAR` for VASP, which are {dargs:argument}`relax_incar <init_bulk_jdata/relax_incar>` and {dargs:argument}`md_incar <init_bulk_jdata/md_incar>` respectively. Use {dargs:argument}`relax_kpt <init_bulk_jdata[ABACUS]/relax_kpt>` and {dargs:argument}`md_kpt <init_bulk_jdata[ABACUS]/md_kpt>` for the relative path for `KPT` files of relaxation and MD simulations. They two can be omitted if `kspacing` (in unit of 1/Bohr) or `gamma_only` has been set in corresponding INPUT files. If {dargs:argument}`from_poscar <init_bulk_jdata/from_poscar>` is set to `false`, you have to specify {dargs:argument}`atom_masses <init_bulk_jdata[ABACUS]/atom_masses>` in the same order as `elements`.
