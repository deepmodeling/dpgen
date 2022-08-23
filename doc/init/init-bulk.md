## Init: Preparing Initial Data

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
`init_bulk` support both VASP and ABACUS for first-principle calculation. You can choose the software by specifying the key `init_fp_style`. If `init_fp_style` is not specified, the default software will be VASP. 

When using ABACUS for `init_fp_style`, the keys of the paths of `INPUT` files for relaxation and MD simulations are the same as `INCAR` for VASP, which are `relax_incar` and `md_incar` respectively. Use `relax_kpt` and `md_kpt` for the relative path for `KPT` files of relaxation and MD simulations. They two can be omitted if `kspacing` (in unit of 1/Bohr) or `gamma_only` has been set in corresponding INPUT files. If `from_poscar` is set to `false`, you have to specify `atom_masses` in the same order as `elements`.

The following table gives explicit descriptions on keys in `PARAM`.

The bold notation of key (such as **Elements**) means that it's a necessary key.

 Key  | Type          | Example                                                      | Description                                                      |
| :---------------- | :--------------------- | :-------------------------------------- | :-------------------------------------------------------------|
| **stages** | List of Integer | [1,2,3,4] | Stages for `init_bulk`
| **Elements** | List of String | ["Mg"] | Atom types
|  cell_type | String  | "hcp" | Specifying which typical structure to be generated. **Options** include fcc, hcp, bcc, sc, diamond.
| latt | Float | 4.479 | Lattice constant for single cell.
| from_poscar | Boolean | True | Deciding whether to use a given poscar as the beginning of relaxation. If it's true, keys (`cell_type`, `latt`) will be aborted. Otherwise, these two keys are **necessary**.
| from_poscar_path | String | "....../C_mp-47_conventional.POSCAR" | Path of POSCAR for VASP or STRU for ABACUS. **Necessary** if `from_poscar` is true.
| relax_incar | String | "....../INCAR" | Path of INCAR for VASP or INPUT for ABACUS for relaxation in VASP. **Necessary** if `stages` include 1.
| md_incar | String |  "....../INCAR_md" | Path of INCAR for VASP or INPUT for ABACUS for MD in VASP. **Necessary** if `stages` include 3.|
| **scale** | List of float | [0.980, 1.000, 1.020] | Scales for transforming cells.
| **skip_relax** | Boolean | False | If it's true, you may directly run stage 2 (perturb and scale) using an unrelaxed POSCAR.
| **pert_numb** | Integer | 30 | Number of perturbations for each POSCAR.
| **pert_box** | Float | 0.03 | Percentage of Perturbation for cells.
| **pert_atom** | Float | 0.01 | Perturbation of each atoms (Angstrom).
| **md_nstep** | Integer | 10 | Steps of AIMD in stage 3. If it's not equal to settings via `NSW` in `md_incar`, DP-GEN will follow `NSW`.
| **coll_ndata** | Integer | 5000 | Maximal number of collected data.
| type_map | List | [ "Mg", "Al"] | The indices of elements in deepmd formats will be set in this order.
| init_fp_style | String | "ABACUS" or "VASP" | First-principle software. If this key is absent, the default value will be "VASP".
| relax_kpt | String | "....../KPT" | Path of `KPT` file for relaxation in stage 1. Only useful if `init_fp_style` is "ABACUS".
| md_kpt | String | "....../KPT" | Path of `KPT` file for MD simulations in stage 3. Only useful if `init_fp_style` is "ABACUS".
| atom_masses | List of float | [24] | List of atomic masses of elements. The order should be the same as `Elements`. Only useful if `init_fp_style` is "ABACUS".
