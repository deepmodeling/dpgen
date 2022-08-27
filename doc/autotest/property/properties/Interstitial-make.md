## Interstitial-make

**Step 1.** For each element in `insert_ele` list, `InterstitialGenerator` module in [pymatgen.analysis.defects.generators](https://pymatgen.org/pymatgen.analysis.defects.generators.html) would help to generate interstitial structure. The structure would be appended into a list if it can meet the requirements in `conf_filters`.

**Step 2.** If `refine` is `True`, we do [refine process](https://github.com/deepmodeling/dpgen/wiki/Refine:-get-started-and-input-examples). If `reprod-opt` is `True` (the default is **False**), we do [reproduce process](https://github.com/deepmodeling/dpgen/wiki/Reproduce:-get-started-and-input-examples). Else, the vacancy structure (`POSCAR`) and supercell information (`supercell.out`) are written in the task directory, for example, in `confs/mp-*/interstitial_00/task.000000` with the check and possible removing of the old input files like before.

**Step 3.** In `interstitial` by VASP, `ISIF = 3`. In `interstitial` by LAMMPS, the same `in.lammps` as that in [EOS (change_box is True)](https://github.com/deepmodeling/dpgen/wiki/EOS:-make) would be generated with `scale` set to one. 
