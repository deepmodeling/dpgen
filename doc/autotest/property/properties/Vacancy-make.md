## Vacancy make

**Step 1.** The `VacancyGenerator` module in [pymatgen.analysis.defects.generators](https://pymatgen.org/pymatgen.analysis.defects.generators.html) is used to generate a set of structures with vacancy.

**Step 2.** If there are `init_from_suffix` and `output_suffix` parameter in the `properties` part, the [refine process](../../refine/Refine-get-started-and-input-examples) follows. If reproduce is evoked, the [reproduce process](../../reproduce/Reproduce-get-started-and-input-examples) follows. Otherwise, the vacancy structure (`POSCAR`) and supercell information (`supercell.out`) are written in the task directory, for example, in `confs/mp-*/vacancy_00/task.000000` with the check and possible removing of the old input files like before.

**Step 3.** When doing `vacancy` by VASP, `ISIF = 3`. When doing `vacancy` by LAMMPS, the same `in.lammps` as that in [EOS (change_box is True)](./EOS-make) would be generated with `scale` set to one.
