## Surface make

**Step 1.** Based on the equilibrium configuration, `generate_all_slabs` module in [pymatgen.core.surface](https://pymatgen.org/pymatgen.core.surface.html) would help to generate surface structure list with using `max_miller`, `min_slab_size`, and `min_vacuum_size` parameters. 

**Step 2.** If `refine` is True, we do [refine process](../../refine/Refine-get-started-and-input-examples). If `reprod-opt` is True (the default is False), we do [reproduce process](../../reproduce/Reproduce-get-started-and-input-examples). Otherwise, the surface structure (`POSCAR`) with perturbations in xz and miller index information (`miller.out`) are written in the task directory, for example, in `confs/mp-*/interstitial_00/task.000000` with the check and possible removing of the old input files like before.
