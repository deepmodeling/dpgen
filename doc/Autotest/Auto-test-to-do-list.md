For the implementation, one should do :
1. Clearly know the input/output of the function/class. How to handle exceptions.
2. Finish coding
3. Provide Unittest
4. Provide Document: what does the user provide in each section of the parameter file (json format)

common.py
- make_*
- run_*
- post_*

Property
- EOS
- Elastic
- Vacancy
- Interstitial
- Surface



Task:
- VASP
- DEEPMD_LMP
- MEAM_LMP


Specific functions:
1. Property.make_confs : Make configurations needed to compute the property. 
        The tasks directory will be named as path_to_work/task.xxxxxx
        IMPORTANT: handel the case when the directory exists.
2. Property.cmpt :  Compute the property. 
3. Task.make_input_file(Property.task_type): Prepare input files for a computational task.
        For example, the VASP prepares INCAR.
        LAMMPS (including DeePMD, MEAM...) prepares in.lammps.
        The parameter of this task will be stored in 'output_dir/task.json'