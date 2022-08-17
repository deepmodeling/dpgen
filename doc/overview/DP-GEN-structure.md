# Overview of the structure of the DP-GEN repository
Let's look at the home page of DP-GEN. https://github.com/deepmodeling/dpgen
````
├── build
├── CITATION.cff
├── conda
├── dist
├── doc
├── dpgen
├── dpgen.egg-info
├── examples
├── LICENSE
├── README.md
├── requirements.txt
├── setup.py
└── tests
````
In `examples`, we prepared templates for PARAM and MACHINE files for different software, versions and tasks.
Most of the code related to DP-GEN functions is in the `dpgen` directory. Open the `dpgen` directory, and we can see
````
├── arginfo.py
├── auto_test
├── collect
├── data
├── database
├── _date.py
├── dispatcher
├── generator
├── __init__.py
├── main.py
├── __pycache__
├── remote
├── simplify
├── tools
├── util.py
└── _version.py
````
`auto_test` corresponds to `dpgen autotest`
`collect` corresponds to `dpgen collect`
`data` corresponds to `dpgen init_bulk` and `dpgen init_surf`
`simplify` corresponds to `dpgen simplify`
A large part of `dispatcher` has been moved to `dpdispatcher`, see https://github.com/deepmodeling/dpdispatcher
`generator` is the core part of DP-GEN, let's open this folder.

````
├── arginfo.py
├── ch4
├── __init__.py
├── lib
└── run.py
````
run.py is the core of DP-GEN, corresponding to `dpgen run`, we can find `make_train`, `run_train`, ... `post_fp`, and other steps related functions here.

## Find how a parameter is used in the code
It is strongly recommended that you use the `find in files` function of Visual Studio software, or similar functions of other software. Type in the name of the parameter you're looking for, and you'll see where it's read in, and the procedure used.
Of course, you can also search for the relevant code according to the above guide.

## Want to modify a function?
If you have special requirements, you can make personalized modifications in the code corresponding to the function. If you think your modification can benefit the public, and it does not conflict with the current DP-GEN function; or if you fix a bug, please make a pull request to contribute the optimization to the DP-GEN repository.

## DP-GEN dependencies
dpdispatcher and dpdata are dependencies of DP-GEN. dpdispatcher is related to task submission, monitoring and recovery, and dpdata is related to data processing. If you encounter an error and want to find the reason, please judge whether the problem comes from DP-GEN, dpdispatcher or dpdata according to the last line of `Traceback`.

## About the update of the parameter file
You may have noticed that there are arginfo.py files in many folders. This is a file used to generate parameter documentation.
If you add or modify a parameter in DP-GEN and intend to export it to the main repository, please sync your changes in arginfo.
