This is an example for `dpgen simplify`. `data` contains a simplistic data set based on MAPbI3-scan case. Since it has been greatly reduced, do not take it seriously. It is just a demo.
`simplify_example` is the work path, which contains `INCAR` and templates for `simplify.json` and `machine.json`. You can use the command `nohup dpgen simplify simplify.json machine.json 1>log 2>err &` here to test if `dpgen simplify` can run normally.

Kindly reminder:
1. `machine.json` is supported by `dpdispatcher 0.4.15`, please check https://docs.deepmodeling.com/projects/dpdispatcher/en/latest/ to update the parameters according to your `dpdispatcher` version.
2. `POTCAR` should be prepared by the user.
3. Please check the path and files name and make sure they are correct.
