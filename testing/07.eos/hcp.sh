#!/bin/bash

targets=`seq 0.95 0.005 1.05`
lmp_cmd=/home/wanghan/Soft/lammps/lammps-16Mar18/src/lmp_mpi

for ii in $targets
do
    sed "s/SCALE/$ii/g" in.hcp > in.tmp
    $lmp_cmd -i in.tmp &> tmp.out
    epa=`grep ENER_PER_ATOM tmp.out | awk '{print $2}'`
    epv=`grep VOLM_PER_ATOM tmp.out | awk '{print $2}'`
    coa=`grep COA tmp.out | awk '{print $2}'`
    echo $epv $epa $coa
done

