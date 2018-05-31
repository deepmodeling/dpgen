#!/bin/bash

targets=`seq 3.78 0.01 4.33`
lmp_cmd=/home/wanghan/Soft/lammps/lammps-16Mar18/src/lmp_mpi

for ii in $targets
do
    sed "s/LATT_A/$ii/g" in.fcc > in.tmp
    $lmp_cmd -i in.tmp &> tmp.out
    epa=`grep ENER_PER_ATOM tmp.out | awk '{print $2}'`
    epv=`grep VOLM_PER_ATOM tmp.out | awk '{print $2}'`
    echo $epv $epa
done

rm -f in.tmp
