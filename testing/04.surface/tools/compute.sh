#!/bin/bash

lmp_cmd=$HOME/local/bin/lmp_mpi_09

if test $# -ne 1; then
    echo usage
    echo $0 Eo
    exit
fi
Epa=$1
# -3.74378767003927

if test ! -f millers.out; then
    echo confs not converted, please run cvt_confs.sh first
    exit
fi

test ! -f lmp; mkdir -p lmp
cd lmp
ln -sf ../*pb .
ln -sf ../potential.mod .
cd ..

cwd=$PWD; 
BASEDIR=$(dirname "$0")
cd $BASEDIR; BASEDIR=$PWD; 
cd $cwd

for ii in `cat millers.out`
do
    cd lmp
    ln -sf ../confs/$ii.lmp conf.lmp
    sed "s/EPA/$Epa/g" $BASEDIR/in.relax > in.lammps
    ener=`$lmp_cmd -i in.lammps | grep Surface_energy | awk '{print $3}'`
    echo $ii $ener
    cd ..
done



