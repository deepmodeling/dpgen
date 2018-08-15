#!/bin/bash

lmp_cmd=$HOME/local/bin/lmp_mpi_010

if test $# -lt 1; then
    echo usage
    echo $0 Eo [TYPE]
    exit
fi
Epa=$1
if test $# -eq 2; then
    atom_type=$2
else 
    atom_type=0
fi
atom_type=$(($atom_type+1))
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
    ln -sf ../confs/$ii.lmp conf.orig.lmp
    rm -f conf.lmp
    cp -L conf.orig.lmp conf.lmp    
    ../tools/set_type.py conf.lmp $atom_type
    sed "s/EPA/$Epa/g" $BASEDIR/in.relax > in.lammps
    ener=`$lmp_cmd -i in.lammps | grep Surface_energy | awk '{print $3}'`
    echo $ii $ener
    cd ..
done



