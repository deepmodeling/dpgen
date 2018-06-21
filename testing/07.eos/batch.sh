#!/bin/bash

lmp_cmd=~/SCR/wanghan/local/bin/lmp_serial_090_gpu

if test $# -ne 1; then
    echo usage
    echo $0 crystal
    exit
fi
crystal=$1

cwd=`pwd`
target_dir=/scratch/gpfs/linfengz/wanghan/deepgen/
cd $target_dir
job_dir=`ls | grep ^[0-9]`

for ii in $job_dir
do
    cd $ii
    job_idx=`echo $ii | cut -d '.' -f 1`
    iter_dir=`ls | grep ^iter `
    if test ! -d $cwd/$job_idx/$crystal; then
	mkdir -p $cwd/$job_idx/$crystal
    fi
    for jj in $iter_dir
    do
	cd $jj
	iter_idx=`echo $jj | cut -d '.' -f 2`
	echo $job_idx/$crystal $iter_idx `pwd`
	tmp_cwd=`pwd`
	cd $cwd/$job_idx/$crystal
	if test ! -f in.$crystal; then
	    ln -sf ../../in.$crystal .
	    ln -sf ../../$crystal.sh .
	fi
	if test ! -f env.sh; then
	    ln -sf ../../env.sh .
	fi
	if test ! -f potential.mod; then
	    ln -sf ../../potential.mod .
	    ln -sf ../../relax.mod .
	fi
	if test -f $tmp_cwd/00.train/graph.000.pb; then
	    if test ! -f eos.$job_idx-$iter_idx.log; then
		ln -sf $tmp_cwd/00.train/*pb .
		./$crystal.sh > eos.$job_idx-$iter_idx.log
	    fi
	fi
	cd $tmp_cwd
	cd .. # cd $jj
    done
    cd .. # cd $ii
done
    

