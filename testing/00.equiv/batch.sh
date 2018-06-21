#!/bin/bash

lmp_cmd=~/SCR/wanghan/local/bin/lmp_serial_090_gpu
cwd=`pwd`
target_dir=/scratch/gpfs/linfengz/wanghan/deepgen/
cd $target_dir
job_dir=`ls | grep ^[0-9]`

for ii in $job_dir
do
    cd $ii
    job_idx=`echo $ii | cut -d '.' -f 1`
    iter_dir=`ls | grep ^iter `
    if test ! -d $cwd/$job_idx; then
	mkdir -p $cwd/$job_idx
    fi
    for jj in $iter_dir
    do
	cd $jj
	iter_idx=`echo $jj | cut -d '.' -f 2`
	echo $job_idx $iter_idx `pwd`
	tmp_cwd=`pwd`
	cd $cwd/$job_idx
	if test ! -f in.equiv; then
	    ln -s ../in.equiv .
	fi
	if test ! -f potential.mod; then
	    ln -s ../potential.mod .
	    ln -s ../init.mod .
	    ln -s ../relax.mod .
	fi
	if test -f $tmp_cwd/00.train/graph.000.pb; then
	    if test ! -f $job_idx-$iter_idx.log; then
		ln -sf $tmp_cwd/00.train/*pb .
		$lmp_cmd -i in.equiv > $job_idx-$iter_idx.log
	    fi
	fi
	cd $tmp_cwd
	cd ..
    done
    cd ..
done
    

