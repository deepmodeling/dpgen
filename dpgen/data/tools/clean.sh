#!/bin/bash

if test $# -ne 1; then
    echo usage
    echo $0 nframe
    exit 1
fi
nframes=$1

for ii in [0-9]*[0-9];
do
    ll=`grep -i 'total-force' $ii/OUTCAR | wc -l`;
    echo $ii $ll;
    if test $ll -ne $nframes ; then
	mv $ii bk.$ii
    fi
done

