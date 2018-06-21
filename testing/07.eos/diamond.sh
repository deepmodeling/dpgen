#!/bin/bash

source env.sh

targets=`seq 4.74 0.02 5.76`

script_file=in.tmp.$$
out_file=tmp.out.$$
for ii in $targets
do
    sed "s/LATT_A/$ii/g" in.diamond > $script_file
    $lmp_cmd -i $script_file &> $out_file
    epa=`grep ENER_PER_ATOM $out_file | awk '{print $2}'`
    epv=`grep VOLM_PER_ATOM $out_file | awk '{print $2}'`
    echo $epv $epa
done

rm -f $script_file $out_file
