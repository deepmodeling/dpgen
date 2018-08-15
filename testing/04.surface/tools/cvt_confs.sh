#!/bin/bash

ovito_cvt=/home/wanghan/study/deep.md/deep.gen/generator/lib/ovito_file_convert.py
if test ! -f $ovito_cvt; then
    echo no file $ovito_cvt
    exit
fi

miller_idx=""
for ii in cifs/*cif
do
    tmp=`echo $ii | sed 's/.*(\(.*\)).cif/\1/g'`
    miller_idx="$miller_idx $tmp"
    element=`echo $ii | cut -d - -f 1| cut -d / -f 2`
done

echo find element: $element
echo find millers: $miller_idx
echo $miller_idx > millers.out

if test ! -d confs; then
    mkdir confs
fi
rm -f confs/POSCAR*

for ii in $miller_idx
do
    echo miller: $ii
    obabel cifs/$element-\($ii\).cif -o POSCAR  -O confs/POSCAR-$ii
    $ovito_cvt confs/POSCAR-$ii confs/$ii.lmp
done






