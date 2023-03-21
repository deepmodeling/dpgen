#!/bin/bash

for ii in *; do ll=`grep -i 'total-force' $ii/OUTCAR | wc -l`; echo $ii $ll; done
