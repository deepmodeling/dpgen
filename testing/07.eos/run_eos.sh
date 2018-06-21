#!/usr/bin/env bash

#EXE=lmp_mpi
EXE=lmp_serial

$EXE -i in.bulk > eos.log
#$EXE -i in.bulk.pv > eos.log

grep "@@@" eos.log | awk -F ":" '{print $NF}' > abc.dat
grep "VVV" eos.log | awk '{print $3,$4}' > ve.dat

rm eos.log
cat ve.dat
