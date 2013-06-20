#!/bin/bash
BIN=./FTLE_CPU_x86_64
DIR=confs/T2
pwd
list=`ls $DIR/*.cfg`
echo $list
for f in $list; do
cp $f Configuration.cfg
$BIN
rm Configuration.cfg
done
