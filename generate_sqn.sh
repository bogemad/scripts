#!/bin/bash

workdir=$(readlink -f $1)
sbt=$(readlink -f $2)
species="$3"
isolate="$4"
strain="$5"

cd $workdir
prefix=$(echo *.tbl | sed 's/.tbl//g')
mkdir -p tmp
cp $prefix.tbl tmp/genome.tbl
cp $prefix.fsa tmp/genome.fsa
tbl2asn -p tmp/ -t $sbt -M n -Z $prefix.discrepency.report.txt -j "[organism=$species] [isolate=$isolate] [strain=$strain]" -V vb -c fx
mv tmp/errorsummary.val $prefix.errorsummary.val
mv tmp/genome.sqn $prefix.sqn
mv tmp/genome.val $prefix.val
mv tmp/genome.gbf $prefix.gbk
rm -r tmp
