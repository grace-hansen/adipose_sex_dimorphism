#!/bin/bash

if [ $# -ne 1 ]
  then
    echo "Usage: collect_eQTL.sh <suffix>"
    exit
fi

suffix=$1


cd /project2/nobrega/grace/QTL_analyses/eQTL/$suffix
mkdir results

for i in $(seq 1 22); do
	cd chr$i
	cat chr$i.chr$i* > chr${i}_eQTL.txt
	mv chr${i}_eQTL.txt ../results
	#rm chr$i.chr$i*
	cd ..
done
