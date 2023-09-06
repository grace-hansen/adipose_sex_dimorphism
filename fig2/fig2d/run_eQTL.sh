#!/bin/bash

if [ $# -ne 4 ]
  then
    echo "Usage: run_eQTL.sh <tissue> <include> <suffix> <data_source>"
    exit
fi

tissue=$1
include=$2
suffix=$3
data_source=$4

cd /project2/nobrega/grace/QTL_analyses/eQTL/$suffix


# Copy expression data over
#if [[ $tissue != "cortex" ]]; then
#	mkdir exp
#	/project2/nobrega/grace/QTL_analyses/eQTL/scripts/prep_exp.py $tissue $suffix
#else
#	exp="To be determined later"
#fi

echo `pwd`

if [[ $include != "" ]]; then
	cut -f1 $include > subjs.txt
	for i in $(seq 1 22); do
		mkdir chr$i
		~/bin/FastQTL/bin/fastQTL.static \
		--vcf /project2/nobrega/grace/genos/${data_source}/vcf/chr${i}_rsid.vcf.gz \
		--include-samples subjs.txt \
		--bed exp/${tissue}_chr${i}_exp.bed.gz \
		--cov /project2/nobrega/grace/QTL_analyses/eQTL/$suffix/covariates/chr${i}_${suffix}_PCs \
		--window 500000 \
		--out chr$i/chr$i \
		--commands 10 chr$i/chr$i
	
		cat ../scripts/eQTL_example.sbatch | sed "s|<chr>|$i|g" > eQTL_${i}.sbatch
	
		cat chr$i/chr$i >> eQTL_${i}.sbatch
		rm chr$i/chr$i
	
		#sbatch eQTL_${i}.sbatch
	done
else 
	for i in $(seq 1 22); do
		mkdir chr$i
		~/bin/FastQTL/bin/fastQTL.static \
		--vcf /project2/nobrega/grace/genos/${data_source}/vcf/chr${i}_rsid.vcf.gz \
		--bed exp/${tissue}_chr${i}_exp.bed.gz \
		--cov /project2/nobrega/grace/QTL_analyses/eQTL/$suffix/covariates/chr${i}_${suffix}_PCs \
		--window 500000 \
		--out chr$i/chr$i \
		--commands 10 chr$i/chr$i
	
		cat ../scripts/eQTL_example.sbatch | sed "s|<chr>|$i|g" > eQTL_${i}.sbatch
	
		cat chr$i/chr$i >> eQTL_${i}.sbatch
		rm chr$i/chr$i
	
		#sbatch eQTL_${i}.sbatch
	done
fi
