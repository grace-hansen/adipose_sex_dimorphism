#!/bin/sh
if [ $# -lt 2 ]
  then
    echo "Usage: run_PCA_genotypes.sh <data_source> <subjlist> <suffix> \n
	Sex can be M, F, or B
	Data source can be GTEx_v8, GTEx_v7, or CMC
    "
    exit
fi

datasource=$1
subjlist=$2
suffix=$3
pcanum=3


cd /project2/nobrega/grace/QTL_analyses/eQTL/$suffix
mkdir covariates
if [ $subjlist != "" ]; then
	paste $subjlist $subjlist > subjs.txt
	for i in $(seq 1 22); do
		plink --bfile ~/midway/genos/${datasource}/plink/chr${i}_rsids \
		--keep subjs.txt \
		--pca $pcanum \
		--out /project2/nobrega/grace/QTL_analyses/eQTL/$suffix/covariates/chr${i}_pca_$suffix
		Rscript ~/midway/QTL_analyses/eQTL/scripts/prep_genotype_PCAs.R $suffix $i
		rm subjs.txt
	done
else 
	pass
	for i in $(seq 1 22); do
		plink --bfile ~/midway/genos/${datasource}/plink/chr${i}_rsids \
		--pca $pcanum \
		--out /project2/nobrega/grace/QTL_analyses/eQTL/$suffix/covariates/chr${i}_pca_${suffix}
		Rscript ~/midway/QTL_analyses/eQTL/scripts/prep_genotype_PCAs.R $suffix $i
	done
fi

