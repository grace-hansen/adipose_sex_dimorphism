cd /project2/yangili1/grace/expression/WHR

python2 ~/midway/expression/scripts/ldsc/munge_sumstats.py \
--out ~/midway/expression/WHR/WHR_GWAS_ \
--N-col N \
--a1  Tested_Allele \
--a2 Other_Allele \
--sumstats /project2/yangili1/grace/expression/WHR/Whradjbmi.giant-ukbb.meta-analysis.combined.23May2018.HapMap2_only.txt.gz \
--p P > ~/midway/expression/WHR/prep_GWAS.log

gunzip WHR_GWAS_.sumstats.gz
mv WHR_GWAS_.sumstats WHR_GWAS_sumstats

cut -f1 WHR_GWAS_sumstats | cut -f1 -d':' > col1
cut -f2-5 WHR_GWAS_sumstats > cols2-5
paste col1 cols2-5 > temp
mv temp WHR_GWAS_sumstats