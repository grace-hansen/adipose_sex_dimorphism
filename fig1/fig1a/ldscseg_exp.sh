#! /bin/bash

cd ~/midway/ldsc_seg

trait=$1

mkdir -p output/$trait/GTEx_exp

N_tests=$(ls annot/Multi_tissue_gene_expr_1000Gv3_ldscores/* | cut -f3 -d'/' | cut -f1-2 -d'.' | sort -u | wc -l)
bf_thresh=$(awk 'BEGIN {print 0.05/"'$N_tests'"}')
echo $bf_thresh > output/$trait/GTEx_exp/bonferroni_threshold


ls annot/Multi_tissue_gene_expr_1000Gv3_ldscores/* | cut -f3 -d'/' | cut -f1-2 -d'.' | sort -u | sed '$ d' | while read line; do
    tissue=$(grep "/${line}[.]" annot/Multi_tissue_gene_expr_GTEx.ldcts | cut -f1)
    input=$(grep "/${line}[.]" annot/Multi_tissue_gene_expr_GTEx.ldcts | cut -f2 | sed "s|Multi|annot/Multi|g")
    python2 scripts/ldsc/ldsc.py \
    --h2 ${trait}_GWAS_sumstats \
    --ref-ld-chr $input \
    --w-ld-chr weights_hm3_no_hla/weights. \
    --overlap-annot \
    --frqfile-chr 1000G_Phase3_frq/1000G.EUR.QC. \
    --out output/$trait/GTEx_exp/${tissue}
done

head=$(ls annot/Multi_tissue_gene_expr_1000Gv3_ldscores/* | cut -f3 -d'/' | cut -f1-2 -d'.' | sort -u | head -1)
tissue=$(grep "/${head}[.]" annot/Multi_tissue_gene_expr_GTEx.ldcts | cut -f1)
echo -e "Tissue\t$(head -1 output/$trait/GTEx_exp/$tissue.results)" > output/$trait/GTEx_exp/all_tissues.results

ls annot/Multi_tissue_gene_expr_1000Gv3_ldscores/* | cut -f3 -d'/' | cut -f1-2 -d'.' | sort -u | sed '$ d' | while read line; do
    tissue=$(grep "/${line}[.]" annot/Multi_tissue_gene_expr_GTEx.ldcts | cut -f1)
    echo -e "${tissue}\t$(head -2 output/$trait/GTEx_exp/$tissue.results | tail -1)"  >> output/$trait/GTEx_exp/all_tissues.results
done

