#!/usr/bin/python
import os, argparse
import pandas as pd
parser = argparse.ArgumentParser()
parser.add_argument("tissue", help="tissue being analyzed")
parser.add_argument("prefix", help="prefix of eQTL study")
args = parser.parse_args()
#############################################################
tissue=args.tissue
prefix=args.prefix

###################### Load gtf #############################
if tissue != "cortex":
    exp=pd.read_csv("/project2/nobrega/grace/expression/GTEx_Analysis_v8_eQTL_expression_matrices/%s.v8.normalized_expression.bed.gz"%tissue,sep='\t',compression='gzip')
    gtf=pd.read_csv("/project2/nobrega/grace/gencode.v29.annotation_hg38.gtf.gz",sep='\t',compression='gzip',header=None,skiprows=range(0,5))
    gtf.columns=['chr','source','type','start','stop','dot','strand','zeroes','info']
    gtf['ENSG']=[f.split('"')[1] for f in [e[0] for e in gtf['info'].str.split(';')]]
    gtf=gtf[gtf['type']=="gene"]
    gtf=gtf[gtf['info'].str.contains('protein_coding')]

    #get gene name from info
    info=gtf['info'].str.split(';')
    genes=list()
    for line in info:
        ind=[line.index(i) for i in line if 'gene_name' in i]
        if len(ind)==0:
            gene_name="NA"
        else:
            gene_name=line[ind[0]].split('"')[1]
        genes.append(gene_name)

    gtf['gene']=genes
    gtf = gtf[gtf.gene != "NA"]

    #Merge GTF data with exp
    gtf=gtf[['gene','ENSG']]
    merged=pd.merge(gtf,exp,left_on="ENSG",right_on="gene_id")
    merged=merged.drop(columns=['ENSG','gene_id'])

    #Rearrange columns
    cols=merged.columns
    gene=cols[0]
    front_cols=list(cols[1:4])
    front_cols.append(gene)
    other_cols=list(cols[5:len(cols)])
    cols=front_cols+other_cols
    merged=merged[cols]

else: 
    exp=pd.read_csv("/project2/nobrega/grace/expression/obesity/cortex/COMMONMIND.RNA.rank_norm.PEER.gexp",sep=' ')

    exp_ids=list(exp.columns) #Add genotype IDs
    geno_ids=pd.read_csv("/project2/nobrega/grace/genos/CMC/COMMONMIND.update_id",sep=' ') 
    geno_ids=geno_ids[['Genotyping_Sample_ID',"Individual_ID"]]
    sort = [s for s in exp_ids if s in geno_ids.Individual_ID.unique()]
    geno_ids = geno_ids.set_index('Individual_ID').loc[sort].reset_index()
    exp.columns=list(geno_ids['Genotyping_Sample_ID'])

    exp['gene']=exp.index.values

    gtf=pd.read_csv("/project2/nobrega/grace/gencode.v19.annotation_hg19.gtf.gz",sep='\t',compression='gzip',header=None,skiprows=range(0,5))
    gtf.columns=['#chr','source','type','start','end','dot','strand','zeroes','info']
    gtf=gtf[gtf['type']=="gene"]
    gtf=gtf[gtf['info'].str.contains('protein_coding')]

    #get gene name from info
    info=gtf['info'].str.split(';')
    genes=list()
    for line in info:
        ind=[line.index(i) for i in line if 'gene_name' in i]
        if len(ind)==0:
            gene_name="NA"
        else:
            gene_name=line[ind[0]].split('"')[1]
        genes.append(gene_name)

    gtf['gene']=genes
    gtf = gtf[gtf.gene != "NA"]

    #Merge GTF data with exp
    gtf=gtf[['#chr','start','end','gene']]
    merged=pd.merge(gtf,exp,on="gene")

for i in range(1,23):
    merged_chr=merged[merged['#chr']=='chr%s'%i]
    merged_chr=merged_chr.sort_values(by=['start'])
    merged_chr=merged_chr.drop_duplicates(subset=['gene'],keep=False)
    merged_chr.to_csv("/project2/nobrega/grace/QTL_analyses/eQTL/%s/exp/%s_chr%s_exp.bed"%(prefix,tissue,i),sep='\t',index=False)
    os.system("bgzip /project2/nobrega/grace/QTL_analyses/eQTL/%s/exp/%s_chr%s_exp.bed"%(prefix,tissue,i))
    os.system("tabix -p bed /project2/nobrega/grace/QTL_analyses/eQTL/%s/exp/%s_chr%s_exp.bed.gz"%(prefix,tissue,i))
