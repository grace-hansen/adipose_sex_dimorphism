import os, gzip
import pandas as pd


##############Make plinks and compute_TWAS_weights#################
#This script combines GTex data with obesity cortex data so that xCell can be run on both together.
###################################################################
CMC=pd.read_csv("~/midway/expression/obesity/cortex/cortex_CMC_tpm.txt.gz",compression="gzip",sep='\t')
CMC=CMC.drop(['coord'], axis=1)
CMC=CMC.drop(['chr'], axis=1)


#tissues=["Adipose_Subcutaneous","Adipose_Visceral_Omentum","Brain_Amygdala",
#           "Brain_Anterior_cingulate_cortex_BA24","Brain_Caudate_basal_ganglia",
#           "Brain_Cerebellar_Hemisphere","Brain_Cerebellum","Brain_Cortex",
#           "Brain_Frontal_Cortex_BA9","Brain_Hippocampus","Brain_Hypothalamus",
#           "Brain_Nucleus_accumbens_basal_ganglia","Brain_Putamen_basal_ganglia","Brain_Spinal_cord_cervical_c-1",
#           "Brain_Substantia_nigra","Breast_Mammary_Tissue","Colon_Transverse","Esophagus_Mucosa",
#           "Heart_Atrial_Appendage","Heart_Left_Ventricle","Kidney_Cortex",
#           "Liver","Lung","Minor_Salivary_Gland","Muscle_Skeletal",
#           "Pancreas","Pituitary","Prostate",
#           "Skin_Sun_Exposed_Lower_leg","Skin_Not_Sun_Exposed_Suprapubic","Small_Intestine_Terminal_Ileum","Stomach",
#           "Thyroid","Uterus","Whole_Blood"]    


#i=1
#for tissue in tissues:
#	dat=pd.read_csv("/project2/nobrega/grace/expression/GTEx_Analysis_v8_eQTL_expression_matrices/%s.v8.normalized_expression.bed.gz"%tissue,sep='\t',compression="gzip")
#	dat=dat.iloc[:,3:]
#	ENSGs=list(dat["gene_id"])
#	ENSGs=[ENSG.split('.')[0] for ENSG in ENSGs]
#	gene_names=[ENSG_lookup.get(ENSG,'.') for ENSG in ENSGs]
#	dat=dat.drop(['gene_id'], axis=1)
#	dat.columns=[ID+'_'+tissue for ID in list(dat.columns)]
#	dat["gene"]=gene_names
#	if i==1:
#		out=pd.merge(cortex,dat,on="gene")
#		i=i+1
#	else:
#		out=pd.merge(out,dat,on="gene")

GTEx=pd.read_csv("/project2/yangili1/GTEx_v8/expression/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz",sep='\t',compression="gzip",skiprows=2)
#Convert gene ID from ENSG to gene symbol
GTEx=GTEx.drop(['Name'], axis=1)
GTEx=GTEx.rename(columns={"Description":"gene"})

#Merge GTEx and cCMC data
out=pd.merge(CMC,GTEx,on="gene")


out.to_csv("~/midway/expression/obesity/cortex/CMC_GTEx_xCell_tpm.txt",sep='\t',index=False)