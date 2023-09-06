#!/usr/bin/Rscript
args=commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("USAGE: prep_genotype_PCAs <prefix> <chr> \n", call.=FALSE)
}
library(tidyverse)

prefix=args[1]
chr=args[2]

dat<-read_table2(paste("~/midway/QTL_analyses/eQTL/",prefix,"/covariates/chr",chr,"_pca_",prefix,".eigenvec",sep=''),col_names=FALSE)
dat<-t(dat)
dat<-dat[-1,]
rownames(dat)<-c("ID",1:(nrow(dat)-1))

write.table(dat,paste("~/midway/QTL_analyses/eQTL/",prefix,"/covariates/chr",chr,"_",prefix,"_PCs",sep=''),sep='\t',
            col.names=FALSE,quote=FALSE)
