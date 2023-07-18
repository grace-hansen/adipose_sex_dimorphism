#!/usr/bin/Rscript
args=commandArgs(trailingOnly=TRUE)
library(tidyverse)
library(xCell)
library(data.table)


dat<-fread(args[1],data.table=FALSE)
rownames(dat)<-make.names(dat[,1],unique=TRUE)
dat<-dat[,-1]

working_dir<-dirname(args[1])
prefix<-strsplit(basename(args[1]),'[.]')[[1]][1]

print(working_dir)
print(prefix)
if (len(args)>1) {
  celltypes<-strsplit(args[2],',')[[1]]
} else {
  celltypes<-c("Adipocytes","Epithelial cells","Hepatocytes","Keratinocytes","Myocytes","Neurons","Neutrophils")
}

out<-xCellAnalysis(dat)
out<-out[which(rownames(out)%in%celltypes),]
write.table(out,paste(working_dir,"/",prefix,"_results_7celltypes.txt",sep=''),quote=FALSE,sep='\t')