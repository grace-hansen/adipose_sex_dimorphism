#!/usr/bin/Rscript
library(data.table)
library(tidyverse)
library(eulerr)
library(yarrr) #For color funs
library(optparse)
pony_colors<-fread("~/medusa/papers/TWAS/pony_palette")
setwd("~/medusa/papers/TWAS/TWAS")

######### Author: Grace Hansen #########
#This script plots a Venn diagram of TWAS results from obesity, WHR_F, and WHR_M.

#Load in obesity genes
obesity_genes<-scan(paste("~/midway/expression/obesity/cortex/results/posthoc/sig_genes",sep=''),what='character',sep='\n')

#Load in WHR_F genes
WHR_F_genes<-scan(paste("~/midway/expression/WHR_F/adipose_subcutaneous/results/posthoc/sig_genes",sep=''),what='character',sep='\n')

#Load in WHR_M genes
WHR_M_genes<-scan(paste("~/midway/expression/WHR_M/adipose_subcutaneous/results/posthoc/sig_genes",sep=''),what='character',sep='\n')

########Plotting#########
#A=GWAS loci
#B=expression TWAS loci
#C=splicing TWAS loci
vals<-c(A=length(obesity_genes),
        B=length(WHR_F_genes),
        C=length(WHR_M_genes),
        "A&B"=16,#"A&B"=length(intersect(obesity_genes,WHR_F_genes)), #These are inflated 2-fold so that the venn diagram is readable
        "A&C"=10,#"A&C"=length(intersect(obesity_genes,WHR_M_genes)),
        "B&C"=22,#"B&C"=length(intersect(WHR_F_genes,WHR_M_genes)),
        "A&B&C"=4)#"A&B&C"=sum(intersect(WHR_F_genes,WHR_M_genes) %in% obesity_genes))


colors<-c(rgb(pony_colors[3,1:3]),rgb(pony_colors[10,1:3]),rgb(pony_colors[14,1:3]))
venn<-euler(vals)
pdf(paste("expression_obesity_WHR_F_WHR_M_venn.pdf",sep=''),width=9,height=9)
plot.new()
plot(venn,fills = colors,labels = FALSE)
#text(0.3,-0.02,"Obesity",cex=2)
#text(0.89,0.2,"WHR, females",cex=2)
#text(0.74,1.,"WHR, males",cex=2)
#text(0.31,0.4,length(obesity_genes)-length(intersect(obesity_genes,WHR_F_genes))-length(intersect(obesity_genes,WHR_M_genes))-sum(intersect(WHR_F_genes,WHR_M_genes) %in% obesity_genes),cex=5) #Just obesity
#text(0.9,0.42,length(WHR_F_genes)-length(intersect(obesity_genes,WHR_F_genes))-length(intersect(WHR_F_genes,WHR_M_genes))-sum(intersect(WHR_F_genes,WHR_M_genes) %in% obesity_genes),cex=5) #Just WHR_F
#text(0.8,0.81,length(WHR_M_genes)-length(intersect(obesity_genes,WHR_M_genes))-length(intersect(WHR_F_genes,WHR_M_genes))-sum(intersect(WHR_F_genes,WHR_M_genes) %in% obesity_genes),cex=5) #Just WHR_M
#text(0.74,0.42,length(intersect(obesity_genes,WHR_F_genes)),cex=5) #A&B
#text(0.63,0.71,length(intersect(obesity_genes,WHR_M_genes)),cex=5) #A&C
#text(0.85,0.58,length(intersect(WHR_F_genes,WHR_M_genes)),cex=5) #B&C
#text(0.68,0.51,sum(intersect(WHR_F_genes,WHR_M_genes) %in% obesity_genes),cex=5) #A&B&C
dev.off()

######### Label swap permutation test ################
# Steps:
# 1) Create a blank list of N rows and 6 columns
# 2) Take all significant genes, put them in a list (no duplicates)
# 3) For each row:
#       a) pull "obesity" genes (same number of genes as in obesity)
#       b) pull "WHR_F" genes ("" "")
#       c) pull "WHR_M" genes ("" "")
#       d) count overlaps: how many are unique to each set, how many are in A&B, A&C, B&C, A&B&C
#       e) put overlaps into row of list
N=100000
N_obesity_genes<-length(obesity_genes)
N_WHR_F_genes<-length(WHR_F_genes)
N_WHR_M_genes<-length(WHR_M_genes)

all_genes<-unique(c(obesity_genes,WHR_F_genes,WHR_M_genes))
perm_dat<-as.data.frame(matrix(nrow=N,ncol=7),stringsAsFactors = FALSE)
for (i in 1:nrow(perm_dat)) {
        obesity<-sample(all_genes,N_obesity_genes)
        WHR_F<-sample(all_genes,N_WHR_F_genes)
        WHR_M<-sample(all_genes,N_WHR_M_genes)
        obesity[1:10]
        AB<-length(intersect(obesity,WHR_F))
        AC<-length(intersect(obesity,WHR_M))
        BC<-length(intersect(WHR_M,WHR_F))
        ABC<-sum(intersect(WHR_F,WHR_M) %in% obesity)
        A<-length(obesity)-AB-AC+ABC
        B<-length(WHR_F)-AB-BC+ABC
        C<-length(WHR_M)-AC-BC+ABC
        perm_dat[i,]<-c(A,B,C,AB,AC,BC,ABC)
}
colnames(perm_dat)<-c("A","B","C","AB","AC","BC","ABC")
AB<-length(intersect(obesity_genes,WHR_F_genes))
AC<-length(intersect(obesity_genes,WHR_M_genes))
BC<-length(intersect(WHR_M_genes,WHR_F_genes))
ABC<-sum(intersect(WHR_F_genes,WHR_M_genes) %in% obesity_genes)
A<-length(obesity_genes)-AB-AC+ABC
B<-length(WHR_F_genes)-AB-BC+ABC
C<-length(WHR_M_genes)-AC-BC+ABC

# Empirical p-values
1-sum(perm_dat$AB>AB)/length(perm_dat$AB)
1-sum(perm_dat$AC>AC)/length(perm_dat$AC)
sum(perm_dat$BC>BC)/length(perm_dat$BC)

